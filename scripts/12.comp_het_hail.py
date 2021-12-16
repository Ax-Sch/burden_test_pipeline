# 12.comp_het_hail.py
# This script performs gene based collapsing for recessive variants - a gene in a person is counted if it contains two qualifying heterozygous variants that are more than 50bp apart or if it contains a homozygous qualifying variant.
# Make sure python, pandas and hail are installed (tested with python 3.6.11, hail 0.2.58-3f304aae6ce2 and pandas 1.1.4 on a local machine).

#paths - please set them. The naming should be compatible with the script 11.regenieAnalysis.sh
#path to QCed plink binaries AND TO THE PLINK FILES THAT UNDERWENT FILTERING WITH THE COMMON LIST (see 11.regenieAnalysis.sh)
#pathPlink='/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/sequence.file.normID.noChrM.bcftools_haileur.vcf.gz'
#path to a temp-folder
pathTmp='/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/tmp/'
#path to regenie inputs, where the mask.def, set.list, anno.file, and aaf.file are located. mask.def can be found on the git
pathReg="/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/"
# path where outputs will be written to (outputs plink files and a file with all variants potentially compatible with recessive inheritance ( ending: _rez_genes.txt) )
pathPlinkCollapsed='/media/axel/Dateien/Arbeit_Gen/covid_resequencing_check_gwas_datafreeze/analysis_pipeline/cp_het/'
# set the reference genome as needed:
ref_genome="GRCh37"

import pandas as pd
import hail as hl
# make sure hail can run properly by allocating enough RAM e.g. here on a local machine (change the 20g as needed/available, set local[n] with n=number of cores). The script is rather memory intensive.
hl.init(spark_conf={'spark.driver.memory': '20g'}, min_block_size=128, master='local[4]',tmp_dir= pathTmp) # initialize hail

#read mask/annotation/AAF file 
mask_def=pd.read_csv(pathReg + "regenie.mask.def.txt", delim_whitespace=True, names=["mask","funs"])
annotation=pd.read_csv(pathReg + "regenie.anno.file.txt", sep=' ', names=["snp","ENSG","fun"])
aaf_file=pd.read_csv(pathReg + "regenie.aaf.file.txt", sep=' ', names=["snp","AF"])

def build_annotation(AFs, masks):
	relevant_snps=pd.DataFrame()
	for mask_var in masks:
		for AF_var in AFs:
			print(AF_var, mask_var)
			relevant_mask=mask_def[mask_def["mask"]==mask_var]
			relevant_funs=relevant_mask.funs.tolist()[0]
			relevant_funs=relevant_funs.split(",")
			relevant_annotation=annotation[annotation.fun.isin(relevant_funs)].copy()
			relevant_aaf=aaf_file[aaf_file["AF"]<AF_var].copy()
			relevant_snps_temp=relevant_annotation.loc[relevant_annotation.snp.isin(relevant_aaf["snp"]),:].copy()
			relevant_snps_temp["ensg_mask_af"]=relevant_snps_temp["ENSG"] + "." + str(mask_var) + "." + str(AF_var)
			relevant_snps=relevant_snps.append(relevant_snps_temp)
			print(relevant_snps_temp)
	relevant_snps_grouped=relevant_snps.groupby(['snp']).ensg_mask_af.apply(list).reset_index()
	t = hl.Table.from_pandas(relevant_snps_grouped)
	t=t.key_by("snp")
	return(t)

def annotate_mt(mt,t,pathTmp): # annotate MatrixTable with gene names ("genecol") of relevant SNPs
	mt = mt.key_rows_by('rsid')
	mt = mt.annotate_rows(genecol=t[mt.rsid].ensg_mask_af)
	mt = mt.explode_rows("genecol")
	mt = mt.key_rows_by(mt.locus, mt.alleles)
	mt = hl.variant_qc(mt, name='variant_qc')
	mt = mt.filter_rows((hl.is_defined(mt.genecol)))	#   & (mt.variant_qc.call_rate>0.9)
	mt.write(pathTmp+'annotate_tmp.mt',overwrite=True)   # execute annotation / filtering
	mt = hl.read_matrix_table(pathTmp+'annotate_tmp.mt')
	return(mt)

def count_carriers(mt,gene_mt,het_count,totalhom,CT_gene,ref_genome_f): # functions based on Liz's code to collapse to genes
	#het_count, totalhom, CT_gene = 'mac','totalhom',genecol
	gene_mt=gene_mt.annotate_rows(alleles=["C","A"])
	#now pick a variant belonging to that gene
	mtr=mt.rows()
	mtr=mtr.key_by() # key removed
	mtr=mtr.key_by(CT_gene)
	mtr=mtr.select('locus')
	gene_mt=gene_mt.annotate_rows(**mtr[gene_mt[CT_gene]])
	gene_mt=gene_mt.key_rows_by()
	#regenie doesn't like chrY or chrM, so I change it to chr1 for now
	gene_mt=gene_mt.annotate_rows(locus=hl.case(missing_false=True).when((gene_mt.locus.contig=="chrY")|(gene_mt.locus.contig=="chrM"),hl.locus("chr1",gene_mt.locus.position, reference_genome=ref_genome_f )).default(gene_mt.locus))
	gene_mt=gene_mt.key_rows_by('locus','alleles')
	gene_mt=gene_mt.annotate_entries(GT=hl.case(missing_false=True)
		                 .when(gene_mt[totalhom]>0,hl.call(1,1))
		                 .when(gene_mt[het_count]==0,hl.call(0,0))
		                 .when(gene_mt[het_count]==1,hl.call(1,0))
		                 .when((gene_mt[het_count]>1)&(gene_mt.diff>=50),hl.call(1,1))
		                 .when(gene_mt[het_count]>1,hl.call(1,0))
		                 .default(hl.call(0,0)))
	gene_mt=gene_mt.annotate_entries(CP=hl.case(missing_false=True)
		                 .when((gene_mt[het_count]>1)&(gene_mt.diff>=50),1)
		                 .default(0))
	gene_mt=gene_mt.key_rows_by(CT_gene)
	return(gene_mt)#here is a full collapsing function:


def collapse_cp_het(mt,location,filename,genecol,ref_genome_f, pathTmp):	#from forum: http://discuss.hail.is/t/logistic-regression-burden-tests/206/3
	gene_mt = mt.group_rows_by(mt[genecol]).aggregate(
		het_count = hl.agg.sum(hl.cond((mt.GT.is_het()) & (mt.GT.is_diploid()),1,0)),
		totalhom=hl.agg.sum( 
		hl.cond(((mt.variant_qc.AF[1] <= 0.5)|(hl.is_missing(mt.variant_qc.AF[1])))&((mt.GT.is_hom_var())|((mt.GT.is_haploid())&(mt.GT.n_alt_alleles()>0))),1,
		hl.cond((mt.variant_qc.AF[1] > 0.5)&((mt.GT.is_hom_ref())|((mt.GT.is_haploid())&(mt.GT.n_alt_alleles()==0))),1,0))),
		diff=hl.agg.max(hl.cond(mt.GT.is_het(),mt.locus.position,hl.null(hl.tint32)))-hl.agg.min(hl.cond(mt.GT.is_het(),mt.locus.position,hl.null(hl.tint32))))
	gene_mt=count_carriers(mt,gene_mt,'het_count','totalhom',genecol,ref_genome_f)
	gene_mt.write(pathTmp+filename+'tmp.gene_mt',overwrite=True)    #export to plink
	gene_mt=hl.read_matrix_table(pathTmp+filename+'tmp.gene_mt') 
	mt2=mt.annotate_entries(REZ=gene_mt[mt.genecol,mt.s].GT)
	comp_gene_sample_ids=mt2.aggregate_entries(hl.agg.filter((mt2.entry.REZ==hl.call(1,1)) & (mt2.entry.GT!=hl.call(0,0)), hl.agg.collect([mt2.genecol, mt2.s, hl.str(mt2.locus), hl.str(mt2.GT)])))
	with open(location+filename +'_rez_genes.txt', 'w') as filehandle:
		filehandle.writelines("%s\n" % variant for variant in comp_gene_sample_ids)
	gene_mt=gene_mt.key_rows_by().select_rows(genecol,'locus','alleles').key_rows_by('locus','alleles')
	gene_mt=gene_mt.select_entries('GT')
	hl.export_plink(gene_mt,location+filename,call=gene_mt.GT,fam_id=gene_mt['s'],ind_id=gene_mt['s'],varid=gene_mt[genecol])



def count_carriers_dom_rez(mt,gene_mt,het_count,totalhom,CT_gene,ref_genome_f): # functions based on Liz's code to collapse to genes
	#het_count, totalhom, CT_gene = 'mac','totalhom',genecol
	gene_mt=gene_mt.annotate_rows(alleles=["C","A"])
	#now pick a variant belonging to that gene
	mtr=mt.rows()
	mtr=mtr.key_by() # key removed
	mtr=mtr.key_by(CT_gene)
	mtr=mtr.select('locus')
	gene_mt=gene_mt.annotate_rows(**mtr[gene_mt[CT_gene]])
	gene_mt=gene_mt.key_rows_by()
	#regenie doesn't like chrY or chrM, so I change it to chr1 for now
	gene_mt=gene_mt.annotate_rows(locus=hl.case(missing_false=True).when((gene_mt.locus.contig=="chrY")|(gene_mt.locus.contig=="chrM"),hl.locus("chr1",gene_mt.locus.position, reference_genome=ref_genome_f )).default(gene_mt.locus))
	gene_mt=gene_mt.key_rows_by('locus','alleles')
	gene_mt=gene_mt.annotate_entries(GT=hl.case(missing_false=True)
		                 .when(gene_mt[totalhom]>0,hl.call(1,1))
		                 .when(gene_mt[het_count]>0,hl.call(1,0))
		                 .default(hl.call(0,0)))
	gene_mt=gene_mt.key_rows_by(CT_gene)
	return(gene_mt)#here is a full collapsing function:

def collapse_dom_rez(mt,location,filename,genecol,ref_genome_f, pathTmp):	#from forum: http://discuss.hail.is/t/logistic-regression-burden-tests/206/3
	gene_mt = mt.group_rows_by(mt[genecol]).aggregate(
		het_count = hl.agg.sum(hl.cond((mt.GT.is_het()) & (mt.GT.is_diploid()),1,0)),
		totalhom=hl.agg.sum( 
		hl.cond(((mt.variant_qc.AF[1] <= 0.5)|(hl.is_missing(mt.variant_qc.AF[1])))&((mt.GT.is_hom_var())|((mt.GT.is_haploid())&(mt.GT.n_alt_alleles()>0))),1,
		hl.cond((mt.variant_qc.AF[1] > 0.5)&((mt.GT.is_hom_ref())|((mt.GT.is_haploid())&(mt.GT.n_alt_alleles()==0))),1,0))))
	gene_mt=count_carriers_dom_rez(mt,gene_mt,'het_count','totalhom',genecol,ref_genome_f)
	gene_mt.write(pathTmp+filename+'tmp.gene_mt',overwrite=True)    #export to plink
	gene_mt=hl.read_matrix_table(pathTmp+filename+'tmp.gene_mt') 
	gene_mt=gene_mt.key_rows_by().select_rows(genecol,'locus','alleles').key_rows_by('locus','alleles')
	gene_mt=gene_mt.select_entries('GT')
	hl.export_plink(gene_mt,location+filename,call=gene_mt.GT,fam_id=gene_mt['s'],ind_id=gene_mt['s'],varid=gene_mt[genecol])





# collapse based on Plink files generated only with gnomAD 
masks=["M1","M3","M4","M0"]
AFs=[0.01]
t=build_annotation(AFs,masks) # build a hail Matrix containing SNP/gene/Mask/AF combinations 
mt = hl.import_plink(bed='filtered_maf001.bed',
	bim='filtered_maf001.bim',
	fam='filtered_maf001.fam',
	reference_genome=ref_genome)
mt=annotate_mt(mt,t,pathTmp) # annotate MatrixTable
collapse_cp_het(mt, pathPlinkCollapsed, "recessive_collapsed_comp_gnomAD001", "genecol",ref_genome, pathTmp) # run collapse-function
collapse_dom_rez(mt, pathPlinkCollapsed, "dom_rez_gnomAD001", "genecol",ref_genome, pathTmp) # run collapse-function


# collapse based on Plink files generated only with gnomAD 
masks=["M1","M3","M4","M0"]
AFs=[0.001]
t=build_annotation(AFs,masks) # build a hail Matrix containing SNP/gene/Mask/AF combinations 
mt = hl.import_plink(bed='filtered_maf0001.bed',
	bim='filtered_maf0001.bim',
	fam='filtered_maf0001.fam',
	reference_genome=ref_genome)
mt=annotate_mt(mt,t,pathTmp) # annotate MatrixTable
collapse_cp_het(mt, pathPlinkCollapsed, "recessive_collapsed_comp_gnomAD0001", "genecol",ref_genome, pathTmp) # run collapse-function
collapse_dom_rez(mt, pathPlinkCollapsed, "dom_rez_gnomAD0001", "genecol",ref_genome, pathTmp) # run collapse-function


# collapse based on Plink files generated only with gnomAD 
masks=["MF"]
AFs=[109]
t=build_annotation(AFs,masks) # build a hail Matrix containing SNP/gene/Mask/AF combinations 
mt = hl.import_plink(bed='filtered_all.bed',
	bim='filtered_all.bim',
	fam='filtered_all.fam',
	reference_genome=ref_genome)
mt=annotate_mt(mt,t,pathTmp) # annotate MatrixTable
collapse_dom_rez(mt, pathPlinkCollapsed, "functional_collapsed", "genecol",ref_genome, pathTmp) # run collapse-function

hl.stop()
