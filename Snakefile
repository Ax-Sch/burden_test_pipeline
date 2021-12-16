#ancestries=["afr","amr","eas","eur","sas"]
ancestries=["eur"]
#contigs=["1"]
contigs=["1","2","3","4","5","6","8","9","10","11","12","14","15","17","19","20","21","X"]
contigs_wo_X=["1","2","3","4","5","6","8","9","10","11","12","14","15","17","19","20","21"]

rule all:
	input:
		#expand("exclusion_list/common.one.perceur.txt", anc=ancestries),
		#expand("regenieInputseur/regenie.aaf.file.txt", anc=ancestries),
		#"split_annot_06_07/sequence.file.normID.noChrM.bcftools_haileur_finalAnnot.Y.txt",
		#expand("generate_PCA_covar_05eur/rarePCA.eigenval",anc = ancestries),
		#expand("regenieInputseur/regenie.set.list.txt", anc=ancestries),
		#expand("regenie_11/regenieResultseur_all/", anc=ancestries),
		#expand("regenie_11/regenieResultseur_male/", anc=ancestries),
		#expand("regenie_11/regenieResultseur_female/", anc=ancestries),
		"split_annot_06_07/sequence.file.normID.noChrM.bcftools_haileur_finalAnnot.vcf",
		"assoc_results/filtered_all.assoc.fisher",
		#"assoc_results/functional.frq.cc"



rule select_ancestries_04:
	input:
		pathID="input_sample_list.txt",
		pathQC="../data/reseqKiel_all_samples_normed_20210528.vcf.gz",
		exclude_vcf="../indels_list_homo5polymer.vcf.gz"
	output:
		pre_indel_filter="sequence.file.normID.noChrM.bcftools_noIndel_filter.vcf.gz",
		post_indel_filter="sequence.file.normID.noChrM.bcftools_haileur.vcf.gz"
	resources: cpus=1, mem_mb=18000, time_job=720
	params:
		partition='batch'
	shell:
		"""
		bcftools view -S {input.pathID} {input.pathQC} -Ob | bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' -Ob | bcftools sort -Oz > {output.pre_indel_filter}
		
		tabix -p vcf {output.pre_indel_filter}
		
		bcftools isec {output.pre_indel_filter} {input.exclude_vcf} -c none -n~10 -w1 -Ob | bcftools view -e'QD<5 | QUAL<120' -Oz -o {output.post_indel_filter}
		
		tabix -p vcf {output.post_indel_filter}
		
		"""








#05
#05
#05 05.PCA.sh




rule split_annot_06_07:
	input:
		pathAncestry="sequence.file.normID.noChrM.bcftools_haileur.vcf.gz"
	output:
		pathAnnot="split_annot_06_07/sequence.file.normID.noChrM.bcftools_haileur_finalAnnot.vcf",
		pathAnnot_tsv="split_annot_06_07/sequence.file.normID.noChrM.bcftools_haileur_finalAnnot.tsv"
	resources: cpus=1, mem_mb=18000, time_job=720
	params:
		partition='batch',
		pathCache="/media/axel/Seagate Backup Plus Drive/annotation_sources/"
		
	conda:
		"env/covid_qc.yml"
	shell:
		"""
		ln -sf "{params.pathCache}dbNSFP/dbNSFP4.0_hg19.gz" dbNSFP4.0_hg19.gz
		ln -sf "{params.pathCache}dbNSFP/dbNSFP4.0_hg19.gz.tbi" dbNSFP4.0_hg19.gz.tbi
		
		  vep -i {input.pathAncestry} \
		    --plugin dbNSFP,dbNSFP4.0_hg19.gz,Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
		    -everything \
		    --buffer_size 3000 \
		    --force_overwrite \
		    --offline \
		    --vcf \
		    --pick_allele_gene \
		    --dir_cache "{params.pathCache}" \
		    --cache -o STDOUT | \
		    cut -f1-8 > {output.pathAnnot}
		
		bcftools +split-vep -l {output.pathAnnot} | cut -f2 | tr '\\n' ';' | awk 'BEGIN {{ FS = ";"}} ;{{ print "ID;"$0}}' > {output.pathAnnot_tsv}
		bcftools +split-vep -d -f'%CHROM:%POS:%REF:%ALT;%CSQ\n' -A ";" {output.pathAnnot} >> {output.pathAnnot_tsv}
		
		"""


rule get_aaf_files_hail:
	input:
		pathAnnot_tsv="split_annot_06_07/sequence.file.normID.noChrM.bcftools_haileur_finalAnnot.tsv",
		pathAncestry="sequence.file.normID.noChrM.bcftools_haileur.vcf.gz"
	output:
		"cp_het/recessive_collapsed_comp_gnomAD001.fam",
		"cp_het/dom_rez_gnomAD001.fam",
		"cp_het/recessive_collapsed_comp_gnomAD0001.fam",
		"cp_het/dom_rez_gnomAD0001.fam",
		"filtered_maf0001.bed"
	resources: cpus=1, mem_mb=18000, time_job=720
	shell:
		"""
		Rscript scripts/make_aaf.R # 1%  // 0.1% ; LoF; Missense + LoF; synonym;   
		plink --vcf {input.pathAncestry} --max-maf 0.01 --make-bed --double-id --out filtered_maf001
		plink --vcf {input.pathAncestry} --max-maf 0.001 --make-bed --double-id --out filtered_maf0001
		plink --vcf {input.pathAncestry} --double-id --max-maf 0.49 --make-bed --out filtered_all
		python scripts/12.comp_het_hail.py
		Rscript scripts/correct_fam_files.R
		plink --bfile cp_het/dom_rez_gnomAD0001 --bmerge cp_het/dom_rez_gnomAD001 --make-bed --out cp_het/dom_rez_comb
		plink --bfile cp_het/dom_rez_gnomAD0001 --bmerge cp_het/dom_rez_gnomAD001 --make-bed --out cp_het/dom_rez_comb_x_to_autosome
		sed -i "s/^23/22/g" cp_het/dom_rez_comb_x_to_autosome.bim
		
		#merge rec. files
		cp cp_het/recessive_collapsed_comp_gnomAD001* .
		rename "s/gnomAD001/gnomAD001_cpHet/g" recessive_collapsed_comp_gnomAD001* -f
		awk 'BEGIN {{ FS = "\t"}} ;{{ print $1,$2"_cpHet",$3,$4,$5,$6}}' cp_het/recessive_collapsed_comp_gnomAD001.bim > recessive_collapsed_comp_gnomAD001_cpHet.bim
		mv recessive_collapsed_comp_gnomAD001_cpHet* cp_het/
		plink --bfile cp_het/dom_rez_gnomAD001 --bmerge cp_het/recessive_collapsed_comp_gnomAD001_cpHet --make-bed --out cp_het/recessive_fileset_001
		sed -i "s/^23/22/g" cp_het/recessive_fileset_001.bim
		
 		"""

rule plink_testing:
	input:
		plink_cp_het="cp_het/recessive_collapsed_comp_gnomAD001.fam"
	output:
		pathOut="assoc_results/functional.frqx"
	resources: cpus=1, mem_mb=18000, time_job=720
	params:
		plink_cp_het_001="cp_het/recessive_collapsed_comp_gnomAD001",
		plink_dom_rez_comb="cp_het/dom_rez_comb",
		plink_dom_rez_comb_xtoauto="cp_het/dom_rez_comb_x_to_autosome",
		plink_dom_rez_gnomAD001="cp_het/dom_rez_gnomAD001",
		plink_cp_het_0001="cp_het/recessive_collapsed_comp_gnomAD0001",
		functional="cp_het/functional_collapsed",
		covariates="covariates.cov",
		plink_rec_set="cp_het/recessive_fileset_001"

	conda:
		"env/covid_qc.yml"
	shell:
		"""
		######
		###### fisher
		# Idea: primary analysis: Fisher's exact test, additive model. 
		# then additional analyses, that compound these analyses: Recessive, Males only, functional, (below 60 years)
		# additive; alpha: 0,05  /  3 * (genes) * 2 (AF) -> ca. 0,0001 (no p value below this threshold); if we say tests are strongly correlated: 0.05 / 100 -> also no value below this threshold
		plink --bfile {params.plink_dom_rez_comb_xtoauto} --assoc fisher counts --out assoc_results/additive_fisher_comb --mperm 1000
		
		# recessive !!! xtoauto was used here
		plink --bfile {params.plink_dom_rez_comb_xtoauto} --model fisher rec --cell 0 --out assoc_results/recessive_fisher_comb --mperm 1000
		
		# functional, do in recessive mode?
		plink --bfile {params.functional} --assoc fisher counts --out assoc_results/functional --mperm 1000
		
		# males
		# functional / males
		plink --bfile {params.functional} --assoc fisher counts --filter-males --cell 0 --out assoc_results/functional_males --mperm 1000
		
		# additive / males
		plink --bfile {params.plink_dom_rez_comb} --assoc fisher counts --filter-males --chr X --out assoc_results/additive_males_comb_fisher --mperm 1000

		
		
		######
		###### firth
		# additive
		plink2 --bfile {params.plink_dom_rez_comb} --glm firth hide-covar --keep italian_participants.dat --covar covariates.cov --covar-variance-standardize --out assoc_results/additive_firth_comb_italian
		
		plink2 --bfile {params.plink_dom_rez_comb} --glm firth hide-covar --remove italian_participants.dat --covar covariates.cov --covar-variance-standardize --out assoc_results/additive_firth_comb_spanish
		# meta analysis:
		tools/METAL-2020-05-05/build/bin/metal scripts/metal_script.txt
		
		# recessive
		plink2 --bfile {params.plink_dom_rez_comb} --glm firth hide-covar recessive --covar covariates.cov --covar-variance-standardize --out assoc_results/recessive_firth_comb
		
		# functional
		plink2 --bfile {params.functional} --glm firth hide-covar --covar covariates.cov --covar-variance-standardize --out assoc_results/functional_firth
		
		# males
		# functional / males
		plink2 --bfile {params.functional} --glm firth hide-covar --covar covariates.cov --covar-variance-standardize --filter-males --out assoc_results/functional_males_firth 
		
		# additive / males		
		plink2 --bfile {params.plink_dom_rez_comb} --glm firth hide-covar --covar covariates.cov --chr X --covar-variance-standardize --filter-males --out assoc_results/additive_males_comb_firth

		
		
		######
		###### other
		plink --bfile {params.plink_dom_rez_comb} --freqx --within cluster_file_cmh_external.dat --out assoc_results/__for_external_cmh_comb 
		
		mkdir -p assoc_results
		mkdir -p check
		plink --bfile filtered_maf001 --recode vcf --out check/filtered_maf001_pre_collapse # for control purposes
		plink --bfile {params.plink_dom_rez_comb} --recode vcf --out check/filtered_maf001_post_collapse # for control purposes
		cat check/filtered_maf001_post_collapse.vcf |grep -v "^##" |  cut -f3- > check/post_collapse.tsv
		cat check/filtered_maf001_pre_collapse.vcf | grep -v "^##" | cut -f3- > check/pre_collapse.tsv
		"""


rule plink_testing_single_SNVs:
	input:
		"filtered_maf0001.bed",
		pathAncestry="sequence.file.normID.noChrM.bcftools_haileur.vcf.gz"
	output:
		pathOut="assoc_results/filtered_all.assoc.fisher"
	resources: cpus=1, mem_mb=18000, time_job=720
	params:
		filtered_all="filtered_all",
		covariates="covariates.cov"
	shell:
		"""
		# single variant tests - logistic?

		plink --bfile {params.filtered_all} --model fisher  --out assoc_results/filtered_all
		plink --bfile {params.filtered_all} --assoc fisher --out assoc_results/filtered_all
		plink --bfile {params.filtered_all} --freqx --out assoc_results/filtered_all
		plink --bfile {params.filtered_all} --freq case-control --out assoc_results/filtered_all
		plink --bfile {params.filtered_all} --freq case-control --filter-males --out assoc_results/filtered_males
				
		plink2 --bfile {params.filtered_all} --glm firth hide-covar  --covar-variance-standardize --out assoc_results/filtered_all_glm_add --covar {params.covariates}
		plink2 --bfile {params.filtered_all} --glm firth hide-covar recessive  --covar-variance-standardize --out assoc_results/filtered_all_glm_rec --covar {params.covariates} 
		
		"""
		
		
		
#https://www.cambridge.org/core/journals/twin-research-and-human-genetics/article/statistical-properties-of-singlemarker-tests-for-rare-variants/AB37EB4C1C8945875E8A41FBD889690E
#https://zzz.bwh.harvard.edu/plink/perm.shtml
#https://www.cog-genomics.org/plink/2.0/assoc#glm_footnote1
#http://www.nealelab.is/blog/2017/9/11/details-and-considerations-of-the-uk-biobank-gwas




