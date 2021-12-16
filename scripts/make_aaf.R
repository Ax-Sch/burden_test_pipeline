library(data.table)
vars<-fread("./split_annot_06_07/sequence.file.normID.noChrM.bcftools_haileur_finalAnnot.tsv", fill=TRUE, header = TRUE, na.strings = ".")
vars<-vars[, c("chrom","pos","ref","alt"):=tstrsplit(ID, ":", fixed=TRUE)]

###### MASKS ######
vars<-vars[, ':=' (pph1=like(Polyphen2_HDIV_pred,"D"), 
                  pph2=like(Polyphen2_HVAR_pred,"D"), 
                  lrt=like(LRT_pred,"D"), 
                  mt=like(MutationTaster_pred,"D") | like(MutationTaster_pred,"A"), 
                  sif=like(SIFT_pred,"D"),
                  missense=like(Consequence, "missense_variant"),
                  syn=like(Consequence, "synonymous_variant"),
                  pLOF=(IMPACT=="HIGH"))]


vars<-vars[, ':=' (five_path=(pph1 & pph2 & lrt & mt & sif & missense), 
             one_path=((pph1 | pph2 | lrt | mt | sif)& missense),
           moderate_non_missense=(!missense & (IMPACT=="MODERATE")))]

vars<-vars[, by=c("ID", "Gene"), mask_anno:=ifelse(sum(pLOF)>0,"pLoF", 
             ifelse(sum(five_path)>0, "missense.5in5",
             ifelse(sum(one_path)>0, "missense.1in5",
             ifelse(sum(moderate_non_missense)>0, "moderate.non.missense",
             ifelse(sum(syn)>0, "synonymous", "not_relevant")
             )  ) ) )
            ]

vars<-vars[, M0:=syn]
vars<-vars[, M1:=pLOF]
vars<-vars[, M3:=(pLOF| moderate_non_missense | five_path)]
vars<-vars[, M4:=(M3 | one_path)]

#synonymous_variant
###### ADD CASANOVA ######
vars_casanova<-fread("./casanova_results/patho_variants.txt", fill=TRUE, header = TRUE, na.strings = ".")
anno_list<-rbindlist(list(vars_casanova, vars[, c("ID", "Gene", "mask_anno")]), use.names = FALSE)

vars<-vars[, functional:=(ID %in% vars_casanova$SNP)]

###### AF ######
vars=vars[, maxAF:=pmax(AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF, AA_AF, EA_AF,
              gnomAD_AFR_AF, gnomAD_AMR_AF, gnomAD_ASJ_AF, gnomAD_EAS_AF, gnomAD_FIN_AF, gnomAD_NFE_AF, gnomAD_OTH_AF, gnomAD_SAS_AF, 
              na.rm = TRUE)]

vars=vars[, minAF:=pmin(AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF, AA_AF, EA_AF,
                        gnomAD_AFR_AF, gnomAD_AMR_AF, gnomAD_ASJ_AF, gnomAD_EAS_AF, gnomAD_FIN_AF, gnomAD_NFE_AF, gnomAD_OTH_AF, gnomAD_SAS_AF, 
                        na.rm = TRUE)]

vars=vars[,AF_category:=ifelse( (maxAF<0.001 | minAF>0.999 | is.na(maxAF)), 0.0005,
                        ifelse( (maxAF<0.01 | minAF>0.99), 0.005, 10)) ]


###### SET LIST #######
collapsed_vars<-vars[,concats:=paste(ID, collapse=","), by=Gene]
collapsed_vars<-unique(collapsed_vars[,c("Gene","chrom","concats")])
collapsed_vars<-collapsed_vars[,pos:=1:nrow(collapsed_vars)]


fwrite(unique(anno_list), file = "regenie.anno.file.txt", col.names = FALSE, sep=" ")
fwrite(unique(vars[, c("ID", "AF_category")]), file = "regenie.aaf.file.txt", col.names = FALSE, sep=" ")
fwrite((collapsed_vars[,c("Gene","chrom","pos","concats")]), file = "regenie.set.list.txt", col.names = FALSE, sep=" ")

fwrite(vars[, -"concats"], file = "variables_categorized.tsv", col.names = TRUE, sep="\t")

     
