library(tidyverse)

gene_masks_w_evidence<-c("ENSG00000115267.M1.0.01","ENSG00000115267.M3.0.01","ENSG00000115267.M4.0.01", "ENSG00000027697.M4.0.01",
                         "ENSG00000027697.M3.0.01","ENSG00000130234.M3.0.001","ENSG00000196664.M4.0.01","TLR7_functional.MF.109")

genes_w_evidence<-str_split(gene_masks_w_evidence,fixed("."), simplify=TRUE)[,1]

single_marker_files<-c("filtered_all.assoc.fisher") # add additive_single_SNV_0001.assoc.logistic
single_marker_AF_files<-c("freq0001_single_SNV.frqx", "freq001_single_SNV.frqx")
single_marker_count_files<-c("freq0001_single_SNV.frq.cc", "freq001_single_SNV.frq.cc")
vep_annotation<-read_delim("split_annot_06_07/sequence.file.normID.noChrM.bcftools_haileur_finalAnnot.tsv", delim=";")
set_list_file<-"regenie.set.list.txt"

AAF_file<-read_delim("regenie.aaf.file.txt", delim=" ", col_names=c("SNP","AF_category"))
Annotation_file<-read_delim("regenie.anno.file.txt", delim=" ", col_names=c("SNP","ENSG","CONS"))

#set_list<-read_delim(set_list_file, delim = " ", col_names=c("ENSG","Chrom","Pseudopos", "variants"))
#set_list$variants<-str_split(set_list$variants, fixed(","))

relevant_sets<-Annotation_file %>% filter(ENSG %in% genes_w_evidence)

relevant_vars<-unlist(relevant_sets$SNP)

vars_single_marker<-data.frame()

for (listx in single_marker_files) {
  interim<-read_delim(paste0("assoc_results/",listx), delim=" ", trim_ws=TRUE) %>% 
    filter(SNP %in% relevant_vars) %>%
    mutate(file=listx)
  vars_single_marker<-rbind(vars_single_marker, interim)
}

vars_single_marker_AF<-data.frame()
for (listx in single_marker_AF_files) {
  interim<-read_tsv(paste0("assoc_results/",listx), trim_ws = TRUE) %>% 
    filter(SNP %in% relevant_vars) %>%
    select(-CHR)
  vars_single_marker_AF<-rbind(vars_single_marker_AF, interim)
}
vars_single_marker_AF<-unique(vars_single_marker_AF)

vars_single_marker_count<-data.frame()
for (listx in single_marker_count_files) {
  interim<-read_delim(paste0("assoc_results/",listx), delim=" ",trim_ws = TRUE) %>% 
    filter(SNP %in% relevant_vars) %>%
    select(-CHR, -A1, -A2)
  vars_single_marker_count<-rbind(vars_single_marker_count, interim)
}
vars_single_marker_count<-unique(vars_single_marker_count)



compiled_AF_stats<-vars_single_marker %>% 
  #left_join(vars_single_marker_AF, by=c("SNP"="SNP")) %>%
  left_join(vars_single_marker_count, by=c("SNP"="SNP"))%>%
  left_join(AAF_file, by=c("SNP"="SNP") ) %>%
  left_join(Annotation_file, by=c("SNP"="SNP"))%>%
  left_join(vep_annotation, by=c("SNP"="ID") ) %>%
  select(-F_A, -F_U)%>%
  mutate(AC_A=round(MAF_A*NCHROBS_A,2), AC_U=round(MAF_U*NCHROBS_U,2)) %>%
  distinct()

first_cols<-c("SNP","AC_A", "AC_U","P", "OR", "SYMBOL", "CONS", "Protein_position",	"Amino_acids")
new_order<-c(first_cols, colnames(compiled_AF_stats %>% select(-first_cols)) )


write_tsv(compiled_AF_stats[,new_order], file = "compiled_AF_stats.tsv")
  




