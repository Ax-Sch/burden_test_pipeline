library(tidyverse)

combined_translated<-read_tsv(file="../all_sample_infos.tsv")%>%
  mutate(cohort=as.integer(cohort=="Spain"))%>%
  mutate(sex_age2=sex*age2)


for (fam_path in c("cp_het/recessive_collapsed_comp_gnomAD001.fam",
		"cp_het/dom_rez_gnomAD001.fam",
		"cp_het/recessive_collapsed_comp_gnomAD0001.fam",
		"cp_het/dom_rez_gnomAD0001.fam",
		"cp_het/functional_collapsed.fam"
		)) {
fam_file<-read_tsv(file=fam_path, col_names=c("FI","II","PA","MA","S","P"))
fam_file<-fam_file %>%
	left_join(combined_translated, by = c("II"="X1")) %>%
	mutate(FID=0)

write_tsv(fam_file[,c("FID","II","FID","FID","sex","phenotype")], file = fam_path, col_names = FALSE)
}


for (fam_path in c("filtered_maf001.fam",
                   "filtered_maf0001.fam",
                   "filtered_all.fam"
)) {
  fam_file<-read_delim(file=fam_path, col_names=c("FI","II","PA","MA","S","P"), delim = " ")
  fam_file<-fam_file %>%
    left_join(combined_translated, by = c("II"="X1")) %>%
    mutate(FID=0)
  write_tsv(fam_file[,c("FID","II","FID","FID","sex","phenotype")], file = fam_path, col_names = FALSE)
}






write_tsv(fam_file[,c("FID","II","age","age2","sex_age","PC1","PC2",
                                 "PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")], file = "covariates.cov", col_names = FALSE)

# generate cluster file for CMH
write_tsv(fam_file[,c("FID","II","cohort")], file="cluster_file.dat")

# generate file of italian participants
write_tsv(fam_file %>% filter(cohort==0) %>% dplyr::select(FID, II), file="italian_participants.dat")

# generate cluster file for external CMH:
fam_file<-fam_file %>% mutate(cmh_external=paste(cohort, sex, phenotype, sep="_"))

write_tsv(fam_file[,c("FID","II","cmh_external")], file="cluster_file_cmh_external.dat")
