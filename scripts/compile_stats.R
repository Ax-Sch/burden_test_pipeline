library(tidyverse)
library(ggrepel)
setwd("..")

stats_add<-read_delim("assoc_results/additive_fisher_xtoauto_comb.model", delim = " ", trim_ws=TRUE) %>% 
  filter(TEST %in% c("REC", "ALLELIC"))

stats_1_cp<-read_delim("assoc_results/cp_het_fisher_001.model", delim = " ", trim_ws=TRUE)%>% 
  filter(TEST %in% c("REC"),  CHR!=23) %>% 
  mutate(TEST="REC_CP_HET")
stats_01_cp<-read_delim("assoc_results/cp_het_fisher_0001.model", delim = " ", trim_ws=TRUE)%>% 
  filter(TEST %in% c("REC"), CHR!=23)%>% 
  mutate(TEST="REC_CP_HET")


stats_combined_autosomes<-rbind(stats_1_add, stats_01_add, stats_1_cp, stats_01_cp) %>% arrange(CHR,SNP)

stats_1_sexmales<-read_delim("assoc_results/sexchr_males_fisher_001.assoc.fisher", delim = " ", trim_ws=TRUE) %>%
  filter(CHR==23) %>%
  mutate(TEST="additive_males")

stats_01_sexmales<-read_delim("assoc_results/sexchr_males_fisher_0001.assoc.fisher", delim = " ", trim_ws=TRUE) %>%
  filter(CHR==23)%>%
  mutate(TEST="additive_males")

stats_1_sex<-read_delim("assoc_results/sexchr_fisher_001.assoc.fisher", delim = " ", trim_ws=TRUE) %>%
  filter(CHR==23) %>%
  mutate(TEST="ALLELIC")

stats_01_sex<-read_delim("assoc_results/sexchr_fisher_0001.assoc.fisher", delim = " ", trim_ws=TRUE) %>%
  filter(CHR==23)%>%
  mutate(TEST="ALLELIC")

stats_TLR7funct_males<-read_delim("assoc_results/functional_males.assoc.fisher", delim = " ", trim_ws=TRUE) %>%
  filter(CHR==23) %>%
  mutate(TEST="additive_males")

stats_01_TLR7funct<-read_delim("assoc_results/functional.assoc.fisher", delim = " ", trim_ws=TRUE) %>%
  filter(CHR==23)%>%
  mutate(TEST="ALLELIC")

stats_combined_sex_chroms<-rbind(stats_1_sexmales, stats_01_sexmales, stats_1_sex, 
                                 stats_01_sex, stats_TLR7funct_males, stats_01_TLR7funct)


glm_stats_01_TLR7funct<-read_tsv("assoc_results/functional_glm.PHENO1.glm.firth", trim_ws=TRUE) %>% mutate(TEST="ALLELIC")
#glm_stats_01_TLR7funct_males<-read_tsv("assoc_results/functional_glm", trim_ws=TRUE)

glm_stats_1_add<-read_tsv("assoc_results/additive_glm_001.PHENO1.glm.firth", trim_ws=TRUE) %>% mutate(TEST="ALLELIC")
glm_stats_01_add<-read_tsv("assoc_results/additive_glm_0001.PHENO1.glm.firth", trim_ws=TRUE)  %>% mutate(TEST="ALLELIC")

glm_stats_1_rec<-read_tsv("assoc_results/recessive_glm_001.PHENO1.glm.firth", trim_ws=TRUE)  %>% mutate(TEST="REC")
glm_stats_01_rec<-read_tsv("assoc_results/recessive_glm_0001.PHENO1.glm.firth", trim_ws=TRUE) %>% mutate(TEST="REC")

glm_stats_1_cprec<-read_tsv("assoc_results/cp_het_glm_001.PHENO1.glm.firth", trim_ws=TRUE) %>% mutate(TEST="REC_CP_HET")
glm_stats_01_cprec<-read_tsv("assoc_results/cp_het_glm_0001.PHENO1.glm.firth", trim_ws=TRUE) %>% mutate(TEST="REC_CP_HET")

#glm_stats_1_males<-read_tsv("assoc_results/sexchr_males_glm_001.PHENO1.glm.firth", trim_ws=TRUE)
#glm_stats_01_males<-read_tsv("assoc_results/sexchr_males_glm_0001.PHENO1.glm.firth", trim_ws=TRUE)


stats_combined_glm<-rbind(glm_stats_01_TLR7funct, glm_stats_1_add, glm_stats_01_add, 
                          glm_stats_1_rec, glm_stats_01_rec, glm_stats_1_cprec, glm_stats_01_cprec)
colnames(stats_combined_glm)<-paste0(colnames(stats_combined_glm), "_glm_firth")

stats_fisher_and_glm<-stats_combined_autosomes %>% left_join(stats_combined_glm, by=c("SNP"="ID_glm_firth", "TEST"="TEST_glm_firth")) %>% filter(P!=1)
write_tsv(stats_fisher_and_glm, file="stats_fisher_and_glm.tsv")


ggplot(stats_fisher_and_glm, aes(x=-log10(P), y=-log10(P_glm_firth))) +
  geom_abline(intercept=c(0,0), slope = 1, color="grey")+
  geom_point(color="lightblue")+
  geom_text_repel(data=stats_fisher_and_glm %>% filter(-log10(P)>2 | -log10(P_glm_firth)>2), 
            aes(x=-log10(P), y=-log10(P_glm_firth), label=paste(SNP, TEST, sep="\n")), size=2.5, min.segment.length=0.1)+
  theme_bw()+
  coord_cartesian(xlim=c(0,5), ylim=c(0,5))


stats_combined_glm_X<-rbind(glm_stats_1_add, glm_stats_01_add, glm_stats_01_TLR7funct) %>% filter(`#CHROM`=="X")
colnames(stats_combined_glm_X)<-paste0(colnames(stats_combined_glm_X), "_glm_firth")

stats_fisher_and_glm_X<-stats_combined_sex_chroms %>% left_join(stats_combined_glm_X, by=c("SNP"="ID_glm_firth", "TEST"="TEST_glm_firth")) %>% filter(P!=1)
write_tsv(stats_fisher_and_glm_X, file="stats_fisher_and_glm_X.tsv")


ggplot(stats_fisher_and_glm_X, aes(x=-log10(P), y=-log10(P_glm_firth))) +
  geom_abline(intercept=c(0,0), slope = 1, color="grey")+
  geom_point(color="lightblue")+
  geom_text_repel(data=stats_fisher_and_glm_X %>% filter(-log10(P)>2 | -log10(P_glm_firth)>2), 
            aes(x=-log10(P), y=-log10(P_glm_firth), label=paste(SNP, TEST, sep="\n")), size=3,
            min.segment.length=0.1)+
  theme_bw()+
  coord_cartesian(xlim=c(0,5), ylim=c(0,5))


single_var_stats_add<-read_delim("assoc_results/filtered_all.assoc.fisher", delim = " ", trim_ws=TRUE,guess_max = 50000) %>% 
  rename(P_fisher_add=P, OR_fisher_add=OR)
single_var_stats<-read_delim("assoc_results/filtered_all.model", delim = " ", trim_ws=TRUE,guess_max = 50000)
single_var_stats_rec_fisher<-single_var_stats %>% 
  filter(TEST=="REC")%>% 
  rename(P_fisher_rec=P)
single_var_stats_add_glm<-read_tsv("assoc_results/filtered_all_glm_add.PHENO1.glm.firth",  trim_ws=TRUE, guess_max = 50000) %>% 
  rename(P_firth_add=P, OR_firth_add=OR)
single_var_stats_rec_glm<-read_tsv("assoc_results/filtered_all_glm_rec.PHENO1.glm.firth",  trim_ws=TRUE ,guess_max = 50000)%>% 
  rename(P_firth_rec=P, OR_firth_rec=OR)
single_var_freqs<-read_delim("assoc_results/filtered_all.frq.cc", delim = " ", trim_ws=TRUE,guess_max = 50000) %>% 
  mutate(AC_cases=round(MAF_A * NCHROBS_A), AC_controls=round(MAF_U * NCHROBS_U))
single_var_freqs_males<-read_delim("assoc_results/filtered_males.frq.cc", delim = " ", trim_ws=TRUE,guess_max = 50000) %>% 
  filter(CHR==23) %>%
  mutate(HEMI_cases=round(MAF_A * NCHROBS_A), HEMI_controls=round(MAF_U * NCHROBS_U))

var_freqs<-single_var_freqs %>% 
  left_join(single_var_freqs_males %>% select(SNP,HEMI_cases, HEMI_controls), by=c("SNP"="SNP")) %>%
  left_join(single_var_stats %>% filter(TEST=="GENO") %>% select(SNP, AFF, UNAFF), by=c("SNP"="SNP"))


single_var_stats_merged<-var_freqs %>% 
  left_join(single_var_stats_add %>% select(SNP, P_fisher_add, OR_fisher_add) , by=c("SNP"="SNP"))%>%
  left_join(single_var_stats_rec_fisher %>% select(SNP, P_fisher_rec) , by=c("SNP"="SNP"))%>%
  left_join(single_var_stats_add_glm %>% select(ID, P_firth_add, OR_firth_add) , by=c("SNP"="ID"))%>%
  left_join(single_var_stats_rec_glm %>% select(ID, P_firth_rec, OR_firth_rec) , by=c("SNP"="ID"))

vars<-read_tsv("variables_categorized.tsv", guess_max = 50000)

vars<-vars %>% left_join(single_var_stats_merged,  by=c("ID"="SNP"))

write_tsv(vars, file="variants_with_p_values_etc.tsv")





#### check against cmh: ######
cmh_ps<-rbind(read_tsv("./assoc_results/for_external_cmh001.frq.stratexact_cmh_p_values.tsv"),
      read_tsv("./assoc_results/for_external_cmh0001.frq.stratexact_cmh_p_values.tsv")) %>%
  rename(p_val_CMH=p_val)
cmh_ps <- cmh_ps %>% left_join(stats_fisher_and_glm %>% filter(TEST=="ALLELIC"), by=c("snps"="SNP"))
ggplot(cmh_ps, aes(x=-log10(P), y=-log10(p_val_CMH))) +
  geom_abline(intercept=c(0,0), slope = 1, color="grey")+
  geom_point(color="lightblue")+
  geom_text_repel(data=cmh_ps %>% filter(-log10(P)>2 | -log10(p_val_CMH)>2), 
                  aes(x=-log10(P), y=-log10(p_val_CMH), label=paste(snps, TEST, sep="\n")), size=3,
                  min.segment.length=0.1)+
  theme_bw()+
  coord_cartesian(xlim=c(0,5), ylim=c(0,5))


cmh_ps<-rbind(read_tsv("./assoc_results/for_external_cmh001.frq.stratexact_cmh_p_values.tsv"),
              read_tsv("./assoc_results/for_external_cmh0001.frq.stratexact_cmh_p_values.tsv")) %>%
  rename(p_val_CMH=p_val)

cmh_ps <- cmh_ps %>% left_join(stats_fisher_and_glm_X %>% filter(TEST=="ALLELIC"), by=c("snps"="SNP"))
ggplot(cmh_ps, aes(x=-log10(P), y=-log10(p_val_CMH))) +
  geom_abline(intercept=c(0,0), slope = 1, color="grey")+
  geom_point(color="lightblue")+
  geom_text_repel(data=cmh_ps %>% filter(-log10(P)>2 | -log10(p_val_CMH)>2), 
                  aes(x=-log10(P), y=-log10(p_val_CMH), label=paste(snps, TEST, sep="\n")), size=3,
                  min.segment.length=0.1)+
  theme_bw()+
  coord_cartesian(xlim=c(0,5), ylim=c(0,5))

