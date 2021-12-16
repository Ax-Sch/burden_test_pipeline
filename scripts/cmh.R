library(stats)
library(tidyverse)
for (filename in c("./assoc_results/for_external_cmh001.frq.strat", "./assoc_results/for_external_cmh0001.frq.strat")){
for_cmh001<-read_delim(filename, delim=" ", trim_ws=TRUE, guess_max = 5000) %>%
  mutate(sex=str_split(CLST, fixed("_"), simplify=TRUE)[,2])%>%
  mutate(phenotype=str_split(CLST, fixed("_"), simplify=TRUE)[,3])%>%
  mutate(cohort=str_split(CLST, fixed("_"), simplify=TRUE)[,1])%>%
  mutate(REF=NCHROBS-MAC)%>%
  rename(ALT=MAC)
# double

for_cmh001<-for_cmh001 %>% 
  gather(key="allele", value="counts", c("ALT","REF")) %>%
  select(cohort, SNP, phenotype, allele, counts)

# needs: center, gene, affection, variant, count per row
get_p_vals<-function(x, data){
  Table = xtabs(counts ~ allele + phenotype + cohort,
                data=data %>% filter(SNP==x))
  CMH<-mantelhaen.test(Table, exact=TRUE)
  return(data.frame(
    p_val=CMH$p.value,
    est=CMH$estimate
  ))
}

snps<-unique(for_cmh001$SNP)
p_vals<-lapply(snps, function(x){get_p_vals(x, for_cmh001)}) 
p_vals<-do.call(rbind, p_vals)
p_vals_cmh<-cbind(snps, p_vals) 
p_vals_cmh %>% arrange(-p_val) %>% filter(p_val<0.05)
write_tsv(p_vals_cmh, file=paste0(filename, "exact_cmh_p_values.tsv"))
}
