## Generate PAV matrix and introgression matrix
## Kaichi Huang 2022 May

library(tidyverse)

my.genomes <- c("HA412", "XRQ", "PSC8", "RHA438", "IR", "HA89", "LR1", "OQP8", "HA300", "PI659440")

# PAV matrix and introgression matrix
my.PAV <- data.frame(Rep_gene=character())
my.introgression <- data.frame(Rep_gene=character())
for (genome in my.genomes) {
  table_file <- paste0(genome,"_synteny_filtered.txt")
  # Read the table
  my.table <- read_delim(table_file, delim=" ", col_types="ccddccddcddcdd")
  my.table <- my.table %>% group_by(Rep_gene) %>% slice_head(n=1) %>% ungroup()
  my.table <- my.table %>% filter(!grepl("Chr00",get(my.name))) # genes on small contigs
  # PAV
  tmp.PAV <- my.table %>% select(Rep_gene, all_of(my.name)) %>% mutate(pav=case_when((is.na(get(my.name))|grepl("MT|CP",get(my.name)))~0, T~1)) %>% select(Rep_gene, pav)
  colnames(tmp.PAV) <- c("Rep_gene", my.name)
  my.PAV <- full_join(my.PAV, tmp.PAV)
  # Introgression
  if (genome != "PI659440") {
    # Introgression in only the cultivars
    my.i1 <- read_tsv(paste0(genome,".i1.tsv"), col_names=F) %>% mutate(ID=paste0(X1,":",X2,"-",X3)) %>% select(ID)
    my.i2 <- read_tsv(paste0(genome,".i2.tsv"), col_names=F) %>% mutate(ID=paste0(X1,":",X2,"-",X3)) %>% select(ID)
    tmp.introgression <- my.table %>%
      mutate(chr=case_when(!is.na(get(my.name))~chr2, T~chr_norep)) %>%
      mutate(start=case_when(!is.na(get(my.name))~start2, T~start_norep)) %>%
      mutate(end=case_when(!is.na(get(my.name))~end2, T~end_norep)) %>%
      select(Rep_gene,chr,start,end) %>%
      mutate(chr=paste0(genome,chr)) %>%
      mutate(ID=paste0(chr,":",start,"-",end))
    tmp.introgression$i <- 0
    tmp.introgression$i[which(tmp.introgression$ID %in% my.i1$ID)] <- 1
    tmp.introgression$i[which(tmp.introgression$ID %in% my.i2$ID)] <- 2
    tmp.introgression <- tmp.introgression %>% select(Rep_gene, i)
    colnames(tmp.introgression) <- c("Rep_gene", my.name)
    my.introgression <- full_join(my.introgression, tmp.introgression)
  }
}
write_tsv(my.PAV, "my.PAV.txt")
write_tsv(my.introgression, "my.introgression.txt")
