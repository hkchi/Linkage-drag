## Perform stepwise merging of SVs from each genome comparison using the methods in Audano et al. (2019)
## Kaichi Huang 2022 Mar

library(tidyverse)

SV_type <- "inversions" # "translocations"

# Function to calculate region overlap
my.overlap <- function(start.1, end.1, start.2, end.2) {
  vec1 <- seq(floor(start.1), floor(end.1))
  vec2 <- seq(floor(start.2), floor(end.2))
  overlap <- length(intersect(vec1, vec2))
  return (overlap)
}

cmp.list <- c("HA412_XRQ", "HA412_PSC8", "HA412_RHA438", "HA412_IR", "HA412_HA89")

# Read in the first file as reference dataset
i <- 1
cmp <- cmp.list[i]
my.data <- read_tsv(paste0(cmp,".b1000.c1000.i90.l1000.m.",SV_type,".tsv"), col_names=F) %>%
  rename(chr1=X1, start1=X2, end1=X3, chr2=X4, start2=X5, end2=X6)
# Integrate alignment info into ANN field
my.data <- my.data %>% mutate(ANN=paste0(chr2,":",start2,"-",end2)) %>% select(-chr2, -start2, -end2)

# Loop and integrate each dataset stepwise
for (i in 2:length(cmp.list)) {
  cmp <- cmp.list[i]
  my.sv <- read_tsv(paste0(cmp,".b1000.c1000.i90.l1000.m.",SV_type,".tsv"), col_names=F) %>%
    rename(chr1=X1, start1=X2, end1=X3, chr2=X4, start2=X5, end2=X6)
  # Copyright: Kaichi Huang
  for (j in 1:nrow(my.sv)) {
    sv.chr1=my.sv$chr1[j]
    sv.start1=my.sv$start1[j]
    sv.end1=my.sv$end1[j]
    sv.chr2=my.sv$chr2[j]
    sv.start2=my.sv$start2[j]
    sv.end2=my.sv$end2[j]
    cmp.data <- my.data %>% filter(chr1==sv.chr1 & !(start1>sv.end1 | end1<sv.start1))
    if (nrow(cmp.data)==0) {
      # a new SV
      tmp.data <- data.frame(chr1=sv.chr1, start1=sv.start1, end1=sv.end1, ANN=paste0(sv.chr2,":",sv.start2,"-",sv.end2))
      my.data <- rbind(my.data, tmp.data)
    } else {
      cmp.data.ol <- cmp.data %>%
        rowwise() %>% mutate(sv.overlap=my.overlap(start1, end1, sv.start1, sv.end1)) %>% ungroup() %>%
        filter(sv.overlap/(end1-start1) >0.5 & sv.overlap/(sv.end1-sv.start1) > 0.5)
      if (nrow(cmp.data.ol)==0) {
        # a new SV
        tmp.data <- data.frame(chr1=sv.chr1, start1=sv.start1, end1=sv.end1, ANN=paste0(sv.chr2,":",sv.start2,"-",sv.end2))
        my.data <- rbind(my.data, tmp.data) # (*O.O*)
      } else {
        # overlap with an existing discovery variant - merge and keep the original reference coordinates
        if (nrow(cmp.data.ol)>1) {print(cmp.data.ol)}
        k <- which(my.data$chr1==cmp.data.ol$chr1[1] & my.data$start1==cmp.data.ol$start1[1] & my.data$end1==cmp.data.ol$end1[1])
        my.data$ANN[k] <- paste(my.data$ANN[k], paste0(sv.chr2,":",sv.start2,"-",sv.end2), sep=",")
      }
    }
  }
}

my.data <- my.data %>% arrange(chr1, start1, end1)
write_tsv(my.data, paste0("all.svs.",SV_type,".tsv"), col_names=F)
