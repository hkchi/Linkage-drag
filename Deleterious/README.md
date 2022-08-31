# Deleterious

Annotate VCF with snpEff

```
mkdir ./snpEff; cd ./snpEff
snpEff="/scratch/hkchi/program/snpEff/snpEff.jar"
configFile="/scratch/hkchi/program/snpEff/snpEff.config" # with customized genome annotation
java -Xmx20G -jar $snpEff -c $configFile -no-utr -no-downstream -no-upstream -no-intergenic HA412 all.snps.vcf > all.snps.anno.vcf
```

Extract synonymous mutations, nonsynonymous mutations and alternative stop codons

```
grep -E '^#|synonymous_variant' all.snps.anno.vcf > all.snps.syn.vcf
grep -E '^#|missense_variant' all.snps.anno.vcf > all.snps.nonsyn.vcf
grep -E '^#|stop_gained|stop_lost' all.snps.anno.vcf > all.snps.stop_codon.vcf
```

Count numbers of sites per window

```
win=500000
vcftools --vcf all.snps.syn.vcf --SNPdensity $win --out all.snps.syn.$win
vcftools --vcf all.snps.nonsyn.vcf --SNPdensity $win --out all.snps.nonsyn.$win
vcftools --vcf all.snps.stop_codon.vcf --SNPdensity $win --out all.snps.stop_codon.$win
```
