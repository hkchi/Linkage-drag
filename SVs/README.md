# SVs

Calculate SNP density

```
win=500000
vcftools --vcf all.snps.vcf --SNPdensity $win --out all.snps.$win
vcftools --vcf all.indels.len50.vcf --SNPdensity $win --out all.indels.$win
```

Windows with > 50% coverage by one-to-one MUMmer alignments (for filtering)

```
win=500000
bedtools makewindows -g HA412.chr_len.txt -w $win > HA412.$win.bed
for cmp in HA412_XRQ HA412_PSC8 HA412_RHA438 HA412_IR HA412_HA89
do
awk '{print $10"\t"$1"\t"$2}' $cmp.b1000.c1000.i90.l1000.1coords > tmp.$cmp.1coords.bed
done
cat tmp.*.1coords.bed | sort -k1,1 -k2,2n | bedtools merge -i - > tmp.all.1coords.bed
bedtools intersect -a HA412.$win.bed -b tmp.all.1coords.bed -wo | \
sort -k1,1 -k2,2n -k5,5n | \
awk 'BEGIN{c="";s="";e="";a=0}{if($1==c&&$2==s&&$3==e){a=a+$7}else{if(a>0&&a/(e-s)>0.5){print c"\t"s"\t"e};c=$1;s=$2;e=$3;a=$7}}END{if(a>0&&a/(e-s)>0.5){print c"\t"s"\t"e}}' \
> HA412.$win.c50.bed
rm tmp.*.bed
```
