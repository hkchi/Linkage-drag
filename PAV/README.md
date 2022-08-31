# PAV

Introgressed regions for each genome

```
for genome in HA412 XRQ PSC8 RHA438 IR HA89 LR1 OQP8 HA300
do
awk 'NR>1 && $4==1{print $1"\t"int($2)"\t"int($3)}' my.regions.$genome.txt > tmp.$genome.1.bed
awk 'NR>1 && $4>1{print $1"\t"int($2)"\t"int($3)}' my.regions.$genome.txt > tmp.$genome.2.bed
done
```

Define introgression value for each gene in each genome

```
for genome in HA412 XRQ PSC8 RHA438 IR HA89 LR1 OQP8 HA300
do
gene_table="${genome}_synteny_filtered.txt"
awk 'NR>1' $gene_table | \
	awk '{if($5=="NA"){print "'$genome'"$9"\t"$10"\t"$11}else{print "'$genome'"$6"\t"$7"\t"$8}}' | \
	awk '{if($3>=$2){print}else{print $1"\t"$3"\t"$2}}' \
	> tmp.${genome}_genes.bed
cat tmp.${genome}_genes.bed | \
	sort -k1,1 -k 2,2n | \
	bedtools intersect -a - -b tmp.$genome.1.bed -wo | \
	sort -k1,1 -k2,2n -k5,5n | \
	awk 'BEGIN{c="";s="";e="";a=0}{if($1==c&&$2==s&&$3==e){a=a+$7}else{if(a>0&&a/(e-s)>0.5){print c"\t"s"\t"e};c=$1;s=$2;e=$3;a=$7}}END{if(a>0&&a/(e-s)>0.5){print c"\t"s"\t"e}}' \
	> $genome.i1.tsv
cat tmp.${genome}_genes.bed | \
	sort -k1,1 -k 2,2n | \
	bedtools intersect -a - -b tmp.$genome.2.bed -wo | \
	sort -k1,1 -k2,2n -k5,5n | \
	awk 'BEGIN{c="";s="";e="";a=0}{if($1==c&&$2==s&&$3==e){a=a+$7}else{if(a>0&&a/(e-s)>0.5){print c"\t"s"\t"e};c=$1;s=$2;e=$3;a=$7}}END{if(a>0&&a/(e-s)>0.5){print c"\t"s"\t"e}}' \
	> $genome.i2.tsv
done
rm tmp.*.bed
```

Generate PAV matrix and introgression matrix

```
Rscript Analyses.PAV.R
```
