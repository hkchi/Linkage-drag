# Variants

Run MUMmer to align genomes

```
cd /scratch/hkchi/linkage_drag/MUMmer/
bash run_mummer.sh /scratch/hkchi/ref/Ha412HO.onlychr.fasta /scratch/hkchi/ref/HanXRQv2.onlychr.fasta HA412_XRQ
bash run_mummer.sh /scratch/hkchi/ref/Ha412HO.onlychr.fasta /scratch/hkchi/ref/HanPSC8.onlychr.fasta HA412_PSC8
bash run_mummer.sh /scratch/hkchi/ref/Ha412HO.onlychr.fasta /scratch/hkchi/ICSG_genomes/genomes/HanRHA438.onlychr.fasta HA412_RHA438
bash run_mummer.sh /scratch/hkchi/ref/Ha412HO.onlychr.fasta /scratch/hkchi/ICSG_genomes/genomes/HanIR.onlychr.fasta HA412_IR
bash run_mummer.sh /scratch/hkchi/ref/Ha412HO.onlychr.fasta /scratch/hkchi/ICSG_genomes/genomes/HanHA89.onlychr.fasta HA412_HA89
bash run_mummer.sh /scratch/hkchi/ref/Ha412HO.onlychr.fasta /scratch/hkchi/ICSG_genomes/genomes/HanLR1.onlychr.fasta HA412_LR1
bash run_mummer.sh /scratch/hkchi/ref/Ha412HO.onlychr.fasta /scratch/hkchi/ICSG_genomes/genomes/HanOQP8.onlychr.fasta HA412_OQP8
bash run_mummer.sh /scratch/hkchi/ref/Ha412HO.onlychr.fasta /scratch/hkchi/ICSG_genomes/genomes/HanHA300.onlychr.fasta HA412_HA300
```

Extract homologous SNPs and small InDels (< 50bp) from MUMmer outputs

```
mkdir ./snps_indels; cd ./snps_indels/
bin="/home/hkchi/program/mummer-4.0.0beta2"
for cmp in HA412_XRQ HA412_PSC8 HA412_RHA438 HA412_IR HA412_HA89 HA412_LR1 HA412_OQP8 HA412_HA300
do
$bin/show-snps -ClrT ../$cmp.b1000.c1000.i90.l1000.1delta > ../$cmp.b1000.c1000.i90.l1000.1snps
awk 'NR>4{print $11"\t"$1"\t"$12"\t"$4"\t"$2"\t"$3}' ../$cmp.b1000.c1000.i90.l1000.1snps > $cmp.mod.1snps
done
# SNP set
for cmp in HA412_XRQ HA412_PSC8 HA412_RHA438 HA412_IR HA412_HA89 HA412_LR1 HA412_OQP8 HA412_HA300
do
echo '##fileformat=VCFv4.2' > $cmp.snps.vcf
echo "##INFO=<ID=QRY_POS,Number=1,Type=String,Description=\"Position on the query genome\">" >> $cmp.snps.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" >> $cmp.snps.vcf
awk '$5!="."&&$6!="."{print $1"\t"$2"\t.\t"$5"\t"$6"\t.\tPASS\tQRY_POS="$3":"$4"\t."}' $cmp.mod.1snps >> $cmp.snps.vcf
bgzip $cmp.snps.vcf; tabix $cmp.snps.vcf.gz
done
# InDels set
for cmp in HA412_XRQ HA412_PSC8 HA412_RHA438 HA412_IR HA412_HA89 HA412_LR1 HA412_OQP8 HA412_HA300
do
echo '##fileformat=VCFv4.2' > $cmp.indels.vcf
echo "##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"Length of the InDel\">" >> $cmp.indels.vcf
echo "##INFO=<ID=V_TYPE,Number=1,Type=String,Description=\"Type of the InDel\">" >> $cmp.indels.vcf
echo "##INFO=<ID=QRY_POS,Number=1,Type=String,Description=\"Position on the query genome\">" >> $cmp.indels.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" >> $cmp.indels.vcf
cat $cmp.mod.1snps | perl mummer_indels.pl | sort -k 1,1 -k 2,2n >> $cmp.indels.vcf
grep '^#' $cmp.indels.vcf > $cmp.indels.len50.vcf # small InDels (<50bp)
grep -v '^#' $cmp.indels.vcf | awk -F '[\t;=]' '$9<=50{print}' >> $cmp.indels.len50.vcf
bgzip $cmp.indels.len50.vcf; tabix $cmp.indels.len50.vcf.gz
done
```

Combine results from each genome

```
bcftools merge *.snps.vcf.gz --info-rules QRY_POS:join -o all.snps.vcf
bcftools merge *.indels.len50.vcf.gz --info-rules LEN:join,V_TYPE:join,QRY_POS:join -o all.indels.len50.vcf
```

Identify large SVs with SyRI

```
mkdir ./syri; cd ./syri/
bash run_syri.sh
```

Merge SVs

```
Rscript merge_SV.R # > all.svs.inversions.tsv
bedtools merge -i all.svs.inversions.tsv -c 1,4 -o count,collapse -delim ";" > all.svs.inversions.count
awk '$3-$2>1000000' all.svs.inversions.count > all.svs.inversions.1Mb.count # 37 large inversions
bedtools merge -i all.svs.translocations.tsv -c 1,4 -o count,collapse -delim ";" > all.svs.translocations.count
```

Identify PAVs/CNVs with SVMU

```
mkdir ./svmu; cd ./svmu/
svmu="/home/hkchi/program/svmu/svmu"
genome1=HA412
for genome2 in XRQ PSC8 RHA438 IR HA89
do
cmp="${genome1}_${genome2}"
prefix="$cmp.b1000.c1000.i90.l1000"
fasta1=$(head -1 ../$prefix.mdelta | cut -d " " -f1)
fasta2=$(head -1 ../$prefix.mdelta | cut -d " " -f2)
$svmu ../$prefix.mdelta $fasta1 $fasta2 l lastz $cmp
done
```
