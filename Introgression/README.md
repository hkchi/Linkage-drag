# Introgression

Convert VCF

```
ref="HA412" # XRQ PSC8 RHA438 IR HA89 LR1 OQP8 HA300
vcf="${ref}.vcf.gz"
while read chr
do
bcftools view -r ${chr} $vcf | perl vcf2str.pl sample_info.3ref.txt ${ref}.${chr}
cat ${ref}.${chr}.txt.tmp | python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" | sed 's/\.//g' > ${ref}.${chr}.txt
rm ${ref}.${chr}.txt.tmp
done < chr_list.${ref}.txt
```

Run the 'site-by-site' linkage admixture model

```
structure="/home/hkchi/program/STRUCTURE/structure_kernel_src/structure"
while read chr
do
n_loci=$(awk 'NR==2{print NF}' ${ref}.${chr}.txt)
n_loci=$(($n_loci-3))
if [ -s ${ref}.str_cmd.txt ]; then rm ${ref}.str_cmd.txt; fi
echo "$structure -m mainparams.txt -e extraparams.txt -K 7 -L $n_loci -N 49 -i ${ref}.${chr}.txt -o ${ref}.${chr}.str.out" >> ${ref}.str_cmd.txt
done < chr_list.${ref}.txt
parallel -j 17 < ${ref}.str_cmd.txt
```
