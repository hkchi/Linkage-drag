#!/bin/bash

## Run SyRI to identify SVs from MUMmer outputs
## Kaichi Huang 2022 Mar

#
echo "Starting run at: `date`" >> syri.log

for cmp in HA412_XRQ HA412_PSC8 HA412_RHA438 HA412_IR HA412_HA89
do
	prefix="$cmp.b1000.c1000.i90.l1000"
	genome1=$(head -1 ../$prefix.mdelta | cut -d " " -f1)
	genome2=$(head -1 ../$prefix.mdelta | cut -d " " -f2)
	python ~/program/syri/syri/bin/syri \
		-r $genome1 -q $genome2 \
		-c ../$prefix.mcoords -d ../$prefix.mdelta \
		--lf $prefix.m.log --log DEBUG \
		-k \
		--nc 17 \
		-s /home/hkchi/program/mummer-4.0.0beta2/show-snps \
		--prefix $prefix.m.
	awk '$11=="INV" {print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8}' $prefix.m.syri.out > $prefix.m.inversions.tsv
	awk '$11=="TRANS" || $11=="INVTR" {print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8}' $prefix.m.syri.out > $prefix.m.translocations.tsv
done

#
echo "Program finished with exit code $? at: `date`" >> syri.log
