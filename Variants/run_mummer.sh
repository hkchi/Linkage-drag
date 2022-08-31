#!/bin/bash

## Run MUMmer to align genomes
## Kaichi Huang 2022 Mar

genome1=$1 # /scratch/hkchi/ref/Ha412HO.onlychr.fasta
genome2=$2
cmp=$3

JOBINFO=$cmp.mummer.log
#
echo "Starting run at: `date`" >> $JOBINFO

bin="/home/hkchi/program/mummer-4.0.0beta2"

$bin/nucmer $genome1 $genome2 -p $cmp.b1000.c1000 -b 1000 -c 1000 -t 48

$bin/delta-filter -1 -i 90 -l 1000 $cmp.b1000.c1000.delta > $cmp.b1000.c1000.i90.l1000.1delta
$bin/show-coords -THrd $cmp.b1000.c1000.i90.l1000.1delta > $cmp.b1000.c1000.i90.l1000.1coords

$bin/delta-filter -m -i 90 -l 1000 $cmp.b1000.c1000.delta > $cmp.b1000.c1000.i90.l1000.mdelta
$bin/show-coords -THrd $cmp.b1000.c1000.i90.l1000.mdelta > $cmp.b1000.c1000.i90.l1000.mcoords

#
echo "Program finished with exit code $? at: `date`" >> $JOBINFO
