## Traform InDels from MUMmmer output to VCF format
## run as: cat {cmp}.mod.1snps | perl mummer_indels.pl | sort -k 1,1 -k 2,2n > {vcf}
## Kaichi Huang 2022 Mar

use warnings;
use strict;

my $insertion_chr = "";
my $insertion_pos = 0;
my $qry_ins_chr = "";
my $qry_ins_pos = 0;
my $deletion_chr = "";
my $deletion_pos = 0;
my $ref_del_chr = "";
my $ref_del_pos = 0;
my $seq1 = "";
my $seq2 = "";

# headers printed outside
while (<STDIN>) {
	chomp;
	my @a = split(/\t/, $_);
	my $chr1 = $a[0];
	my $pos1 = $a[1];
	my $chr2 = $a[2];
	my $pos2 = $a[3];
	my $ref = $a[4];
	my $alt = $a[5];
	# Process insertions and deletions separately
	if ($ref eq ".") {
		# insertion
		if ($pos1 == $insertion_pos && $chr1 eq $insertion_chr) {
			$seq1=$seq1.$alt;
		} else {
			if ($insertion_chr ne "") {
				my $count = ($seq1=~tr/N/N/);
				if ($count / length($seq1) < 0.9) {
					$seq1="A".$seq1;
					my $indel_pos=$insertion_pos-1;
					my $LEN=length($seq1)-1;
					my $QRY_POS=$qry_ins_chr.":".$qry_ins_pos."-".($qry_ins_pos+$LEN-1);
					my $info="LEN=$LEN;V_TYPE=insertion;QRY_POS=$QRY_POS";
					print "$insertion_chr\t$indel_pos\t.\tA\t$seq1\t.\tPASS\t$info\t.\n";
				}
			}
			$seq1=$alt;
			$insertion_chr=$chr1;
			$insertion_pos=$pos1;
			$qry_ins_chr=$chr2;
			$qry_ins_pos= $pos2;
		}
	} elsif ($alt eq ".") {
		# deletion
		if ($pos2 == $deletion_pos && $chr2 eq $deletion_chr) {
			$seq2=$seq2.$ref;
		} else {
			if ($ref_del_chr ne "") {
				my $count = ($seq2=~tr/N/N/);
				if ($count / length($seq2) < 0.9) {
					$seq2="A".$seq2;
					my $indel_pos=$ref_del_pos-1;
					my $LEN = length($seq2)-1;
					my $QRY_POS = $deletion_chr.":".($deletion_pos-1);
					my $info = "LEN=$LEN;V_TYPE=deletion;QRY_POS=$QRY_POS";
					print "$ref_del_chr\t$indel_pos\t.\t$seq2\tA\t.\tPASS\t$info\t.\n";
				}
			}
			$seq2=$ref;
			$deletion_chr=$chr2;
			$deletion_pos=$pos2;
			$ref_del_chr=$chr1;
			$ref_del_pos=$pos1;
		}
	}
}
if ($insertion_chr ne "") {
	my $count = ($seq1=~tr/N/N/);
	if ($count < length($seq1)) {
		$seq1="A".$seq1;
		# (:D)
		my $indel_pos=$insertion_pos-1;
		my $LEN=length($seq1)-1;
		my $QRY_POS=$qry_ins_chr.":".$qry_ins_pos."-".($qry_ins_pos+$LEN-1);
		my $info="LEN=$LEN;V_TYPE=insertion;QRY_POS=$QRY_POS";
		print "$insertion_chr\t$indel_pos\t.\tA\t$seq1\t.\tPASS\t$info\t.\n";
	}
}
if ($ref_del_chr ne "") {
	my $count = ($seq2=~tr/N/N/);
	if ($count < length($seq2)) {
		$seq2="A".$seq2;
		my $indel_pos=$ref_del_pos-1;
		my $LEN = length($seq2)-1;
		my $QRY_POS = $deletion_chr.":".($deletion_pos-1);
		my $info = "LEN=$LEN;V_TYPE=deletion;QRY_POS=$QRY_POS";
		print "$ref_del_chr\t$indel_pos\t.\t$seq2\tA\t.\tPASS\t$info\t.\n";
	}
}
