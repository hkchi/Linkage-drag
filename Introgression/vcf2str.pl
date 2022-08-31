## Transforms VCF to STRUCTURE input format
## Kaichi Huang 2021 Jun

use warnings;
use strict;
use List::MoreUtils qw(uniq);

my $out_prefix = $ARGV[1];
my %pop_info;
my @fields;

my %pop_num = ('Landrace'=>'1',
			'ANN'=>'2',
			'ARG'=>'3', #(*o_o*)#
			'PetPet'=>'4',
			'PetFal'=>'5',
			'NIV'=>'6',
			'DEB'=>'7');

open POP, $ARGV[0];
while (<POP>){
	chomp;
	my @a = split(/\t/, $_);
	$pop_info{$a[0]} = $pop_num{$a[1]};
}
close POP;

my $current_pos;
open OUT, ">$out_prefix.txt.tmp";
open OUT_SNP, ">$out_prefix.snp";
while (<STDIN>) {
	if(/^##/){next;}
	chomp;
	if(/^#/){
		my @a = split(/\t/, $_);
		my @label = ('Ref','Ref');
		my @pop = ('0','0');
		my @flag = ('0','0'); #(-^-)#
		foreach my $i (9..$#a) {
			push @label, ($a[$i],$a[$i]);
			push @pop, ($pop_info{$a[$i]},$pop_info{$a[$i]});
			push @flag, ('1','1');
		}
		print OUT ". ", join(" ",@label),"\n";
		print OUT ". ", join(" ",@pop),"\n";
		print OUT ". ", join(" ",@flag),"\n";
		next;
	}
	my @a = split(/\t/, $_);
	my $chr = $a[0];
	my $pos = $a[1];
	my $dis;
	if ($current_pos) {
		$dis = $pos - $current_pos;
	} else {
		$dis = -1;
	}
	my @gt_line = ('0','0');
	foreach my $i (9..$#a) {
		my @tmp = split(/:/,$a[$i]);
		my $gt = $tmp[0];
		# Copyright: Kaichi Huang
		if ($gt eq "0/0" || $gt eq "0|0") {
			push @gt_line, ('0','0');
		} elsif ($gt eq "1/1" || $gt eq "1|1") {
			push @gt_line, ('1','1');
		} elsif ($gt eq "0/1" || $gt eq "0|1") {
			push @gt_line, ('0','1');
		} else {
			push @gt_line, ('9','9');
		}
	}
	print OUT "$dis ",join(" ",@gt_line),"\n";
	print OUT_SNP "$chr\t$pos\n";
	$current_pos = $pos;
}
close OUT;
close OUT_SNP;
