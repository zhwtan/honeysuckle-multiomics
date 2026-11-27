#!/usr/bin/perl -w 
#"Chr01   889246  895565  Chr08   38303407        38297084        color=chr1"
open IN,$ARGV[0];
while(<IN>){
	chomp;
	@l=split;
	$len1=abs($l[2]-$l[1]);
	$len2=abs($l[5]-$l[4]);
	if ($len1 >=$ARGV[1] && $len2 >=$ARGV[1]){print "$_\n";}
}
