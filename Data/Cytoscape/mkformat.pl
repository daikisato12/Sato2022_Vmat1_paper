use warnings;
use strict;

#perl mkformat.pl CytoScape_result.txt CytoScape_result_format.tsv
open(IN,$ARGV[0]);
open(OUT,">$ARGV[1]");
print OUT "Cluster\tSymbol\n";
while(<IN>){
chomp;
	if($_=~/^[1-2]/){
		my @li = split(/\t/, $_);
		my @genes = split(/, /, $li[4]); 
		foreach my $gene (@genes){
			print OUT "$li[0]\t$gene\n";
		}
	}
}
