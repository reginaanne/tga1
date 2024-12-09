use strict;
use warnings;

my $file=$ARGV[0];
my $cutoff=$ARGV[1];

open(FILE, '<', $file) or die "no file dufus\n";
open(OUT, '>', "$file.consensus.vcf");
open(BAD, '>', "$file.consensus.bad.txt");
my $space=0;
my $mult_allele=0;
my $neither=0;
my $snps=0;
while(<FILE>){
	if($_=~m/\#\#/){ print OUT $_; next; }
	chomp;
	my @line=split("\t",$_);
	if($line[0] eq "#CHROM"){ print OUT join("\t", @line[0..8]), "\t", "CON\n"; next}
	$snps++;
	my $chrom=$line[0];
	my $start=$line[1];
	my $ref=$line[3];
	my $alt=$line[4];
	my @alt=split(',',$alt);
	$alt=(split(",",$alt))[0];
	my $ref_count=0; 
	my $alt_count=0;
	my $consensus="";
	for my $nam (9 .. 28){
		my @geno=split(/\|/,$line[$nam]);
		if($geno[0]==0){ $ref_count++; }
		else{ $alt_count++; }	
		if($geno[1]==0){ $ref_count++; }
		else{ $alt_count++; }	
	}
	if($alt_count/($alt_count+$ref_count)>$cutoff){
		$consensus="1|1";
	}
	elsif($ref_count/($alt_count+$ref_count)>$cutoff){
		$consensus="0|0";
	}
	else{$consensus=".|."; $neither++; print BAD "$_\n";}
	if($#alt>0){$consensus="-999|-999"; $mult_allele++; print BAD "$_\n";}
	print OUT join("\t", @line[0..8]),"\t",$consensus,"\n";	
}
print "DONE! Processed $snps SNPs. Of these, $mult_allele had multiple alleles and were dropped, and $neither had no allele that matched your threshold of $cutoff. These SNPs were coded as missing .|. and the original data were printed to $file.consensus.bad.txt.\n";
