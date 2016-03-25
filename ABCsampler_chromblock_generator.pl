#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;

######## To Do ########
# add option for file with specfic mutation rates for particular contigs
# add option for file with specfic recombination rates for particular contigs
######################

my $version = '0.010';

die (qq/
Version $version

Description: 
Use to format a .par input file for DNA sequence simulation with fastsimcoal 
so that each contig represents a block on a single chromosome. 
This script only works for a single chromosome.

Usage: ABCsampler_chromblock_generator.pl [options]

Options: 
-p FILE  ABCsampler.par file with chromosome block info to be appended
-f FILE  fasta file of contigs (i.e. chromosome blocks)
-o CHAR  outfile name preffix
-r FLOAT recombination rate between adjacent contigs and within contigs
-s FILE  list of recombination rates for particular contigs (Contig\tRate)
-u FLOAT mutation rate
-v FILE  list of mutation rates for particular contigs (Contig\tRate)
-t FLOAT transition rate (0.33 implies no transition bias) [0.33]
-w FILE  list of transition rates for particular contigs (Contig\tRate)\n
\n/) if !@ARGV;

my %opts = (
	p => undef,
	f => undef,
	r => 0.5, #need to find
	s => undef,
	u => 0.00001, #need to find
	v => undef,
	o => undef,
	t => 0.33,
	w => undef
	);

getopts ('p:f:r:s:u:v:o:t:w:', \%opts);

open(PAR, '<', $opts{p});
open(OUT, '>', "$opts{o}.par");
# write original .par file to output file
while (<PAR>) {
	print OUT $_;
}
close PAR;

#define global variables
my (
$name, 
$cont_number, 
$prev_number, 
$prev_name
);

my $seq_length = 0;

#initialize values
open(FASTA, '<', $opts{f});my $first_line = <FASTA>;
$prev_name = $1 if ($first_line =~ />(\w+)\b/);
$prev_number = $1 if ($first_line =~ />[A-Za-z]*(\d+)/);

while (<FASTA>) {
	if ($_ =~ />(\w+)\b/) {
		$name = $1;
		$cont_number = $1 if ($_ =~ />[A-Za-z]*(\d+)/);
		if ($cont_number != $prev_number) {
			print OUT "DNA $seq_length $opts{r} $opts{u} $opts{t} //$prev_name\n";
			$seq_length = 0;
		}
		$prev_number = $cont_number;
		$prev_name = $name;
	} else {
		my $seq = $_;
		chomp $seq;
		$seq =~ s/N//gi;
		$seq_length += length($seq);
		print OUT "DNA $seq_length $opts{r} $opts{u} $opts{t} //$prev_name" if (eof);
	} 
}
close FASTA;
close OUT;
