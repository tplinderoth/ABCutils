#!/usr/bin/perl

#fastsimcoal_sampler.pl

use warnings;
use strict;
use Getopt::Long;
#use File::Path qw(remove_tree);

#########################################################################
# TO DO:
# add more prior sampling distributions (can only do unif and logunif now)
##########################################################################

my $version = '1.1.4'; # 10/30/2015

die (qq/
Version $version

Usage: 
fastsimcoal_sampler.pl [options]

fastsimcoal must be installed and in path

Options:\n 
--outfile	CHAR		Outfile name preffix (NO PATH)\n
--estfile	FILE		Fastsimcoal .est format file ("rules" section last)\n
--tplfile	FILE		Fastsimcoal .tpl format file\n
--recomb	FLOAT		Recombination rate between adjacent nucleotides within sequence blocks [0]\n
--mut		FLOAT\/CHAR 	Mutation rate as float or parameter name from .est file [2.2e-09]\n 
--trans		FLOAT		Transition rate (0.33 implies no mutation bias) [0.33]\n
--seqlist	FILE		Fasta format file of sequence blocks (list of contigs)\n
--numseq	INT		Number of sequence blocks to simulate at a time (NA when meta2DSFS = 1) [500]\n
--numsim	INT		Number of simulations (each sim uses different parameter values) [1000]\n
--meta2DSFS	0|1		Calculate 2D-SFS for two metapopulations if INT = 1 (each comprising >= 1 demes)\n	
--pop1		INT 		First population ID number in joint SFS (from top of tpl file, 0-based)
				If meta2DSFS = 1, these are the >= 1 demes comprising metapopulation 1 (1-based)\n 
--pop2		INT		Second population ID number in joint SFS (from top of tpl file, 0-based)
				If meta2DSFS = 1, these are the >= 1 demes comprising metapopulation 2 (1-based)\n
--norm		0|1		Normalize SFS if INT = 1, else keep absolute counts if INT = 0 [0]\n
--rmvfixed	0|1		Remove fixed categories in SFS if INT = 1 (prior to standardization if --norm 1), 
				else keep fixed categories if INT = 0 [0]\n
--rmvMutation	CHAR		Remove SNP sites of certain mutation type(s); only works when meta2DSFS = 1\n
--p1missing	FLOAT		Discard site if at least FLOAT percent of individuals in pop1 are missing data at the site. 
				Otherwise, missing data will be treated as major\/ancestral allele [50]\n
--p2missing 	FLOAT		Discard site if at least FLOAT percent of individuals in pop2 are missing data at the site.
				Otherwise, missing data will be treated as major\/ancestral allele [50]\n
Notes:\n 
The .tpl input file:
-Should not include fastsimcoal "genetic info" section; file should end after "historical event" section\n
The .est input file:
-should list sections in the order [parameters], [complex parameters], [conditional], [rules]
-All parameter values drawn from priors or calculated must only contain uppercase letters and
 '_' in the name
-Parameters must be defined in .est file before other parameters that depend on them are defined
-Loguniform priors should not have zero as a boundary\n
The output file:
-Names for unfolded SFS indicate element in the 2k+1 SFS: E<pop1_sfs_category><pop2_sfs_category>\n
Calculating 2DSFS with metapopulations (meta2DSFS = 1):
-pop1 and pop2 arguments each take a list of >= 1 INT values for the demes comprisiong each
 metapopulation, respectively, in the arlequin project file (ex: --pop1 1 2 3 --pop2 4 5 6).
 The deme numbers correspond to the ascending order in which they are listed in the .tpl file (1-indexed).
-The type of mutation to remove is denoted by <major_allele><minor_allele>; the major allele is considered
 ancestral and the minor allele is derived. Example to remove C->T and G->A mutations: --rmvMutation CT GA
\n/
) if !@ARGV;

my @command_log = @ARGV;

## check input

my $options = '--outfile --estfile --tplfile --recomb --mut --trans --seqlist --numseq --numsim --pop1 --pop2 --norm --rmvfixed --meta2DSFS --rmvMutation --p1missing --p2missing ';
foreach (@ARGV) {
	if ($_ =~ /^\-\-/) {
		die("ERROR: Unknown command $_\n") if ($options !~ /$_\s+/);
	}
}


my $outfile = undef;
my $estfile = undef;
my $tplfile = undef;
my $recomb = 0;
my $mut = 2.2e-09;
my $trans = 0.33;
my $seqlist = undef;
my $numseq = 500;
my $numsim = 1000;
my @pop1;
my @pop2;
my @rmvMutation;
my $norm = 0;
my $rmvfixed = 0;
my $meta2DSFS = 0;
my $p1missing = 50;
my $p2missing = 50;

GetOptions
(
'outfile=s' => \$outfile,
'estfile=s' => \$estfile,
'tplfile=s' => \$tplfile,
'recomb=f' => \$recomb,
'mut=s' => \$mut,
'trans=f' => \$trans,
'seqlist=s' => \$seqlist,
'numseq=i' => \$numseq,
'numsim=i' => \$numsim,
'pop1=i{1,}' => \@pop1,
'pop2=i{1,}' => \@pop2,
'norm=i' => \$norm, 'rmvfixed=i' => \$rmvfixed,
'rmvMutation=s{0,}' => \@rmvMutation,
'meta2DSFS=i{1,1}' => \$meta2DSFS,
'p1missing=f' => \$p1missing,
'p2missing=f' => \$p2missing
); 

my ($pop1, $pop2);
if ($meta2DSFS == 0) {
	$pop1 = $pop1[0];
	$pop2 = $pop2[0];
}

die("ERROR: --recomb must be numeric\n") if ($recomb !~ /\d+/);
die("ERROR: --trans must be numeric\n") if ($trans !~ /\d+/);
die("ERROR: --numsim must be numeric\n") if ($numsim !~ /\d+/);
die("ERROR: --numseq must be numeric\n") if ($numseq !~ /\d+/);
if (!@pop1 || !@pop2) {
	die("ERROR: Need to define what populations in tpl file to examine for constructing the 2D-SFS\n");
} else {
	my $p1arg = join('', @pop1);
	my $p2arg = join('', @pop2);
	die("ERROR: --pop2 takes only integer values\n") if ($p2arg !~ /\d+/);
	die("ERROR: --pop1 takes only integer values\n") if ($p1arg !~ /\d+/);
	undef $p1arg;
	undef $p2arg;	
}
if ($meta2DSFS == 0) {
	die("ERROR: --pop1 and --pop2 can only take multiple values if meta2DSFS = 1") if (scalar(@pop1) > 1 || scalar(@pop2) > 1);
}
die ("ERROR: --norm must be 0 or 1\n") unless ($norm == 0 || $norm == 1);
die("ERROR: --rmvfixed must be 0 or 1\n") unless ($rmvfixed == 0 || $rmvfixed == 1);
die("ERROR: --mut must be a character string without whitespace or float\n") if ($mut =~ /\s+/);
die("ERROR: --meta2DSFS must be 0 or 1\n") unless ($meta2DSFS == 0 || $meta2DSFS == 1);
die("ERROR: --p1missing must be a float in (0, 100]\n") if ($p1missing <= 0 || $p1missing > 100);
die("ERROR: --p2missing must be a float in (0, 100]\n") if ($p2missing <= 0 || $p2missing > 100);

## set up global variables, hashes, arrays
my %params;
my @param_names;
my %priors;
my %complex;
my %conditional;
my @rules;
my $end_fasta = 0;
my $mu;
my @sfs;
my $simnum = 0;
my (@complex_ord, @conditional_ord);
my ($pop1n, $pop2n);

# set up outfile prefixes
my $filename = $1 if ($tplfile =~ /(.+).tpl/);
my ($tplname, $directory) = ($1, $`) if ($tplfile =~ /\/([^\/]+).tpl$/);

# start log

my $time_specs = localtime;
my $logaddress = "$directory/$outfile.log";
open(LOG, '>', $logaddress) or die("ERROR: Could not open LOG file $logaddress: $!\n-->exiting\n"); 
print LOG "fastsimcoal_sampler version $version\n\n--outfile $outfile --estfile $estfile --tplfile $tplfile --recomb $recomb --mut $mut --trans $trans --seqlist $seqlist --numseq $numseq --numsim $numsim --pop1 @pop1 --pop2 @pop2 --norm $norm --rmvfixed $rmvfixed --meta2DSFS $meta2DSFS --p1missing $p1missing --p2missing $p2missing ";
print LOG @rmvMutation ? "--rmvMutation @rmvMutation\n\n" : "--rmvMutation NONE\n\n";
print LOG "START: $time_specs\n";
open(my $fastafh, '<', $seqlist) or die("ERROR: Could not open FASTA file $seqlist: $!\n-->exiting\n");

### main ###
if ($meta2DSFS == 0) {
	# get sample sizes for populations in the joint SFS
	($pop1n, $pop2n) = sample_n($pop1, $pop2);
	# get rules
	getrules();
	# open output file
	my $sampaddress = "$directory/$outfile.samp";
	open(my $samplefh, '>', $sampaddress) or die("ERROR: Could not open SAMPLE file $sampaddress: $!\n-->exiting\n");
	# main loop
	while ($simnum < $numsim) {
		# get parameter values
		good_params();
		# replace parameter names with values in .tpl file to create .val file
		get_par_input();
		# get mutation rate, etc, for genetic info
		genetic_rates();
		seek $fastafh, 0, 0;
		while ($end_fasta == 0) {
			# add genetic info to .val file to create a .prepar and .par file
			genetic_info();
			# run simulations and update sfs
			my $return = simulate();
			die ("ERROR: Failed simulation\n") if ($return != 0); # TO DO: Make better exception handling
			# delete files
			delete_files();
		}
		modify_sfs(); # perform operations on sfs prior to output
		# write header for output if first iteration
		if ($simnum == 0) {
			write_head($samplefh, \@param_names);  
			# print joint SFS categories
			unfolded_header($pop1n, $pop2n, $samplefh);
		}
		write_output($samplefh, \@param_names); # write parameters and summary stats
		$end_fasta = 0; #reset fasta processing loop
		undef @sfs; 
		$simnum++;
		# track progress of simulations
		print "$simnum simulations done\n" if ($simnum % 1000 == 0);
	}
	close $samplefh;
} elsif ($meta2DSFS == 1) {
	# get rules
	getrules();
	# prepare .par file template 
	genetic_rates();
	my $simtpl = prepParFile($directory, $fastafh, $recomb, $mut, $trans, $tplfile);
	#perform simulation iterations
	open(my $par_template, '<', $simtpl) or die("ERROR: Could not open .par template file $simtpl: $!\n-->exiting\n");
	open(my $sample, '>', "$directory/$outfile.samp") or die("ERROR: Could not open SAMPLE file $directory/$outfile.samp: $!\n-->exiting\n");
	my $parvals;
	my ($p1demes, $p2demes) = (':', ':');
	foreach (@pop1) {
		$p1demes .= "$_:";
	}
	foreach (@pop2) {
		$p2demes .= "$_:";
	}
	my $rmv_str;
	if (@rmvMutation) {
		$rmv_str = join(':', @rmvMutation);
		die("ERROR: Invalid mutation type for removal!\n") if $rmv_str =~ /[^ACGT:]/i;
	}	
	printMetaHeader($estfile, $tplfile, $sample, $p1demes, $p2demes, \@param_names);
	while ($simnum < $numsim) {
		# replace parameter variables with values
		good_params();
		$parvals = addParVal($directory, $par_template);
		# simulate
		my $simreturn;
		if ($rmv_str) {
			$simreturn = system("fastsimcoal -i $parvals -n 1 -q");
		} else {
			$simreturn = system("fastsimcoal -i $parvals -n 1 -s -q");
		}
		if ($simreturn != 0)
		{
			print STDERR "WARNING: fastsimcoal run failed -> discarding run\n";
			#eraseOldSim($directory);
			delete_files();
			next;
		}
		# calculate 2D-SFS
		if (!-e "$directory/modelpar_1_1.arp")
		{
			print STDERR "Something went wrong with this simulation --> discarding it ...\n";
			delete_files();
			next;
                }
		open(my $arpfile, '<', "$directory/modelpar_1_1.arp") or die ("ERROR: Couldn't open arlequin project file $directory/modelpar_1_1.arp: $!\n-->exiting\n");
		# make sure arlequin file actually has something in it
		if (-z "$directory/modelpar_1_1.arp") {
			#eraseOldSim($directory);
			delete_files();
			next;
		}
		my ($p1mat, $p2mat) = getMatrix($arpfile, $p1demes, $p2demes, $rmv_str, $p1missing, $p2missing);
		close $arpfile;
		if ($p1mat eq 'N') {
			#eraseOldSim($directory);
			delete_files();
			next;
		} elsif ($p1mat == 1) {
			print STDERR "WARNING: Discarding bad simulation\n";
			#eraseOldSim($directory);
			delete_files();
			next;
		}
		my $ncol = scalar(@$p1mat);
		my $nrow = scalar(@$p2mat);
		my $sfs2d = calc2dsfs($p1mat, $p2mat, $ncol, $nrow);
		undef $p1mat;
		undef $p2mat;
		# output simulation results
		printMetaSim($sfs2d, $sample, $simnum, \@param_names);
		# delete this iteration's files
		#eraseOldSim($directory);
		delete_files();
		# track simulation progress
		$simnum++;
	}
	close $par_template;
	unlink "$directory/modelpar.template" or die ("ERROR: Couldn't delete the modelpar.template file $directory/modelpar.template\n-->exiting\n");
} else
	{
		print STDERR "Unrecognized argument for --meta2DSFS: $meta2DSFS\n-->exiting...\n";
	}
# end log

$time_specs = localtime;
print LOG "FINISHED $simnum simulations: $time_specs\n";
close LOG;
close $fastafh;

print "Finished $simnum simulations!\n";
exit;

### subroutines ###

sub good_params {
	my $rule_fail;
	do {
		$rule_fail = 0;
		getpriors();
		getparams();
		if (@rules) {
			foreach (@rules) {
				$rule_fail++ if ( eval($_) == 0);
			}
		}
	} until ($rule_fail == 0);
}

sub getrules {
	open(EST, '<', $estfile) or die("ERROR: Could not open EST file $estfile: $!\n-->exiting\n");
	while (<EST>) {
		if ($_ !~ /^\[COMPLEX PARAMETERS\]/i) {
			next;
		} else {
			last;
		}
	} 
	while (<EST>) {
		if ($_ !~ /^\[CONDITIONAL\]/i) {
			next if ($_ !~ /^\d/);
			my @complex_line = split (/\s+/, $_);
			push @complex_ord, $complex_line[1];
		} else {
			last;
		}
	}
	while (<EST>) {
		if ($_ !~ /^\[RULES\]/i) {
			next if ($_ !~ /^\d/);
			my @condition_line = split (/\s+/, $_);
			push @conditional_ord, $condition_line[1];
		} else {
			last;
		}
	}
	while (<EST>) {
		if ($_ =~ /\w+/) {
			my $rules = $_;
			chomp $rules;
			$rules =~ s/([A-Z|_|0-9]+)/\$params{$1}/ig;
			push @rules, $rules;
		}
	}
close EST;
}

sub getpriors {
	open(EST, '<', $estfile) or die("ERROR: Could not open EST file $estfile: $!\n-->exiting");
	my @est_line;
	# get 'Parameters' section
	while (<EST>) {
		if ($_ !~ /\[COMPLEX PARAMETERS\]/i) {
			if ($_ =~ /^[1|0]\s+/) {
				@est_line = split(/\s+/, $_);
				$priors{$est_line[1]} = [$est_line[0], @est_line[2..4]];
			}	
		} else {
			last;
		}
	} 
	# get 'Complex Parameters' section
	while (<EST>) { 
		if ($_ !~ /\[CONDITIONAL\]/i) { 
			if ($_ =~ /^[1|0]\s+/) {
				 my ($type, $comp_name, $equation) = ($1, $2, $3) if ($_ =~ /(\d+)\s+(\w+)\s*=\s*(.+)/);
				 $equation =~ s/([A-Z|_]+)/\$params{$1}/g;
				 $complex{$comp_name} = [$type, $equation];
			}	
		} else {
			last;
		}
	}
	# get 'Conditional' section
	while (<EST>) {
		if ($_ !~ /\[RULES\]/i) {
			if ($_ =~ /^[1|0]/) {
					my $statement = $_;
					my ($con1, $test, $con2, $val, $integer, $equation);
					if ($statement =~ /(\d)\s+(\S+)\s*=\s*(\S+)\s+if\s+(\S+\s*\S+\s*\S+)\s+else\s+(\S+)/) {
						$integer = $1;
						$val = $2;		
						$con1 = $3;
						$test = $4;
						$con2 = $5;
					}
				$test =~ s/([A-Z|_]+)/\$params{$1}/g;
				$conditional{$val} = [$test, $con1, $con2, $integer];
			}
		} else {
			last;
		}
	}
close EST;
}			 		


sub getparams {
	my $error = 0;
	# draw parameter values from priors
	foreach ( keys(%priors) ) 
	{
		$params{$_} = unif_sample(@{$priors{$_}}[0 .. 3], \$error);
		die ("-->exiting\n") if ($error > 0);
	}

	# calculate complex parameters
	if (%complex) {
		foreach (@complex_ord) {
			$params{$_} = eval(${$complex{$_}}[1]);
			$params{$_} = int($params{$_}) if ( ${$complex{$_}}[0] == 1);
		}
	}
	
	# determine conditional paramters
	if (%conditional) {
		foreach(@conditional_ord) {
			$params{$_} = eval(${$conditional{$_}}[0]) ? ${$conditional{$_}}[1]:${$conditional{$_}}[2];
			if ($params{$_} =~ /[A-Z_]+/) { # parameter depends on other parameters
				if ($params{$_} =~ /log|exp|\*|\/|\+|-/) { # parameter is a function of other parameters
						$params{$_} =~ s/([A-Z_]+)/\$params{$1}/g;
						$params{$_} = ${$conditional{$_}}[3] == 0 ? eval($params{$_}) : int(eval($params{$_}));		
				} else { # parameter equals the value of a different parameter's value
					$params{$_} = $params{$params{$_}};	
				}
			}
		}
	}

		sub unif_sample {
		my ($type, $log, $min, $max, $error) = @_;
		my $priortype = 0;
		my $absdist = 0;
		my $range = 0;
		if ($log =~ /logunif/i) 
		{
			if ($min > $max)
			{
				$error = 1;
				print STDERR "minimum prior value $min > maximum prior value $max in 'unif_sample' subroutine\n";
				return 0;
			}
			elsif ($min == 0 || $max == 0)
			{
				$error = 1;
				print STDERR "Zero set as parameter boundary for logunif distribution\n";
				return 0;
			}
			elsif ($min > 0 && $max > 0)
			{
				$priortype = 1;
				$range = log($max) - log($min);
			}
			elsif ($min < 0 && $max < 0)
			{
				$priortype = 2;
				$range = log(-1*$min) - log(-1*$max);
				$min = -1*$max;
			}
			elsif ($min < 0 && $max > 0)
			{
				$priortype = 3;
				$absdist = -1*$min;
				$min = 1;
				$range = log($max + $absdist + $min);
			}
		}
		else
		{
			$range = $max - $min;
		}
		my $u = 0;
		do { $u = rand() } while ($u == 0);
		my $draw;
		if ($priortype == 0)
		{
			$draw = $min + $u * $range; # draw from uniform
		}
		else
		{
			$draw = $min * exp($u * $range); # draw from log uniform
			if ($priortype == 2)
			{
				$draw *= -1;
			}
			elsif ($priortype == 3)
			{
				$draw -= $absdist + $min;
			}
		}
		$draw = int($draw) if ($type == 1);
		return ($draw);
	}
}

sub get_par_input {
	open(TPL, '<', $tplfile) or die("ERROR: Could not open TPL file $tplfile: $!\n-->exiting\n");
	open(VAL, '>', "$filename.val") or die("ERROR: Could not open VAL file $filename.val: $!\n-->exiting\n");
	# get values for parameters in the .tpl file
	while (<TPL>) {
		my $l = $_;
		if ($l =~ /^\/\//) {
			print VAL $l;
			next;
		}
		foreach ( keys(%params) ) {
			$l =~ s/$_/$params{$_}/g;	
		}
		print VAL $l;
	}
	close TPL;
	close VAL;
}

sub genetic_info {
	my $fafh;
	# add genetic information to the .val file created by get_par_input
	open(VAL, '<', "$filename.val") or die("ERROR: Could not open VAL file $filename.val: $!\n-->exiting\n");
	open(PREPAR, '>', "$filename.prepar") or die("ERROR: Could not open PREPAR file $filename.prepar: $!\n-->exiting\n");
	# copy contents of .val to .prepar file
	while (<VAL>) {
		print PREPAR $_;
	}
	close VAL;	
	# add genetic information to .par file
	my $processed = 0;
	my $seq_length = 0;
	print PREPAR '//Number of indepdent loci [contigs]',"\n",'CHUNK_SIZE_TMP 1', "\n";
	while ($processed < $numseq && !eof($fafh)) {	 
		my $fasta = <$fafh>;
		if ($fasta =~ /^>/ && $seq_length != 0) {
			print PREPAR '//Per contig: Number of linkage blocks',"\n","1\n";
			print PREPAR '//Per block: data type, num sites, rec. rate, mut. rate, tran. rate',"\n";
			print PREPAR "DNA $seq_length $recomb $mu $trans\n";
			$processed++;
			$seq_length = 0;
			next;
		} elsif ($fasta !~ /^>/) {
			chomp(my $seq = $fasta);
			$seq =~ s/N//gi;
			$seq_length += length($seq);
			if (eof($fafh)) {
				print PREPAR '//Per contig: Number of linkage blocks',"\n","1\n";
				print PREPAR '//Per block: data type, num sites, rec. rate, mut. rate, tran. rate',"\n";
				print PREPAR "DNA $seq_length $recomb $mu $trans\n";
				$processed++;
			} 
		}
	}
	close PREPAR;
	# finalize .par file
	open(PAR, '>', "$filename.par") or die("ERROR: Could not open PAR file $filename.par: $!\n-->exiting\n");
	open(PREPAR, '<', "$filename.prepar") or die("ERROR: Could not open PREPAR file $filename.prepar: $!-->exiting\n");
	while (<PREPAR>) {
			$_ =~ s/CHUNK_SIZE_TMP/$processed/;
			print PAR $_;
	}
	close PREPAR;
	close PAR;
	# check for the end of the fasta
	if (eof($fafh)) {
		$end_fasta = 1;	
	}
}

sub genetic_rates {
	# check parameter hash for mutation rate
	if ($mut =~ /\d+/) {
		$mu = $mut;
	} else {
		$mu = $params{$mut}
	}
}

sub simulate {
	# run fastsimcoal
	my $returnval;
	$returnval = system("fastsimcoal -i $filename.par -x -s -d -n 1 -q");
	return 1 if ($returnval != 0);
	# read joint SFS file
	open(SFS, '<', "$filename\_jointDAFpop$pop2\_$pop1.obs") or die("ERROR: Could not open SFS file $filename\_jointDAFpop$pop2\_$pop1.obs: $!\n-->exiting\n");
	if (@sfs) {
		# add new sfs to existing sfs
		my @sfstemp;
		while (<SFS>) {
			next if (1 ..2);
			chomp ($_);
			my @rowtemp = split(/\s+/, $_);
			push @sfstemp, @rowtemp[1 .. $#rowtemp];
		}
		@sfs = map{$sfs[$_] + $sfstemp[$_]} (0 .. $#sfs);
	} else {
		# construct new sfs vector
		while (<SFS>) {
			# concatenate rows of 2dsfs to vectorize, leaving out row/column headers
			next if (1 .. 2);
			chomp($_);
			my @row = split(/\s+/, $_);
			push @sfs, @row[1 .. $#row];
		} 	
	}
	close (SFS);
}

sub modify_sfs {
	
	# remove fixed sfs categories
	@sfs = (@sfs[1 .. ($pop1n-1)], @sfs[($pop1n+1) .. ($#sfs-($pop1n+1))], @sfs[($#sfs-($pop1n-1)) .. ($#sfs-1)]) if ($rmvfixed == 1);
	
	# fold sfs - requires the same sample sizes for the 2 populations - DISABLED
	#@sfs = ((map{$sfs[$_] + $sfs[$#sfs - $_]} (0 .. (($#sfs/2)-1))),$sfs[$#sfs/2]) if ($opts{d} eq 'm'); # error with folding - use ABCUtils.pl
	
	# normalize sfs
	if ($norm == 1) {
		my $total_counts = 0;
		foreach (@sfs) {
			$total_counts += $_;
		}
		@sfs = map{$sfs[$_] / $total_counts} (0 .. $#sfs);
	}
}

sub delete_files {
	# delete fastsimcoal output files
	opendir(SIMDIR, $directory) or die("ERROR: Could not open SIMDIR directory $directory: $!\n-->exiting\n");
	while (readdir(SIMDIR)) {
		# check that file is simulation output file
		if ($_ =~ /$tplname|modelpar/i) {
			my $file = $directory . "/$_"; #complete filename w/ directory
			if (-d $file) { # check for directory that fastsimcoal creates 
				my $del_success = rmdir $file;
				if ($del_success) {
					next;
				} else {
					opendir(NEWDIR, $file) or die("ERROR: Could not open NEWDIR directory $file: $!\n-->exiting\n");
					# remove all files in NEWDIR
					while (readdir(NEWDIR)) {
						my $newfile = $file . "/$_";
						unlink $newfile;	
					}
					# remove emptied directory created by fastsimcoal
					close NEWDIR;
					rmdir($file);
					next;
				}
			} elsif (-f $file) {
				unlink $file if ($file !~ /\.tpl$|\.est$|\.log$|\.template$|\.samp$/);
			} 
		}
	}
	close SIMDIR;
}

sub write_head {
	my ($outfile, $parnames) = @_;
	print $outfile "Sim ";
	# print parameter names
	foreach ( @$parnames ) 
	{
		print $outfile "$_ ";
	}
}

sub unfolded_header {
	my ($p1, $p2, $outfh) = @_;
	for (my $i = 0; $i <= $p2; $i++) {
		for (my $j = 0; $j <= $p1; $j++) {
			if ($rmvfixed == 1) { # fixed removed
				(print $outfh $j == $p1-1 && $i == $p2 ? ("E$j$i\n") : ("E$j$i ")) unless (($i == 0 && $j == 0) || ($i == $p2 && $j == $p1) || ($i == 0 && $j == $p1) || ($i == $p2 && $j == 0));	
			} else { # fixed kept
			print $outfh $j == $p1 && $i == $p2 ? ("E$j$i\n") : ("E$j$i ");
			}		
		}
	}
}

sub printMetaHeader {
	my ($est, $tpl, $sampfh, $deme1pops, $deme2pops, $parnames) = @_;
	print $sampfh "Sim ";
	my @param_names;
	open(ESTFILE, '<', $est) or die ("ERROR: Can't open est file $est: $!\n-->exiting\n");
	while (<ESTFILE>) 
	{
		last if ($_ =~ /^\[RULES\]/i);
		if ($_ =~ /^[0|1]\s+(\w+)\s+\S+/) 
		{
			push @$parnames, $1;
		}
	}
	close ESTFILE;
	foreach (@$parnames) 
	{
		print $sampfh "$_ ";
	}
	my ($meta1n, $meta2n) = (0, 0);
	my $section = 0;
	open(TPLFILE, '<', $tpl) or die ("ERROR: Can't open tpl file $tpl: $!\n-->exiting");
	while (<TPLFILE>) {
		$section++ if $_ =~ /^\/\/\S*/;
		if ($section == 3) {
			my $subpop = 1;
			my $estline;
			do { $estline = <TPLFILE>;
				if ($estline =~ /^(\d+)\/*/) {
					my $demen = $1;
					if ($demen != 0) {
						if ($deme1pops =~ /:$subpop:/) {
							$meta1n += $demen;
						} elsif ($deme2pops =~ /:$subpop:/) {
							$meta2n += $demen;
						}
					}
					$subpop++;
				}
			} until ($estline =~ /^\/\/\S*/);	
			$section++;
		}
	last if $section == 4;
	}
	close TPLFILE;
	unfolded_header($meta1n, $meta2n, $sampfh);
}

sub sample_n {
	my ($p1, $p2) = @_;
	open (TPL, '<', $tplfile) or die("ERROR: Could not open TPL file: $!");
	my ($p1n, $p2n);
	my $section = 0;
	while (<TPL>) {
		$section++ if ($_ =~ /^\/\//);
		if ($section == 3) {
			my $size_num = 0;
			my $est_line;
			do {
				$est_line = <TPL>;
				my $est_pop = $1 if ($est_line =~ /^(\d+)/);
				$p1n = $est_pop if ($size_num == $p1);
				$p2n = $est_pop if ($size_num == $p2);
				$size_num++;
			} until ($est_line =~ /^\/\//);
			last;
		} else {
			next;
		}
	}
	close TPL;
	return ($p1n, $p2n);
}

sub write_output {
	my ($outfile, $parnames) = @_;
	# write sim number
	print $outfile "$simnum ";
	# write parameters to output
	foreach ( @$parnames ) 
	{
		print $outfile "$params{$_} "
	}
	# write SFS
	$" = " ";
	print $outfile "@sfs\n";
}

sub prepParFile {
	my ($outdir, $fastafh, $recomb, $mutation, $transbias, $tpl) = @_;
	
	my %genetic_info;
	my $contig_num = 0;
	my $fasta_line;
	
	my $siminfo = $outdir . '/modelpar.template';
	open(PARTEMP, '>', $siminfo);
	open(TEMPLATE, '<', $tpl);
	while (<TEMPLATE>) {
		print PARTEMP $_;
	}
	close(TEMPLATE);
	
	my $fasta_name;
	do {chomp($fasta_line = <$fastafh>); $fasta_name = $1 if $fasta_line =~ /^>(\S+)/} until ($fasta_line =~ /^>\S+/ || eof($fastafh));
	die ("ERROR: Sequences have no names in fasta file\n") unless $fasta_name;
	undef $fasta_name;
	while (<$fastafh>) {
		if ($_ !~ /^>\S+/) {
			chomp(my $fasta_line = $_);
			my $seq = $fasta_line;
			if (!eof($fastafh)) {
				do {chomp($fasta_line = <$fastafh>); $seq .= $fasta_line if $fasta_line !~ /^>\S+/} until ($fasta_line =~ /^>\S+/ || eof($fastafh));
			}
			$contig_num++;
			$genetic_info{$contig_num} = length($seq);	
		}
	}
	
	print PARTEMP "// begin genetic info\n$contig_num 1\n";
	for (my $i = 1; $i <= $contig_num; $i++) {
		print PARTEMP "//\n1\n//\nDNA $genetic_info{$i} $recomb $mutation $transbias\n";
	}

	close PARTEMP;
	return ($siminfo);
}

sub addParVal {
	my ($outdir, $prepar) = @_;
	my $par_complete_address = $outdir . '/modelpar.val';
	open(PARFINAL, '>', $par_complete_address);
	seek $prepar, 0, 0;
	while (<$prepar>) {
		if ($_ =~ /^\/\/ begin genetic info/) {
			print PARFINAL $_;
			while (<$prepar>) {
				print PARFINAL $_;
			}
		} else {
			my $parline = $_;
			if ($parline =~ /^\/\//) {
				print PARFINAL $parline;
				next;
			}
			map { $parline =~ s/$_/$params{$_}/g } keys(%params);
			print PARFINAL $parline;	
		}
	}
	close PARFINAL;
	return($par_complete_address);
}

sub getMatrix {
	my ($arp, $deme1, $deme2, $remove, $miss_cutoff1, $miss_cutoff2) = @_;
	my $meta1 = [];
	my $meta2 = [];
	my ($arpline, $snp);
	my ($row1, $row2) = (0, 0);
	# make sure there are polymorphic sites
	while (<$arp>) {
		next unless $_ =~ /SampleData=/i;
		$arpline = <$arp>;
		$snp = $1 if $arpline =~ /^\d+_\d+\s+\d+\s+(.+)$/;
		return('N', 'N') if (!$snp);
		if ($remove) {
			$snp =~ s/[^ACGT]/N/ig;
		} else { 
			$snp =~ s/[^01]/N/g;
		} 
		return('N', 'N') unless $snp =~ /\w/;
		last;	
	}
	# make sure arlequin project file is complete
	seek $arp, -3, 2;
	return (1, 1) if (<$arp> !~ /^\s+\}$/);
	# now okay to proceed
	seek $arp, 0, 0;
	do {$arpline = <$arp>} until ($arpline =~ /SampleName=/);
	my $demenum = $1 if $arpline =~ /SampleName="Sample\s+(\d+)"$/;
	while (<$arp>) {
		if ($deme1 =~ /:$demenum:/) {
			do {$arpline = <$arp>;
				if ($arpline =~ /^\d+_\d+\s+\d+\s+(.+)$/) {			
					$snp = $1;
					return('N', 'N') if (!$snp);
					if ($remove) {
						$snp =~ s/[^ACGT]/N/ig;
					} else { 
						$snp =~ s/[^01]/N/g;
					}
					$$meta1[$row1] = [split(//, $snp)];
					$row1++;
					return (1, 1) if (eof($arp)); # premature end of file; bad simulation
				}		
			} until ($arpline =~ /^}/);		
			$demenum = 'NA';
		} elsif ($deme2 =~ /:$demenum:/) {
			do {$arpline = <$arp>;
				if ($arpline =~ /^\d+_\d+\s+\d+\s+(.+)$/) {			
					$snp = $1;
					if ($remove) {
						$snp =~ s/[^ACGT]/N/ig;
					} else { 
						$snp =~ s/[^01]/N/g;
					}
					$$meta2[$row2] = [split(//, $snp)];
					$row2++;
					return (1, 1) if (eof($arp)); # premature end of file; bad simulation 				
				}		
			} until ($arpline =~ /^}/);
			$demenum = 'NA';
		}
	$demenum = $1 if $_ =~ /SampleName="Sample\s+(\d+)"$/;
	last if $_ =~ /^\[\[structure\]\]/i;
	return(1, 1) if (eof($arp)); # premature end of file; bad simulation
	}
	my $pop1missing = 0;
	my $pop2missing = 0;
	my %allele_count;
	my @sorted_counts;
	my @rmv_index;
	my ($mutation, $major);
	# consider major allele ancestral (0) and minor allele derived (1) - mispolarization won't matter for folded spectrum
	my $siteidx = scalar(@{$$meta1[0]}) - 1;
	for (my $i = 0; $i <= $siteidx ; $i++) 
	{
		if ($remove)
		{
			%allele_count = ('A' => 0, 'C' => 0, 'G' => 0, 'T' => 0, 'N' => 0);
		}
		else
		{
			%allele_count = ('1' => 0, '2' => 0, 'N' => 0);
		}
		map { $allele_count{$$meta1[$_]->[$i]}++} (0 .. $#$meta1);
		$pop1missing = $allele_count{'N'};
		map { $allele_count{$$meta2[$_]->[$i]}++} (0 .. $#$meta2);
		$pop2missing = $allele_count{'N'} - $pop1missing;
		# check for excess missing data
		if ( 100 * ($pop1missing/($#$meta1 + 1)) >= $miss_cutoff1 || 100 * ($pop2missing/($#$meta2 + 1)) >= $miss_cutoff2)
		{
			push @rmv_index, $i;
			next;
		}
		if ($remove)
		{
			@sorted_counts = sort { $allele_count{$a} <=> $allele_count{$b} } keys(%allele_count);
			$mutation = $sorted_counts[4] . $sorted_counts[3];
			# remove sites with particular substitutions
			if ($remove =~ /$mutation/i && $allele_count{$sorted_counts[4]} != $allele_count{$sorted_counts[3]}) 
			{
				push @rmv_index, $i;
				next;
			}
			# check if site is diallelic
			if ($allele_count{$sorted_counts[2]} > 0) 
			{
				push @rmv_index, $i;
				next;
			}
			# check for erroneous fixed site
			if ($allele_count{$sorted_counts[3]} == 0)
			{
				push @rmv_index, $i;
				next; 
			}
			$major = $sorted_counts[4];
			map { $$meta1[$_]->[$i] = ($$meta1[$_]->[$i] eq $major || $$meta1[$_]->[$i] eq 'N') ? 0 : 1} (0 .. $#$meta1);
			map { $$meta2[$_]->[$i] = ($$meta2[$_]->[$i] eq $major || $$meta2[$_]->[$i] eq 'N') ? 0 : 1} (0 .. $#$meta2);
		}
	}
	# remove sites
	if (@rmv_index)
	{
		my ($k, $index);
		for (my $i = 0; $i <= $#$meta1; $i++) 
		{
			$k = 0;
			map { $index = $_ - $k; splice(@{$$meta1[$i]}, $index, 1); $k++ } @rmv_index; 
		}
		for (my $j = 0; $j <= $#$meta2; $j++) 
		{
			$k = 0;
			map { $index = $_ - $k; splice(@{$$meta2[$j]}, $index, 1); $k++ } @rmv_index;
		}     
	}
	return($meta1, $meta2);
}	

sub calc2dsfs {
	my ($meta1, $meta2, $meta1n, $meta2n) = @_;
	my $maxn = $meta2n < $meta1n ? $meta1n : $meta2n; 
	# initiate 2dsfs matrix
	my $sfsmat = [];
	for (my $i = 0; $i <= $meta2n; $i++) {
		map {$$sfsmat[$i]->[$_] = 0} (0 .. $meta1n);
	}
	#calculate 2dsfs matrix by tabulation
	my $nsites = scalar(@{$$meta1[0]});
	my ($p1cat, $p2cat) = (0, 0);
	for (my $i = 0; $i <= $nsites - 1; $i++) {
		map { $p1cat += $$meta1[$_]->[$i] } (0 .. $#$meta1);
		map { $p2cat += $$meta2[$_]->[$i] } (0 .. $#$meta2);
		$$sfsmat[$p2cat]->[$p1cat]++;
		$p1cat = 0;
		$p2cat = 0;
	}
	return($sfsmat);
}

sub printMetaSim {
	my ($sfs, $output, $simn, $parnames) = @_;
	print $output "$simn ";
	foreach ( @$parnames ) 
	{
		print $output "$params{$_} ";
	}
	my $lastindex = $#$sfs;
	my $iter = 0;
	foreach my $row (@$sfs) {
        print $output $iter == $lastindex ? "@$row\n" : "@$row ";
        $iter++;
	}
	
}

sub eraseOldSim {
	my $simdir = $_[0];
	if (-e "$simdir/modelpar.val")
	{
		unlink "$simdir/modelpar.val" or die ("ERROR: Couldn't delete the .val file\n");
	}
	if (-e "$simdir/modelpar_1.arb")
	{
		unlink "$simdir/modelpar_1.arb" or die ("ERROR: Couldn't delete the .arb file\n");
	}
	if (-e "$simdir/modelpar_1.simparam")
	{
		unlink "$simdir/modelpar_1.simparam" or die ("ERROR: Couldn't delete the .simparam file\n");     
	}
	if (-e "$simdir/modelpar_1_1.arp")
	{
		unlink "$simdir/modelpar_1_1.arp" or die ("ERROR: Couldn't delete arp file\n ");
	}
	if (-e "$simdir/modelpar" && -d "$simdir/modelpar")
	{
		#remove_tree("$simdir/modelpar") or die ("ERROR: Couldn't delete directory dumped by fastsimcoal\n");
	}
}
