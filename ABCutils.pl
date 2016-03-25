#!/usr/bin/perl

#ABCutils.pl

use warnings;
use strict;
use Getopt::Std;
use Getopt::Long;
use Cwd qw(getcwd abs_path);
use File::Basename;
use Scalar::Util qw(openhandle);
use File::Path qw(rmtree);

my $version = "1.4.1"; # as of 6/15/2015

### define variables to avoid 'will not stayed shared' warnings ###
my $fixed;
my @oldfiles;
my @newfiles;

&main;
exit;

sub main {
    &usage if (@ARGV<1);
    my $command = shift(@ARGV);
    my %fun = 
    	(
    	catsims=>\&catsims,
    	MaskCats=>\&MaskCats,
    	rmvFixedZero=>\&rmvFixedZero,
    	SetSpace=>\&SetSpace,
    	DistReject=>\&DistReject,
    	format2Dsim=>\&format2Dsim,
    	format2Dobs=>\&format2Dobs,
    	gofReject=>\&gofReject,
    	best2Dsfs=>\&best2Dsfs,
    	fold2D=>\&fold2D,
    	StatDistr=>\&StatDistr,
    	);
    die ("Unknown command \"$command\"\n") if (!defined($fun{$command}));
    &{$fun{$command}};
}


########### main functions ############

sub catsims {

die (qq/
Usage: ABCutils.pl catsims <folder_containing_sims> <simulation_file>\n
This script loops through a directory of folders or files containing
simulations and appends the simulation output to the "simulation_file"\n
Note: -The folder_containing_sims can be a parent directory
      -If simulation_file does not exist, it will be created
\n/) if @ARGV == 0;
die "folder containing sims is not a directory!\n" if (!-d $ARGV[0]);
$ARGV[0] =~ s/\/$//;
die "file for appending is not a file!\n" if (-e $ARGV[1] && !-f $ARGV[1]);
open(MASTERFILE, '>>', $ARGV[1]) or die("ERROR: Can't open output file: $!"); # file to append sims to
opendir(MAIN, $ARGV[0]) or die("ERROR: Can't open directory containing simulations: $!"); # directory containing sims
my $simcount = 0;
while (readdir(MAIN)) {   
    my $dir = "$ARGV[0]/$_";
    next if ($dir =~ /$ARGV[0]\/\.+$/);
    if (-d $dir) {   
        $simcount = readSubs($dir, $simcount);
    } else {
        $simcount = writeSims($dir, *MASTERFILE, $simcount) if ($_ =~ /\.(samp){1}(fold)?$/i && $dir !~ /$ARGV[1]/); 
    }   
}
close MAIN;
close MASTERFILE;

### catsims subroutines ###

sub readSubs {
    my ($dirpath, $simnum) = @_;
    my @subs;
    my @subs2 = "$dirpath/*";
    my $sub2size = @subs2;
    while ($sub2size > 0) {
        foreach (@subs2) {
            @subs = shift @subs2;
            foreach (<@subs>) {
                if (-d $_) {
                    push @subs2, "$_/*";
                    shift @subs if $subs[0];
                } else {
                    $simnum = writeSims($_, *MASTERFILE, $simnum) if ($_ =~ /\.(samp){1}(fold)?$/i && $_ !~ /$ARGV[1]/);
                }
            }
        }
        $sub2size = @subs2;
    }
    return ($simnum);
}

sub writeSims {
    my ($simfile, $master, $simnum) = @_;
    open(SIMS, '<', $simfile) or die("ERROR: Can't open simulation file: $!");
    if (-z $master) {
    	select((select($master), $|=1)[0]); # make MASTERFILE hot
    	my $head = <SIMS>;
    	print $master $head;
    	select((select($master), $|=0)[0]); # make MASTERFILE not hot
    } else {
        <SIMS>; # skip header if masterfile already has it
    }
    while (<SIMS>) {
        my $simline = $' if ($_ =~ /^\d+\s+/);
        print $master "$simnum $simline";
        $simnum++;
    }
    close SIMS;
    return ($simnum);
}

}


sub format2Dobs {
die (qq/
Description: Formats the 2D-SFS for use as the observed SFS input file for DistReject\n
Usage: ABCutils.pl format2Dobs [options]\n
Options:
--insfs       FILE    2D-SFS input file
--outfile     FILE    output file name including path
--fold        0|1     fold joint SFS if 1 [0]
--rmvFixed    0|1     remove invariable categories in joint spectrum if 1 [0]
--maskFixed   0|1     mask invariable categories in joint spectrum if 1 [0]
--norm	      0|1     normalize 2D-SFS if 1 [0]
--bin	      INT     calculate diagonal bins of 2D-SFS categories; INT is 2 x the number of categories on
                      either side of each bin's diagonal (~ bin width, which must be a factor of 2)
--do1dsfs     INT     calculate folded 1D-SFS for pop1 (INT = 0), pop2 (INT = 1), or both populations (INT = 2) [null]\n  
NOTES:
-format2Dobs assumes that the input agrees with fastsimcoal joint sfs format: 
 pop1 categories are COLUMNS and pop2 categories are ROWS
-the input 2D-SFS should be the same format as 2dsfs.pl output
-if no outfile name is provided the output is dumped as a .fmt file to the same directory
 as the input
\n/) if !@ARGV;   
   	# UNFOLDED PROCESSING NOT TESTED!!!
    # check input
    my $options = '--insfs --outfile --fold --rmvFixed --maskFixed --norm --bin --do1dsfs ';
    foreach (@ARGV) 
    {
        if ($_ =~ /^\-\-/) 
        {
            die("ERROR: Unknown command $_\n") if ($options !~ /$_\s+/);
        }
    }
    
    # get options
    my $insfs = undef;
    my $outfile = undef;
    my $fold = 0;
    my $norm = 0;
    my $rmvFixed = 0;
    my $maskFixed = 0;
    my $bin = undef;
    my $do1dsfs = undef;
   
    GetOptions
    (
    'insfs=s'=>\$insfs,
    'outfile=s'=>\$outfile,
    'fold=i'=>\$fold,
    'rmvFixed=i'=>\$rmvFixed,
    'maskFixed=i'=>\$maskFixed,
    'norm=i'=>\$norm,
    'bin=i' =>\$bin,
    'do1dsfs=i'=>\$do1dsfs
    );

    # open files and check inputs
    die ("ERROR: --fold can only take an argument of 0 or 1\n") unless ($fold == 0 || $fold == 1);
    die ("ERROR: --rmvFixed can only take an argument of 0 or 1\n") unless ($rmvFixed == 0 || $rmvFixed == 1);
    die ("ERROR: --maskFixed can only take an argument of 0 or 1\n") unless ($maskFixed == 0 || $maskFixed == 1);
    die ("ERROR: --rmvFixed and --maskFixed cannot both be set to 1\n") if ($rmvFixed == 1 && $maskFixed == 1);
    die ("ERROR: --norm can only take an argument of 0 or 1\n") unless ($norm == 0 || $norm == 1);
    if (defined $bin)
    {
    	die ("ERROR: If binning, --bin must be zero or a positive factor of 2\n") if ($bin > 0 && $bin % 2 != 0 || $bin < 0);
    }
    if (defined $do1dsfs)
	{
		die("ERROR: --do1dsfs takes only arguments of 0, 1, or 2 for calculating the 1D-SFS\n") if ($do1dsfs < 0 || $do1dsfs > 2 || $do1dsfs =~ /[^0-2]/);
	}
    open(my $sfs, '<', $insfs) or die("ERROR: Can't open input sfs file: $!");
    if ($outfile) {
    	open(OUT, '>', $outfile) or die("ERROR: Can't open outfile: $!");
    } else {
    	open(OUT, '>', "$insfs.fmt") or die("ERROR: Can't open default output: $!");
    }

	# read in joint SFS and start folding
	my $sfsmat = [];
	my @sfsvec;
	my $pop2n = -1;
	while (<$sfs>) 
	{
		chomp;
		my @row = split(/\s+/, $_);
		$pop2n++;
		if ($fold == 0) 
		{
			push @sfsvec, \@row;
		} 
		else 
		{
			@row = ((map { $row[$_] + $row[$#row - $_] } ( 0 .. (($#row/2)-1))), $row[$#row/2]);
			@{$$sfsmat[$pop2n]} = @row;
		}
	}
	close $sfs;	
	#finish folding
	if ($fold == 1) 
	{
		map { do{ my $j = $_; push @sfsvec, [map { $$sfsmat[$j][$_] + $$sfsmat[$#$sfsmat - $j][$_] } 0 .. $#{$$sfsmat[0]}] }; } 0 .. ($#$sfsmat/2)-1;
		push @sfsvec, \@{$$sfsmat[$#$sfsmat/2]}; 
	}
	my $pop1n = $#{$sfsvec[0]};
	undef $sfsmat;
	
	# remove invariable categories in joint SFS -- untested!
	if ($rmvFixed == 1) 
	{
		if ($fold == 1) 
		{
			shift @{$sfsvec[0]};
		} 
		else 
		{
			@{$sfsvec[0]} = @{$sfsvec[0]}[1 .. $pop1n-1];
			@{$sfsvec[$pop2n]} = @{$sfsvec[$pop2n]}[1 .. $pop1n-1];
		}
	}
	
	# mask invariable categories in joint SFS as 0 -- untested for unfolded!
	if ($maskFixed == 1) 
	{
		if ($fold == 1) 
		{
			${$sfsvec[0]}[0] = 0;
		} 
		else 
		{
			(${$sfsvec[0]}[0], ${$sfsvec[0]}[$pop1n], ${$sfsvec[$pop2n]}[0], ${$sfsvec[$pop2n]}[$pop1n]) = (0, 0, 0, 0);
		}
	}
	
	# normalize joint SFS -- untested!
    if ($norm == 1)
    {
        my $total = 0;
        map { do { my $i = $_; map {$total += ${$sfsvec[$i]}[$_] } 0 .. $#{$sfsvec[$i]}  };  } 0 .. $#sfsvec;
        map { do { my $i = $_; map {${$sfsvec[$i]}[$_] /= $total} 0 .. $#{$sfsvec[$i]} }; } 0 .. $#sfsvec;
    }
    
    ## calculate 1D-SFS ##
    my @sfsout = ();
    if (defined $do1dsfs)
    {
    	if ($do1dsfs == 0)
    	{
    		calc1DSFS(\@sfsvec, \@sfsout, 0);
    	}
    	elsif ($do1dsfs == 1)
    	{
    		calc1DSFS(\@sfsvec, \@sfsout, 1);
    	}
    	elsif ($do1dsfs == 2)
    	{
    		calc1DSFS(\@sfsvec, \@sfsout, 0);
    		calc1DSFS(\@sfsvec, \@sfsout, 1);
    	}
    	else
    	{
    		die("ERROR: '$do1dsfs' is an unknown argument for calculating 1D-SFS\n");
    	}
    }
    
    my ($diagonal, $offdiagonal) = (0, 0);
    if (defined $bin)
    {
    	$offdiagonal = binJointSFS(\@sfsvec, \@sfsout, $bin); # bins along off-diagonal
    	if (!defined $do1dsfs)
    	{
    		my @rev = reverse @sfsvec;
    		$diagonal = binJointSFS(\@rev, \@sfsout, $bin); # bins along diagonal
    	}
    }
    else
    {
    	map { push(@sfsout, @{$sfsvec[$_]}) } 0 .. $#sfsvec; 
    }
   
    
	# print header
	my $k = 0;
	if (defined $do1dsfs)
	{
		if ($do1dsfs == 0)
		{
			for ($k = 0; $k <= $pop1n; $k++)
			{
				print OUT "FX$k ";
			}
		}
		elsif ($do1dsfs == 1)
		{
			for ($k = 0; $k <= $pop2n; $k++)
			{
				print OUT "FY$k ";
			}
		}
		elsif ($do1dsfs == 2)
		{
			for ($k = 0; $k <= $pop1n; $k++)
			{
				print OUT "FX$k ";
			}
			for ($k = 0; $k <= $pop2n; $k++)
			{
				print OUT "FY$k ";
			}
		}
		else
		{
			die("ERROR: '$do1dsfs' is an unknown argument for calculating 1D-SFS\n");
		} 
	}
	
	if (defined $bin)
	{
		map {print OUT "OB$_ "} 1..$offdiagonal;
		map {print OUT "DB$_ "} 1..$diagonal if (!defined $do1dsfs);
	}
	
	#print joint SFS or bins
	print OUT "\n@sfsout\n";
	close OUT;
}


sub MaskCats {
die (qq/
Usage: ABCutils.pl MaskCats [options]\n
Force nonvariable 2D-SFS categories to be variable by changing the value of the nonvariable
category in a randomly chosen simulation to 1e-100. This is in order to run ABCestimator on 2DSFS
simulations for which all values in a category among simluations is 0.\n
Options:
-i    FILE    simulation input file; should be .samp format outputted by fastsimcoal_sampler
-o    CHAR    outfile name including path
-n    INT    iterate through INT simulations to determine if a 2DSFS category is fixed for 0,
        otherwise, if INT = 0, iterate through all simulations [0]
\n/) if @ARGV == 0;

    my %opts = (i=>undef, o=>undef, n=>0);
    getopts('i:o:n:', \%opts);

    open(SIMS, '<', $opts{i});
   
    # find first field containing SFS info
   
    my $field = split(/\s+/, $`) if (<SIMS> =~ /\s+E\d+E*\d*\s+/);
   
    # find the total number of fields
   
    my $ncol = split(/\s+/, <SIMS> );
   
    #find the number of simulations
   
    seek (SIMS, 0, 0);
    my $nsims = 0;
    while (<SIMS>) {
        $nsims++;
    }
    $nsims--; # account for zero index and header
   
    # find which categories are nonvariable
   
    my @nonvar;
    for (my $i = $field; $i < $ncol; $i++) {
        my $randsim = int(rand()*($nsims-1) + 1 );
        push(@nonvar, "$i\_$randsim");
        my $iter = 0;   
        seek(SIMS, 0, 0);
        <SIMS>;
        while (<SIMS>) {
            my @sim = split(/\s+/, $_);
            if ($sim[$i] != 0) {
                pop @nonvar;
                last;
            }
       
            if ($opts{n} != 0) {
                $iter++;
                last if ($iter = $opts{n});   
            }
        }
        print STDERR "$i categories done...\n";
    }
    my $novar = join(" ", @nonvar);
    undef @nonvar;
   
    # change 0 values when necessary and print
   
    open(OUT, '>', $opts{o});
    seek(SIMS, 0, 0);
    $" = "\t";
    $. = -1; # reset and account for header
    while (<SIMS>) {
        if ($novar =~ /_$.\b/) {
            my @match = ($novar =~ /\b(\d+)_$.\b/g);
            chomp;
            my @line = split(/\s+/, $_);
            foreach (@match) {
                $line[$_] = 1e-100;
            }
            $novar =~ s/\b\d+_$.\b//g;
            print OUT "@line\n";
        } else {
            print OUT $_;
        }
    }   
       
    close SIMS;
    close OUT;
}

sub rmvFixedZero {
die (qq/
The rmvFixedZero command removes summary statistics fixed for 0 from simulation files and the .obs file
used for ABCestimator. If a summary statistic is fixed for 0 in at least one model, this statistic
is removed from all simulation files and the .obs file.\n
Usage: ABCutils.pl rmvFixedZero -n INT -c INT <.obs_file> <simulation_file(s)>\n
Notes: - New .obs and simulation file(s) are dumped in the directory of the original files
         with 'estimator' appended to the name
       - The INT for option -n specifies how many simulations in each sim file to read
         through in order to determine which summary statistics are fixed. The default is
     to read all simulations.
       - The INT for option -c specifies that if a summary statistic has INT or fewer simulations with
     a nonzero value, the statistic is considered fixed and removed.  
\n/) if !@ARGV;

my %opts = (n=>0, c=>7);
getopts('n:c:', \%opts);

# define variables
# some variables used for this major subroutine are defined as global variables at top of script
my $fileref = \@newfiles;
@oldfiles = @ARGV;

# make array of new file names
filename($fileref);

# make array of fixed summary stats names
print STDERR "Finding summary statistics fixed for zero...\n";
$fixed = getfixed($opts{n}, $opts{c});
print STDERR "Done finding fixed summary statistics!\n";

# remove fixed statistics and output new files;
print STDERR "Removing fixed summary statistics from .obs and simulation files...\n";
rmvStats($fileref);

print STDERR "Finished!\n";

sub filename {
    # get working directory
    my $fileref = $_[0];
    my $working = getcwd;
    $working = $working . '/';
    my $order = 1;
    foreach (@ARGV) {
        my ($path, $name);
        if ($_ =~ /([^\/]+)$/) {
            $path = $` ne '' ? $` : $working;
            $name = $path . $1;
            $name =~ s/\.\w+$//; # remove redundant suffix
            $name = $order == 1 ? $name . '_estimator.obs' : $name . '_estimator.txt';
        }
        push @{$fileref}, $name;
        $order++;
    }
}

sub getfixed {
    my $limit = $_[0];
    my %sumstats;
    my $entry = 1;
    foreach (@ARGV[1 .. $#ARGV]) {
        my $file = $_;
        open(SIM, '<', $file);
        my $param = -1;
        $. = 0;
        while (<SIM>) {
            last if ($. > $limit && $limit != 0);
            my @sim;
            # get header information
            if ($. == 1) {
                my $pos = 0;
                @sim = split(/\s+/, $_);
                foreach (@sim) {
                    if ($_ =~ /^E\d+E*\d*$/) {
                        ${$sumstats{$pos}}[0] = $_;
                        push @{$sumstats{$pos}}, 0;
                        $pos++;
                    } else {
                        $param++;
                    }           
                }
            } else {
                # find fixed sum stats   
                @sim = split(/\s+/, $_);
                my ($i, $k) = (0, 0);
                foreach (@sim) {
                    if ($i <= $param) {
                        $i++;
                        next;
                    }
                    ${$sumstats{$k}}[$entry]++ if ($_ != 0);
                    $k++;
                }
            }
        }
        print STDERR "finished searching $file...\n";
        close SIM;
        $entry++;
    }
    my $cutoff = $_[1];
    my $fix = '';
    foreach my $a (keys %sumstats) {
        foreach ( @{$sumstats{$a}}[1..$#{$sumstats{$a}}] ) {
            if ($_ <= $cutoff) {
                $fix = $fix . "${sumstats{$a}[0]} ";
                last;
            }
        }
    }
    return($fix);
}

sub rmvStats {
    my $file;
    my $j = 0;
    $" = " ";
    foreach (@oldfiles) {
        open(INFILE, '<', $_);
        # create array of fixed sum stats fields
        my @zeros;
        my $zeroref = \@zeros;
        my $i = 0;
        foreach ( split(/\s+/, <INFILE>) ) {
            push @{$zeroref}, $i if $fixed =~ /$_ /;
            $i++;
        }
        my $arrlen = @zeros;
        my $zerocut = setsplice($zeroref, $arrlen);
        # remove fixed sum stats
        seek INFILE, 0, 0;
        $file = $newfiles[$j];
        open(OUT, '>', $file);
        print STDERR "Writing file $file...\n";
        while (<INFILE>) {
            chomp;
            my $prevsize = 0;
            my @line = split(/\s+/, $_);
            foreach (@{$zerocut}) {
                my ($index, $size) = (($1 - $prevsize), $2) if ($_ =~ /(\d+)_(\d+)/);
                splice @line, $index, $size;
                $prevsize += $size;
            }
            print OUT "@line\n";
        }
    close OUT;
    close INFILE;
    $j++;
    }   
}

sub setsplice {
    my $field = $_[0];
    my $len = $_[1];
    my %numcut;
    my $group = ${$field}[0];
    for (my $i=0; $i < $len; $i++) {
        my $j = $i + 1;
        if ($j == $len) {
            $numcut{$group}++;
            last;
        } elsif (${$field}[$j] == ${$field}[$i] + 1) {
            $numcut{$group} = 0 if !exists $numcut{$group};
            $numcut{$group}++;
        } else {
            $numcut{$group}++;
            $group = ${$field}[$j];
            $numcut{$group} = 1;
        }
    }
    my @z;
    my $zref = \@z;
    foreach ( sort { $a <=> $b } keys %numcut ) {
        push @{$zref}, "$_\_$numcut{$_}";
    }
    return $zref;
}
}

sub SetSpace {
die (qq/
Sets all whitespace in file to a single space\n
Usage: ABCutils.pl SetSpace <input_file> <outfile_name>\n
\n/) if !@ARGV;
    open(INPUT, '<', $ARGV[0]);
    open(OUT, '>', $ARGV[1]);

    while (<INPUT>) {
        chomp;
        my $line = $_;
        $line =~ s/\s+/ /g;
        print OUT "$line\n";
    }

    close INPUT;
    close OUT;
}

sub DistReject {
### TO DO: ###
### make the format more flexible (i.e. blank line compatible) ###

my $numsim = 5000;
my $numkeep = 1000;
my $percent = 10;
my $stat_start = undef;
my $sim = undef;
my $obs_sfs = undef;
my $outfile = undef;
my $multiobs = undef;

die (qq/
Description: Finds closest simulations to observed data in terms of Euclidean distance 
             using a specified tolerance\n
Usage: ABCutils.pl DistReject [options]\n
Options:
    --numsim     INT     number of simulations to examine [$numsim]
    --numkeep    INT     number of closest simulations to retain [$numkeep]
    --percent    FLOAT   percentage of closest simulations to retain [$percent]
    --stat_start INT     starting column containing summary statistics in simulation file
                         (the first column is 1)
    --sim        FILE    file containing simulations
    --obs_sfs    FILE    file containing observed data
    --outfile    CHAR    outfile name including path
    --multiobs           compare nth simulation to nth observation (requires multiple observations)\n
Note: options --numkeep and --percent indirectly set the tolerance in an ABC rejection framework.
      When both options are specified, whichever specifies the fewest simulations to retain is used.           
\n/) if !@ARGV;

# check input

my $options = '--numsim --numkeep --percent --stat_start --sim --obs_sfs --outfile --multiobs ';
foreach (@ARGV) {
    if ($_ =~ /^\-\-/) {
        die("ERROR: Unknown command $_\n") if ($options !~ /$_\s+/);
    }
}

# get options and check input

GetOptions('numsim=i{1}' => \$numsim,
'numkeep=i{1}' => \$numkeep,
'percent=f{1}' => \$percent,
'stat_start=i{1}' => \$stat_start,
'sim=s{1}' => \$sim,
'obs_sfs=s{1}' => \$obs_sfs,
'outfile=s{1}' => \$outfile,
'multiobs' => \$multiobs);

die ("ERROR: --numsim must be a positive integer\n-->exiting\n") if ($numsim <= 0);
die ("ERORR: --numkeep must be a positive integer\n-->exiting\n") if ($numkeep <= 0);
die ("ERROR: --percent must be a number between 0 and 100\n-->exiting\n") if ($percent <= 0 || $percent > 100);
die ("ERROR: --stat_start must be a positive integer") if ($stat_start <= 0);
open(SIMULATED, '<', $sim) or die("ERROR: Can't open file of simulations: $!\n-->exiting\n");
open(OBSERVED, '<', $obs_sfs) or die("ERROR: Can't open observed SFS file: $!\n-->exiting\n");
open(OUT, '>', $outfile) or die("ERROR: Can't open outfile: $!\n-->exiting"); 

# determine which summary statistics in observed to use

my $obs = <OBSERVED>;
my @obs = split(/\s+/, $obs);
my @obs2 = @obs;

my @sim = split(/\s+/, <SIMULATED>);
my $simlen = @sim;
$simlen--;

my @statcol;
my %bounds;

my $shiftnum = 0;
my $n = 0;
my ($obscat, $simcat);
for(my $i = $stat_start-1; $i <= $simlen; $i++) {
    foreach (@obs) {
    	$obscat = $1 if $obs[$shiftnum] =~ /(\S+)/i;
    	$simcat = $1 if $sim[$i] =~ /(\S+)/i; 
        if ($obscat eq $simcat) 
        {
            push @statcol, $n;
            $bounds{$i} = [0, 0, $obs2[$n], $n];
            shift @obs;
            $n++;
            last;
        } else {
            $shiftnum++;
            $n++;
        }
           
    }
}

@obs = split(/\s+/, <OBSERVED>);
undef @obs2;
close OBSERVED if (!defined $multiobs);

# get euclidean distance between observed and sim, and statistic range

print STDERR "Finding Euclidean distances ...\n";
my %distance;
my $simnum = 1;
while (<SIMULATED>) 
{
	if (!@obs)
	{
		print STDERR "Number of simulations and observations differ: --multiobs invalid\n-->exiting\n";
		close OBSERVED; close SIMULATED;
		return -1;
	}
    next if $_ =~ /^\s*\n/;
    chomp;
    @sim = split(/\s+/, $_);
    my $sum = 0;
    my $k = 0;
    for(my $i = $stat_start-1; $i <= $simlen; $i++) {
        $sum += ($obs[$statcol[$k]] - $sim[$i])**2;
        # get range for statistic
        ${$bounds{$i}}[0] = $sim[$i] if ( ${$bounds{$i}}[0] < $sim[$i] || $. == 2 ); # max
        ${$bounds{$i}}[1] = $sim[$i] if ( ${$bounds{$i}}[1] > $sim[$i] || $. == 2); # min
        $k++;   
    }
    $distance{$simnum} = sqrt($sum);
    $simnum++;
    if (defined $multiobs)
    {
    	do { @obs = split(/\s+/, <OBSERVED>) } until (@obs || eof(OBSERVED));
    }
}
close OBSERVED if (defined $multiobs);

# find stats out of range of simulatons -- only if one set of observations is used
if (!defined $multiobs)
{
	foreach (sort { $a <=> $b } keys %bounds) 
	{
    	if ( $obs[${$bounds{$_}}[3]] > ${$bounds{$_}}[0] || $obs[${$bounds{$_}}[3]] < ${$bounds{$_}}[1] ) 
    	{
        	print STDERR "WARNING: \"${$bounds{$_}}[2]\" in observed data is outside simulated range!\n";
    	}
	}
}

# find closest simulations

$simnum--;
my $frac = sprintf("%.0f", ($simnum*($percent/100)));
my $bestnum = $numkeep <= $frac ? $numkeep : $frac;
$bestnum = $simnum + 1 if ($bestnum > $simnum + 1); # make sure number sims to retain does not exceed total number sims
print STDERR "Finding $bestnum closest simulations to observed data ...\n";
$bestnum--;

my %smalldist;
my @best = (sort { $a <=> $b} values(%distance))[0 .. $bestnum];
my @bestsims;
foreach (@best) {
    my $dist = $_;
    my @keys = grep { $distance{$_} == $dist } keys %distance;
    foreach (@keys) {
        if (exists $distance{$_}) {
            push @bestsims, $_;
            $smalldist{$_} = $distance{$_};
            delete $distance{$_};
        }
    }
}
undef %distance; # free memory

# adjust number of retained simulations if necessary

my $bestlen = @bestsims;
if ($bestlen > $bestnum + 1) {
    do {
        pop @bestsims;
        $bestlen = @bestsims;
    } until ($bestlen == $bestnum + 1);
}

undef %bounds; # free memory

# output simulations and distances

print STDERR "Writing file $outfile ...\n";

@bestsims = sort { $a <=> $b } @bestsims;

$" = " ";

seek(SIMULATED, 0, 0);

@sim = split(/\s+/, <SIMULATED>);
@sim = ($sim[0], "Dist", @sim[1 .. $#sim]);
print OUT "@sim\n";

$simnum = 1;
my $p = 0;
my $targetline = $bestsims[$p];
while ($simnum <= $bestsims[$#bestsims]) {
    my $simline = <SIMULATED>;
    next if $simline =~ /^\s*\n/;
    if ($simnum == $targetline) {
        chomp($simline);
        @sim = split(/\s+/, $simline);
        @sim = ($sim[0], $smalldist{$targetline}, @sim[1 .. $#sim]);
        print OUT "@sim\n";
        $p++;
        $targetline = $bestsims[$p];
    }
    $simnum++;
}
close SIMULATED;
close OUT;

print STDERR "Finished!\n";
}

sub format2Dsim {

die (qq/
Description: Manipulates the set of simulated 2D-SFS outputted by fastsimcoal_sampler.pl
	     format2Dsim automatically folds the joint SFS.\n 
Usage: ABCutils.pl format2Dsim [options]\n
Options:\n
--pop1n     INT    number of diploid individuals in first (x-axis) population in joint SFS
--pop2n	    INT    number of diploid individuals in second (y-axis) population in joint SFS
--sfsfile   FILE   .samp format file of 2D-SFS simulations outputted by fastsimcoal_sampler.pl
--outfile   CHAR   outfile name (including path). If not specified the output will be dumped to the same
                   directory as the .samp input file
--norm      0|1    normalize the 2D-SFS if 1 [0]
--maskFixed 0|1    mask fixed categories in 2D-SFS if 1 [0]
--bin	    INT	   calculate diagonal bins of 2D-SFS categories; INT is 2 x the number of categories on
		   either side of each bin's diagonal (~ bin width, which must be a factor of 2)
--do1dsfs   INT    calculate folded 1D-SFS for pop1 (INT = 0), pop2 (INT = 1), or both populations (INT = 2) [null]\n
Note:
--maskFixed 1 is automatically induced if fastsimcaol_sampler was run with --rmvfixed 1
\n/) if !@ARGV;

	### check input ###

	my $options = '--pop1n --pop2n --sfsfile --outfile --norm --maskFixed --bin --do1dsfs ';
	foreach (@ARGV) 
	{
    	if ($_ =~ /^\-\-/) 
    	{
        	die("ERROR: Unknown command $_\n") if ($options !~ /$_\s+/);
    	}
	}

	### get options ###

	my $pop1n;
	my $pop2n;
	my $sfsfile;
	my $outfile;
	my $norm = 0;
	my $maskFixed = 0;
	my $bin = undef;
	my $do1dsfs = undef;

	GetOptions(
		'pop1n=i' => \$pop1n, 
		'pop2n=i' => \$pop2n, 
		'sfsfile=s' => \$sfsfile, 
		'outfile=s' => \$outfile, 
		'norm=i' => \$norm, 
		'maskFixed=i' => \$maskFixed,
		'bin=i' => \$bin,
		'do1dsfs=i' => \$do1dsfs
	);

	die("ERROR: --norm must be either 0 or 1\n") unless ($norm == 0 || $norm ==1);
	die("ERROR: --maskFixed must be either 0 or 1\n") unless ($maskFixed == 0 || $maskFixed == 1);
	die("ERROR: --pop1n and --pop2n can only take integer values\n") if ($pop1n =~ /\D+/ || $pop2n =~ /\D+/);
	if (defined $bin)
	{
		die("ERROR: If binning, --bin must be zero or a positive factor of 2\n") if ($bin > 0 && $bin % 2 != 0 || $bin < 0);
	}
	if (defined $do1dsfs)
	{
		die("ERROR: --do1dsfs takes only arguments of 0, 1, or 2 for calculating the 1D-SFS\n") if ($do1dsfs < 0 || $do1dsfs > 2 || $do1dsfs =~ /[^0-2]/);
	}

	### open file of simulations ###
	open(SFS, '<', $sfsfile) or die("ERROR! Couldn't open file SFS: $!");
	my $out_default = $sfsfile . 'fold';
	if ($outfile) 
	{
    	open(OUT, '>', $outfile) or die("ERROR! Couldn't open file OUT: $!");
	} 
	else 
	{
    	open(OUT, '>', $out_default) or die("ERROR! Couldn't open file OUT: $!");
	}

	### determine if fixed categories are removed from sims and get param fields

	my $fix_keep = 0;
	my $params_end = 0;
	foreach ( split(/\s+/, <SFS>) ) 
	{
    	if ($_ =~ /E\d+/) 
    	{
        	$fix_keep++ if ($_ =~ /E00/);
        	last;
    	} 
    	else 
    	{
        	$params_end++;
        	print OUT "$_ "; # print the parameter names in header
    	}
	}

	### print header ###
	my $numcats;
	if (defined $do1dsfs)
	{
		my $i = 0;
		if ($do1dsfs == 0)
		{
			for ($i = 0; $i <= $pop1n; $i++)
			{
				print OUT "FX$i ";
			}
		}
		elsif ($do1dsfs == 1)
		{
			for ($i = 0; $i <= $pop2n; $i++)
			{
				print OUT "FY$i ";
			}
		}
		elsif ($do1dsfs == 2)
		{
			for ($i = 0; $i <= $pop1n; $i++)
			{
				print OUT "FX$i ";
			}
			for ($i = 0; $i <= $pop2n; $i++)
			{
				print OUT "FY$i ";
			}
		}
		else
		{
			die("ERROR: '$do1dsfs' is an unknown argument for calculating 1D-SFS\n");
		}
	}
	if (!defined $bin)
	{
		$numcats = (($pop1n + 1) * ($pop2n + 1)) - 1; # number folded categories
		for (my $i = 0; $i <= $numcats; $i++) 
		{
    		print OUT $i == $numcats ? "F$i\n" : "F$i ";
		}
	}
	else # number of bins
	{
		$numcats = roundup(($pop1n + 1 - ($bin+1))/(2*$bin+1)) + 1 + roundup(($pop2n + 1 - ($bin + 1))/(2*$bin+1));
		my $i = 1;
		for ($i = 1; $i <= $numcats; ++$i) 
		{
			if (!defined $do1dsfs)
			{
    			print OUT "OB$i ";
			}
			else
			{
				print OUT $i < $numcats ? "OB$i " : "OB$i\n"; 
			}
		}
		if (!defined $do1dsfs)
		{
			for ($i = 1; $i <= $numcats; ++$i)
			{
				print OUT $i < $numcats ? "DB$i " : "DB$i\n";
			}
		}
	}

	### process simulations ###
	my ($pop1chr, $pop2chr) = ($pop1n * 2, $pop2n * 2); # get number of chromosomes in pops
	my @param_vals;
	while (<SFS>)
	{   
    	chomp;  
    	my @sfsvec = split(/\s+/, $_);
    	@param_vals = splice @sfsvec, 0, $params_end;
    	print OUT "@param_vals "; 
    	my @sfsfold;

    	## insert fixed categories into 2dsfs vector if they were already removed - needs testing !!!
    	if ($fix_keep == 0)
    	{
        	foreach (0, $pop1chr, $#sfsvec + 3, -$pop1chr) 
        	{
            	splice @sfsvec, $_, 0, 0;   
        	} 
    	}
    	
    	## fold
		for (my $i = 0; $i <= $#sfsvec; $i += $pop1chr + 1) 
    	{
        	my @row = @sfsvec[$i .. $i + $pop1chr];
        	@row = ((map { $row[$_] + $row[$#row - $_] } ( 0 .. (($#row/2)-1))), $row[$#row/2]);
        	push @sfsfold, \@row;
    	}
		my $lastindx= $#sfsfold;
		map { do{ my $j = $_; map { $sfsfold[$j][$_] += $sfsfold[$lastindx - $j][$_] } 0 .. $#{$sfsfold[0]}; pop @sfsfold; }} 0 .. ($lastindx/2)-1;
		undef @sfsvec;

    	## mask fixed categories (optional)
   		${$sfsfold[0]}[0] = 0 if ($maskFixed == 1 && $fix_keep);

		## normalize 2D-SFS (optional) - untested !!!
   		if ($norm == 1)
    	{
        	my $total = 0;
        	map { do { my $i = $_; map {$total += ${$sfsfold[$i]}[$_] } 0 .. $#{$sfsfold[$i]}  };  } 0 .. $#sfsfold;
        	map { do { my $i = $_; map {${$sfsfold[$i]}[$_] /= $total} 0 .. $#{$sfsfold[$i]} }; } 0 .. $#sfsfold;
    	}
    
    	## calculate 1D-SFS ##
		my @sfsout = ();
    	if (defined $do1dsfs)
    	{
    		if ($do1dsfs == 0)
    		{
    			calc1DSFS(\@sfsfold, \@sfsout, 0);
    		}
    		elsif ($do1dsfs == 1)
    		{
    			calc1DSFS(\@sfsfold, \@sfsout, 1);
    		}
    		elsif ($do1dsfs == 2)
    		{
    			calc1DSFS(\@sfsfold, \@sfsout, 0);
    			calc1DSFS(\@sfsfold, \@sfsout, 1);
    		}
    		else
    		{
    			die("ERROR: '$do1dsfs' is an unknown argument for calculating 1D-SFS\n");
    		}
    	}
    
    	## bin 2D-SFS ##
    	if (defined $bin)
    	{
    		binJointSFS(\@sfsfold, \@sfsout, $bin);
    		if (!defined $do1dsfs)
    		{
    			my @rev = reverse @sfsfold;
    			binJointSFS(\@rev, \@sfsout, $bin);
    		}
    	}
    	else
    	{
   			map { push(@sfsout, @{$sfsfold[$_]}) } 0 .. $#sfsfold;
    	}

		## output SFS ##
		print OUT "@sfsout\n";
	}
	close SFS;
	close OUT;	
}
	

sub binJointSFS # calculates diagonal bins from 2DSFS 
{
	# width = 2 X number of categories on each side of diagonal (~ bin width)
	# bins are ordered from upper left to lower right bin in returned array reference
	my ($jointsfs, $binarr, $width) = ($_[0], $_[1], int($_[2]));
	if ($width > 0 && $width % 2 != 0 || $width < 0)
	{
		print STDERR ("ERROR: Bin width must be zero or a positive factor of 2\n");
		return(-1);
	}
	my $p1n = @{$$jointsfs[0]}; # number of categories for population 1
	my $p2n = @$jointsfs; # number of categories for population 2

	# bin the joints SFS
	my $numxbin = roundup(($p1n - ($width+1))/(2*$width+1)) + 1;
	my $numybin = roundup(($p2n - ($width + 1))/(2*$width+1)) + 1;
	my $binsum = 0;
	my @bins;
	
	 #print STDERR "p1 N: $p1n\n"; # debug
	 #print STDERR "p2 N: $p2n\n"; # debug
	 #print STDERR "num x bins: $numxbin\n"; # debug
	 #print STDERR "num y bins: $numybin\n"; # debug

	# bin along the y-axis (pop2)
	my ($i, $x, $k, $y) = (0, 0, 0, 0);
	@bins = map { 0 } 1..($numybin-1);
	for ($i = 1; $i < $numybin; ++$i)
	{
		$y = $i * 2 * $width + $i;
		$binsum = 0;
		for ($x = 0; ($x - $width) < $p1n; ++$x)
		{
			 #print "X: $x\n"; # debug
			 #my @index = ($y - $width .. $y + $width); # debug
             #print STDERR "@index\n"; # debug
			for ($k = $y-$width; $k <= $y+$width; ++$k)
			{
				$binsum += ${$$jointsfs[$k]}[$x] if ($k < $p2n && $x < $p1n && $k >= 0);
			}
			++$y;
		}
		$bins[$numybin-(1+$i)] = $binsum;
	}
	
	# bin along the x-axis (pop1)
	for ($i = 0; $i < $numxbin; ++$i)
	{
		$y = 0;
		$binsum = 0;
		for ($x = $i * 2 * $width + $i; ($x - $width) < $p1n; ++$x)
		{
			#print "X: $x\n"; # debug
			#my @index = ($x - $width .. $x + $width); # debug
			#print STDOUT "@index\n"; # debug
			for ($k = $x-$width; $k <= $x+$width; ++$k)
			{
				$binsum += ${$$jointsfs[$y]}[$k] if ($y < $p2n && $k < $p1n && $k >= 0);
			}
			++$y;
		}
		push @bins, $binsum;
	}
	push @$binarr, @bins;
	return (scalar(@bins));
}

sub roundup # works as a ceiling function
{
	my $n = shift;
	return( ($n == int($n)) ? $n : int($n + 1) );
}

sub calc1DSFS # extracts 1D-SFS from 2D-SFS
{
	my ($jointsfs, $sfsout, $popID) = @_;
	my $p1n = $#{$$jointsfs[0]}; # pop1 n
	my $p2n = $#$jointsfs; # pop2 n
	my @sfs1d;
	
	# extract 1D-SFS from joint SFS
	my $i = 0;
	if ($popID == 0) # calculate 1D-SFS for population 1
	{
		map { do { $i = $_; $sfs1d[$i] = 0; map {$sfs1d[$i] += $$jointsfs[$_]->[$i]} 0 .. $p2n;  }; } 0 .. $p1n;
	} 
	elsif ($popID == 1) # calculate 1D-SFS for population 2
	{
		map { do { $i = $_; $sfs1d[$i] = 0; map {$sfs1d[$i] += $$jointsfs[$i]->[$_]} 0 .. $p1n; }; } 0 .. $p2n;
	}
	else
	{
		print STDERR "Invalid population ID for calculating 1D-SFS in calc1DSFS function\n";
		return -1;
	}
	push @$sfsout, @sfs1d;
	return 0;
}

sub gofReject
{
die (qq/
Description:\n
Calculates goodness of fit statistic as (observed - expected)^2\/expected between simulated and observed
summary statistics and retains simulations closest to the observed data\n
Usage:\n
ABCutils.pl gofReject [options]\n
Options:\n
--obsfile	FILE	Observed summary stats (e.g. SFS bins as outputted by format2Dobs)
--simfile	FILE	Simulated summary stats (e.g. SFS bins as outputted by format2Dsim)
--outfile	CHAR	Output file name where results will be dumped
--stat_start	INT	Field number of the first summary statistic in simfile (1-based) to use in fitting.
			All summary statistics in fields >= INT will be used.
--num_keep	INT	Find the closest INT simulations to the observed joint spectrum [0]
--percent_keep	FLOAT	Find the closest FLOAT percent of simulations to the observed joint spectra [10] 
--numsim	INT	Number of simulations (from beginning of simfile) to calculate goodness of fit for\n
Notes:\n
The smaller number of simulations to find between --num_keep and --percent_keep will be used 
\n/) if (!@ARGV);
	
	### check input ###
	my $options = '--obsfile --simfile --outfile --stat_start --num_keep --percent_keep --numsim ';
    	foreach (@ARGV) 
	{
        	if ($_ =~ /^\-\-/) 
		{
        		die("ERROR: Unknown command $_\n") if ($options !~ /$_\s+/);
        	}
    	}

	### get options ###
	my $obsfile = undef;
	my $simfile = undef;
	my $outfile = undef;
	my $stat_start = undef;
	my $num_keep = 0;
	my $percent_keep = 10;
	my $numsim = 500;

	GetOptions (
		'obsfile=s{1}' => \$obsfile,
		'simfile=s{1}' => \$simfile,
		'outfile=s{1}' => \$outfile,
		'stat_start=i{1}' => \$stat_start,
		'num_keep=i{1}' => \$num_keep,
		'percent_keep=f{1}' => \$percent_keep,
		'numsim=i{1}' => \$numsim
	);

	if ($num_keep =~ /\./)
	{
		print STDERR "WARNING: --num_keep is not an integer and so will be rounded down\n";
		$num_keep = int($num_keep);
	}
 	if ($numsim =~ /\./)
        {
                print STDERR "WARNING: --numsim is not an integer and so will be rounded down\n";
                $numsim = int($numsim);
        }
	die("ERROR: --percent_keep must be in the interval (0,100] --> exiting\n") if ($percent_keep <= 0 || $percent_keep > 100);
	die("ERROR: --num_keep must be a positive integer < number of simulations --> exiting\n") if ($num_keep <= 0);
	die("ERROR: --numsim must be a positive integer --> exiting\n") if ($numsim <= 0);
	die("ERROR: --stat_start must be a positive integer --> exiting\n") if ($stat_start <= 0 || $stat_start =~ /\.\d*[^0]+/);
	open(OBS, '<', $obsfile) or die("ERROR: Can't open observed file $obsfile --> exiting\n");
	open(SIM, '<', $simfile) or die("ERROR: Can't open file of simulations $simfile --> exiting\n");
	open(OUT, '>', $outfile) or die("ERROR: Can't open outfile $outfile --> exiting\n");

	# determine how many simulations to keep
	my $percent_num = int(($percent_keep/100) * $numsim);
	my $keep = $percent_num < $num_keep ? $percent_num : $num_keep;
	if ($keep < 1)
	{
		unlink $outfile;
		die ("ERROR: Finding 0 closest simulations. Increase --percent_keep or --numsim\n--> exiting\n");
	}
	# get observed stat vector	
	<OBS>;
	my @obstat = split(/\s+/, <OBS>);
	die ("ERROR: Problem storing observed stats...check observed file\n") unless @obstat;
	close OBS;	

	# go through simulations and calculate goodness of fit
	my $simpos;
	my @header = split(/\s+/, <SIM>);
	my @simstat;
	my $lastidx = $#header;
	my $gof;
	my %fit;
	my $retain = 0;
	use bytes;
	while (<SIM>)
	{
		$simpos = tell(SIM) - length($_);;
		chomp;
		@simstat = (split(/\s+/, $_))[$stat_start - 1 .. $lastidx];
		$gof = pearsonChiSq(\@obstat, \@simstat);
		if ($gof == -1)
		{
			die ("ERROR: Number of statistics in observed and simulated differ --> exiting\n");
		}
		elsif ($gof == -2) # observed too different than expected
		{
			next;
		}
		else
		{
			push @{$fit{$gof}}, $simpos;
			$retain++;
		}
	}

	# find the smallest gof statistics
	my $subpercent;
	my @best;
	my $fitkey;
	my $counter = 0;
	if ($retain > $keep)
	{
     		$subpercent = ($keep / $retain) * 100;
                print STDERR "Best $keep simulations represent $subpercent% of retained simulations\n";

		foreach (sort {$a <=> $b} keys  %fit)
		{
			$fitkey = $_;
			foreach ( @{$fit{$fitkey}} )
			{
				push @best, [$_, $fitkey];
				$counter++;
				last if $counter == $keep;
			}
			last if $counter == $keep;
		}
	}
	else
	{
		die("Number of best fitting sims to keep exceeds overall number of retained simulations --> exiting\n");
	}
	
	# print best fitting lines
	print OUT "$header[0] ", "GOF " , "@header[1 .. $#header]\n";
	my @goodline;
	foreach ( sort {$a->[0] <=> $b->[0]} @best )
	{
		seek SIM, $_->[0], 0;
		@goodline = split(/\s+/, <SIM>);
		print OUT "$goodline[0] ", "$_->[1] ", "@goodline[1 .. $#goodline]\n";
	}
	
	close SIM;
	close OUT; 
	
}

sub pearsonChiSq
{
	my ($obs, $exp) = @_;
	return -1 if (scalar(@$obs) != scalar(@$exp));
	my $stat = 0;
	my $idx = 0;
	foreach (@$exp)
	{
		if ($$obs[$idx] != 0 && $_ == 0)
		{
			# simulation should be discarded
			return (-2);
		}
		elsif ($$obs[$idx] == 0 && $_ == 0)
		{
			# goes to zero in the limit when exp and obs approach zero
			# numerator goes to zero faster than denominator
			$stat += 0;
		}
		else
		{
			$stat += ($$obs[$idx] - $_)**2 / $_;
		}
		$idx++;
	} 
	return $stat;
}

sub flexibleHeader
{
	# header codes for $type INT
	# 0 = antidiagonal bins
	# 1 = antidiagonal bins + diagonal bins
	# 2 = antidiagonal bins + pop1 1D-SFS
	# 3 = antidiagonal bins + pop2 1D-SFS
	# 4 = antidiagonal bins + pop1 1D-SFS + pop2 1D-SFS 
	# 5 = categories in folded 2D-SFS
	
	my ($pop1n, $pop2n, $type, $outfile, $bin) = @_;
	my $numcats;
	my $k = 0;
	
	# determine number of 2D-SFS categories or bins
	if ($type == 5)
	{
		$numcats = (($pop1n + 1) * ($pop2n + 1)) - 1; # number folded 2D-SFS categories	
	}
	else
	{
		if (defined $bin)
		{
			$numcats = ( int(($pop1n - $bin + 2) / (2 * $bin + 1)) 
			+ roundup(($pop1n + 1 - (int(($pop1n - $bin + 2) / (2* $bin + 1))) * (2 * $bin + 1)) / ($bin + 1)) 
			+ int(($pop2n - $bin + 2) / (2 * $bin + 1)) 
			+ roundup(($pop2n + 1 - (int(($pop2n - $bin + 2) / (2* $bin + 1))) * (2 * $bin + 1)) / ($bin + 1)) ) - 2;
		}
	}
	
	# print headers
	if ($type == 2)
	{
		map {print $outfile "0F$_ "} 0 .. $pop1n;
	}
	elsif ($type == 3)
	{
		map {print $outfile "1F$_ "} 0 .. $pop2n;
	}
	elsif ($type == 4)
	{
		map {print $outfile "0F$_ "} 0 .. $pop1n;
		map {print $outfile "1F$_ "} 0 .. $pop2n;
	} 
	elsif ($type == 5)
	{
		map {print $outfile $_ == $numcats ? "F$_\n" : "F$_ "} 0 .. $numcats;
		return 0;	
	}
	elsif ($type == 1)
	{
		map {print $outfile "B$_ "} 0 .. $numcats;
	}
	else
	{
		return 1;
		print STDERR "Invalid header type in call to flexibleHeader\n";
	}
	map {print $outfile $_ == $numcats ? "AB$_\n" : "AB$_ "} 0 .. $numcats;
	return 0;
}

sub best2Dsfs {
	# Doesn't work for removed fixed categories -- need to implement
	# parse arguments
	my ($pop1n, $pop2n, $simfile, $distfile, $outfile, $fit, $dofold) = (undef, undef, undef, undef, undef, 2, undef);
	if (@ARGV < 10 || $ARGV[0] =~ /help/i)
	{
		best2DsfsMessage($fit);
		return 0;
	}
	my $args = '--pop1n --pop2n --simfile --distfile --outfile --fit --dofold ';
	foreach (@ARGV)
	{
		if ($_ =~ /^\-\-[^\-]+/)
		{
			if ($args !~ /$_/)
			{
				print STDERR "Unknown argument: $_\n-->exiting\n";
				return 0;
			}
		}
	}
	GetOptions(
		'pop1n=i{1}'=>\$pop1n, 
		'pop2n=i{1}'=>\$pop2n, 
		'simfile=s{1}'=>\$simfile, 
		'distfile=s{1}'=>\$distfile, 
		'outfile=s{1}'=>\$outfile, 
		'fit=i{1}'=>\$fit,
		'dofold'=>\$dofold);
	if ($fit-1 < 0)
	{
		print STDERR "--fit must be and integer >= 1\n-->exiting\n";
		return 0;
	}
	--$fit;
	if ($pop1n <= 0)
	{
		print STDERR "--pop1n needs to be a positive integer\n-->exiting\n";
		return 0;
	}
	if ($pop2n <= 0)
	{
		print STDERR "--pop2n needs to be a positive integer\n-->exiting\n";
		return 0;
	}
	open(my $simfh, '<', $simfile) or die("Couldn't open file of simulations $simfile: $!\n-->exiting\n");
	open(my $distfh, '<', $distfile) or die("Couldn't open file containing distances $distfile: $!\n-->exiting\n");
	# get information from header
	my @header;
	do { @header = split(/\s+/, <$simfh>) } until ($header[0] || eof($simfile));
	my $sfsfield = 0;
	foreach (@header)
	{
		if ($_ =~ /(^F\d+$|^E\d+$)/)
		{
			if ($1 =~ /E/)
			{
				print STDERR "\nUnfolded spectrum detected\n\n";
				$pop1n = 2*$pop1n + 1;
				$pop2n = 2*$pop2n + 1;
			}
			else
			{
				print STDERR "\nFolded spectrum detected\n\n";
				undef $dofold if $dofold;
				++$pop1n;
				++$pop2n;
			}
			last;
		}
		++$sfsfield;
	}
	
	if ($sfsfield == $#header)
	{
		print STDERR "Couldn't find SFS start field: Check simfile header.\n-->exiting\n";
		return 0;
	}
	# check that population sizes match joint spectra dimensions
	my $specsize = $#header - $sfsfield + 1;
	if ($specsize/$pop1n != $pop2n)
	{
		print STDERR "Population sizes provided do not match simfile SFS dimensions\n-->exiting\n";
		return 0;
	}
	# open output file
	open(my $outfh, '>', $outfile) or die("Couldn't open output file $outfile: $!\n-->exiting\n");
	# find simulation with minimum distance from observed
	my $bestdist = undef;
	my $simnum = MinDist($fit, $distfh, \$bestdist);
	return -1 if ($simnum < 0);
	print STDERR "\nBest simulation:\n\n"; print STDOUT $bestdist; print STDERR "\n";
	close $distfh;
	# find matching joint SFS simulation
	my @jointsfs;
	ExtractSFS($simnum, $simfh, $sfsfield, $pop1n, \@jointsfs);
	close $simfh;
	# print output
	PrintBestSFS($outfh, \@jointsfs, $dofold);
	close $outfh;
	return 0;
}

sub FoldJoint {
	# takes unfolded 2D-SFS matrix as input - only works if sample size is an even number for both pops
	my ($unfolded, $sfsfold) = @_;
	foreach my $row (@$unfolded) 
    {
    	@$row = ((map { $$row[$_] + $$row[$#$row - $_] } ( 0 .. (($#$row/2)-1))), $$row[$#$row/2]);
    	push @{$sfsfold}, $row;
    }
		my $lastindx= $#$sfsfold;
		map { do{ my $j = $_; map { $$sfsfold[$j][$_] += $$sfsfold[$lastindx - $j][$_] } 0 .. $#{$$sfsfold[0]}; pop @$sfsfold; }} 0 .. ($lastindx/2)-1;
}

sub PrintBestSFS {
	my ($ofh, $sfs, $fold) = @_;
	my $outsfs = [];
	if (defined $fold)
	{
		FoldJoint($sfs, $outsfs);
	}
	else
	{
		$outsfs = $sfs;
	}
	foreach (@$outsfs)
	{
		print $ofh "@$_\n";
	}
}

sub ExtractSFS {
	my ($sim, $fh, $start, $xdim, $sfs) = @_;
	my @line;
	while (<$fh>)
	{
		if ($_ =~ /^(\d+)\s/)
		{
			if ($sim == $1)
			{
				@line = split(/\s+/, $_);
				@line = @line[$start..$#line];
			}
		}
		else
		{
			next;
		}
	}
	my $i = 0;
	while (@line)
	{
		for (my $j = 0; $j < $xdim; ++$j)
		{
			$$sfs[$i][$j] = shift @line;
		}
		++$i;
	}
	return 0;
}

sub MinDist {
	my ($distcol, $fh, $bestline) = @_;
	my @line;
	do { my $l = <$fh>; @line = split(/\s+/, $l); $$bestline = $l } until ($line[0] =~ /^\d+/ || eof($fh));
	if (eof($fh))
	{
		print STDERR "No data found in distfile\n";
		return -1;
	}
	my ($sim, $min) = ($line[0], $line[$distcol]);
	while (<$fh>)
	{
		@line = split(/\s+/, $_);
		if ($line[$distcol] < $min)
		{
			$sim = $line[0];
			$min = $line[$distcol];
			$$bestline = $_;
		}
	}
	return $sim;	
}

sub best2DsfsMessage {
my $fit = shift;
die(qq/
best2Dsfs extracts the simulated joint spectrum that most closely fits the observed joint spectrum.\n
Usage: ABCutils.pl best2Dsfs <sim_file> [arguments]\n
Arguments:
simfile    FILE  File of simulations containing joint SFS
distfile   FILE  File containing distances of observed from expected joint spectra
outfile    CHAR  Name of output file
pop1n      INT   Number of diploid individuals in population 1
pop2n      INT   Number of diploid individuals in population 2
fit        INT   Field number (1-based) of column specifying distance of observed from expected [$fit]
dofold           If specified, folds the joints spectrum\n
Notes:
--simfile should be output from fastsimcoal_sampler.pl or ABCutils 'format2Dsim' 
--distfile should be output from ABCutils 'DistReject' or 'gofReject'\n
\n/)
}

sub fold2D {
	if (!@ARGV || $ARGV[0] =~ /help/i)
	{
		fold2DMessage();
		return 0;
	}
	open(my $infh, '<', $ARGV[0]) or die("Couldn't open input joint SFS $ARGV[0]: $!\n-->exiting\n");
	my @sfsvec;
	while (<$infh>) { push @sfsvec, [split(/\s+/, $_)]; }
	close $infh;
	my @folded;
	FoldJoint(\@sfsvec, \@folded);
	foreach (@folded) { print STDOUT "@$_\n"; }
}

sub fold2DMessage {
die (qq/
fold2D folds a 2-dimensional site frequency spectrum\n
Usage: ABCutils.pl fold2D <unfolded_SFS_file>\n
Note: The input 2D-SFS should not contain any headers.
\n/)
}

sub StatDistr {
	# parse arguments
	my $error = 0;
	my ($parfile, $outfile, $numsim, $arlexec, $infsites, $metasamp, $sitefst) = (undef, undef, 100, undef, undef, undef, undef);
	my (@pop1id, @pop2id);
	if (@ARGV < 6 || $ARGV[0] =~ /help/i)
	{
		StatDistrMessage($numsim);
		return 0;
	}
	my $args = '--parfile --outfile --numsim --arlexec --infsites --metasamp --pop1id --pop2id --sitefst ';
	foreach (@ARGV)
	{
		if ($_ =~ /^\-\-[^\-]+/)
		{
			if ($args !~ /$_/)
			{
				print STDERR "Unknown argument: $_\n-->exiting\n";
				return 0;
			}
		}
	}
	
	GetOptions(
	'parfile=s{1}'=>\$parfile,
	'outfile=s{1}'=>\$outfile,
	'numsim=i{1}'=>\$numsim,
	'arlexec=s{1}'=>\$arlexec,
	'infsites'=>\$infsites,
	'metasamp'=>\$metasamp,
	'pop1id=i{1,}'=>\@pop1id,
	'pop2id=i{1,}'=>\@pop2id,
	'sitefst'=>\$sitefst);
	
	$arlexec = abs_path($arlexec);
	if (!-x $arlexec)
	{
		print STDERR "Couldn't find arlsumstat exectuable $arlexec: $!\n-->exiting";
		return 0;
	}
	my $arldir = dirname($arlexec);
	if (!-d $arldir)
	{
		print STDERR "Couldn't find directory with arlsumstat files $arldir: $!\n-->exiting";
		return 0;
	}
	my $exname = basename($arlexec);
	$parfile = abs_path( $parfile );
	if (!-f $parfile)
	{
		print STDERR "Couldn't find parameter file $parfile ... Check if it exists\n-->exiting\n";
		return 0;
	}
	if (defined $metasamp)
	{
		foreach (@pop1id)
		{
			++$error if ($_ < 0);
		}
		foreach (@pop2id)
		{
			++$error if ($_ < 0);
		}
		if ($error)
		{
			print STDERR "All subpopulation IDs must be >= 0\n-->exiting\n";
			return 0;
		}
	}
	my $outfh;
	open($outfh, '>', $outfile) or die ("Couldn't open output file $outfile: $!\n-->exiting");
	close $outfh;
	my $siteout = $outfile . ".sitefst";
	open($$outfh, '>', $siteout) or die("Coulnd't open output file of Fst per site $siteout: $!\n-->exiting") if (defined $sitefst);
	$outfile = abs_path($outfile);
	# get into correct directory
	#my $dir = $1 if ($parfile =~ /(.+)\/.+$/);
	my $dir = dirname($parfile);
	chdir($dir);
	# run simulations
	$parfile = $1 if ($parfile =~ /.+\/(.+)$/); # this is because fastsimcoal does weird things about where it dumps output based on path
	$error = runfastsimcoal($parfile, $numsim, $infsites);
	if ($error)
	{
		print STDERR "fastsimcoal failed\n-->exiting;";
		return -1;
	}
	# process simulations
	my ($inarpfh, $outarpfh);
	my $simdir = $dir . '/' . basename($parfile);
	$simdir =~ s/\.\S*$//;
	#my $simdir = $dir . '/modelpar';
	my $targetdir = $simdir;
	my ($dirh, $file) = (undef, undef);
	if ($metasamp)
	{
		my $pooldir = $simdir . '/meta';
		if (mkdir($pooldir) == 0)
		{
			print STDERR "Couldn't create directory $pooldir to hold modified .arp files: $!\n-->exiting\n";
			StatDistrErrorExit($outfh, $simdir);
			return -1;
		}
		my @subpops = (\@pop1id, \@pop2id);
		my $metanum = 2;
		opendir($dirh, $simdir) or die("Couldn't open directory containing .arp files $simdir: $!\n-->exiting\n");
		while (readdir($dirh))
		{
			$file = $simdir . "/$_";
			if (-f $file && $_ =~ /\.arp$/)
			{
				open($inarpfh, '<', $file) or die ("Couldn't open arlequin project file $file: $!\n-->exiting");
				open($outarpfh, '>', "$pooldir/$_") or die ("Couldn't open alequin project file $pooldir/$_: $!\n-->exiting");
				poolSubPops($inarpfh, $outarpfh, \@subpops, $metanum);
				close $inarpfh;
				close $outarpfh;
			}
		}
		close $dirh;
		$targetdir = $pooldir;
	}
	# compute summary statistcs
	print STDERR "\nCalculating Summary Statistics with arlsumstat ...\n\n";
	runArlSumStat($arldir, $exname, $targetdir, $outfile);
	# compute per site Fst
	if (defined $sitefst)
	{
		print STDERR "Calculating per site Fst for all simulations ...\n\n";
		opendir ($dirh, $targetdir) or die("Couldn't open directory containing .arp files $simdir: $!\n-->exiting\n");
		while (readdir $dirh)
		{
			$file = $targetdir . "/$_";
			my @fstarr = ();
			if (-f $file && $_ =~ /\.arp$/)
			{
				open(my $arpfh, '<', $file);
				siteFst($arpfh, \@fstarr);
				print $outfh "@fstarr\n";
				close $arpfh;
			}
	 	
		}
	}
	close $outfh if (defined $sitefst);
	#delete temp files
	if (!rmtree $simdir)
	{
		print STDERR "Couldn't remove modelpar directory in 'StatDistr' subroutine\n-->exiting";
		return -1;
	}
	return 0;
}

sub siteFst {
	my ($arpfh, $fstvec) = @_;
	my ($popn, $maxn) = (0, 2);
	my %seq;
	# collect sequences
	while (<$arpfh>)
	{
		if ($popn < $maxn)
		{
			my $line = $_;
			if ($line =~ /SampleData= {/)
			{
				chomp(my $line = <$arpfh>);
				while ($line =~ /\S+\s+\d+\s+(.+)$/)
				{
					my $seqonly = $1;
					$seqonly =~ s/\s+//g;
					push @{$seq{$popn}}, [split(//, $seqonly)];
					chomp($line = <$arpfh>);
				}
				++$popn;
			}
		}
		else
		{
			last;
		}	
	}
	# calculate Reynold's distance
	ReynoldDist(\%seq, $fstvec, $maxn);
	return 0;
}

sub ReynoldDist
{
	my ($seqdat, $distvec, $maxn) = @_;
	my @hapn = (scalar @{$$seqdat{0}}, scalar @{$$seqdat{1}});
	my @dipn = ($hapn[0]/2, $hapn[1]/2);
	my $poolhapn = 0;
	foreach (@hapn)
	{
		$poolhapn += $_;
	}
	my $pooldipn = 0;
	foreach (@dipn)
	{
		$pooldipn += $_;
	}
	my $nsites =  scalar @{${$$seqdat{0}}[0]};
	my $i = 0;
	my %counts;
	my @alleles = ('A', 'C', 'G', 'T');
	for (my $site = 0; $site < $nsites; ++$site)
	{
		for($i = 0; $i < $maxn; ++$i)
		{
			foreach (@alleles)
			{
				$counts{$i}{$_} = 0;
			}
			foreach my $ind (@{$$seqdat{$i}})
			{
				++$counts{$i}{$$ind[$site]};
			}
			foreach(@alleles)
			{
				$counts{$i}{$_} /= $hapn[$i];
			}
		}
		# for 2 pops
		my $diffsum = 0;
		my $alpha1 = 0;
		my $alpha2 = 0;
		foreach (keys %{$counts{0}})
		{
			$diffsum += ($counts{0}{$_} - $counts{1}{$_})**2;
			$alpha1 = $counts{0}{$_}**2;
			$alpha2 = $counts{1}{$_}**2;
		}
		$diffsum *= 0.5;
		$alpha1 = 1 - $alpha1;
		$alpha2 = 1 - $alpha2;
		my $within = $diffsum - ($pooldipn*($dipn[0]*$alpha1 + $dipn[1]*$alpha2))/(4*$dipn[0]*$dipn[1]*($pooldipn-1));
		my $total = $diffsum + ((4*$dipn[0]*$dipn[1] - $pooldipn)*($dipn[0]*$alpha1 + $dipn[1]*$alpha2))/(4*$dipn[0]*$dipn[1]*($pooldipn-1));
		my $fst = $within/$total;
		push @$distvec, $fst < 0 ? 0 : $fst;
	}
	return 0;
}

sub StatDistrErrorExit {
	foreach (@_)
	{
		if (openhandle($_))
		{
			close $_ or die("Couldn't close filehandle $_: $!\n-->exiting\n");
		}
		if (-d $_)
		{
			rmtree($_) or die("Couldn't remove directory $_: $!\n-->exiting\n");
		}
	}
}

sub runArlSumStat {
	my ($arldir, $exec, $arpdir, $outfile) = @_;
	my $currentdir = getcwd();
	chdir($arldir);
	my $returnval = 0;
	my $header = 1;
	opendir(my $arpdirh, $arpdir) or die("Couldn't open directory with .arp files $arpdir in call to 'runArlSumStat': $!\n-->exiting");
	my $simnum = 0;
	while (readdir($arpdirh))
	{
		my $arpfile = $arpdir . "/$_";
		if (-f $arpfile && $arpfile =~ /\.arp$/)
		{
			$header = $simnum == 0 ? 1 : 0;
			my $returnval = system("./$exec $arpfile $outfile 1 $header");
			if ($returnval)
			{
				print STDERR "arlsumstat run failed\n";
				last;
			}
			++$simnum;
		}
	}
	chdir($currentdir);
	return $returnval;
}

sub runfastsimcoal {
	my ($parfile, $numsim, $inf) = @_;
	my $input = "-i $parfile -n $numsim -q";
	$input .= ' -I' if (defined $inf);
	my $returnval = system("fastsimcoal $input");
	return $returnval;
}

sub poolSubPops {
	my ($arpin, $arpout, $popids, $metan) = @_;
	my $i = 1;
	my %pops;
	for ($i = 1; $i <= $metan; ++$i)
	{
		map { $pops{$i}{$_} = 1 } @{$$popids[$i-1]};	
	}
	# print header info
	my $line = <$arpin>;
	while ($line !~ /SampleName=/)
	{
		$line =~ s/NbSamples=\d+/NbSamples=2/;
		print $arpout $line;
		$line = <$arpin>;
	}
	# get sequence information
	my %popdata;
	my $subpop;
	for ($i = 1; $i <= $metan; ++$i)
	{
		$popdata{$i}{size} = 0;
	}
	while ($line !~ /\[\[Structure\]\]/)
	{
		if ($line =~ /SampleName="Sample\s(\d+)"/)
		{
			$subpop = $1; 
			while ($line !~ /^\d+_\d+/)
			{
				$line = <$arpin>;
			}
			while ($line =~ /^\d+_\d+/)
			{
				for ($i = 1; $i <= $metan; ++$i)
				{	
					if (exists $pops{$i}{$subpop})
					{
						push @{$popdata{$i}{seq}}, $1 if ($line =~ /\S+\s+\d+\s+(.+)$/);
						++$popdata{$i}{size};
						last;
					}
				}
				$line = <$arpin>;
			}
			++$subpop;	
		}
		$line = <$arpin>;
	}
	# print sequence data
	for ($i = 1; $i <= $metan; ++$i)
	{
		print $arpout "\t\tSampleName=\"Sample $i\"\n\t\tSampleSize=$popdata{$i}{size}\n\t\tSampleData= {\n";
		my $j = 1;
		if (!@{$popdata{$i}{seq}})
		{
			print STDERR "No sequences for metapopulation $i\n";
			return -1;
		}
		foreach(@{$popdata{$i}{seq}})
		{
			print $arpout "${i}_$j\t1\t$_\n";
			++$j;
		}
		print $arpout "\n}\n";
	}
	# print structure info
	print $arpout "\n[[Structure]]\n\n\tStructureName=\"Simulated data\"\n\tNbGroups=1\n\tGroup={\n";
	for($i = 1; $i <= $metan; ++$i)
	{
		print $arpout "\t  \"Sample $i\"\n";
	}
	print $arpout "\t}\n";
	return 0;
}

sub StatDistrMessage {
my ($numsim) = @_;
die(qq/
Distr generates summary statistic distributions for a given demographic model\n
Usage: ABCutils.pl StatDistr [arguments]\n
Requirements: 
(1)fastsimcoal executable must be in the user's path
(2)arlsumstat executable, arl_run.ars, and ssdefs.txt file must all be in the same directory
- arl_run.ars is the settings file generated by the Windows version of Arlequin
- ssdefs.txt file lists the summary statistics to be computed\n
Arguments:
parfile    FILE   fastsimcoal format file of model parameters
outfile    FILE   File to write results to
numsim     INT    Number of simulations
arlexec    FILE   arlsumstat executable file
infsites          Use infinite sites model
metasamp          Pool subpopulations into 2 metapopulations
pop1id     INT    Subpopulations comprising metapop 1 (from top of parfile, 1-based)
pop2id     INT    Subpopulations comprising metapop 2 (from top of parfile, 1-based)
sitefst           Calculate Fst per site for each simulation\n
\nOutput format of .sitefst file is rows = simulations and columns = sites (in linear order)
\n/)
}

sub usage {
die (qq/
### version $version ###\n
Usage: ABCutils.pl [command] [arguments]\n
Commands:\n
catsims         Concatenate simulations from multiple .samp files outputted by fastsimcoal_sampler.pl\n
MaskCats    	Make fixed zero categories among multiple 2DSFS simulations variable\n   
rmvFixedZero    Remove summary statistics fixed for zero in at least one demographic model
                from .obs and simulation files\n
SetSpace    	Change all whitespace in a file to a single space\n
DistReject    	Find simulations with the closest Euclidean distance to the observed data\n
format2Dsim   	Tools for manipulating the output of fastsimcoal_sampler prior to using DistReject\n
format2Dobs    	Formats the observed 2D-SFS for use with DistReject\n
gofReject	Find simulations that fit most closely to the observed data using goodnesss-of-fit\n
best2Dsfs	Extracts the simulated joint spectrum that most closely fits the observed joint spectrum\n
fold2D          Folds 2D-SFS\n
StatDistr       Generates distribution of summary statistics for a given demographic model\n
\n/)
}