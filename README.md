ABCutils
========

ABCutils is a collection of scripts for inferring population demography by fitting the 2D-site frequency spectrum (2D-SFS) using Approximate Bayesian Computation (ABC).

## Installation

Download ABCutils by either cloning it from github

% git clone https://github.com/tplinderoth/ABCutils.git

or download the ZIP file.

## fastsimcoal_sampler.pl

fastsimcaol_sampler is simply a wrapper around the coalescent simulator, fastsimcoal (Excoffier & Foll 2011), and should also be compatible with [fastsimcoal2](http://cmpg.unibe.ch/software/fastsimcoal2/) (Excoffier et al. 2013). As such, the documentation for those programs is a good resource for setting up simulations. For a description of how to run the script, simply run fastsimcoal\_sampler.pl without arguments:

% ./fastsimcoal_sampler.pl

Version 1.1.4

Usage: 
fastsimcoal_sampler.pl [options]

fastsimcoal must be installed and in path

Options:
 
--outfile	CHAR		Outfile name preffix (without path)

--estfile	FILE		Fastsimcoal estimation (*.est) format file ("rules" section last)

--tplfile	FILE		Fastsimcoal template (*.tpl) format file

--recomb	FLOAT		Recombination rate between adjacent nucleotides within sequence blocks [0]

--mut		FLOAT/CHAR 	Mutation rate as float or parameter name from estimation file [2.2e-09]
 
--trans		FLOAT		Transition rate (0.33 implies no mutation bias) [0.33]

--seqlist	FILE		Fasta format file of sequence blocks (list of chromosomes or contigs)

--numseq	INT		Number of sequence blocks to simulate at a time (NA when meta2DSFS = 1) [500]

--numsim	INT		Number of simulations (each sim uses different parameter values) [1000]

--meta2DSFS	0|1		Calculate 2D-SFS for two metapopulations if (1) (each comprising >= 1 demes)
	
--pop1		INT 		First population ID number in joint SFS (from top of template file, 0-based)
				If meta2DSFS = 1, these are the >= 1 demes comprising metapopulation 1 (1-based)
 
--pop2		INT		Second population ID number in joint SFS (from top of template file, 0-based)
				If meta2DSFS = 1, these are the >= 1 demes comprising metapopulation 2 (1-based)

--norm		0|1		Normalize SFS if 1, else keep absolute counts if (0) [0]

--rmvfixed	0|1		Remove fixed categories in SFS if (1) (prior to standardization if --norm 1), 
				else keep fixed categories if (0) [0]

--rmvMutation	CHAR		Remove SNP sites of certain mutation type(s) (only works when meta2DSFS = 1)

--p1missing	FLOAT		Discard site if at least FLOAT percent of individuals in pop1 are missing data for the site. 
				Otherwise, missing data will be treated as major/ancestral allele [50]

--p2missing 	FLOAT		Discard site if at least FLOAT percent of individuals in pop2 are missing data for the site.
				Otherwise, missing data will be treated as major/ancestral allele [50]


The template (*.tpl) input file:
-Should not include fastsimcoal "genetic info" section
-The file should end after "historical event" section

The estimation (*.est) input file:
-Should list sections in the order [parameters], [complex parameters], [conditional], [rules]
-All parameter values drawn from priors or calculated must only contain uppercase letters and
 '_' in the name
-Parameters must be defined in estimation file before other parameters that depend on them are defined
-Loguniform priors should not have zero as a boundary

-The fasta format file supplied to --seqlist is used to determine how many independent linkage blocks to simulate.
 Each sequence in the fasta file is treated as a completely independent, non-recombining locus (linkage block), which has
 length equal to the sequence length. For example, each fasta sequence could be a chromosome. The recombination rate
 between adjacent sites within a linkage block is specified with the --recomb argument.

The output file:
-Names for unfolded SFS indicate element in the 2k+1 SFS: E<pop1_sfs_category><pop2_sfs_category>

Calculating 2DSFS with metapopulations (meta2DSFS = 1):
-pop1 and pop2 arguments each take a list of >= 1 INT values for the demes comprisiong each
 (meta)population, respectively, in the arlequin project file (ex: --pop1 1 2 3 --pop2 4 5 6).
 The deme numbers correspond to the ascending order in which they are listed in the template file (1-indexed).
-The type of mutation to remove is denoted by <ancestral allele><derived allele>
 Example to remove C->T and G->A mutations: --rmvMutation CT GA

## ABCutils.pl

ABCutils.pl is for manipulating the simulations generated with fastsimcoal\_sampler.pl and estimating demographic parameters. Some of the functions make the output of fastsimcoal_sampler.pl compatible with [ABCtoolbox](http://cmpg.unibe.ch/software/ABCtoolbox/) (Wegmann et al. 2010). For an overview of the subroutines run ABCutils.pl without arguments, and for help on each subroutine run 'ABCutils.pl <subroutine name>' without arguments.

% ./ABCutils.pl

Usage: ABCutils.pl [command] [arguments]

Commands:

catsims         Concatenate simulations from multiple .samp files outputted by fastsimcoal_sampler.pl

MaskCats    	Make fixed zero categories among multiple 2D-SFS simulations variable
   
rmvFixedZero    Remove summary statistics fixed for zero in at least one demographic model
                from observed (*.obs) and simulation files

SetSpace    	Change all whitespace in a file to a single space

DistReject    	Find simulations with the closest Euclidean distance to the observed data

format2Dsim   	Tools for manipulating the output of fastsimcoal_sampler prior to using DistReject

format2Dobs    	Formats the observed 2D-SFS for use with DistReject

gofReject	Find simulations that fit most closely to the observed data based on goodnesss-of-fit

best2Dsfs	Extracts the simulated 2D-SFS that most closely fits the observed 2D-SFS

fold2D          Folds 2D-SFS

StatDistr       Generates distribution of summary statistics for a given demographic model

### ABCutils.pl subroutines

#### catsims

Usage: ABCutils.pl catsims <folder_containing_sims> <simulation_file>

Loops through a directory (and subdirectories) containing simulation files and appends the simulation output to the "simulation_file"

Notes: -The 'folder\_containing_sims' can be a parent directory
      -If simulation_file does not exist, it will be created

./ABCutils.pl MaskCats

#### MaskCats

Usage: ABCutils.pl MaskCats [options]

Force nonvariable 2D-SFS categories to be variable by changing the value of the nonvariable category in a randomly chosen simulation to 1e-100. This is in order to run ABCtoolsbox ABCestimator function on 2D-SFS simulations for which all values in a category among simluations is 0.

Options:
-i    FILE    Simulation input file; should be *.samp format outputted by fastsimcoal_sampler.pl)
-o    CHAR    Outfile name (including path)
-n    INT     Iterate through INT simulations to determine if a 2DSFS category is fixed for 0,
              otherwise, if (0) iterate through all simulations [0]
#### rmvFixedZero

Removes summary statistics fixed for 0 from simulation files and the *.obs file for use with ABCtoolbox ABCestimator. If a summary statistic is fixed for 0 in at least one model, this statistic
is removed from all simulation files and the *.obs file.

Usage: ABCutils.pl rmvFixedZero -n INT -c INT <*.obs\_file> <simulation_file(s)>

-n specifies how many simulations in each sim file to read through in order to determine which summary statistics are fixed. The default is to read all simulations.

-c specifies that if a summary statistic has INT or fewer simulations with a nonzero value, the statistic is considered fixed and removed.

Notes: New *.obs and simulation file(s) are dumped in the directory of the original files with 'estimator' appended to the name

#### DistReject

Finds closest simulations to observed data in terms of Euclidean distance using a specified tolerance

Usage: ABCutils.pl DistReject [options]

Options:
--numsim     INT     number of simulations to examine [5000]
--numkeep    INT     number of closest simulations to retain [1000]
--percent    FLOAT   percent of closest simulations to retain [10]
--stat_start INT     starting column containing summary statistics in simulation file (the first column is 1)
--sim        FILE    file containing simulations
--obs_sfs    FILE    file containing observed data
--outfile    CHAR    outfile name (including path)
--multiobs           compare nth simulation to nth observation (requires multiple observations)

Note: options --numkeep and --percent indirectly set the tolerance in an ABC rejection framework.
      When both options are specified, whichever specifies the fewest simulations to retain is used.

#### format2Dsim

Manipulates the set of simulated 2D-SFS outputted by fastsimcoal_sampler.pl
 
Usage: ABCutils.pl format2Dsim [options]

Options:
--pop1n     INT    number of diploid individuals representing the horizontal (x-axis) of 2D-SFS
--pop2n	    INT    number of diploid individuals representing the vertical (y-axis) of 2D-SFS
--sfsfile   FILE   *.samp format file containing 2D-SFS simulations outputted by fastsimcoal_sampler.pl
--outfile   CHAR   outfile name (including path). If not specified the output will be dumped to the same directory as the *.samp input file
--norm      0|1    normalize the 2D-SFS if (1) [0]
--maskFixed 0|1    mask fixed categories in 2D-SFS if (1) [0]
--bin	    INT	   calculate diagonal bins of 2D-SFS categories; INT is 2 x the number of categories on
		   either side of each bin's diagonal (~bin width, which must be a factor of 2)
--do1dsfs   INT    calculate folded 1D-SFS for pop1 if (0), pop2 if (1), or both populations (2) [null]

Notes:
- Format2Dsim automically folds the 2D-SFS.
'--maskFixed 1' is automatically induced if fastsimcaol_sampler was run with --rmvfixed 1

#### format2Dobs

Formats the 2D-SFS for use as the observed input file for DistReject

Usage: ABCutils.pl format2Dobs [options]

Options:
--insfs       FILE    2D-SFS input file
--outfile     FILE    output file name (including path)
--fold        0|1     fold 2D-SFS if (1) [0]
--rmvFixed    0|1     remove invariable categories of 2D-SFS if (1) [0]
--maskFixed   0|1     mask invariable categories of 2D-SFS if (1) [0]
--norm	      0|1     normalize 2D-SFS if 1 [0]
--bin	      INT     calculate diagonal bins of 2D-SFS categories; INT is 2 x the number of categories on
                      either side of each bin's diagonal (~bin width, which must be a factor of 2)
--do1dsfs     INT     calculate folded 1D-SFS for pop1 (0), pop2 (1), or both populations (2) [null]
  
NOTES:
-format2Dobs assumes that the input agrees with fastsimcoal joint sfs format: 
 a matrix where pop1 allele frequency categories are COLUMNS and pop2 categories are ROWS
-if no outfile name is provided the output is dumped as a *.fmt file to the same directory
 as the input

#### gofReject

Calculates goodness-of-fit statistic as (observed - expected)^2/expected between simulated and observed
summary statistics and retains simulations closest to the observed data

ABCutils.pl gofReject [options]

Options:
--obsfile	FILE	Observed summary stats (e.g. SFS bins as outputted by format2Dobs)
--simfile	FILE	Simulated summary stats (e.g. SFS bins as outputted by format2Dsim)
--outfile	CHAR	Output file name where results will be dumped
--stat_start	INT	Field number of the first summary statistic in simfile (1-based) to use in fitting.
			All summary statistics in fields >= INT will be used.
--num_keep	INT	Find the closest INT simulations to the observed 2D-SFS [0]
--percent_keep	FLOAT	Find the closest FLOAT percent of simulations to the observed 2D-SFS [10] 
--numsim	INT	Number of simulations (from beginning of simfile) to calculate goodness-of-fit for

Notes:
The smaller number of simulations to find between --num\_keep and --percent_keep will be used

#### best2Dsfs

Extracts the simulated 2D-SFS that most closely fits the observed 2D-SFS.

Usage: ABCutils.pl best2Dsfs <sim_file> [arguments]

Arguments:
--simfile    FILE  File of 2D-SFS simulations
--distfile   FILE  File containing distances of observed from expected 2D-SFS
--outfile    CHAR  Name of output file
--pop1n      INT   Number of diploid individuals in population 1
--pop2n      INT   Number of diploid individuals in population 2
--fit        INT   Field number (1-based) of column specifying distance of observed from expected [2]
--dofold           folds the joints spectrum

Notes:
--simfile should be output from fastsimcoal_sampler.pl or ABCutils 'format2Dsim' 
--distfile should be output from ABCutils 'DistReject' or 'gofReject'

#### fold2D

Folds a 2-dimensional site frequency spectrum

Usage: ABCutils.pl fold2D <unfolded_SFS_file>

Note: The input 2D-SFS should not contain any headers.

#### StatDistr

Generates summary statistic distributions under a given demographic model. The program arlsumstat (Excoffier & Lischer 2010)
is used to calculate some summary statistics so it's README is a useful resource, http://cmpg.unibe.ch/software/arlequin35/man/arlsumstat_readme.txt.

Usage: ABCutils.pl StatDistr [arguments]

Requirements:
(1) fastsimcoal executable must be in the user's path

(2) arlsumstat executable, arl_run.ars, and ssdefs.txt file must all be in the same directory
    - arl_run.ars can be generated by the Windows version of Arlequin
    - ssdefs.txt file lists the summary statistics to be computed

Arguments:
--parfile    FILE   fastsimcoal parameter format (*.par) file of model parameters
--outfile    FILE   File to write results to
--numsim     INT    Number of simulations
--arlexec    FILE   arlsumstat executable file
--infsites          Use infinite sites model
--metasamp          Pool subpopulations into 2 metapopulations
--pop1id     INT    Subpopulations comprising metapop 1 (from top of parfile, 1-based)
--pop2id     INT    Subpopulations comprising metapop 2 (from top of parfile, 1-based)
--sitefst           Calculate per site Fst for each simulation

Notes:
Output format of .sitefst file is rows = simulations and columns = sites (in increasing order)

## ABCsampler_chromblock_generator.pl

Use to format a *.par input file for DNA sequence simulation with fastsimcoal so that each contig represents a block on a single chromosome. 
This script only works for a single chromosome.

Usage: ABCsampler\_chromblock_generator.pl [options]

Options: 
-p FILE  fastsimcoal_sampler parameter (*.par) file to which chromosome block info will be appended
-f FILE  fasta file of contigs (i.e. chromosome blocks)
-o CHAR  outfile name preffix
-r FLOAT recombination rate between adjacent contigs and within contigs
-s FILE  list of recombination rates for particular contigs (contig rate)
-u FLOAT mutation rate
-v FILE  list of mutation rates for particular contigs (contig rate)
-t FLOAT transition rate (0.33 implies no transition bias) [0.33]
-w FILE  list of transition rates for particular contigs (contig rate)

## Example commands to obtain demographic parameter posterior distributions by fitting 2D-SFS bins

### Simulate 2D-SFS between pooled metapopulations for historic and modern sampling using parameter values drawn from prior distributions

% ./fastsimcoal\_sampler.pl --outfile alpinusG\_out --estfile ./ABCutils/example\_files/alpyosG.est --tplfile ./ABCutils/example\_files/alpyosG.tpl --recomb 0 --mut 2.2e-9 --trans 0.725 --seqlist ./ABCutils/example\_files/alpinus\_test\_chr.fa --numsim 25000 --pop1 1 2 3 4 8 10 --pop2 11 12 13 14 15 16 17 19 --norm 0 --rmvfixed 0 --meta2DSFS 1 --p1missing 90 --p2missing 90 --rmvMutation CT GA

### Calculate simulated 2D-SFS bins

% ./ABCutils.pl format2Dsim --pop1n 48 --pop2n 56 --sfsfile ./ABCutils/example\_files/alpinusG\_out.samp --outfile alpinusG\_fold\_bin2.txt --norm 0 --maskFixed 1 --bin 2

### Calculate observed 2D-SFS bins

% ./ABCutils.pl format2Dobs --insfs ./ABCutils/example\_files/T\_alpinus\_Yosemite_unfolded.2dsfs --outfile ./ABCutils/example\_files/ynp\_alpinus\_fold\_bin2\_mask.obs --fold 1 --maskFixed 1 --bin 2

### Use rejection sampling to obtain parameter posterior distributions

% ./ABCutils.pl DistReject --numsim 25000 --numkeep 200 --stat_start 12 --sim ./ABCutils/example\_files/alpinusG\_fold\_bin2.txt --obs\_sfs ./ABCutils/example\_files/ynp\_alpinus\_fold\_bin2\_mask.obs --outfile ./ABCutils/example\_files/alpinusG\_fold\_bin2.dist

##### Notes on testing and depencies

Versions of programs that fastsimcoal_sampler.pl and ABCutils.pl have been tested with:
fastsimcoal 1.1.2
arlsumstat 3.5.2.1 (30.04.15)

##### Contact
Tyler Linderoth: tylerp.linderoth@gmail.com
