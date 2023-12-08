#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: sam_uncollapse.pl
### Created: Nov 21, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Requirements: GetOpt::Long
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> PRAGMAS / PACKAGES <---#
use strict;
use warnings;
use Getopt::Long;

#---> USAGE <---#
my $usage_info = &usage;

#---> OPTIONS / ARGUMENTS <---#
my $usage = '';
my $quiet = '';
my $suffix = '';
my $in = '';
my $out = '';
my $options_result = GetOptions (
        'usage|help' => \$usage,
        'quiet' => \$quiet,
        'suffix' => \$suffix,
        #-----------------------#
        'i|in=s' => \$in,
        'o|out=s' => \$out
);
die $usage_info if $usage || !$options_result;
die $usage_info if !$in || !$out;

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting '$0'...\n" unless $quiet;

#---> BODY <---#

        #---> Read & re-write file <---#
        &sam_uncollapse($in, $out, $suffix);

#---> STATUS MESSAGE <---#
print "Done.\n" unless $quiet;

#---> PROGRAM EXIT <---#
exit 0;
#================#
#    MAIN END    #
#================#


#=======================#
#   SUBROUTINES START   #
#=======================#
sub usage {
### Function: Returns usage information for current script
### Accepts: n/a
### Returns: String with usage information
### Type: Specialized
'Usage: perl ./sam_uncollapse.pl [OPTIONS] --in [SAM] --out [SAM]

Description: Reverses the collapsing of reads with identical sequences as done with "fastx_collapser" ("FASTX Toolkit") or similar. Reads and writes files in SAM format. Each line is printed n times, where n is the suffix appended to the read/query name via a dash.

==================================================
Required arguments:
--in    Input SAM file.
--out   Output SAM file.
==================================================
Optional arguments:
--suffix	Add serial number suffix to each QNAME during uncollapsing (separated by a ".") to allow distinction of multimappers by QNAME
--usage|help    Show this information and die
--quiet Shut up!

Comments:
CAUTION: Only marginal validation of the input file type/format performed!

Version 1.2.1 (2023-12-08)
Written by Alexander Kanitz on 2013-11-21
';
}
#-----------------------#
sub sam_uncollapse {
### Function: For each line of a SAM file, parses the identifier QNAME for the presence of a number n appended to its end via a dash ('-') and re-writes the line n times. Header lines are reproduced as they are.
### Accepts: 1. Filename [SAM]
### Returns: n/a
### Dependencies: n/a
### Type: Generic

        #---> PASS ARGUMENTS ---#
        my $in = shift;
        my $out = shift;
	my $suffix = shift;

        #---> STATUS MESSAGE <---#
        print STDERR "Processing file '$in'..." . "\n" unless $quiet;

        #---> SUBROUTINE VARIABLES <---#
        my $line;

        #---> BODY <---#

                #---> Open files <---#
                open IN, "<", $in;
                open OUT, ">", $out;

                #---> Push line to array <---#
                while ($line = <IN>) {
                        # Print header line
                        print OUT $line and next if $line =~ m/\A\@\w\w\t/;
                        # Get QNAME
                        my ($id, $rest) = split /\t/, $line, 2;
                        # Find and remove appended copy number n
                        $id =~ /^([^_-]+)-(\d+).?\d+/;
                        # Write appended copy number n to variable
                        my $read = $1;
                        my $repeat = $2;
			# If --suffix option is set...
			if ($suffix) {
				# Iterate over number of identical reads/alignments
				for my $suffix (1..$repeat) {
		                        # Recreate line with suffix
                		        $line = join "\t", "$read.$suffix", $rest;
		                        # Print line
		                        print OUT $line;
				}

			}
			# Else...
			else {
                	        # Recreate line
        	                $line = join "\t", $read, $rest;
                        	# Print line n times
	                        print OUT $line x $repeat;
			}
                }

                #---> Close file <---#
                close OUT;
                close IN;

        #---> STATUS MESSAGE <---#
        print STDERR "File '$out' written." . "\n" unless $quiet;

}
#=======================#
#    SUBROUTINES END    #
#=======================#
