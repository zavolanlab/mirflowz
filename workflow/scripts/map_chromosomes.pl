#! /usr/bin/perl -w

# MAY 2019, Paula Iborra
# University of Basel

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

my @in = ();
my $column_delimiters_href_split = {
	'TAB' => q{\t},
	'COMMA' => ",",
	'DASH' => "-",
	'UNDERSCORE' => "_",
	'PIPE' => q{\|},
	'DOT' => q{\.},
	'SPACE' => " "
};

my $column_delimiters_href_join = {
        'TAB' => qq{\t},
        'COMMA' => ",",
        'DASH' => "-",
        'UNDERSCORE' => "_",
        'PIPE' => "|",
        'DOT' => ".",
        'SPACE' => " "
};

# a wrapper for converting between UCSC and ensembl chromosome representations from within galaxy
# convert_UCSC_ensembl.pl [input] [col] [delimiter] [genome] [out_file1]

die "Check arguments: $0 [input] [col] [delimiter] [map] [out_file1]\n" unless @ARGV == 5;
die "No columns specified: $ARGV[1]\n" if looks_like_number($ARGV[1]) == 0;
die "Delimeter must be one of TAB, COMMA, DASH, UNDERSCORE, PIPE, DOT, SPACE\n" unless defined $column_delimiters_href_split->{$ARGV[2]};

# process input
my $input = $ARGV[0];
$ARGV[1] =~ s/\s+//g;
my $col = --$ARGV[1];
my $delim = $ARGV[2];
my $map_file = $ARGV[3];
my $output = $ARGV[4];
my $delim_split = $column_delimiters_href_split->{$delim};
my $delim_join = $column_delimiters_href_join->{$delim};

open (MAP, "<$map_file") or die "Cannot open map file $map_file:$!\n";
my %chr_map;
while(my $line = <MAP>) {
	chop $line;
	next if grep /^#/, $line;
	my @map = split /\t/, $line;
	$map[1] = "remove" unless $#map;
	$chr_map{$map[0]} = $map[1];
}
close MAP;

open (IN,  "<$input") or die "Cannot open $input:$!\n";
open (OUT, ">$output") or die "Cannot create $output:$!\n";
while (my $line = <IN>) {
	chop $line;
	@in = split /$delim_split/, $line; 
	if(defined $in[$col] && defined $chr_map{$in[$col]}) {
		$in[$col] = $chr_map{$in[$col]};
		if($in[$col] eq "remove") {
			print "Removed line \"$line\" as chromosome does not have a proper mapping\n";
		} else {
			print OUT join($delim_join, @in), "\n";
		}
	} elsif(grep /^#/, $in[0]) {
		print OUT join($delim_join, @in), "\n";
	} else {
		print "Removed line \"$line\" as \"$in[$col]\" is not a valid chromosome name\n";
	}
}
close IN;
close OUT;


