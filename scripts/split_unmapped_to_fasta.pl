#!/usr/bin/perl -w

#  split_unmapped_to_fasta.pl
#  (c) 2012 - Ryan M. Layer
#  Hall Laboratory
#  Quinlan Laboratory
#  Department of Computer Science
#  Department of Biochemistry and Molecular Genetics
#  Department of Public Health Sciences and Center for Public Health Genomics,
#  University of Virginia
#  rl6sf@virginia.edu
# 
# Licenced under the GNU General Public License 2.0 license.


use strict;
use Getopt::Long;
use File::Basename;
my $prog = basename($0);

sub soft_clip_len {
	my ($cigar) = @_;

	my $len = 0;
	#my $s = "$cigar\t";
	while ($cigar =~ m/(\d+)(M|I|D|N|S|H|P|X|=)/g) {
		#$s = $s .  "$1,$2\t";
		if ($2 ne 'M') {
			$len += $1;
		}
	}

	return $len;
}

sub print_usage {
	my ($msg) = @_;
    warn <<"EOF";

$msg

USAGE
  $prog [options]

DESCRIPTION

OPTIONS
	-h	Print this help message
	-b	Minimum number of unmapped baspairs 

EOF
	exit;
}

my $b;

my $help = 0;
GetOptions ("b=i"		=> \$b,
			"h"			=> \$help) or print_usage(); 

print_usage() if $help;

print_usage("No minimum basepairs") if not($b);

my $unmapped_flag = 4;


my $id = 0;
while ( my $l = <STDIN>) {
	chomp($l);
	if ($l =~ /^\@SQ/) {
		print $l . "\n";
	} else {
		my @a = split(/\t/, $l);
		my $flag = $a[1];
		my $cigar = $a[5];
		my $print_it = 0;

		if ( $flag & $unmapped_flag ) {
			$print_it = 1;
		} elsif (soft_clip_len($cigar) > $b) {
			$print_it = 1;
		}

		if ($print_it == 1) {
			my $name = $a[0];
			my $seq = $a[9];
			my $qual = $a[10];

			print "\@$name\_$id\n$seq\n+\n$qual\n";

			++$id;
		}
	}
}
