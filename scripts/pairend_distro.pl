#!/usr/bin/perl -w
#  pairend_distro.pl
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
use Statistics::Descriptive;
use Getopt::Long;
use File::Basename;
my $prog = basename($0);

sub print_usage {
	my ($msg) = @_;
    warn <<"EOF";

$msg

USAGE
  $prog [options]

DESCRIPTION

OPTIONS
	-h	Print this help message
	-rl	Read length
	-X	Number of stdevs from mean to extend
	-N	Number to sample 
	-o	Out file

EOF
	exit;
}

my $read_length;
my $X;
my $N;
my $o_file;

my $help = 0;
GetOptions ("rl=i"		=> \$read_length,
			"X=i"		=> \$X,
			"N=i"		=> \$N,
			"o=s"		=> \$o_file,
			"h"			=> \$help) or print_usage(); 

print_usage() if $help;

print_usage("No rl") if not($read_length);
print_usage("No X") if not($X);
print_usage("No N") if not($N);
print_usage("No o") if not($o_file);

my $required = 67; # read paired, propper pair, first in pair
my $restricted = 384; # read unmapped, mate unmapped

my @L;
my $c = 0;

while ( defined((my $l = <STDIN>)) && ($c <= $N) ) {
	chomp($l);
	my @a = split(/\t/, $l);
	
	if ( ( ($a[1] & $required) == $required ) &&
		 ( ($a[1] & $restricted) == 0) &&
		 ( $a[8] >= 0 ) ) {

		 push(@L, $a[8]);
		 $c++;
	}
}

my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@L);

my $mean = $stat->mean();
my $stdev = $stat->standard_deviation();

my $start = $read_length;
my $end = int($mean + $X*$stdev);


my @H = (0) x ($end - $start + 1);
my $sum = 0;
for (my $i = 0; $i < $c; $i++) {
	if (($L[$i] >= $start) && ($L[$i] <= $end)) {
		my $j = $L[$i] - $start;
		$H[ $j ] = $H[ $L[$i] - $start ] + 1;
		$sum = $sum + 1;
	}
}

open(FILE, "> $o_file");
for (my $i = 0; $i < ($end - $start); $i++) {
    #print FILE "$i\t" . ($H[$i]/$sum) . "\n";
    if ($H[$i]/$sum == 0) {
        last;
    } else {
	print FILE "$i\t" . ($H[$i]/$sum) . "\n";
    }
}

close(FILE);
print "mean:$mean\tstdev:$stdev\n";
