#!perl
#

use PDL::IO::Dcm qw/load_dcm_dir parse_dcms/;
use PDL;
use 5.10.0;
use strict;
use PDL::IO::Sereal qw/wsereal/;


my $dir=shift; 
my $pre=shift;

my $dcms=load_dcm_dir($dir);
die "no data!" unless (keys %$dcms);
my $data=parse_dcms($dcms);

for my $pid (keys %$data) {
	$$data{$pid}->wsereal("$pre\_$pid.srl");
}
