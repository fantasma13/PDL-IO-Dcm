#!perl
#

use PDL::IO::Dcm qw/load_dcm_dir parse_dcms/;
use strict "vars";
use PDL::IO::Sereal qw/wsereal/;
use Getopt::Tabular;

my ($nifti,$sereal,$usage);

my @opts=(
	['-n','boolean',0,\$nifti, 'Create .nii and .txt files'],
	['-s','boolean',1,\$sereal, 'Create sereal (defautl)'],
	['-h', 'boolean',1,\$usage, 'print this help'],
);
&GetOptions (\@opts, \@ARGV) || exit 1;

if ($usage or $#ARGV ==0 ) {
	print "read_dcm.pl [options] <directory> <output_prefix>\n"; 
	print "Reads and sorts Siemens dicom files and converts them into either serealised piddles or nifti/text.\n";
	for my $arg (@opts) { 
		print "$$arg[0] [ $$arg[1] ] : -- $$arg[-1]\n"; 
	} 
	exit; 
}

my $dir=shift; 
my $pre=shift;

# loads all dicom files in this directory
my $dcms=load_dcm_dir($dir);
die "no data!" unless (keys %$dcms);
print "Read data; ProtIDs: ",join ', ',keys %$dcms,"\n";
# sort all individual dicoms into a hash of piddles.
my $data=parse_dcms($dcms);

# save all data to disk
for my $pid (keys %$data) {
	print "Processing $pid \n";
	if ($nifti) {
		require (PDL::IO::Nifti) || die "Make sure PDL::IO::Nifti is installed!";
		print "Creating Nifti $pid ",$$data{$pid}->info,"\n";
		my $ni=PDL::IO::Nifti->new;
		$ni->img($$data{$pid}->double); #->reorder(0,1,9,4,2,3,5,6,7,8,10));
		$ni->write_nii($pre."_$pid.nii");
		open F,">",$pre."_$pid.txt";
		for my $k (sort keys %{$$data{$pid}->hdr->{ascconv}} ) 
			{print F "$k = ",$$data{$pid}->hdr->{ascconv}->{$k},"\n" }
		for my $k (sort keys %{$$data{$pid}->hdr->{dicom}} ) 
			{print F "$k: ",$$data{$pid}->hdr->{dicom}->{$k},"\n" }	
		close F;
	} 
	print "Nifti? $nifti Sereal? $sereal write? ",((! $sereal) and $nifti),"\n";
	unless ((! $sereal) and $nifti) { 
		print "Writing file $pre\_$pid.srl\n";
		$$data{$pid}->wsereal("$pre\_$pid.srl"); 
	}
}
