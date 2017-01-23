#!perl
#

use PDL::IO::Dcm qw/printStruct load_dcm_dir parse_dcms/;
use strict "vars";
use PDL::NiceSlice;
use PDL;
use PDL::IO::Sereal qw/wsereal/;
use Getopt::Tabular;
use Data::Dumper;
my %opt;
my $plugin="Primitive";
my ($d,$nifti,$sereal,$usage,$t,$sp);
my @opts=(
	#['-d','boolean',0,\$d, 'split not into lProtID but dicom series number'],
	['-u', 'string' ,1,\$plugin,'specify the plugin to process data'],
	['-h', 'boolean',1,\$usage, 'print this help'],
	['-p', 'boolean',0,\$sp, 'split slice groups'],
	['-s','boolean',1,\$sereal, 'Create sereal (defautl)'],
	['-t','boolean',0,\$t, 'serialise cardiac phases and time'],
	['-i','boolean',0,\$nifti, 'Create .nii and .txt files'],
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
$opt{plugin}=$plugin;
$plugin = "PDL/IO/Dcm/Plugins/$opt{plugin}"; # set plugin

require ($plugin.'.pm') || die "Plugin $plugin could not be loaded!\n"; 
print "Plugin: PDL::IO::Dcm::Plugins::$opt{plugin}\n";
eval("PDL::IO::Dcm::Plugins::$opt{plugin}")->import( qw/setup_dcm/);


setup_dcm(\%opt);
$opt{split}=$sp;
# loads all dicom files in this directory
my $dcms=load_dcm_dir($dir,\%opt);
die "no data!" unless (keys %$dcms);
print "Read data; ProtIDs: ",join ', ',keys %$dcms,"\n";
$opt{Nifti}=$t;
# sort all individual dicoms into a hash of piddles.
my $data=parse_dcms($dcms,\%opt);

# save all data to disk
for my $pid (keys %$data) {
	print "Processing $pid.\n";
	if ($t) { # use Dicom series number
		#print "-t: ",$$data{$pid}->info," \n";
		#$$data{$pid}->hdrcpy(1);
		#$$data{$pid}=$$data{$pid}->clump(6,3);
		#pop @{$$data{$pid}->hdr->{Dimensions}};
		#print "-t: ",$$data{$pid}->info," \n";
	}
	if ($nifti) {
		require (PDL::IO::Nifti) || die "Make sure PDL::IO::Nifti is installed!";
		print "Creating Nifti $pid ",$$data{$pid}->info,"\n";
		my $ni=PDL::IO::Nifti->new;
		$ni->img($$data{$pid}->double->(,-1:0,)); #->reorder(0,1,9,4,2,3,5,6,7,8,10));
		$ni->write_nii($pre."_$pid.nii");
		open F,">",$pre."_$pid.txt";
		print F "## Generated using PDL::IO::Dcm\n\n";
		print F "dimensions: ",join ' ',(join ' ',@{$$data{$pid}->hdr->{Dimensions}})."\n\n";
		print F "### ASCCONV BEGIN ###\n";
		for my $k (sort keys %{$$data{$pid}->hdr->{ascconv}} ) 
			{ print F "$k = ",$$data{$pid}->hdr->{ascconv}->{$k},"\n" ;}
		print F "### ASCCONV END ###\n\n";
		print F "*** Parameters extracted from dicom fields \n";
		for my $k (sort keys %{$$data{$pid}->hdr->{dicom}} ) {
			print F printStruct( $$data{$pid}->hdr->{dicom}->{$k},$k,) ;	
			#print F Dumper( $$data{$pid}->hdr->{dicom}->{$k}); 
		}	
		print F "*** End of Parameters \n\n";
		close F;
	} 
	print "Nifti? $nifti Sereal? $sereal write? ",((! $sereal) and $nifti),"\n";
	unless ((! $sereal) and $nifti) { 
		print "Writing file $pre\_$pid.srl\n", $$data{$pid}->info;
		$$data{$pid}->wsereal("$pre\_$pid.srl"); 
	}
	#$$data{$pid}->(,,0,0,0,0,0,0,0;-)->double->wpic("$pre\_$pid.png");
}
