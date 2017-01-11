#!perl
#

use PDL::IO::Dcm qw/load_dcm_dir parse_dcms/;
use strict "vars";
use PDL::NiceSlice;
use PDL::IO::Sereal qw/wsereal/;
use Getopt::Tabular;
use Data::Dumper;


# copied and modified from stackoverflow or perlmonks thread (can't remember atm)
sub printStruct {
	my ($struct,$structName,$pre)=@_;
#    print "-----------------\n" unless (defined($pre));
#   
	my $res;
	#if (!ref($struct)){ # $struct is a scalar.
	if (ref($struct) eq "ARRAY") { # Struct is an array reference
#return ("ARRAY(".scalar(@$struct).")") if (@$struct>100);
		for(my$i=0;$i<@$struct;$i++) {
			if (ref($struct->[$i]) eq "HASH") {
				$res.=printStruct($struct->[$i],$structName."->[$i]",$pre." ");
			} elsif (ref($struct->[$i]) eq "ARRAY") { # contents of struct is array ref
				$res.= "$structName->"."[$i]: ()\n" if (@{$struct->[$i]}==0);
				my $string = printStruct($struct->[$i],$structName."->[$i]",$pre." ");
				$res.= "$structName->"."[$i]: $string\n" if ($string);
			} else { # contents of struct is a scalar, just print it.
				$res.= "$structName->"."[$i]: $struct->[$i]\n";
			}
		}
		#return($res);
	} elsif (ref($struct) eq "HASH"){ # $struct is a hash reference or a scalar
		foreach (sort keys %{$struct}) {
			if (ref($struct->{$_}) eq "HASH") {
				$res.=printStruct($struct->{$_},$structName."->{$_}",$pre." ");
			} elsif (ref($struct->{$_}) eq "ARRAY") { # contents of struct is array ref
				my $string = printStruct($struct->{$_},$structName."->{$_}",$pre." ");
				$res.= "$structName->"."{$_}: $string\n" if ($string);
			} else { # contents of struct is a scalar, just print it.
				$res.= "$structName->"."{$_}: $struct->{$_}\n";
			}
		}
		#return($res);
	} else {
		$res.= "$structName: $struct\n";
	}
#print "------------------\n" unless (defined($pre));
	return($res);
}

my ($d,$nifti,$sereal,$usage,$t);

my @opts=(
	['-n','boolean',0,\$nifti, 'Create .nii and .txt files'],
	['-t','boolean',0,\$t, 'serialise cardiac phases and time'],
	['-d','boolean',0,\$d, 'split not into lProtID but dicom series number'],
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

# how should we split series?
my $id=sub {$_[0]->hdr->{ascconv}->{"lProtID"};};
$id=sub {my $ret=$_[0]->hdr->{dicom}->{"Series Number"}; 
	$ret=~ s/^\s+|\s+$//g; $ret;} if $d;
#$id=~ s/^\s+|\s+$//g;
# loads all dicom files in this directory
my $dcms=load_dcm_dir($dir,$id);
die "no data!" unless (keys %$dcms);
print "Read data; ProtIDs: ",join ', ',keys %$dcms,"\n";
# sort all individual dicoms into a hash of piddles.
my $data=parse_dcms($dcms);

# save all data to disk
for my $pid (keys %$data) {
	print "Processing $pid.\n";
	if ($t) {
		print "-t: ",$$data{$pid}->info," \n";
		$$data{$pid}=$$data{$pid}->clump(6,3);
		print "-t: ",$$data{$pid}->info," \n";
	}
	if ($nifti) {
		require (PDL::IO::Nifti) || die "Make sure PDL::IO::Nifti is installed!";
		print "Creating Nifti $pid ",$$data{$pid}->info,"\n";
		my $ni=PDL::IO::Nifti->new;
		$ni->img($$data{$pid}->double->(,-1:0,)); #->reorder(0,1,9,4,2,3,5,6,7,8,10));
		$ni->write_nii($pre."_$pid.nii");
		open F,">",$pre."_$pid.txt";
		print F "### ASCCONV BEGIN ###\n";
		for my $k (sort keys %{$$data{$pid}->hdr->{ascconv}} ) 
			{print F "$k = ",$$data{$pid}->hdr->{ascconv}->{$k},"\n" }
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
		print "Writing file $pre\_$pid.srl\n";
		$$data{$pid}->wsereal("$pre\_$pid.srl"); 
	}
	#$$data{$pid}->(,,0,0,0,0,0,0,0;-)->double->wpic("$pre\_$pid.png");
}
