#!/usr/bin/perl
#
package PDL::IO::Dcm;

=head1 NAME

PDL::IO::Dcm - Reads dicom files, sorts them and stores the result into piddles with headers 

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS


This is inteded to read and sort dicom images created by medical imaging devices.

	# loads all dicom files in this directory
	my $dcms=load_dcm_dir($dir);
	die "no data!" unless (keys %$dcms);
	print "Read data; ProtIDs: ",join ', ',keys %$dcms,"\n";
	# sort all individual dicoms into a hash of piddles.
	my $data=parse_dcms($dcms);


This software is tailored to Siemens MRI data based on the author's need. For 
general usage, the read_dcm function should probably be moved to vendor specific
plugins. 

=head1 Some notes on Dicom fields and how they are stored/treated

The image data field is stored as the piddle, the other dicom
elements are stored in the header under the raw_dicom key. 

The Siemens protocol ASCCONV part is stored in the ascconv key. 

Key 0029,1010 is the Siemens specific field that contains the ICE
miniheaders with dimension information - and position in matrix
0029,1020 is deleted from the header, it is big, containing the whole
protocol. The important part is parsed into ascconv.

Keys are parsed into a hash under the dicom key using DicomPack
methods to unpack. For more complex structures this does not seem to
work well.

The header fields IceDims and IcePos are used for sorting datasets. 

Piddles are created for each lProtID value. 


=head1 SUBROUTINES/METHODS

=head2 read_text_hdr

parses the ASCCONV part of Siemens data header into the ascconv field of the
piddle header. All special characters except [a-z0-9]i are converted to _ -- no
quoting of hash keys required!

=head2 load_dcm_dir

reads all dicom files in a dicrectory and returns a hash of piddles containing
sorted N-D data sets. Uses ascconv field lProtID as keys. 

=head2 parse_dcms 

Parses and sorts a hash of hashes of dicoms (such as returned by load_dcm_dir)
based on lProtID and the ICE_Dims field in 0029_1010. Returns a hash of piddles 
(lProtID).

=head2 unpack_field

unpacks dicom fields and Walks subfield structures recursively.


=head2 read_dcm

reads a dicom file and creates a piddle-with-header structure.

=cut

use PDL;
use PDL::NiceSlice;
use List::MoreUtils; # qw{any};
use Data::Dumper;
use DicomPack::IO::DicomReader;
use Storable qw/dclone/;
use DicomPack::DB::DicomTagDict qw/getTag getTagDesc/;
use DicomPack::DB::DicomVRDict qw/getVR/;
use Exporter;
#use PDL::IO::Nifti;
use strict;
#use PDL::IO::Sereal;
use 5.10.0;

our @ISA=qw/Exporter/;
our @EXPORT_OK=qw/read_text_hdr read_dcm parse_dcms load_dcm_dir/;

my @key_list=("Instance Number",,'Window Center','Content Time',
	'Nominal Interval','Instance Creation Time','Largest Image Pixel Value',
	'Trigger Time','Window Width','Acquisition Time','Smallest Image Pixel Value',
);

sub read_text_hdr {
    my $f=shift; # File
    my $self=shift;
    open (HDR,'<',\$f) || die "no header !";
    my $l;
    #say "file $f line $l";
    do {$l=<HDR>; } until ($l=~/ASCCONV BEGIN/);
    while (($l=<HDR>)!~/ASCCONV END/) {
        chomp $l;
        if ( $l) {
            chomp (my ($key,$val)=split /\s*=\s*/,$l);
            chomp($key);
            $key=~s/[\[\].]/_/g;
            $self->hdr->{ascconv}->{$key}=$val;
        }
    }
    close HDR;

}

sub unpack_field{ 
	my $id=shift;
	my $tag=shift;
	my $packstring;
	my $value=shift;
	my $return; #=shift;
	#say "id $id value? ", ref( $value) unless ($id=~/0029/);
	if (ref($value) eq 'ARRAY') {
		my @vs=();
		for my $n ($#$value) {
			push @vs,unpack_field ("$id/$n",getTag("$id/$n"),$$value[$n],$return); 
		}
		$return=\@vs;
	} elsif (ref ($value) eq 'HASH') {
		my %vh=();
		for my $v (keys %$value) {
			#(my $w=$v)=~s/([0-9a-fA-F]{4}),([0-9a-fA-F]{4})/$1_$2/ ;
			#say "key $v ";
			$vh{$v}=unpack_field("$id/$v",getTag("$id/$v"),$$value{$v},$return);		
			#say "hash $id/$v:  $vh{$v} " unless $id=~/0029/;
		}
		$return=\%vh;
	} else {
	my $vr=substr($value,0,2);
	$packstring=join ('',(eval {getVR($vr)->{type}}||'a').'*'); 
	#say "ID $id tag $tag ps $packstring ";
	if ($vr eq 'XX' and defined $tag) {
		#say "ID $id, Tag $$tag{desc} ", %{$$tag{vr}};
		$packstring=join '',map {((getVR($_))->{type}||'a').'*'} 
			keys %{$$tag{vr}};
	}
	$return=unpack ($packstring,substr($value,3,));
	# split vector string into list
	#($value=[split /\\/,$value]) if $id =~ /0020,003[27]/; 
	# position and orientation
	#($value=[split /\\/,$value]) if $id =~ /0028,0030/; # pixel size
	#say "Image Rows $vr $packstring $value " if ($$tag{desc} eq 'Rows');
	}
	#say "ID $id value $value";
	$return;
}

sub read_dcm {
	my $file=shift;
	my $dcm=DicomPack::IO::DicomReader->new($file) || return; 
	my $h=unpack('S',substr ($dcm->getValue('Rows','native'),3,2));
	my $w=unpack('S',substr ($dcm->getValue('Columns','native'),3,2));
	my $data=$dcm->getValue('PixelData','native');
	my $datatype= (substr($data,0,2));
	say "datatype $datatype";
	my $pdl=zeroes(ushort,$w,$h) if ($datatype =~/OW|XX/); 
	$pdl->make_physical;
	${$pdl->get_dataref}=substr($data,3);
	$pdl->upd_data;
	#say "Rows $h";
	read_text_hdr($dcm->getValue ('0029,1020','native'),$pdl); # The protocol is in here
	$pdl->hdr->{raw_dicom}=$dcm->getDicomField;
	delete $pdl->hdr->{raw_dicom}->{'0029,1020'}; # Protocol
	#delete $pdl->hdr->{raw_dicom}->{'0029,1010'}; # also big structure, not sure what
		# Ice Dimension string: 13 numbers or X=undef
	#say "Ice Dims ",$dcm->getValue('0029,1010','native')=~/ICE_Dims/;
	my @d=$dcm->getValue('0029,1010','native')=~/ICE_Dims.{92}((_?(X|\d+)){13})/s; 
	my $dims=shift @d;
	#say "Raw Dims @d";
	$pdl->hdr->{IceDims}=$dims || die "No Ice Dims ",$file; #pdl->hdr->{raw_dicom}->{'0029,1010'}; #[split '_',$dims{$pid}=~s/X/0/r];
	delete $pdl->hdr->{raw_dicom}->{'7fe0,0010'}; # Pixel data
	for my $id (keys %{$pdl->hdr->{raw_dicom}}) {
		my $tag=getTag($id);
		my $value=unpack_field($id,$tag,$dcm->getValue($id,'native')); 
		#say "id $id v $value";
		if (defined $tag) {
			#say "ID $id, Tag $$tag{desc} ", %{$$tag{vr}};
			#$packstring=join '',map {((getVR($_))->{type}||'a').'*'} 
			#keys %{$$tag{vr}};
			#say "Packstring $packstring";
			#$value=unpack($packstring,$dcm->getValue($id),'native');
			#say "Vaule $value";
			$pdl->hdr->{dicom}->{$tag->{desc}}=$value;
		} else { 
		}
		$pdl->hdr->{dicom}->{$id} #=~s/([0-9a-fA-F]{4}),([0-9a-fA-F]{4})/$1_$2/r}
			=$value;
	} # for loop over dicom ids
	return $pdl;
}

sub load_dcm_dir {
	my %dcms; #([]);
	my @pid;
	my $dname=shift;
	my %dims; 
	my $id=shift; # field by which to split
	my $n=0;
	opendir (my $dir, $dname) ||die "cannot open directory!";
	for my $file (readdir ($dir)) {
		next unless (-f "$dname/$file"); # =~m/\.dcm$|\.IMA$/;
		#say "file $file";
		my $p=read_dcm("$dname/$file");
		eval{$p->isa('PDL')} ||next;
		$n++;
		#my $pid=$p->hdr->{ascconv}->{lProtID};
		no PDL::NiceSlice;
		my $pid=$id->($p); 
		use PDL::NiceSlice;
		#say "PID: $pid @pid";
		unless (grep (/$pid/,@pid)) {
			$dims{$pid}=zeroes(short,13);
			push @pid,$pid;
		}
		$dcms{$pid}={} unless ref $dcms{$pid};
		say "pos ",$p->hdr->{dicom}->{'0020,0032'}=~s/\\/ /r;
		say "orientation ",$p->hdr->{dicom}->{'0020,0037'}=~s/\\/ /r;
		say "Spacing ",$p->hdr->{dicom}->{'0028,0030'}=~s/\\/ /r;
#say "IceDims ",$p->hdr->{IceDims};
		say "$n Series $pid IceDims ",$p->hdr->{IceDims};
		my $iced=pdl(short,[split ('_',$p->hdr->{IceDims}=~s/X/1/er)]); #badvalue(short)/er)]);
		#say "$file IceDims ",$iced;
#$iced->badflag(1);
		$dims{$pid}.=$dims{$pid}*($dims{$pid}>=$iced)+$iced*($iced>$dims{$pid});
		$p->hdr->{IcePos}=$iced--;
		$dcms{$pid}->{$p->hdr->{IceDims}}=$p; # if ($p->isa('PDL')); # and $pid == $p->hdr->{ascconv}->{lProtID});
	}
	for my $id (@pid) {
		$dcms{$id}->{dims}=dclone $dims{$id};
	}
	\%dcms;
}

sub parse_dcms {
	my %dcms=%{shift()}; # reference to hash of 
	my %data;
	#my (%tes,);
	for my $pid (keys %dcms) {
		my %stack=%{$dcms{$pid}};
		my $dims =dclone $stack{dims};
		delete $stack{dims};
		my $ref=$stack{(keys %stack)[0]};
		#my $z=$ref->hdr->{ascconv}->{sSliceArray_lSize};
		my $x=$ref->hdr->{dicom}->{Columns} ; 
		die "No $x ",$ref->info unless $x;
		my $y=$ref->hdr->{dicom}->{Rows};
		#my $echoes=$ref->hdr->{ascconv}->{lContrasts};
		#my $t=$ref->hdr->{ascconv}->{lRepetitions}+1;
		#for my $i (0..$echoes-1) {
		#	$tes{$ref->hdr->{ascconv}->{"alTE\_$i\_"}}=$i;
		#}
		#my ($ncoils,@coils);
		#if (defined ($ref->hdr->{dicom}->{"0051,100f"}) 
		#		and $ref->hdr->{dicom}->{"0051,100f"}=~/^C:/ ) {
		#	@coils=map {$ref->hdr->{ascconv}->$_ } 
		#	grep {/asCoilSelectMeas_0__asList_\d+__lElementSelected/}
		#	keys %{$ref->hdr->{ascconv}}; #{ } @{$conv{$pid}};
		#	$ncoils=@coils;
		#} else {
		#	$ncoils=1;
		#}
#my $other=(@stack/$z,$t,$echoes/$ncoils);
#my $data{$pid}=zeroes(ushort,$x,$y,$z,$echoes,$ncoils,$t,$other);
		say "Creating piddle $x,$y, $dims";
		# dims: coil echo phase set t ? partition? slice? ? slice ? some_id
		my $order=pdl[6,7,4,1,0,2,3];
		$data{$pid}=zeroes(ushort,$x,$y,$dims($order));
		$data{$pid}->sethdr(dclone($ref->gethdr)); # populate the header
		$data{$pid}->hdr->{diff}={};
		for my $key (@key_list) {
			$data{$pid}->hdr->{dicom}->{$key}=zeroes(list $dims($order));
			#$data{$pid}->hdr->{dicom}->{$key}.=$ref->hdr->{dicom}->{$key};
			say "$key ",$data{$pid}->hdr->{dicom}->{$key}->info;
		}
		$data{$pid}->hdr->{dicom}->{'Image Orientation (Patient)'}=zeroes(6,list $dims($order));
		$data{$pid}->hdr->{dicom}->{'Image Position (Patient)'}=zeroes(3,list $dims($order));
		$data{$pid}->hdr->{dicom}->{'Pixel Spacing'}=zeroes(2,list $dims($order));
		for my $dcm (values %stack) {
			#say $data{$pid}->info,list( $dcm->hdr->{IcePos}->($order));
			$data{$pid}->(,,list $dcm->hdr->{IcePos}->($order)).=$dcm;
			for my $key (@key_list) {
				#say "setting $key ",$dcm->hdr->{IcePos}->($order) ;
				$data{$pid}->hdr->{dicom}->{$key}->(list $dcm->hdr->{IcePos}->($order))
					.=$dcm->hdr->{dicom}->{$key};
			}
			$data{$pid}->hdr->{dicom}->{'Image Orientation (Patient)'}->(,list $dcm->hdr->{IcePos}->($order))
					.=pdl (split /\\/,$dcm->hdr->{dicom}->{'Image Orientation (Patient)'});
			$data{$pid}->hdr->{dicom}->{'Pixel Spacing'}->(,list $dcm->hdr->{IcePos}->($order))
					.=pdl (split /\\/,$dcm->hdr->{dicom}->{'Pixel Spacing'});
			$data{$pid}->hdr->{dicom}->{'Image Position (Patient)'}->(,list $dcm->hdr->{IcePos}->($order))
					.=pdl (split /\\/,$dcm->hdr->{dicom}->{'Image Position (Patient)'});
			for my $field (keys %{$dcm->hdr->{dicom}}) {
				if ($dcm->hdr->{dicom}->{$field} ne $ref->hdr->{dicom}->{$field}) {
					$data{$pid}->hdr->{diff}->{$field}={}
						unless ref ($data{$pid}->hdr->{diff}->{$field});
					$data{$pid}->hdr->{diff}->{$field}->{$dcm->hdr->{IceDims}}=
						$dcm->hdr->{dicom}->{$field};
				}
			}
		}
		#say "The following keys differ from ref:\n\t",join ", ",sort keys (%{$data{$pid}->hdr->{diff}});
		#say Dumper $data{$pid}->hdr->{diff};
		#say $data{$pid}->info;
		#say $data{$pid}->hdr->{dicom}->{Rows};
		$data{$pid}->hdrcpy(1);
		# serialise partitions/slices and phases/sets
		$data{$pid}=$data{$pid}->clump(2,3)->clump(6,7);
		#say $data{$pid}->info;
		#say $data{$pid}->hdr->{dicom}->{Rows};
	} # for my $pid ...
	\%data;
}



BEGIN {
        if ($_[0] eq q/-d/) {
                require Carp;
                $SIG{__DIE__} = sub {print Carp::longmess(@_); die;};
        }
}
1;

=head1 LICENSE AND COPYRIGHT

Copyright 2016 Albrecht Ingo Schmid.

This program is free software; you can redistribute it and/or modify it
under the terms of the the Artistic License (2.0). You may obtain a
copy of the full license at:

L<http://www.perlfoundation.org/artistic_license_2_0>

Any use, modification, and distribution of the Standard or Modified
Versions is governed by this Artistic License. By using, modifying or
distributing the Package, you accept this license. Do not use, modify,
or distribute the Package, if you do not accept this license.

If your Modified Version has been derived from a Modified Version made
by someone other than you, you are nevertheless required to ensure that
your Modified Version complies with the requirements of this license.

This license does not grant you the right to use any trademark, service
mark, tradename, or logo of the Copyright Holder.

This license includes the non-exclusive, worldwide, free-of-charge
patent license to make, have made, use, offer to sell, sell, import and
otherwise transfer the Package with respect to any patent claims
licensable by the Copyright Holder that are necessarily infringed by the
Package. If you institute patent litigation (including a cross-claim or
counterclaim) against any party alleging that the Package constitutes
direct or contributory patent infringement, then this Artistic License
to you shall terminate on the date that such litigation is filed.

Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER
AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY
YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR
CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


=cut

__END__

0020,000d Study Instance UID
0020,000e Series Instance UID
0020,0010 [Study ID]->
0020,0011 [Series Number]->
0020,0012 [Acquisition Number]->
0020,0013 [Instance Number]-> 
0020,0032 Image Position (Patient) - in mm? 
0020,0037 Image Orientation (Patient) - 
0020,0052 Frame of Reference UID
0020,1040 [Position Reference Indicator]->
0020,1041 [Slice Location]->
