#!/usr/bin/perl
#
package PDL::IO::Dcm;

=head1 NAME

PDL::IO::Dcm - Reads dicom files, sorts them and stores the result into piddles with headers 

=head1 VERSION

Version 0.9

=cut

our $VERSION = '0.9002';


=head1 SYNOPSIS


This is inteded to read and sort dicom images created by medical imaging devices.

Either use something like the following from within your module/application

	# loads all dicom files in this directory
	# $id is a code reference returning e.g. the Series Number
	my $dcms=load_dcm_dir($dir,$id); 
	die "no data!" unless (keys %$dcms);
	print "Read data; ProtIDs: ",join ', ',keys %$dcms,"\n";
	# sort all individual dicoms into a hash of piddles.
	my $data=parse_dcms($dcms);

or use the read_dcm.pl script to convert dicom files in a directory to serealised
piddles (PDL::IO::Sereal) or NIFTI files with separate text headers (PDL::IO::Nifti).

This software is tailored to Siemens MRI data based on the author's needs. For 
general usage, the read_dcm function should and probably will be moved to
vendor/modality specific plugin modules in future releases.

=head1 Some notes on Dicom fields and how they are stored/treated

The image data field is stored as the piddle, the other dicom
elements are stored in the header under the raw_dicom key. 

The Siemens protocol ASCCONV part is stored in the ascconv key. 

Key 0029,1010 is the Siemens specific field that contains the ICE
miniheaders with dimension information - and position in matrix
0029,1020 is deleted from the header, it is big, containing the whole
protocol. The important part is parsed into ascconv.

Keys are parsed into a hash under the dicom key using the DicomPack module(s)
to unpack. 

The header fields IceDims and IcePos are used for sorting datasets. The field
Dimensions lists the sorted names of dimensions. 

Piddles are created for each lProtID or Series Number value. 


=head1 SUBROUTINES/METHODS

=head2 is_equal ($dcm1,$dcm2,$pattern)

This is used to check if two dicoms can be stacked based on matrix size, orientation and pixel spacing.

If $pattern matches /d/, only dims are checked

=head2 read_text_hdr

parses the ASCCONV part of Siemens data header into the ascconv field of the
piddle header. All special characters except [a-z0-9]i are converted to _ -- no
quoting of hash keys required!

=head2 load_dcm_dir ( $dir,\%options)

reads all dicom files in a dicrectory and returns a hash of piddles containing
sorted N-D data sets. 

Fields in the options hash include:

=over 

=item id: 

Uses a code reference to access the field by which to split. See read_dcm.pl for details. Currently, the ascconv field lProtID or dicom Series Number are used as keys. 

=item sp: 

Split slice groups, otherwise they are stacked together if xy-dims match, even transposed.

=back 


=head2 parse_dcms 

Parses and sorts a hash of hashes of dicoms (such as returned by load_dcm_dir)
based on lProtID and the ICE_Dims field in 0029_1010. Returns a hash of piddles 
(lProtID).

=head2 unpack_field

unpacks dicom fields and walks subfield structures recursively.


=head2 read_dcm

reads a dicom file and creates a piddle-with-header structure.

=head1 TODO 

write tests! 

Generalise to other modalities. This will be done based on request or as needed.

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
#use 5.10.0;

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
	#say "datatype $datatype";
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
sub is_equal {
	my $a=shift;
	my $b=shift;
	my $opt=shift;
	#say "dims ", any $a->shape-$b->shape,', ', # they have equal dimensions
		#"spacing ",$a->hdr->{dicom}->{'Pixel Spacing'} ,' ne ', $b->hdr->{dicom}->{'Pixel Spacing'},', ',
		#"orientation ",$a->hdr->{dicom}->{'Image Orientation (Patient)'} ,' ne ', $b->hdr->{dicom}->{'Image Orientation (Patient)'};
	return if (any ($a->shape-$b->shape)); # they have equal dimensions
	#say "is_equal: dims ok";
	return 1 if ($opt =~/d/);
	return if $a->hdr->{dicom}->{'Pixel Spacing'} ne $b->hdr->{dicom}->{'Pixel Spacing'};
	#say "is_equal: spacing ok";
	return if $a->hdr->{dicom}->{'Image Orientation (Patient)'} ne $b->hdr->{dicom}->{'Image Orientation (Patient)'};
	#say "is_equal: orientation ok";
	1;
}
sub load_dcm_dir {
	my %dcms; #([]);
	my @pid;
	my $dname=shift;
	my %dims; 
	my $opt=shift; # field by which to split
	my $id=$$opt{id};
	my $sp=$$opt{sp};
	my $n=0;
	my %refs; # reference images for each stack
	opendir (my $dir, $dname) ||die "cannot open directory!";
	for my $file (readdir ($dir)) {
		next unless (-f "$dname/$file"); # =~m/\.dcm$|\.IMA$/;
		#say "file $file";
		my $p=read_dcm("$dname/$file");
		eval{$p->isa('PDL')} ||next;
		$n++;
		#my $pid=$p->hdr->{ascconv}->{lProtID};
		no PDL::NiceSlice;
		my $pid=$id->($p); # Call to subroutine reference 
		#say "PID: $pid @pid";
		#say $p->info;
		#say "ID $pid Instance number ",$p->hdr->{dicom}->{'Instance Number'};
		$dcms{$pid}={} unless ref $dcms{$pid};
		my $ref =$refs{$pid}; 
		#say "ref ID: ",$id->($ref)," ref Instance number ",$ref->hdr->{dicom}->{'Instance Number'},
			#$ref->info,$p->info, "diff? ",is_equal($ref,$p);
		if (defined $ref) {
		unless ( is_equal($ref,$p )) {
			if ( !$sp and is_equal($ref,$p->transpose,'d')) {
				#say "Split? $sp; ", ( !$sp and is_equal($ref,$p->transpose,'d'));
				$p->hdr->{tp}=1;
			} else {
				#say "Difference! $sp";
				my $flag=0;
				my $n='a';
				my $nid;
				do {
					$nid=$id->($p).$n;
					if (ref $dcms{$nid} eq 'HASH'){ # group
						for my $r2 (values %{$dcms{$nid}}){
							$flag=is_equal($r2,$p);
					#		say "r2 vs. p id $nid ",$r2->info,$p->info;
							last unless $flag;
						}
					} else {
						$dcms{$nid}={};
						$pid=$nid;
						$flag=1;
					}
					$n++;
				} until $flag;
				$pid=$nid;
			}
		}
		} # defined $ref
		use PDL::NiceSlice;
		unless (grep (/^$pid$/,@pid)) {
			#say "PID: $pid";
			$dims{$pid}=zeroes(short,13);
			push @pid,$pid;
			$refs{$pid}=$p;
		}
		#say "pos ",$p->hdr->{dicom}->{'0020,0032'}=~s/\\/ /r;
		#say "orientation ",$p->hdr->{dicom}->{'0020,0037'}=~s/\\/ /r;
		#say "Spacing ",$p->hdr->{dicom}->{'0028,0030'}=~s/\\/ /r;
#say "IceDims ",$p->hdr->{IceDims};
		#say "$n Series $pid IceDims ",$p->hdr->{IceDims};
		(my $str=$p->hdr->{IceDims})=~s/X/1/e;
		my @d=split ('_',$str);
		my $iced=pdl(short,@d); #badvalue(short)/er)]);
		#say "$file IceDims ",$iced;
#$iced->badflag(1);
		$dims{$pid}.=$dims{$pid}*($dims{$pid}>=$iced)+$iced*($iced>$dims{$pid});
		$p->hdr->{IcePos}=$iced--;
		$dcms{$pid}->{$p->hdr->{IceDims}}=$p; # if ($p->isa('PDL')); # and $pid == $p->hdr->{ascconv}->{lProtID});
	}
	for my $id (@pid) {
		$dcms{$id}->{dims}=$dims{$id}->copy;
		print "Set dims: id $id, $dims{$id}\n";
	}
	\%dcms;
}



sub parse_dcms {
	my %dcms=%{shift()}; # reference to hash of 
	my %data;
	#my (%tes,);
	for my $pid (keys %dcms) {
		my %stack=%{$dcms{$pid}};
		#next unless (ref $stack{dims} eq 'HASH');
		#say keys %stack;
		my $dims =$stack{dims};
		die "No dims $pid " unless eval {$dims->isa('PDL')};
		delete $stack{dims};
		my $ref=$stack{(keys %stack)[0]};
		#my $z=$ref->hdr->{ascconv}->{sSliceArray_lSize};
		my $x=$ref->hdr->{dicom}->{Columns} ; 
		die "No $x ",$ref->info unless $x;
		my $y=$ref->hdr->{dicom}->{Rows};
		#print "ID: $pid dims $dims transpose? ",$ref->hdr->{tp},"\n";
		# dims: coil echo phase set t ? partition? slice? ? slice ? some_id
		my $order=pdl[6,7,4,1,0,2,3];
		if ($ref->hdr->{tp}) { $data{$pid}=zeroes(ushort,$y,$x,$dims($order));}
		else { $data{$pid}=zeroes(ushort,$x,$y,$dims($order));}
		$data{$pid}->sethdr(dclone($ref->gethdr)); # populate the header
		$data{$pid}->hdr->{diff}={};
		$data{$pid}->hdr->{Dimensions}=[qw/x y z t echo channel set/];
		for my $key (@key_list) {
			$data{$pid}->hdr->{dicom}->{$key}=zeroes(list $dims($order));
			#$data{$pid}->hdr->{dicom}->{$key}.=$ref->hdr->{dicom}->{$key};
			#say "$key ",$data{$pid}->hdr->{dicom}->{$key}->info;
		}
		$data{$pid}->hdr->{dicom}->{'Image Orientation (Patient)'}=zeroes(6,list $dims($order));
		$data{$pid}->hdr->{dicom}->{'Image Position (Patient)'}=zeroes(3,list $dims($order));
		$data{$pid}->hdr->{dicom}->{'Pixel Spacing'}=zeroes(2,list $dims($order));
		for my $dcm (values %stack) {
			#say $data{$pid}->info,list( $dcm->hdr->{IcePos}->($order));
			#say "$x $y ",$ref->info;
			#say "ID $pid: Transpose? ",$dcm->hdr->{tp},$dcm->info;
			if ($dcm->hdr->{tp}) {
				$data{$pid}->(,,list $dcm->hdr->{IcePos}->($order)).=$dcm->transpose;}
			else {$data{$pid}->(,,list $dcm->hdr->{IcePos}->($order)).=$dcm;}
			for my $key (@key_list) {
				#say "setting $key ",$dcm->hdr->{IcePos}->($order) ;
				$data{$pid}->hdr->{dicom}->{$key}->(list $dcm->hdr->{IcePos}->($order))
					.=$dcm->hdr->{dicom}->{$key};
			}
			$data{$pid}->hdr->{dicom}->{'Image Orientation (Patient)'}
				->(,list $dcm->hdr->{IcePos}->($order))
				.=pdl (split /\\/,$dcm->hdr->{dicom}->{'Image Orientation (Patient)'});
			$data{$pid}->hdr->{dicom}->{'Pixel Spacing'}
				->(,list $dcm->hdr->{IcePos}->($order))
				.=pdl (split /\\/,$dcm->hdr->{dicom}->{'Pixel Spacing'});
			$data{$pid}->hdr->{dicom}->{'Image Position (Patient)'}
				->(,list $dcm->hdr->{IcePos}->($order))
				.=pdl (split /\\/,$dcm->hdr->{dicom}->{'Image Position (Patient)'});
			for my $field (keys %{$dcm->hdr->{dicom}}) {
				if ($dcm->hdr->{dicom}->{$field} ne $ref->hdr->{dicom}->{$field}) {
					$data{$pid}->hdr->{diff}->{$field}={}
						unless ref ($data{$pid}->hdr->{diff}->{$field});
					$data{$pid}->hdr->{diff}->{$field}->{$dcm->hdr->{IceDims}}=
						$dcm->hdr->{dicom}->{$field};
				}
			}
		} # for ... values %stack
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
