#!/usr/bin/perl
#
package PDL::IO::Dcm;

=head1 NAME

PDL::IO::Dcm - Reads Siemens dicom files, sorts them and stores the result into piddles with meta-data

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS


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

=head2 read_dcm

reads a dicom file and creates a piddle-with-header structure.

=cut

use PDL;
use PDL::NiceSlice;
use List::MoreUtils; # qw{any};
#use Data::Dumper;
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


sub read_dcm {
	my $file=shift;
	#say "file $file";
	my $dcm=DicomPack::IO::DicomReader->new($file);
	my $h=unpack('s',$dcm->getValue('Rows'));
	my $w=unpack('s',$dcm->getValue('Columns'));
	my $pdl=zeroes(ushort,$w,$h);
	$pdl->make_physical;
	${$pdl->get_dataref}=$dcm->getValue('PixelData');
	$pdl->upd_data;
	read_text_hdr($dcm->getValue ('0029,1020'),$pdl); # The protocol is in here
	$pdl->hdr->{raw_dicom}=$dcm->getDicomField;
	delete $pdl->hdr->{raw_dicom}->{'0029,1020'}; # Protocol
	#delete $pdl->hdr->{raw_dicom}->{'0029,1010'}; # also big structure, not sure what
		# Ice Dimension string: 13 numbers or X=undef
	#say "Ice Dims ",$dcm->getValue('0029,1010')=~/ICE_Dims/;
	my ($dims)=$dcm->getValue('0029,1010')=~/((_?(X|\d+)){13})/; 
	#say "Raw Dims $dims";
	$pdl->hdr->{IceDims}=$dims || die "No Ice Dims ",$pdl->hdr->{raw_dicom}->{'0029,1010'}; #[split '_',$dims{$pid}=~s/X/0/r];
	delete $pdl->hdr->{raw_dicom}->{'7fe0,0010'}; # Pixel data
	for my $id (keys %{$pdl->hdr->{raw_dicom}}) {
		my $tag=getTag($id);
		my $packstring;
		my $value;
		if (defined $tag) {
			#say "ID $id, Tag $$tag{desc} ", %{$$tag{vr}};
			$packstring=join '',map {((getVR($_))->{type}||'a').'*'} 
			keys %{$$tag{vr}};
			#say "Packstring $packstring";
			$value=unpack($packstring,$dcm->getValue($id));
			#say "Vaule $value";
			$pdl->hdr->{dicom}->{$tag->{desc}}=$value;
		} else { 
			$value=$dcm->getValue($id); 
		}
		$pdl->hdr->{dicom}->{$id=~s/([0-9a-fA-F]{4}),([0-9a-fA-F]{4})/Private_$1_$2/r}
			=$value;
	}
	return $pdl;
}

sub load_dcm_dir {
	my %dcms; #([]);
	my @pid;
	my $dname=shift;
	my %dims; 
	my $n=0;
	opendir (my $dir, $dname) ||die "cannot open directory!";
	for my $file (readdir ($dir)) {
		next unless $file =~m/\.dcm$|\.IMA$/;
		defined (my $p=read_dcm("$dname/$file")) ||die "Could not read $dname/$file!";
		$n++;
		my $pid=$p->hdr->{ascconv}->{lProtID};
#say "PID: $pid @pid";
		unless (grep (/$pid/,@pid)) {
			$dims{$pid}=zeroes(short,13);
			push @pid,$pid;
		}
		$dcms{$pid}={} unless ref $dcms{$pid};
#say "pos ",$p->hdr->{raw_dicom}->{'0020,0032'};
#say "IceDims ",$p->hdr->{IceDims};
		say "$n Pid $pid IceDims ",$p->hdr->{IceDims};
		my $iced=pdl(short,[split ('_',$p->hdr->{IceDims}=~s/X/1/er)]); #badvalue(short)/er)]);
#say "IceDims ",$iced;
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
		$data{$pid}=zeroes(ushort,$x,$y,$dims(:8));
		say $data{$pid}->info;
		$data{$pid}->sethdr(dclone($ref->gethdr)); # populate the header
		for my $dcm (values %stack) {
			$data{$pid}->(,,list $dcm->hdr->{IcePos}->(:8)).=$dcm;
		}
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
