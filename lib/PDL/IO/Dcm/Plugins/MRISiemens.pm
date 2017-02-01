#!/usr/bin/perl

package PDL::IO::Dcm::Plugins::MRISiemens;
#use base 'PDL::IO::Dcm';
use Exporter;
use PDL;
use PDL::NiceSlice;
use strict;
use 5.10.0;


our @ISA=qw/Exporter/;
our @EXPORT_OK=qw/populate_header setup_dcm/;

sub setup_dcm {
	my $opt=shift;
	$opt={} unless (ref($opt) eq 'HASH'); # ensure hash context
	# split on series number by default
	$$opt{id}=\&PDL::IO::Dcm::sort_series;
	$$opt{sort}=\&populate_header;
	$$opt{duplicates}=\&handle_duplicates;
	$$opt{delete_raw}=1; # deletes the raw_dicom structure after parsing
	push @PDL::IO::Dcm::key_list,'0051,100f';
	#say join ' ',%{$opt};
	if ($$opt{s_phase_t} ) { 
		$$opt{Dimensions}=[qw/x y z=partitions*slices t*sets*phases echo channel /];
	}else {
		$$opt{Dimensions}=[qw/x y z=partitions*slices t echo channel set phase/];
	}
	# part,sl,t,echo,coil,phase,set
	$$opt{dim_order}=[6,7,4,1,0,2,3];
	$$opt{internal_dims}=[
	#c e p s t ? n l ? i ? ? id
	#2 3 4 5 6 7 8 9 a b c d e
	#
		qw/x y coil echo phase set t ? partition? slice? ? slice ? some_id/];
	# note the order since dims change by clump!
	# partitions and slices
	$$opt{clump_dims}=[[0,1],];
	# phases and set, phases*set and t 
	push (@{$$opt{clump_dims}},[4,5]) if $$opt{s_phase_set};
	push (@{$$opt{clump_dims}},[1,4]) if $$opt{s_phase_t};
	$opt;
}



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

sub sort_protid {
	$_[0]->hdr->{ascconv}->{"lProtID"};
}

sub populate_header {
	my $dicom =shift;
	my $piddle=shift;
	# The protocol is in here:
	#say "populate_header ",$_[1]->info,$_[0]->getValue('0020,0032');
	read_text_hdr($dicom->getValue ('0029,1020','native'),$piddle); 
	delete $piddle->hdr->{raw_dicom}->{'0029,1020'}; # Protocol
	my @ret=$dicom->getValue('0029,1010','native')=~/ICE_Dims.{92}((_?(X|\d+)){13})/s; 
	(my $str=$ret[0])=~s/X/1/e;
	# to make this unique
	$piddle->hdr->{dcm_key}=$dicom->getValue('InstanceNumber').'_'.($dicom->getValue('0051,100f')||0);
	my @d=split ('_',$str);
	my $iced=pdl(short,@d); #badvalue(short)/er)]);
	$iced--;
	$piddle->hdr->{dim_idx}=$iced;
	#say "Dims $str pos $iced";
	return $str;
}

sub handle_duplicates {
	my $stack=shift;
	my $dcm=shift;
	my $opt=shift;
	"This entry (". $dcm->hdr->{dim_idx}->($$opt{dim_order}).
		max ($stack->(,,list $dcm->hdr->{dim_idx}->($$opt{dim_order}))).
		") is already set! This should not happen, please file a bug report!\n";
}

sub init_dims {
	use PDL::NiceSlice;
	my $self=shift;
	my $opt=shift;
	require PDL::Dims || return;
	say diminfo $self,$self->hdr->{ascconv}->{sGroupArray_asGroup_1__nSize};
# center=inner(rot*scale,dim/2)+Image Position (Patient)
# xf=t_linear(matrix=>Pixel Scale*$rot->transpose,post=>Image Position (Patient))
	if (hpar($self,'ascconv','sGroupArray_asGroup_1__nSize')){
		warn "Multiple slice groups, not yet covered!";
	}
	use PDL::NiceSlice;
	my $v=zeroes(3);
	$v(0,).=hpar($self,'ascconv','sSliceArray_asSlice_0__sPosition_dSag') ||0; #x
	$v(1,).=hpar($self,'ascconv','sSliceArray_asSlice_0__sPosition_dCor') ||0; #y
	$v(2,).=hpar($self,'ascconv','sSliceArray_asSlice_0__sPosition_dTra') ||0; #z
	my $ir=hpar($self,'ascconv','sSliceArray_asSlice_0__dInPlaneRot') ||0; #radiant
	my $or=pdl(hpar($self,'dicom','Image Orientation (Patient)'))->(:5,0;-)
		->reshape(3,2)->transpose; #
	my $p=pdl(hpar($self,'dicom','Image Position (Patient)'))->(:2,;-);
	my $pe_dir=hpar($self,'dicom','In-plane Phase Encoding Direction');
# Scaling
	my $sp=zeroes(3);
	$sp(:1;-).=pdl(hpar($self,'dicom','Pixel Spacing'))->(:1,0);
	$sp(2;-).=sqrt sumover(($p(,1;-)-$p(,0;-))**2);
	my $scale=stretcher $sp;
# Rotation
	my $srot=zeroes(3,3);
	$srot(:1,).=$or;
	$srot(2,;-).=norm($p(,1;-)-$p(,0;-));
# Calculate and initialise the transformation
	my $mat=$srot x $scale;
	$mat=$mat->transpose ;#if ($pe_dir =~ /COL/);
	my $xf=t_linear(matrix=>$mat,post=>$p(,0;-));
	my $s=zeroes(3); # matrix size 
	my $fov=zeroes(3); # FOV
	$s(0).=$self->hdr->{dicom}->{Width}||$self->hdr->{dicom}->{Columns};
	$s(1).=$self->hdr->{dicom}->{Height}||$self->hdr->{dicom}->{Rows};
	$fov(0).=$self->hdr->{ascconv}->{sSliceArray_asSlice_0__dReadoutFOV};
	$fov(1).=$self->hdr->{ascconv}->{sSliceArray_asSlice_0__dPhaseFOV};
	$self->hdr->{'3d'}=1 if (($self->hdr->{dicom}->{MRAcquisitionType}||
		hpar($self,'dicom','MR Acquisition Type')) eq '3D'); # 3D
	if ($self->hdr->{'3d'}) {
		$s(2).=$self->hdr->{ascconv}->{'sKSpace_lImagesPerSlab'} ;
		$fov(2).=$self->hdr->{ascconv}->{sSliceArray_asSlice_0__dThickness};
	} else {
		$s(2).=$self->hdr->{ascconv}->{'sSliceArray_lSize'};
		$fov(2).=$self->hdr->{dicom}->{SpacingBetweenSlices}*$s(2);
	}
	my $rot=identity($self->ndims);
	my $inc_d=zeroes(3);
	#say "Pixel Spacing", hpar($self,'dicom','Pixel Spacing');
	$inc_d(:1).=hpar($self,'dicom','Pixel Spacing')->(:1;-);
	$inc_d(2).=$fov(2,0)/$s(2,0);
	my $pos_d=hpar($self,'dicom','Image Position (Patient)')->(:2,0;-);
	$rot(:2,:2).=$srot;
	say "Rot: $rot";
	initdim($self,'x',size=>$s(0),min=>sclr$pos_d(0),inc=>sclr$inc_d(0),unit=>'mm');
	initdim($self,'y',size=>$s(1),min=>sclr$pos_d(1),inc=>sclr$inc_d(1),unit=>'mm');
	initdim($self,'z',size=>$s(2),rot=>$rot,min=>sclr$pos_d(2),inc=>sclr$inc_d(2),unit=>'mm',);
#my $xf=t_linear(matrix=>$rot,pre=>$pos_d,scale=>$inc,post=>$post,dims=>3);
# this seems to wowrk. center as in slice ... position is apply(dim/2)!
	say "size $s min $pos_d inc $inc_d rot $rot";
	$self->hdr->{transform}=$xf;
	$self->hdr->{rotation}=$srot;
	$self->hdr->{matrix}=$mat;
	$self->hdr->{offset}=$p(,0;-);
	say $xf,diminfo $self,$rot;
	my @ors=qw/Cor Sag Tra/;
	$self->hdr->{orientation}=$ors[maximum_ind($rot(2;-))];  # normal vector lookup
	my $pe=$self->hdr->{dicom}->{"In-plane Phase Encoding Direction"} 
	||$self->hdr->{dicom}->{"InPlanePhaseEncodingDirection"}; 
	if ($self->hdr->{orientation} eq 'Tra'){ # transversal slice
		say $self->hdr->{orientation}. " Orientation";
		$self->hdr->{sl}='z';
		if ($pe eq 'COL'){
			$self->hdr->{ro}='x'; 
			$self->hdr->{pe}='y';
		} else { 
			$self->hdr->{ro}='y'; 
			$self->hdr->{pe}='x';
		}
	}
	if ($self->hdr->{orientation} eq 'Cor'){ # coronal slice
		$self->hdr->{sl}='y';
		if ($pe eq 'COL') {
			$self->hdr->{ro}='x'; 
			$self->hdr->{pe}='z';
		} else { 
			$self->hdr->{ro}='z'; 
			$self->hdr->{pe}='x';
		}
	}
	if ($self->hdr->{orientation} eq 'Sag'){ # sagittal slice
		$self->hdr->{sl}='x';
		if ($pe eq 'COL') {
			$self->hdr->{ro}='y'; 
			$self->hdr->{pe}='z';
		} else { 
			$self->hdr->{ro}='z'; 
			$self->hdr->{pe}='y';
		}
	}
	idx($self,'x',dimsize($self,'x')/2);
	idx($self,'y',dimsize($self,'y')/2);
	idx($self,'z',dimsize($self,'z')/2);
	say "orientation : ",hpar($self,'orientation'),diminfo $self;
	# other dimensions
	for my $n (3..$#{$$opt{Dimensionsa}}) { # x,y,z are handled above
		my $dim=$$opt{Dimensions}->[$n];
		if ($dim eq 'echo') {
			initdim ($self,'echo',vals=>hpar($self,'dicom','Echo Time'));
		}
		elsif ($dim eq 't') {
			my $str=('(0),' x ($n-2)).','.('(0),' x ($#{$$opt{Dimensions}}-$n));
			initdim ($self,'t',vals=>hpar($self,'dicom','Achisition Time')->($str));
		} elsif ($dim eq 'channel') {
			my $coil=hpar($self,'dicom','0051,100f');
			if ($self->dim($n)>1) {
				initdim ($self,'t',vals=>hpar($self,'dicom','0051,100f')->flat->(:2));
			} else {
				initdim ($self,'channel',vals=>$coil, size=>1);
			}
		} elsif ($dim eq 'set') {
			initdim ($self,'Set'); # This can be anything, no further info easily available
		} elsif ($dim eq 'phase') {
			my $str=('(0),' x ($n-2)).','.('(0),' x ($#{$$opt{Dimensions}}-$n));
			initdim ($self,'phase',vals=>hpar($self,'dicom','Trigger Time')->($str));
		}
	}
#barf "initdim fails!" unless ($#{dimname($self)}>2);
}


=head1 Specific handling of Simenes MRI data

Key 0029,1010 is the Siemens specific field that contains the ICE
miniheaders with dimension information - and position in matrix
0029,1020 is deleted from the header, it is big, containing the whole
protocol. The important part, the Siemens protocol ASCCONV part, is stored in
the ascconv key, see read_text_hdr.


=head1 FUNCTIONS

=head2 handle_duplicates

What to do if two images with the same position in the stack arrive. Throws an 
error, atm. Should handle duplicate exports

=head2 init_dims

provides support for PDL::Dims. Useful in combination with PDL::IO::Sereal to
have a fully qualified data set.

=head2 populate_header

Here happens the vendor/modallity specific stuff like parsing private fields.
It is required to return a position vector in the series' piddle.

=head2 read_text_hdr

parses the ASCCONV part of Siemens data header into the ascconv field of the
piddle header. All special characters except [a-z0-9]i are converted to _ -- no
quoting of hash keys required! You don't need to load this yourself.

=head2 setup_dcm

sets useful options for this modality. 

=head2 sort_protid

alternative to split based on lProtID (matches raw data key). To activate,
after running setup_dcm, point option id to \&sort_protid.

=head1 Specific options

=over 

=item Nifti

Do we want Nifti output? May be used by your plugin to apply additional steps,
eg. more clumps, reorders, setting header fields ...

=item s_phase_t

Serialize phase and t dimensions

=back

=cut

1;
