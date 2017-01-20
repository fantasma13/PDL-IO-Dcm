#!/usr/bin/perl

package PDL::IO::Dcm::Plugins::MRISiemens;
#use base 'PDL::IO::Dcm';
use Exporter;
use PDL;
use strict;
use 5.10.0;


our @ISA=qw/Exporter/;
our @EXPORT_OK=qw/populate_header setup_dcm/;

sub setup_dcm {
	my $opt=shift;
	$opt={} unless (ref($opt) eq 'HASH'); # ensure hash context
	# split on series number by default
	$$opt{id}=\&PDL::IO::Dcm::sort_series;
	$$opt{dims}=\&populate_header;
	$$opt{delete_raw}=1; # deletes the raw_dicom structure after parsing
	#say join ' ',%{$opt};
	$$opt{Dimension}=[qw/x y z t echo channel set/];
	$$opt{dim_order}=[6,7,4,1,0,2,3];
	$$opt{internal_dims}=[
	#c e p s t ? n l ? i ? ? id
	#2 3 4 5 6 7 8 9 a b c d e
	#
		qw/x y coil echo phase set t ? partition? slice? ? slice ? some_id/];
	# note the order since dims change by clump!
	$$opt{Clump_dims}=[[0,1],[4,5]];
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
	$piddle->hdr->{dcm_key}=$str;
	my @d=split ('_',$str);
	my $iced=pdl(short,@d); #badvalue(short)/er)]);
	$iced--;
	$piddle->hdr->{dim_idx}=$iced;
	#say "Dims $str pos $iced";
	return $str;
}

=head1 Specific handling of Simenes MRI data

Key 0029,1010 is the Siemens specific field that contains the ICE
miniheaders with dimension information - and position in matrix
0029,1020 is deleted from the header, it is big, containing the whole
protocol. The important part, the Siemens protocol ASCCONV part, is stored in
the ascconv key, see read_text_hdr.


=head1 FUNCTIONS

=head2 read_text_hdr

parses the ASCCONV part of Siemens data header into the ascconv field of the
piddle header. All special characters except [a-z0-9]i are converted to _ -- no
quoting of hash keys required! You don't need to load this yourself.

=head2 populate_header

Here happens the vendor/modallity specific stuff like parsing private fields.
It is required to return a position vector in the series' piddle.

=head2 setup_dcm

sets useful options for this modality. 

=head2 sort_protid

alternative to split based on lProtID (matches raw data key)

=cut

1;
