#!/usr/bin/perl

package PDL::IO::Dcm::Plugins::MRISiemens;
#use base 'PDL::IO::Dcm';
use Exporter;
#use PDL::Lite;
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
	say join ' ',%{$opt};
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
	# dicom, piddle
	# The protocol is in here:
	#say "populate_header ",$_[1]->info,$_[0]->getValue('0020,0032');
	read_text_hdr($_[0]->getValue ('0029,1020','native'),$_[1]); 
	delete $_[1]->hdr->{raw_dicom}->{'0029,1020'}; # Protocol
	my @ret=$_[0]->getValue('0029,1010','native')=~/ICE_Dims.{92}((_?(X|\d+)){13})/s; 
	say "Ret: @ret";
	return shift @ret;
}


=head1 FUNCTIONS

=head2 read_text_hdr

parses the ASCCONV part of Siemens data header into the ascconv field of the
piddle header. All special characters except [a-z0-9]i are converted to _ -- no
quoting of hash keys required! You don't need to load this yourself.

=head2 populate_header

here happens the vendor/modallity specific stuff like parsing private fields

=head2 setup_dcm

sets useful options for this modality. 

=head2 sort_protid

alternative to split based on lProtID (matches raw data key)

=cut

1;
