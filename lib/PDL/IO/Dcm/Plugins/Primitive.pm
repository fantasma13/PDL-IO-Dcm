#!/usr/bin/perl

package PDL::IO::Dcm::Plugins::Primitive;
use Exporter;
use PDL;
use strict;
#use 5.10.0;


our @ISA=qw/Exporter/;
our @EXPORT_OK=qw/populate_header setup_dcm/;

sub setup_dcm {
	my $opt=shift;
	$opt={} unless (ref($opt) eq 'HASH'); # ensure hash context
	# split on series number by default
	$$opt{id}=\&PDL::IO::Dcm::sort_series;
	$$opt{sort}=\&populate_header;
	$$opt{delete_raw}=1; # deletes the raw_dicom structure after parsing
	$$opt{Dimensions}=[qw/x y InstanceNumber/];
	$opt;
}

sub populate_header {
	my $dicom =shift;
	my $piddle=shift;
	my $in=$piddle->hdr->{dicom}->{'Instance Number'};
	$piddle->hdr->{dcm_key}=$in; 
	my $pos=pdl(short,$in); 
	$pos--;
	$piddle->hdr->{dim_idx}=$pos;
	return $in;
}

=head1 General

This module provides simple splitting based on intance number and should be used
as template when writing more specific plugin modules.

The setup_dcm creates a template options hash. 

=head1 FUNCTIONS

=head2 populate_header

Here happens the vendor/modallity specific stuff like parsing private fields.
It is required to set the IcePos and dcm_key fields in the piddle header. dcm_key
serves mainly as a unique identifier, IcePos is an index piddle. 

=head2 setup_dcm

sets useful options for this modality. Should accept a hash ref and return one.

=head2 sort_protid

alternative to split based on lProtID (matches raw data key)

=cut

1;
