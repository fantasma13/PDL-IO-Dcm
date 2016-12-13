#!perl 
use 5.006;
use strict;
use warnings FATAL => 'all';
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'PDL::IO::Dcm' ) || print "Bail out!\n";
}

diag( "Testing PDL::IO::Dcm $PDL::IO::Dcm::VERSION, Perl $], $^X" );
