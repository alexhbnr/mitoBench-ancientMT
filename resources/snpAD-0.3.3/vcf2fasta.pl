#!/usr/bin/perl

use strict ;

# output seq (length fix according to reference used for mapping)
my $seq = "N" x 16569 ;
# for PP field
my %base2num = ( 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3 ) ;

while (<STDIN>) {
	next if ( /^#/ ) ;
	chomp ;
	my ( $chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, $sample ) = split /\t/ ;

	my @names = split /:/, $format ;
	my @vals = split /:/, $sample ;
	my %info ;
	for ( my $i = 0 ; $i < @names ; $i++ ) {
		$info{$names[$i]}=$vals[$i] ;
	}

	my @alleles = ( $ref ) ;
	@alleles = ( $ref, split ",", $alt ) if ( $alt ne "." ) ;
#	print "$chrom $pos $alleles[0] $alleles[1] $info{GT}\n" ;
	
	if ($info{DP} ge 3 & $qual ge 50) {
        if ( $info{GT} eq "0/0" ) { substr( $seq, $pos-1, 1 ) = $alleles[0] ; }
        elsif ( $info{GT} eq "1/1" ) { substr( $seq, $pos-1, 1 ) = $alleles[1] ; }
        elsif ( $info{GT} eq "0/1" ) { 
            my @PP = split ",", $info{PP} ;
            if ( $PP[$base2num{$alleles[0]}] < $PP[$base2num{$alleles[1]}] ) {
                substr( $seq, $pos-1, 1 ) = $alleles[0] ;
            } elsif ( $PP[$base2num{$alleles[0]}] > $PP[$base2num{$alleles[1]}] ){ 
                substr( $seq, $pos-1, 1 ) = $alleles[1] ;
            } # else: stick with N
        } elsif ( $info{GT} eq "1/2" ) {
            my @PP = split ",", $info{PP} ;
            if ( $PP[$base2num{$alleles[1]}] < $PP[$base2num{$alleles[2]}] ) {
                substr( $seq, $pos-1, 1 ) = $alleles[1] ;
            } elsif ( $PP[$base2num{$alleles[1]}] > $PP[$base2num{$alleles[2]}] ){ 
                substr( $seq, $pos-1, 1 ) = $alleles[2] ;
            } # else: stick with N
        }
	}
}

print ">MT\n$seq\n" ;

