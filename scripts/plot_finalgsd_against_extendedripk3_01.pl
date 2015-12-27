#! /usr/bin/perl -w

use strict;

# takes the final grain size distribution obtained from a SCORE simulation as CDF = f((Volume/4pi/3)^(1/3)) 
#and plots it against two difference extended ripk3 defined as ripk3(extract descriptive statistics

my $jobid = shift; 			$jobid = &trim($jobid);
#my $fn_gsd = shift; 		$fn_gsd = &trim($fn_gsd);
my $fn_ripk3_1 = shift; 	$fn_ripk3_1 = &trim($fn_ripk3_1);
my $fn_ripk3_64 = shift; 	$fn_ripk3_64 = &trim($fn_ripk3_64);
my $fn_ripk3_512 = shift;	$fn_ripk3_512 = &trim($fn_ripk3_512);
my $nbor = shift; 			$nbor = &trim($nbor);

#if (!-e $fn_gsd) {print("\t\t$fn_gsd does not exist\n"); exit;}
if (!-e $fn_ripk3_1) {print("\t\t$fn_ripk3_1 does not exist\n"); exit;}
if (!-e $fn_ripk3_64) {print("\t\t$fn_ripk3_64 does not exist\n"); exit;}
if (!-e $fn_ripk3_512) {print("\t\t$fn_ripk3_512 does not exist\n"); exit;}
if ( $nbor <= 0 ) {print("\t\tInvalid value for nbor!"); exit; }

#import data
#open ( GSD, "<$fn_gsd") || die ( "ERROR.open.file=($fn_gsd).mode=(r)");
#my @fgsd = <GSD>;
#close ( GSD );

#import RipK3
open ( RIPK3, "<$fn_ripk3_1") || die ( "ERROR.open.file=($fn_ripk3_1).mode=(r)");
my @r1 = <RIPK3>;
close ( RIPK3 );

open ( RIPK3, "<$fn_ripk3_64") || die ( "ERROR.open.file=($fn_ripk3_64).mode=(r)");
my @r64 = <RIPK3>;
close ( RIPK3 );

open ( RIPK3, "<$fn_ripk3_512") || die ( "ERROR.open.file=($fn_ripk3_512).mode=(r)");
my @r512 = <RIPK3>;
close ( RIPK3 );

my $sout = "Vol;r;CDF;Diff64;Diff512\n";

#how many neighbors?
my @data = split(/;/, &trim($r1[0]));
my $nr = @data - 1;
print("Some many r $nr\n");

my $rip1 = 0.0;
my $rip64 = 0.0;
my $rip512 = 0.0;
my $v64 = 0.0;
my $v512 = 0.0;
my $rvalue = 0.0;

for ( my $r = 0; $r < $nr; $r++ ) {

	@data = split(/;/, &trim($r1[0]));
	$rvalue = $data[1+$r];

	#extract the desired line with exactly N neighbors in the box

	$rip1 = 0.0;
	$rip64 = 0.0;
	$rip512 = 0.0;
	
	#HACK!
	$nbor = 100;
	for ( my $i = 0; $i < $nbor; $i++ ) {
		#1
		@data = split(/;/, &trim($r1[1+$i]));
		$rip1 += ($i * $data[1+$r]);
		
		#64
		@data = split(/;/, &trim($r64[1+$i]));
		$rip64 += ($i * $data[1+$r]);
		
		#512
		@data = split(/;/, &trim($r512[1+$i]));
		$rip512 += ($i * $data[1+$r]);
	}

	
	$v64 = -1.0;
	$v512 = -1.0;
	if ( $rip1 > 0.0 ) 	{
		$v64 = (abs($rip64 - $rip1)); #/$rip1;
		$v512 = (abs($rip512 - $rip1)); #/$rip1;
	}
	print ("$r\t\t$rip1;$rip64;$rip512\t\t$v64\t\t$v512\n");
	
	$sout .= ";$rvalue;;$v64;$v512\n";	
}


open ( RESULTS, ">ExtendedRipK3.$jobid.csv") || die ( "ERROR.open.file=(ExtendedRipK3.$jobid.csv).mode=(r)");
print RESULTS $sout;
close ( RESULTS );
$sout = ""; 


###################################################################################################
sub trim {
	my $str = shift;
	$str =~ s/^\s+//;	
	$str =~ s/\s+$//;	
	return ( $str );
}

