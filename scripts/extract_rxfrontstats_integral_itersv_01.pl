#! /usr/bin/perl -w

use strict;

# extract descriptive statistics

#my $fname = shift; $fname = &trim($fname);
my $jobid = shift; $jobid = &trim($jobid);
my $si = shift; $si = &trim($si);
my $sf = shift; $sf = &trim($sf);

#assuming valid files on ID interval [$si, $sf] with increment +1

print("Extract Integral Sv(Iter) dIter\n");
my $afname = "";
my $nlines = 0;
my @data = ();
my @tmp = ();
my @rxf = ();
my $sout = "Filename;Int(Sv(Iter)dIter)\n";
my $integralvalue = 0.0;

for ( my $f = $si; $f <= $sf; $f++ ) {
	$afname = "SCORE.$jobid.MyRank.$f.JobID.$f.RXFrontStats.csv";
	
	if (!-e $afname) {print("\t\t$afname does not exist\n"); exit;}
	
	open ( RXFRONT, "<$afname") || die ( "ERROR.open.file=($afname).mode=(r)");
	@rxf = <RXFRONT>;
	close ( RXFRONT );
	
	shift(@rxf);
	$nlines = @rxf;
	$integralvalue = 0.0;
	
	#integrate
	for (my $i = 0; $i < $nlines; $i++) {
		#$data =~ s/,/;/g;
		@data = split(/;/, &trim($rxf[$i]) );
		
		$integralvalue += $data[3];
	}
	
	print ( "$afname\t\t$integralvalue\n");
	$sout .= "$afname;$integralvalue\n";
}

open ( RESULTS, ">SCORE.IntegralSvIter.$jobid.csv") || die ( "ERROR.open.file=(SCORE.IntegralSvIter.$jobid.csv).mode=(r)");
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

