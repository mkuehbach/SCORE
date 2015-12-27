#! /usr/bin/perl -w

use strict;

# extract descriptive statistics

my $fname = shift; $fname = &trim($fname);

#check input
if (!-e $fname) {print("$fname does not exist\n"); exit;}

#open log.xx.csv
open ( LOGCSV, "<$fname") || die ( "ERROR.open.file=($fname).mode=(r)");
my @alogcsv = <LOGCSV>;
close ( LOGCSV );

print("Extract descriptive statistics\n");
my $nlines = @alogcsv;
my $afname = "";
my $data = "";
my @data = ();
my @tmp = ();
my @gsd = ();
my $sout = "Filename;Time;NGrains;Vref;Mean;Median\n";

for ( my $line = 0; $line < $nlines; $line++ ) {
	$afname = &trim($alogcsv[$line]);
	
	if (!-e $afname) {print("\t\t$afname does not exist\n"); exit;}
	
	open ( GSDCSV, "<$afname") || die ( "ERROR.open.file=($afname).mode=(r)");
	@gsd = <GSDCSV>;
	close ( GSDCSV );
	
	$data = &trim($gsd[0]);
	$data =~ s/,/;/g;
	@data = split(/;/,$data);
	
	$sout .= "$afname;";
	
	@tmp = split(/\=/,$data[0]); #time
	$sout .= substr($tmp[1], 0, -1).";";
	
	@tmp = split(/\=/,$data[1]); #ngrains
	$sout .= "$tmp[1];";

	@tmp = split(/\=/,$data[2]); #vref
	$sout .= "$tmp[1];";

	@tmp = split(/\=/,$data[3]); #mean
	$sout .= "$tmp[1];";

	@tmp = split(/\=/,$data[4]); #median
	$sout .= "$tmp[1];";
	
	$sout .= "\n";
}

open ( RESULTS, ">$fname.DescrStats.csv") || die ( "ERROR.open.file=($fname.DescrStats.csv).mode=(r)");
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

