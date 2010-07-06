#!/bin/perl
#application of functions for particular file. Add to and change as needed

use Getopt::Long;

my $opt_ifile = "";
GetOptions("input=s" => \$opt_ifile);

open IFILE, $opt_ifile;

%flist = {};
$prev_hour = "";
$prev_date = "";

#skip header
if(<IFIlE>) {;}

for($i = 0; <IFILE>; $i++) {
    $_ =~ m/^([0-9\/]+)\t([0-9:]+)\t([0-9.]+)\t([0-9.]+)/;
    $date = $1;
    $hour = $2;
    $x = $3;
    $y = $4;

    #check dup
    if( ($hour eq $prev_hour) && ($date eq $prev_date)) {
	print "DUP FOUND at line " . $i . ":";
	print $_ . "\n";
    }

    #set prevs for next iter dup check
    $prev_hour = $hour;
    $prev_date = $date;

    #print $date . "\t" . $hour . "\t" . $x . "\t" . $y . "\n";

}


sub lendian_to_mendian($date) {
    $day = $date;
    $day =~ m/^([0-9]+)\//;
    $day = $1;

    $month = $date;
    $month =~ m/^[0-9]+\/([0-9]+)\//;
    $month = $1;
    #if leading 0, get rid of it
    if($month =~ m/^0[1-9]+/) {
	$month =~ /^0([1-9]+)/;
	$month = $1;
    }

    $year = $date;
    $year =~ m/^[0-9]+\/[0-9]+\/([0-9]+)/;
    $year = $1;
    
    $mendian = $month . '/' . $day . '/' . $year;

    return($mendian);    
}

sub int_to_hour($hour) {
    return($hour . ":00");
}
