#!/bin/perl
#application of functions for particular file. Add to and change as needed

open IFILE, "../data/Colby_test_data_all_loc8-cluster.txt";

%flist = {};

#skip header
if(<IFIlE>) {;}

for($i = 0; <IFILE>; $i++) {
    $_ =~ m/^([a-zA-Z0-9]+)\t([0-9]+)\t([0-9\/]+)\t([0-9]+)\t([0-9.]+)\t([0-9.]+)/;
    $goat = $1;
    $fixnum = $2;
    $date = $3;
    $date = lendian_to_mendian($date);
    $hour = $4;
    $hour = int_to_hour($hour);
    $x = $5;
    $x = sprintf("%.3f", $x);
    $y = $6;
    $y = sprintf("%.2f", $y);

    $file = "../data/" . $goat . ".dat";
    #little redundant, but nicer than arrray dups
    $flist{$file} = $file;
    open OFILE, ">>" . $file;

    #print $date . "\t" . $hour . "\t" . $x . "\t" . $y . "\n";
    print OFILE $date . "\t" . $hour . "\t" . $x . "\t" . $y . "\n";

}
#append header back in files
while (($key, $val) = each(%flist)){
    if($val =~ /[a-zA-Z0-9]+[.]dat/) {
	#open file and append header to the front of it.
	$header = "Date\tHour\tX\tY\n";
	$file = "../data/" . $val;
	open IFILE, $file;
	@content = <IFILE>;
	open OFILE, ">" . $file;
	print OFILE $header;
	print OFILE @content;
    }
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
