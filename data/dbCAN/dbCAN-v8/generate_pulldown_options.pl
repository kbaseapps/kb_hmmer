#!/usr/bin/perl

if ($#ARGV < 2) {
    print STDERR "Usage: $0 <id_list_file> <category> <category_name>\n";
    exit (-1);
}


$id_list_file = shift @ARGV;
$cat = shift @ARGV;
$cat_name = shift @ARGV;


open (ID_LIST_FILE, $id_list_file);
while (<ID_LIST_FILE>) {
    chomp;
    ($id, @rest) = split (/\t/, $_);

    $ind_1 = "\t\t\t\t";
    $ind = "\t\t\t\t\t";
    print "$ind_1\{\n";
    print "$ind\"value\": \"$id\",\n";
    print "$ind\"display\": \"$cat_name - $id\",\n";
    print "$ind\"id\": \"dbCAN_".$cat."_".$id."\",\n";
    print "$ind\"ui-name\": \"dbCAN_".$cat."_".$id."\"\n";
    print "$ind_1\},\n";
}
