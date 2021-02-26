#!/usr/bin/perl

if ($#ARGV < 2) {
    print STDERR "Usage: $0 <id_list_file> <category> <category_name>\n";
    exit (-1);
}

$base_url = 'http://www.cazy.org';


$id_list_file = shift @ARGV;
$cat = shift @ARGV;
$cat_name = shift @ARGV;

$cat_name_url = $cat_name;
$cat_name_url =~ s/\s/-/g;
$cat_url = $base_url.'/'.$cat_name_url.'.html';

@ids = ();
open (ID_LIST_FILE, $id_list_file);
while (<ID_LIST_FILE>) {
    chomp;
    ($id, @rest) = split (/\t/, $_);
    push (@ids, $id);
}
close (ID_LIST_FILE);


@out = ();
$target_pattern = 'Activities in Family</th><td class="tdsum">';
$alt_target_pattern = 'Note</th><td class="tdsum">';
foreach $id (@ids) {
    $this_url = $base_url.'/'.$id.'.html';
    $page_source = `curl $this_url`;

    $found_desc = '';
    for $line (split (/\n/,$page_source)) {
	if ($line =~ /$target_pattern/) {
	    $desc = $line;
	    $desc =~ s/^.*$target_pattern//;
	    $desc =~ s/<\/td><\/tr>.*$//; 
	    if ($desc !~ /^\s*$/) {
		$found_desc = $desc;
		last;
	    }
	}
	elsif ($line =~ /$alt_target_pattern/) {
	    $desc = $line;
	    $desc =~ s/^.*$alt_target_pattern//;
	    $desc =~ s/<\/td><\/tr>.*$//; 
	    if ($desc !~ /^\s*$/) {
		$found_desc = $desc;
		last;
	    }
	}
    }
    push (@out, join("\t", $id, join (" - ", '<a href="'.$cat_url.'" target="_cazy_cat_link">'.$cat_name.'</a>', '<a href="'.$this_url.'" target="_cazy_fam_link">'.$id.'</a>', $found_desc)));
}

# output
if ($out_file) {
    open (OUT, '>'.$out_file);
    select (OUT);
}
print join ("\n", @out)."\n"  if (@out);
if ($out_file) {
    close (OUT);
    select (STDOUT);
}

# done
exit 0;
