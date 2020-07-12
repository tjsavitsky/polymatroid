#!/usr/bin/perl -w

# convert the flats (as bitmasks) and their ranks representing
# a k-polymatroid into to a more human-readable format
#
#example:
#input 
#1,0 3,1 5,1 9,1 7,2 b,2 d,2 f,3
#
#output
#a,0  ab,1  ac,1  ad,1  abc,2  abd,2  acd,2  abcd,3

my @alphabet = ("a" .. "z");
my $nullset = "{}";
my $maxsize = 26;

while (<>) {
	my $line = $_;
	my @linestr = split / /, $line;
	#print $line;
	foreach $word (@linestr) {
		next if ($word !~ m/,/);
		my ($bitmask, $rank) = split /,/, $word;
		#print $bitmask, ",", $rank, " ";
		my $outset = &bitmask_to_string($bitmask);
		print $outset, ",", $rank, "  ";
	}
	#print "\n\n\n\n\n\n\n";
	print "\n";
}

sub bitmask_to_string {
	my $bitmask = pop;
	my $outset = "";
	if (hex($bitmask) == 0) {
		$outset = $nullset;
	} else {
		for (my $i=0; $i<$maxsize; $i++) {
			if (hex($bitmask) & (1 << $i)) {
				$outset .= $alphabet[$i];
			}
		}
	}
	return $outset;
}
