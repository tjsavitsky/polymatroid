#!/usr/bin/perl -w

my %hash = ();
my $count = 0;
my $intermed = 100000000;

while (<>) {
	my $line = $_;
	$count++;
	$line =~ /(\d+)\s*$/;
	my $rank = $1;
	$hash{$rank}++;
	if ($count % $intermed == 0) {
		print "\nAfter $count lines:\n";
		&printstuff;
	}
}

print "\nFINAL TOTALS:\n";
&printstuff;

sub printstuff {

print "Count\t  Rank\n";
my $sum = 0;
for my $key (sort {$a <=> $b} (keys %hash)) {
	my $value = $hash{$key};
	print "$value\t- rank $key\n";	
	$sum += $value;
}
print "------------------------------\n";
print "$count\t- total\n";
$other = $count - $sum;
print "$other\t- other\n";
}

