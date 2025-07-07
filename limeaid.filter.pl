use strict;
use JSON;

# header:
# ID	Sequence	Element_Hits	Sequence_Length	Element_Annotation	Element_Divergence	Element_Proportion	Element_Percentage	Element_Designation	Orientation	FILTER_RESULTS	Tail_Begins	Tail_Type	Tail_Length	Tail_Seed_Hits	Unique_Element_Count	Twin_Priming_Flag

sub process_rmsk {
	my $rmsk = shift(@_);
	my @rmsk=split ", ", $rmsk;
	for $i(@rmsk) {
		my @tmp=split /\s/, $i;
		my $strand = $tmp[8];
		my $te=$tmp[9];
		my $class=$tmp[10];
		my $start=$strand eq "+"?$tmp[11]:$tmp[13]; 
		my $end = $tmp[12];
		my $left=$strand eq "+"?$tmp[13]:$tmp[11]; 
		my $left=~s/\((.*)\)/$1/;
		my $sv_length = $tmp[6] -$tmp[5]+1;
		my $te_length=$end-$start+1;
		my $te_total=$tmp[12]+$1;
		my %hash;
		$hash{$te}{"class"}=$class;
		$hash{$te}{"te_length"}+=$te_length;
		$hash{$te}{"sv_length"}+=$sv_length;
		$hash{$te}{"frac"}+= $te_length/$te_total;
		$hash{$te}{"left"}=$1 if (not exists $hash{$te}{"left"}) || $1 < $hash{$te}{"left"}; 
		$hash{$te}{"start"}=$start if (not exists $hash{$te}{"start"}) || $start < $hash{$te}{"start"}; 
		$hash{$te}{"end"}=$end if (not exists $hash{$te}{"end"}) || $end > $hash{$te}{"end"};
	};
	return \%hash;
}

sub filter_ltr {
	my @F = @_;
	my $hashref = process_rmsk($F[2]);
	my $tes_proportion_ref = decode_json $F[6];
	my $total_length; # total length of LTR elements matching column 9.
	my $ltr_frac;
	my $int_frac;
	for my $te(keys %$hashref){
		if ($hashref->{$te} eq $F[8]){
			$total_length += $hashref->{$te}->{"sv_length"};
			if ($hashref->{$te} =~ /-int/) {
				$int_frac += $hashref->{$te}->{"frac"}
			}
			else {$ltr_frac += $hashref->{$te}->{"frac"}};
			
		}
	}
	if $total_length / $F[3] > 0.8;
	return 
}


sub filter_alu {
	my @F = @_;
	$rmsk_aoa = process_rmsk($F[2]);

	return 
}


open IN, $ARGV[0];

while (<IN>) {
	chomp;
	my @F= split "\t", $_;
	next if $. ==1;
	next if $F[8] eq "";
	if ($F[8] =~ /^LTR/) {
		my $good_ltr = filter_ltr(@F);
	}elsif ($F[8] eq "LINE/L1") {. # this filter is human specific
		my $good_l1 = filter_l1(@F);
	}elsif ($F[8] eq "SINE/Alu") { # human specific
		my $good_alu = filter_alu(@F);
	}elsif ($F[8] eq "Retroposon/SVA") {
		my $good_sva = filter_sva(@F);
	}
}

close IN;
# print all.



my %hash;



$str = encode_json \%hash;
%lookup = map {(uc $_, $hash{$_})} keys %hash;
$perc=$lookup{uc($F[4])}{"perc"};$start=$lookup{uc($F[4])}{"start"}; 
$end=$lookup{uc($F[4])}{"end"};$left=$lookup{uc($F[4])}{"left"};
print "$_\t$str\t$perc\t$start\t$end\t$left"}






