{
	if (/^>(.+)\s/){
		$id = $1;
		push @ids, $id;
		$h{$id}="";
	}
	else{
		$h{$id}.=$_;
	}
}
END{
	$seq1_gap = $h{$ids[0]};
	$h{$ids[0]} =~ s/-//g;
	$seq2_gap= $h{$ids[1]};
	$h{$ids[1]} =~ s/-//g;
	print join("\t",$ids[0],"1",$seq1_len,$ids[1],"1",$seq2_len,"+","0\n");
	print "$seq1_gap\n";
	print "$seq2_gap\n";
}