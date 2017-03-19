use strict;
use warnings;
use Getopt::Long;

my($file1,$file2,$maf1,$maf2,$mincov1,$mincov2,$noindel,$outname,$bed);

my $usage = "USAGE:\nperl $0 --in1 <varstat file 1> --in2 <varstat file 2> --out <output prefix>\n";
$usage .= "<varstat file> is produced by ~/scripts/vcf_var_stat.pl [Required]\n";
$usage .= "<out> is the output prefix. [Required]\n";
$usage .= "--maf1 <maf1>: minor allele frequency for <varstat file 1>. [Optional]\n";
$usage .= "--maf2 <maf2>: minor allele frequency for <varstat file 2>. [Optional]\n";
$usage .= "--cov1 <cov1>: minor coverage of a variant to be kept in <varstat file 1>. [Optional]\n";
$usage .= "--cov1 <cov2>: minor coverage of a variant to be kept in <varstat file 2>. [Optional]\n";
$usage .= "--noindel: to ignore indels in both files. [Optional]\n";
$usage .= "--bed <bed>: only keep variants in the bed file. [Optional]\n";

GetOptions(
	"in1=s" => \$file1,
	"in2=s" => \$file2,
	"out=s" => \$outname,
	"maf1=s" => \$maf1,
	"maf2=s" => \$maf2,
	"cov1=s" => \$mincov1,
	"cov2=s" => \$mincov2,
	"noindel!" => \$noindel,
	"bed=s" => \$bed,
) or die $usage;

die $usage unless(defined $file1 and defined $file2 and defined $outname);


unless(defined $maf1){
	$maf1 = 0;
}
unless(defined $maf2){
	$maf2 = 0;
}
unless(defined $mincov1){
	$mincov1 = 0;
}
unless(defined $mincov2){
	$mincov2 = 0;
}


my %hash_bed;
if(defined $bed){
	open(IN,"<$bed") or die $!;
	while(<IN>){
		chomp;
		next if($_ =~ /^#/);
		my($chr,$start,$end,$other) = split/\t/;
		push @{$hash_bed{$chr}{regs}}, "$start:$end";
	}
	close IN;
}


my %hash_var;

open(IN,"<$file1") or die $!;

=format
##samplenum=1257
#CHROM	POS	REF	ALT	FILTER	ALLELEnum	HETnum	NAnum	COVfreq	ALLELEfreq
Chr1	1030	TA	T	PASS	1453,13	31	1533	0.484	0.480,0.004
Chr1	1044	TA	T	SnpCluster	951,271	575	1233	0.403	0.314,0.089
Chr1	1054	A	C	SnpCluster	1928,7	29	1066	0.639	0.636,0.002
Chr1	1057	C	T	SnpCluster	1924,13	47	1046	0.639	0.635,0.004
Chr1	1058	T	A	SnpCluster	1875,27	98	1030	0.628	0.619,0.009
Chr1	1067	A	G	SnpCluster	2089,5	2	934	0.691	0.689,0.002
Chr1	1071	TAACCCTA	T,TAAACCCTA	SnpCluster	1676,49,118	274	913	0.608	0.553,0.016,0.039
Chr1	1074	C	A,T	SnpCluster	2034,8,8	86	894	0.677	0.671,0.003,0.003
=cut

while(<IN>){
	chomp;
	next if($_ =~ /^#/);
	my($chr,$pos,$ref,$alts_join,$filter,$alleleNums_join,$hetNum,$naNum,$covFreq,$alleleFreqs_join) = split/\t/;
	my @alts = split/,/,$alts_join;
	my @alleles = ($ref,@alts);
	my @alleleNums = split/,/, $alleleNums_join;
	my @alleleFreqs = split/,/, $alleleFreqs_join;
	next if($alleleNums[0] < $maf1);
	next if($covFreq < $mincov1);
	if(defined $noindel){
		next if(length($ref) > 1);
	}
	for(my $i = 1; $i < @alleles; $i++){
		my $allele = $alleles[$i];
		if(defined $noindel){
			next if(length($allele) > 1);
		}
		my $alleleFreq = $alleleFreqs[$i];
		next if($alleleFreq < $maf1);
		$hash_var{$chr}{$pos}{$ref}{$allele}{freq1} = $alleleFreq;
	}
}
close IN;


open(IN,"<$file2") or die $!;

while(<IN>){
	chomp;
	next if($_ =~ /^#/);
	my($chr,$pos,$ref,$alts_join,$filter,$alleleNums_join,$hetNum,$naNum,$covFreq,$alleleFreqs_join) = split/\t/;
	my @alts = split/,/,$alts_join;
	my @alleles = ($ref,@alts);
	my @alleleNums = split/,/, $alleleNums_join;
	my @alleleFreqs = split/,/, $alleleFreqs_join;
	next if($alleleNums[0] < $maf2);
	next if($covFreq < $mincov2);
	if(defined $noindel){
		next if(length($ref) > 1);
	}
	for(my $i = 1; $i < @alleles; $i++){
		my $allele = $alleles[$i];
		if(defined $noindel){
			next if(length($allele) > 1);
		}
		my $alleleFreq = $alleleFreqs[$i];
		next if($alleleFreq < $maf2);
		$hash_var{$chr}{$pos}{$ref}{$allele}{freq2} = $alleleFreq;
	}
}
close IN;

open(OUT1,">$outname.uniq1");
open(OUT2,">$outname.uniq2");
open(BOTH,">$outname.both");

foreach my $chr(sort keys %hash_var){
	my @positions = sort keys %{$hash_var{$chr}};
	if(defined $bed){
		@positions = pos_in_reg(\@positions,\@{$hash_bed{$chr}{regs}});
	}
	foreach my $pos(sort {$a <=> $b} @positions){
	foreach my $ref(sort keys %{$hash_var{$chr}{$pos}}){
	foreach my $alt(sort keys %{$hash_var{$chr}{$pos}{$ref}}){
		if(exists $hash_var{$chr}{$pos}{$ref}{$alt}{freq1} and exists $hash_var{$chr}{$pos}{$ref}{$alt}{freq2}){
			my $freq1 = $hash_var{$chr}{$pos}{$ref}{$alt}{freq1};
			my $freq2 = $hash_var{$chr}{$pos}{$ref}{$alt}{freq2};
			print BOTH "$chr\t$pos\t$ref\t$alt\t$freq1\t$freq2\n";
		}elsif(exists $hash_var{$chr}{$pos}{$ref}{$alt}{freq1}){
			my $freq1 = $hash_var{$chr}{$pos}{$ref}{$alt}{freq1};
			print OUT1 "$chr\t$pos\t$ref\t$alt\t$freq1\n";
		}else{
			my $freq2 = $hash_var{$chr}{$pos}{$ref}{$alt}{freq2};
			print OUT2 "$chr\t$pos\t$ref\t$alt\t$freq2\n";
		}
	}
	}
	}		
}
close OUT1;
close OUT2;
close BOTH;

sub pos_in_reg{
	my ($arr1,$arr2) = @_;
	my @positions = @{$arr1};
	my @regions = @{$arr2};
	my @outs = ();
	my %hash = ();
	foreach(@regions){
		my($start,$end) = split/\:/;
		if(exists $hash{$start}){
			if($hash{$start}{end} < $end){
				$hash{$start}{end} = $end;
			}
		}else{
			$hash{$start}{end} = $end;
		}
	}
	my @starts = sort {$a <=> $b} keys %hash;
	foreach my $pos(sort {$a <=> $b} @positions){
		for(my $j = 0; $j < @starts; $j++){
			my $start = $starts[$j];
			my $end = $hash{$start}{end};
			last if($start > $pos);
			if($end < $pos){
				splice(@starts,$j,1);
				$j--;
				next;
			}
			push @outs, $pos;
			last;
		}
	}
	return(@outs);
}
