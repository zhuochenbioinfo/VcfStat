use strict;
use warnings;
use Getopt::Long;

my($vcf,$out,$keepList,$all,$CHROM);

my $usage = "USAGE:\nperl $0 --vcf <vcf> --out <out> --keep <keep list> --all --chr <keep chr>\n";
$usage .= "<vcf> is the input vcf fle. [Necessary].\n";
$usage .= "<out> is the output file. [Necessary]\n";
$usage .= "<keep list> is a list of samples to be keep. [Optional]\n";
$usage .= "Use --all to output all variants. Default will output SNPs and InDels with only PASS or SnpCluster tag. [Optional]\n";

GetOptions(
	"vcf=s" => \$vcf,
	"out=s" => \$out,
	"keep=s" => \$keepList,
	"all=s" => \$all,
	"chr=s" => \$CHROM,
) or die $usage;

die $usage unless(defined $vcf and defined $out);

my $keepChrNum;
if(defined $CHROM){
	($keepChrNum) = $CHROM =~ /Chr(\d+)/;
}

my %hash_keep;
if(defined $keepList){
	open(IN,"<$keepList") or die $!;
	while(<IN>){
		next if($_ =~ /^#/);
		my($sample,@others) = split/\t/;
		$hash_keep{$sample} = "";
	}
	close IN;
}

open(IN,"<$vcf") or die $!;
open(OUT,">$out");

my @samples = ();
my @keepRanks = ();

while(<IN>){
	chomp;
	next if($_ =~ /^##/);
	my($chr,$pos,$id,$ref,$alts_join,$qual,$filter,$info,$format,@datas) = split/\t/;
	my @alts = split/,/,$alts_join;
	my @alleles = ($ref,@alts);
	if($_ =~ /^#CHROM/){
		@samples = @datas;
		my $num = 0;
		if(defined $keepList){
			for(my $i = 0; $i < @datas; $i++){
				if(exists $hash_keep{$datas[$i]}){
					push @keepRanks, $i;
					$num++;
				}
			}
		}else{
			my $allnum = @samples;
			$num = $allnum;
			$allnum--;
			@keepRanks = (0..$allnum);
		}
		print OUT "##samplenum=$num\n";
		print OUT "#CHROM\tPOS\tREF\tALT\tFILTER\tALLELEnum\tHETnum\tNAnum\tCOVfreq\tALLELEfreq\n";
		next;
	}
	unless(defined $all){
		next unless($filter eq "PASS" or $filter eq "SnpCluster");
	}
	if(defined $keepChrNum){
		my($chrNum) = $chr =~ /Chr(\d+)/;
		last if($chrNum > $keepChrNum);
		next if($chrNum < $keepChrNum);
	}
	my %hash_tmp;
	my $covNum = 0;
	foreach my $rank(@keepRanks){
		my $spot = $datas[$rank];
		my $gt = "NA";
		if($spot =~ /^(\d+)\/(\d+)/){
			$gt = $1;
			if($1 ne $2){
				$gt = "HET";
			}else{
				$covNum ++;
			}
		}
		$hash_tmp{$gt} ++;
	}
	my @alleleNums = ();
	my @alleleFreqs = ();
	my $naNum = 0;
	my $hetNum = 0;
	if(exists $hash_tmp{HET}){
		$hetNum = $hash_tmp{HET};
		delete($hash_tmp{HET});
	}
	if(exists $hash_tmp{NA}){
		$naNum = $hash_tmp{NA};
		delete($hash_tmp{NA});
	}
	for(my $i = 0; $i < @alleles; $i++){
		my $count = 0;
		if(exists $hash_tmp{$i}){
			$count = $hash_tmp{$i};
		}
		push @alleleNums, $count;
		my $alleleFreq = $count/@keepRanks;
		$alleleFreq = sprintf("%.3f",$alleleFreq);
		push @alleleFreqs, $alleleFreq;
	}
	my $covFreq = $covNum/@keepRanks;
	$covFreq = sprintf("%.3f",$covFreq);
	print OUT "$chr\t$pos\t$ref\t$alts_join\t$filter\t".join(",",@alleleNums)."\t$hetNum\t$naNum\t$covFreq\t".join(",",@alleleFreqs)."\n";
}
close OUT;
