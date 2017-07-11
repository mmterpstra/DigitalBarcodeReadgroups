#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Getopt::Std;

#$ARGV[0]="/home/terpstramm/workspace/FqToBcFq/src/head_400_120830_SN163_0474_BD0WDYACXX_L5_ACTGAT_1.merged.bam";
#$ARGV[1]="/home/terpstramm/workspace/FqToBcFq/src/head_400_120830_SN163_0474_BD0WDYACXX_L5_ACTGAT_1.merged.nugene.bam";

&main();


sub main {
	#use List::Util qw/max/; my $in;
	
	warn INIT_MESSAGE();

	my $opts = {};
	getopts('hp:l:',$opts);
	
	if( $opts -> {'h'} ){
		warn HELP_MESSAGE();
		exit();
	}

	my $bam = "";
	if(defined($ARGV[0]) && -e $ARGV[0]){
		$bam = $ARGV[0];
	}else{
		die "## ".localtime(time())." ## ERROR: file '$bam' does not exist!\n".HELP_MESSAGE();
	}
	


	die "No output bam specified".HELP_MESSAGE() if(not(defined($ARGV[1])) || length($ARGV[1])==0);
	my $bamout=$ARGV[1];
	
	#warn "$use";
	if(not($opts -> {'p'})){
		my $rgdata=CollectNugeneRgData($bam);
		ApplyNugeneRgData($bam,$bamout,$rgdata);
	}else{
		ApplyPicardTag($bam,$bamout,$opts);
	}
}

sub INIT_MESSAGE {
	
	die "## ".localtime(time())." ## ERROR Missing samtools please add samtools to your \$PATH\n" if(`which samtools` eq "");
	my $perlver= `(perl --version | head -n 2| tail -n 1 | perl -wpe 'chomp;s/This is //g') 2>/dev/null`;
	my $samtoolsver .= `samtools 2>/dev/stdout | grep 'Version' | cut -d\" \" -f2-`;
	chomp $samtoolsver;
	my $ver = VERSION_MESSAGE();
	my $time = localtime(time());
	
	my $append;
	my $host = `echo -n \$HOST `;
	$append .= " on host '$host'" if($host);
	my $user = `echo -n \$USER`;
	$append .= " as user '$user'" if($user);

	
	my $init =<<"EOF";
## $time ## INIT running '$0 @ARGV' version '$ver' with perl '$perlver' and samtools '$samtoolsver'$append.
EOF
	
}

sub HELP_MESSAGE {
	my $version = VERSION_MESSAGE();
	return <<"EOF";
use: $0 [-h|-p TAG  -l INT] input.bam output.bam

default mode
 Add readsgroups for each unique UMI. The output.bam will contain with the header and body containg the updated Readgroup/\@RG info.

options
 -p TAG	picard UmiAwareMarkDuplicatesWithMateCigar mode. Adds the umis to TAG instead of creating new readgroups. Recommended value is "RX" for TAG.
 -l INT shorten the UMI to INT length only functions on conjunction with -p

Needs Readgroup info to be added about the sample see: PicardAddOrReplaceReadgroups. 

Notes
	- using PicardAddOrReplaceReadgroups: be practical if it accepts your input then it is ok. 
	- required fields for default mode: LB. 
	- no streaming into this program => it will fail"
	- needs samtools to read/write bam files.
EOF
}

sub VERSION_MESSAGE {
	my $md5cmd ="(md5sum $0 | cut -d \" \" -f1,1)2>/dev/null";
	chomp $md5cmd;
	return `$md5cmd`;
}

sub ApplyPicardTag {
	my $bam = shift @_;
	my $bamout = shift @_;
	my $opts = shift @_;
	my $tag	= $opts -> {'p'};

	my $cmdin = "samtools view -h $bam";
	open(my $in,'-|',$cmdin) or die "## ".localtime(time())." ## ERROR invalid read from command $cmdin";
	
	my $cmdout = "samtools view  -S -b - > $bamout";
	#my $cmdout = "cat - > $bamout";
	open(my $out,'|-',$cmdout) or die "## ".localtime(time())." ## ERROR invalid write to command $cmdout";
	
	while(<$in>){
		chomp;
		my $line = $_;
		#warn $line;
		if($line =~ m/^\@/){
			print $out $line."\n";
		}else{
			my $sam;
			@{$sam} = split("\t",$line);
			#my $rgID = SamGetReadGroupID($sam);
			my $randombc = SamGetRandombc($sam);
			$randombc = substr($randombc,0,$opts -> {'l'}) if(defined($opts -> {'l'}) && $opts -> {'l'} > 0);
			$sam = SamRemoveRandombc($sam);
			
			#a
			$line .= "\t$tag:Z:$randombc";
			#die "$line";
			print $out $line."\n";
		}
	}
	close $out;
	warn "## ".localtime(time())." ## DONE\n"
}

sub CollectNugeneRgData {
	my $bam = shift @_;
	
	my $originalRgs=GetReadGroupData($bam);
	
	my $rgdata = CollectBarcodes($bam,$originalRgs);
	
	return $rgdata;
}
sub GetReadGroupData {
	my $bam = shift @_;
	
	my $cmd = "samtools view -H $bam";
	open(my $in,'-|',$cmd) or die "## ".localtime(time())." ## invalid read from command $cmd";

	my $rgs;

	while(<$in>){
		chomp;
		if(m/\@RG/){
			my $rg = ReadGroupParser($_);
			$rgs->{ $rg->{"ID"} }=$rg;
		}
	}
	
	close($in);
	
	return $rgs;
}

sub ReadGroupParser {
	my $line = shift @_;
	my @t = split("\t",$line);
	shift @t;
	#warn join(',', @t);
	my $rg;
	map{my ($key,$val) = split(":",$_);$rg->{$key}=${val}}(@t);
	#die Dumper($rg);
	
	#minimal required fields
	for my $k ("LB", "SM", "ID"){

		die $k if(not(defined($rg->{$k})));
		
	}
	
	return $rg;
}

sub CollectBarcodes {
	my $bam = shift @_;
	my $rgsOld = shift @_;
	
	my $cmd = "samtools view $bam";
	open(my $in,'-|',$cmd) or die "## ".localtime(time())." ## invalid read from command $cmd";

	my $rgs;
	my $bcs;
	while(<$in>){
		chomp;
		my $line=$_;
		my $sam;
		@{$sam} = split("\t",$line);
		#warn Dumper($sam);
		
		my $randombc = SamGetRandombc($sam);
		my $rgID = SamGetReadGroupID($sam);
		if(not(defined($bcs->{$rgID."_".$randombc}))){
			$bcs->{$rgID."_".$randombc}->{'count'}++;
			$bcs->{$rgID."_".$randombc}->{'ID'}=scalar(keys(%{$bcs}));
			#die "$. test".Dumper(scalar(keys(%{$bcs}))) if($. > 20);
			#warn $randombc;
			#warn $bcs->{$randombc}->{'ID'};
			#die Dumper($bcs);
			#die Dumper($rgsOld);
			%{ $rgs->{ $bcs->{$rgID."_".$randombc}->{'ID'} } } = %{$rgsOld->{$rgID}};
			$rgs->{ $bcs->{$rgID."_".$randombc}->{'ID'} }->{"ID"}=$bcs->{$rgID."_".$randombc}->{'ID'};
			$rgs->{ $bcs->{$rgID."_".$randombc}->{'ID'} }->{"LB"}.="_".$randombc;
			#die "$. test".Dumper($rgs) if($. > 20);
		}
		
	}
	
	my $rgdata;
	$rgdata->{"rgs"} = $rgs;
	$rgdata->{"bcs"} = $bcs;
	
	#die Dumper($rgdata);
	
	return ($rgdata);
}

sub SamGetReadGroupID {
	my $sam = shift @_;
	my $rgId;
	#warn join "#",@{$sam};
	for my $field (@{$sam}){
		if($field =~ m/^RG\:Z\:(.*$)/){
			#warn $1;
			$rgId = $1;
		}
	}
	die "## ".localtime(time())." ## ERROR No RGID present in sam record! Dump of sam record:".Dumper($sam) if(not(defined($rgId)));
	
	return $rgId;
}

sub SamGetRandombc {
	
	my $sam = shift @_;
	my @readname = split(':',$sam->[0]);
	#here
	$readname[2] =~ /([ATCGN]+$)/;
	my $randombc=$1;
	die "## ".localtime(time())." ## ERROR No randombc present in sam record! Dump of sam record:".Dumper($sam) if(not(defined($randombc)));
	
	return $randombc;

}

sub SamRemoveRandombc {
		
	my $sam = shift @_;
	my @readname = split(':',$sam->[0]);
	#here
	$readname[2] =~ s/[ATCGN]+$//g;
	#my $randombc=$1;
	#die "No randombc present in sam record! Dump of sam record:".Dumper($sam) if(not(defined($randombc)));
	$sam->[0]=join("\t",@readname);
	
	return $sam;
}

sub ApplyNugeneRgData{
	my $bam = shift @_;
	my $bamout = shift @_;
	my $rgdata = shift @_;
	
	my $cmdin = "samtools view -h $bam";
	open(my $in,'-|',$cmdin) or die "## ".localtime(time())." ## ERROR invalid read from command $cmdin";
	
	my $cmdout = "samtools view  -S -b - > $bamout";
	#my $cmdout = "cat - > $bamout";
	open(my $out,'|-',$cmdout) or die "## ".localtime(time())." ## ERROR invalid read from command $cmdout";
	
	my $readsGroupsPrinted=0;
	
	while(<$in>){
		chomp;
		my $line = $_;
		#warn $line;
		if($line =~ m/^\@RG/){
			if ($readsGroupsPrinted==0){
				print $out ReadGroupDataToString($rgdata);
				$readsGroupsPrinted++;
			}
		}elsif($line =~ m/^\@/){
			print $out $line."\n";
		}else{
			my $sam;
			@{$sam} = split("\t",$line);
			my $rgID = SamGetReadGroupID($sam);
			my $randombc = SamGetRandombc($sam);
			$sam = SamRemoveRandombc($sam);
			
			#add new random barcode
			my $newRgId=$rgdata->{'bcs'}->{$rgID.'_'.$randombc}->{'ID'};
			$line =~ s/RG:Z:[^\t\n]+/RG:Z:$newRgId/;
			
			print $out $line."\n";
		}
	}
	close $out;
	
}


sub ReadGroupDataToString {
	my $rgdata = shift @_;
	#die "rgdata:".Dumper($rgdata)."line $.";
	my $string = '';
	for my $rgid (keys(%{$rgdata->{"rgs"}})){
		$string.= "\@RG";
		for my $rgfield ("ID","LB","PL","PU","SM"){
			die "Required field '$rgfield' is not present!" if(not(defined($rgdata->{"rgs"}->{$rgid}->{$rgfield})));
			$string.="\t".$rgfield.":".$rgdata->{"rgs"}->{$rgid}->{$rgfield};
			delete($rgdata->{"rgs"}->{$rgid}->{$rgfield});
		}
		for my $rgfield (keys(%{$rgdata->{"rgs"}->{$rgid}})){
			$string.="\t".$rgfield.":".$rgdata->{"rgs"}->{$rgid}->{$rgfield};
		}
		$string.="\n";
	}
	#die $string;
	return $string;
}
