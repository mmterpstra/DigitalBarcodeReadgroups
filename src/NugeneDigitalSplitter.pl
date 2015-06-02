#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
#$ARGV[0]="/home/terpstramm/workspace/FqToBcFq/src/head_400_120830_SN163_0474_BD0WDYACXX_L5_ACTGAT_1.merged.bam";
#$ARGV[1]="/home/terpstramm/workspace/FqToBcFq/src/head_400_120830_SN163_0474_BD0WDYACXX_L5_ACTGAT_1.merged.nugene.bam";

&main();


sub main {
	#use List::Util qw/max/; my $in;
	my $use="$0 input.bam output.bam with the header and body containg the updated Readgroup/\@RG info. Needs Readgroup info to be added about the sample see: PicardAddOrReplaceReadgroups. Note using PicardAddOrReplaceReadgroups: be practical if it accepts your input then it is okay important fields:LB. no streaming into this program => it will fail";
	
	my $bam;
	$bam = $ARGV[0] if(-e $ARGV[0]);
	die "$ARGV[0]: file does not exist!$use" if(! -e $ARGV[0]);

	die "No output bam specified" if(not(defined($ARGV[1])) || length($ARGV[1])==0);
	my $bamout=$ARGV[1];
	
	#warn "$use";
	my $rgdata=CollectNugeneRgData($bam);
	ApplyNugeneRgData($bam,$bamout,$rgdata);
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
	open(my $in,'-|',$cmd) or die "invalid read from command $cmd";

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
	open(my $in,'-|',$cmd) or die "invalid read from command $cmd";

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
	die "No RGID present in sam record! Dump of sam record:".Dumper($sam) if(not(defined($rgId)));
	
	return $rgId;
}

sub SamGetRandombc {
	
	my $sam = shift @_;
	my @readname = split(':',$sam->[0]);
	#here
	$readname[2] =~ /([ATCGN]+$)/;
	my $randombc=$1;
	die "No randombc present in sam record! Dump of sam record:".Dumper($sam) if(not(defined($randombc)));
	
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
	open(my $in,'-|',$cmdin) or die "invalid read from command $cmdin";
	
	my $cmdout = "samtools view  -S -b - > $bamout";
	#my $cmdout = "cat - > $bamout";
	open(my $out,'|-',$cmdout) or die "invalid read from command $cmdout";
	
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
			
			my $newRgId=$rgdata->{'bcs'}->{$rgID.'_'.$randombc}->{'ID'};
			$line =~ s/RG:Z:[^\t\n]+/RG:Z:$newRgId/;
			
			print $out $line."\n";
		}
	}
	close $out;warn $bamout
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