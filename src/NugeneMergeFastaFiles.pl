#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
&main();

sub main {
	#use List::Util qw/max/; my $in;
	my $use="$0 randombarcodes.fq.gz reads_1.fq.gz reads_2.fq.gz";
	
	open(my $randomBcHandle,"-|",'gzip -dc '.$ARGV[0]) 
	or die "Cannot read open randombarcode file ".$ARGV[0];
	
	open(my $fastq1Handle,"-|",'gzip -dc '.$ARGV[1]) 
	or die "Cannot read open fq1 file ".$ARGV[1];
	
	my $fastq2Handle;
	my $fastq2OutHandle;
	
	if(-e $ARGV[2]){
	
		open($fastq2Handle,"-|",'gzip -dc '.$ARGV[2]) 
		or die "Cannot read open fq2 file ".$ARGV[2];
		
		open($fastq2OutHandle,"|-",'gzip -c >'.$ARGV[2].".nug.fq.gz") 
		or die "Cannot read open fq2 file ".$ARGV[2];
	
	}
	
	open(my $fastq1OutHandle,"|-",'gzip -c > '.$ARGV[1].".nug.fq.gz") 
	or die "Cannot read open fq1 file ".$ARGV[1];
	
	
	
	#warn;
	
	while( (((( ! -e $ARGV[2] )|| (my $fq2 = ReadFastq(\$fastq2Handle) ))) && (my $fq1 = ReadFastq(\$fastq1Handle))) && (my $rfq = ReadFastq(\$randomBcHandle)) ){
		
		warn $.if $. =~ /00000$/;
		my $rbc = $rfq->[1];
		my $fcid= getFCID($fq1);
		
		$fq1 = setFCID($fq1,$fcid."_".$rbc);
		WriteFastq(\$fastq1OutHandle,$fq1);
		
		if(-e $ARGV[2]){
			$fq2 = setFCID($fq2,$fcid."_".$rbc);
			WriteFastq(\$fastq2OutHandle,$fq2);
		}
		
	}
}
sub getFCID{
	my $fq = shift @_;
	#die $fq->[0];
	my @h = split(":",$fq->[0]);
	return $h[2];
}
sub setFCID{
	my $fq = shift @_;
	my $fcid = shift @_;
	my @h = split(":",$fq->[0]);
	$h[2]=$fcid;
	$fq->[0]=join(":",@h);
	return $fq;
}

sub ReadFastq {
	my $fqHandle = ${shift(@_)};
	return if(eof($fqHandle));
	my $seqHeader = <$fqHandle>;
	chomp $seqHeader;
	
	$seqHeader = substr($seqHeader,1);
	my $seq = <$fqHandle>;
	
	chomp $seq;
	my $qualHeader = <$fqHandle>;
		
	my $qual = <$fqHandle>;
	chomp $qual;
	
	my $fastq;
	$fastq->[1]=$seq;
	$fastq->[0]=$seqHeader;
	$fastq->[2]=$qual;
	
	#warn $fastq->[0];
	
	return $fastq;	
}

sub WriteFastq {
	my $fqHandle = ${shift(@_)};
	my $fq = shift @_;
	
	print $fqHandle "\@".$fq->[0]."\n";
	print $fqHandle $fq->[1]."\n";
	print $fqHandle "\+\n";
	print $fqHandle $fq->[2]."\n";
	
	#warn "\@".$fq->[0]."\n";
	#print $fqHandle $fq->[1]."\n";
	#print $fqHandle "\+\n";
	#print $fqHandle $fq->[2]."\n";
}
