#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper;
use Getopt::Std;
use Scalar::Util qw(looks_like_number);

&main();

sub main {
	#use List::Util qw/max/; my $in;
	
	my $use="$0 [-q INT] randombarcodes.fq.gz outdir reads_1.fq.gz reads_2.fq.gz.

Note:
now filters for illumina 1.8 base quality < -q INT barcode reads and N bases in randombarcodes.fq.gz";
	
	#open/check some resuired files
	
	my $opts;
        getopts('q:', \%{$opts});
        die "[FATAL] if '-q INT' is specified use integer for INT".$use if(defined($opts -> {'q'}) && not(looks_like_number($opts -> {'q'})));
	
	die $use."\n" if(scalar(@ARGV) < 2);

	open(my $randomBcHandle,"-|",'gzip -dc '.$ARGV[0]) 
	or die "Cannot read open randombarcode file ".$ARGV[0]."\n".$use;
	
	if(! -e $ARGV[1]){
		die "Outdir does not exist! ".$ARGV[1]."\n".$use;;
	}
	
	#open handles fastq files
	
	open(my $fastq1Handle,"-|",'gzip -dc '.$ARGV[2]) 
	or die "[FATAL] Cannot read open fq1 file ".$ARGV[2];
	
	open(my $fastq1OutHandle,"|-",'gzip -c > '.GetOutFileName($ARGV[1],$ARGV[2])) 
	or die "[FATAL] Cannot write fq1 file ".$ARGV[2];
	
	my $fastq2Handle;
	my $fastq2OutHandle;
	
	if(defined($ARGV[3]) && -e $ARGV[3]){
	
		open($fastq2Handle,"-|",'gzip -dc '.$ARGV[3]) 
		or die "Cannot read open fq2 file ".$ARGV[3];
		
		open($fastq2OutHandle,"|-",'gzip -c >'.GetOutFileName($ARGV[1],$ARGV[3])) 
		or die "Cannot write fq2 file ".$ARGV[3];
	
	}
	
	
	
	warn "## ". localtime(time()). " ## INFO ## Starting iteration.\n";
	my $timelast = time();

	#if reading record fails from both fq files then end
	my $stats;
	
	while( ((((not(defined($ARGV[3])))|| (my $fq2 = ReadFastq(\$fastq2Handle) ))) && (my $fq1 = ReadFastq(\$fastq1Handle))) && (my $rfq = ReadFastq(\$randomBcHandle)) ){
		if($. % 1000 == 0){
			my $time = time();
			if($time % 20 == 0 && $time != $timelast){
				$timelast = $time;
	        	        warn "## ". localtime($time). " ## INFO ## Currently at line ". $. ."\n";
			}
		}
		my $rbc = $rfq->[1];
		my $fcid= getFCID($fq1);
		
		$stats -> {'recordcount'}++;
		if(not(defined($opts -> {'q'})) || TestRandomBarcodeQual($opts -> {'q'},$rfq)){
			$stats -> {'passcount'}++;
			
			$fq1 = setFCID($fq1,$fcid."_".$rbc);
			WriteFastq(\$fastq1OutHandle,$fq1);
			
			if(defined($ARGV[3]) && -e $ARGV[3]){
				$fq2 = setFCID($fq2,$fcid."_".$rbc);
				WriteFastq(\$fastq2OutHandle,$fq2);
			}
		}
	}
	if( ((((defined($ARGV[3]) && -e $ARGV[3] ) && (my $fq2 = ReadFastq(\$fastq2Handle) ))) || (my $fq1 = ReadFastq(\$fastq1Handle))) || (my $rfq = ReadFastq(\$randomBcHandle)) ){
		
		unlink(GetOutFileName($ARGV[1],$ARGV[3]))if(-e GetOutFileName($ARGV[1],$ARGV[2]));
		unlink(GetOutFileName($ARGV[1],$ARGV[3]))if(-e GetOutFileName($ARGV[1],$ARGV[3]));
		die "Some fastq handles are still readable:This program produced invalid output.Check your fastq files!";
	}
	 warn "## ". localtime(time()). " ## INFO ## Done with " . $stats -> {'recordcount'} . " records processed and " .$stats -> {'passcount'} . " records passed at line ". $. ."\n";
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

sub GetOutFileName {
	my $outdir = shift @_;
	my $inputFile = shift @_;
	my $outFile = $inputFile;
	
	$outFile =~ s/^.*\/|\.fastq\.gz|\.fq\.gz//g;
	$outFile = $outdir . $outFile . '.fq.gz';
	
	return $outFile;
}

sub TestRandomBarcodeQual {
	my $q = shift @_;
	my $fq = shift @_;
	if(not( $fq->[1] =~ /^[ATCGatcg]*$/)){
		#warn $fq->[1] . $&;
		#die $fq->[1] if($. > 20);
		return 0;
	}
	#warn "something goes wrong here";
	my @quals = split('',$fq->[2]);
	#die scalar(@quals), Dumper(\@quals);
	my @qualsOrd = map{ord()}(@quals);
	for my $qual (@qualsOrd){
		#die Dumper($fq,\@quals);
		if($qual < (33+$q)){
			#die Dumper($fq,\@quals,\@qualsOrd);
			return 0 ;#if($qual < (33+$q));
		}
	}
	return 1;
}
