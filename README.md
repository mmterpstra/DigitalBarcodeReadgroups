

# DigitalBarcodeReadgroups

----

Create readgroups from barcodes with random nucleotides incorporates into them for improving sequencing depth. The analysis for marking / removing duplicates does not handle random barcodes and should be modified to include random barcodes in the analysis. This is a git for noting down the steps.

This is primairly used for the Nugene Ovation® Target Enrichment System.  (link)[http://www.nugen.com/nugen/index.cfm/products/ovation/ovation-target-enrichment-system/]


# basecall to fastq

----

Use bcl to fastq from the illumina website 

install bcltofastq

```
#read http://support.illumina.com/downloads/bcl2fastq_conversion_software.html
mkdir -p path/to/bclToFastq/
cd path/to/bclToFastq/; wget  ftp://webdata:webdata@ussd-ftp.illumina.com/Downloads/Software/bcl2fastq/bcl2fastq-1.8.4.tar.bz2
tar -jxf bcl2fastq-1.8.4.tar.bz2
mv bcl2fastq bcl2fastq-1.8.4
cd bcl2fastq-1.8.4
src/configure --prefix=/full/path/to/bclToFastq/bcl2fastq-1.8.4
make
make install
```

run bcltofastq on data

```
##make makefiles #note first run without commands or without --use-bases-mask y151,i8y6,y151 to infer experment
perl /full/path/to/bclToFastq/bcl2fastq-1.8.4/bin/configureBclToFastq.pl --output-dir /full/path/to/fastq --input-dir /full/path/to/basecalldir/Data/Intensities/BaseCalls --force --no-eamss --use-bases-mask y151,i8y6,y151 --sample-sheet /full/path/to/SampleSheetBclToFastq.csv --mismatches 1 &
#run makefiles with 12 cores in background
nohup make -C path/to/fastq -j 12&
```

# integrate the barcode readgroup into reads

---

For the in house generated data from nugene use:

```
$perl -wpe 'if($. % 4 == 1 ){chomp; my @t =split(" "); my @hd= split(":",$t[0]); $hd[2] =  $hd[2] ."_" .$t[2]; $_=join(" ",(join(":", @hd), $t[1], $t[2]))."\n";}
```

or when you have used bcltofastq:

```
perl NugeneMergeFastaFiles.pl randombc.fq.gz read1.fq.gz read2.fq.gz
```

# trim off landing probes

trim of landing probes using the bbduk.sh programme from the bbmap package. [ http://sourceforge.net/projects/bbmap/ ]. This also needs a probe sequence file 'probes_ET1262F_2_182.fasta' you can generate this using the supplied 'probes_ET1262F_2_182.bed' file and [bedtools getfasta](http://bedtools.readthedocs.org/en/latest/content/tools/getfasta.html) (versus ucsc genome). This is the advised config:
```
bbduk.sh in=in.fq.gz out=out.fq.gz ref=probes_ET1262F_2_182.fasta hdist=1 ktrim=r rcomp=f k=31 mink=11 qtrim=r trimq=20 minlen=20 
```

# align

---

```
bwa-mem ....
```

# add read group info

---

```
java -Xmx6g -jar AddOrReplaceReadGroups.jar ...
```

#add read group info for each barcode / reindex

----

```
perl NugeneMergeFastaFiles.pl in.bam out.bam
java -Xmx6g -jar BuildBamIndex.jar INPUT=out.bam

```
