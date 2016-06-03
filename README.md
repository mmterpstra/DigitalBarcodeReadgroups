
# DigitalBarcodeReadgroups
----

Create readgroups from barcodes with random nucleotides incorporates into them for improving sequencing depth. The analysis for marking / removing duplicates does not handle random barcodes/UMIs (**U**nique **M**olecular **I**ndentifiers) and should be modified to include random barcodes/UMIs in the analysis.

This is a git for noting down the steps using simple `perl` tools.

**The following steps are performed**

1. Basecall to fastq
2. Integrate the barcode into reads
3. Trim off landing probes
4. Align
5. Add read group info
6. Add read group info for each barcode / reindex
7. Mark duplicates

**Requirements**

- perl
- an aligner  eg. bwa-7.12/hisat2-2.0.4
- picard-1.140 and up
- this repo
- [pipeline-util](https://github.com/mmterpstra/pipeline-util)
- bedtools
- samtools

**Installing**

Add the `scr/` dir to your path.

>**notes**
>
>This is primairly used for the Nugene Ovation® Target Enrichment System.  [link](http://www.nugen.com/nugen/index.cfm/products/ovation/ovation-target-enrichment-system/)




## Basecall to fastq
----

Use bcl to fastq from illumina. This is an example for questions and documentation see [here](http://support.illumina.com/downloads/bcl2fastq_conversion_software.html). This is not a complete guide for more detailed install instructions visit the illumina website.

###  bcl2fastq 2 software install / run

----

Install bcltofastq(2)

```
mkdir -p source/bclToFastq2/
cd source/bclToFastq/; wget  ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/bcl2fastq/bcl2fastq2-v2.17.1.14.tar.zip
unzip bcl2fastq2-v2.17.1.14.tar.zip
tar -zxf bcl2fastq2-v2.17.1.14.tar.gz
mv bcl2fastq bcl2fastq-v2.17.1.14
cd bcl2fastq-v2.17.1.14
mkdir -p install/bclToFastq/bcl2fastq-v2.17.1.14
./src/configure --prefix=install/bclToFastq/bcl2fastq-v2.17.1.14
make
make install
```

Run bcl2fastq. This is ran on a paired end run starting from the main run folder. Note the `--use-bases-mask y*,i8y*,y*` and the `--minimum-trimmed-read-length 0` options. 

 * `--use-bases-mask y*,i8y*,y*` will output the N6 random/UMI barcodes a second read. For Single End run only specify `--use-bases-mask y*,i8y*`.
 * `--minimum-trimmed-read-length 0` this will prevent the bcl2fastq software from truncating your N6 barcodes/UMI into sequences of NNNNNN because the software is set by default to mask the cycles with N these when they are less then 22nt long.
 * In the Samplesheet.csv it's also advised to turn on the Illumina adapter trimming when running PE sequencing.

```
bcl2fastq \
 --sample-sheet SampleSheet.csv \
 -i Data/Intensities/BaseCalls/ \
 -R ./ \
 --intensities-dir Data/Intensities/ \
 -o Data/Intensities/BaseCalls/run2 \
 --use-bases-mask y*,i8y*,y* \
 --barcode-mismatches 1 \
 --minimum-trimmed-read-length 0 \
 --create-fastq-for-index-reads
```

----

install bcltofastq(1) depreciated

```
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
the i8y6 part says to treat the random part as read2 and the paired ends as read1/3 

```
##make makefiles #note first run without commands or without --use-bases-mask y151,i8y6,y151 to infer experiment
perl /full/path/to/bclToFastq/bcl2fastq-1.8.4/bin/configureBclToFastq.pl --output-dir /full/path/to/fastq --input-dir /full/path/to/basecalldir/Data/Intensities/BaseCalls --force --no-eamss --use-bases-mask y151,i8y6,y151 --sample-sheet /full/path/to/SampleSheetBclToFastq.csv --mismatches 1 &
#run makefiles with 12 cores in background
nohup make -C path/to/fastq -j 12&
```
## Integrate the barcode into reads
---

For the in-house generated data from nugene use:

```
$perl -wpe 'if($. % 4 == 1 ){chomp; my @t =split(" "); my @hd= split(":",$t[0]); $hd[2] =  $hd[2] ."_" .$t[2]; $_=join(" ",(join(":", @hd), $t[1], $t[2]))."\n";}
```

When you have used bcltofastq. you can use these helper scripts from this github. Here are some examples.

```
perl NugeneMergeFastqFiles.pl  -q 20  randombc_R2.fq.gz outdir reads_R1.fastq.gz reads_R3.fastq.gz
perl NugeneMergeFastqFiles.pl   randombc_R2.fq.gz outdir reads_R1.fastq.gz
```

This will add the barcode to the flowcell ID(fcid) like "${fcid}_${barcode}"

## Trim off landing probes
---

trim of landing probes using the bbduk.sh programme from the bbmap package. This is installed by downloading/extracting the data from the [sourceforge project page]( http://sourceforge.net/projects/bbmap/ ). This also needs a probe sequence file 'probes_ET1262F_2_182.fasta' you can generate this using the supplied 'probes_ET1262F_2_182.bed' file and [bedtools getfasta](http://bedtools.readthedocs.org/en/latest/content/tools/getfasta.html) (versus ucsc genome). This is the advised config:

```
bbduk.sh in=in.fq.gz out=out.fq.gz ref=probes_ET1262F_2_182.fasta hdist=1 ktrim=r rcomp=f k=31 mink=11 qtrim=r trimq=20 minlen=20 
```

A better alternative is to use the probe mapping locations (probes.bed) with the read alignment to filter the landing probes. This step will need to have bedtools,samtools installed also see:[pipeline-util](https://github.com/mmterpstra/pipeline-util) and the [trimByBed.pl script](https://github.com/mmterpstra/pipeline-util/blob/master/bin/trimByBed.pl).

```
bwa mem ... reads_R1.fastq.gz > align.sam || hisat2 ... -U reads_R1.fastq.gz -S align.sam
trimByBed.pl -s align.sam  -b probe.bed  -o reads_R1.cleaned.fastq.gz
bbduk.sh in=in.fq.gz out=out.fq.gz qtrim=r trimq=20 minlen=20 
```
## Align
---

Align with your favourite aligner storing the data in sam format for output data. 

```
bwa-mem .... || hisat2 ...
```

## Add read group info
---

Use AddOrReplaceReadGroups from the picard toolkit to add your readgroups. 

```
java -Xmx6g -jar AddOrReplaceReadGroups.jar ...
```

## Add read group info for each barcode / reindex
----
This wil restore your fcid tag + create a readgroup for each random barcode with the correct sample.

```
perl NugeneDigitalSplitter.pl in.bam out.bam
java -Xmx6g -jar picard.jar BuildBamIndex INPUT=out.bam

```

## Mark duplicates
---

This should work. Although the duplicate filteing might not be completely perfect because It does't take into account the Alignment distances of the different UMIs.

```
java -Xmx6g -jar picard.jar MarkDuplicates...
```

and now you can run the rest of your variant calling pipeline!


