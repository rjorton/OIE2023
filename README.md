# Reference Alignment Practical
## [6th Viral Bioinformatics and Genomics Training Course](https://github.com/centre-for-virus-research/CVR-VBG-2023)
* Monday 21st - Friday 25th August 2023
* Glasgow, UK
* [Medical Research Council - University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/)

## Contact

[Richard Orton](https://www.gla.ac.uk/schools/infectionimmunity/staff/richardorton/)   
[MRC-University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/)  
464 Bearsden Road  
Glasgow  
G61 1QH  
UK  
E-mail: Richard.Orton@glasgow.ac.uk  

## Contents

This practical is associated with a lecture on Reference Alignment of High-Throughoput Sequencing (HTS) reads to a reference sequence.

* [0: Overview](#0-overview)
* [1: Setup](#1-setup)
	+ [1.1: Basic read statistics](#11-basic-read-statistics)
* [2: Read Alignment](#2-read-alignment)
	+ [2.1: Indexing the reference sequence](#21-indexing-the-reference-sequence)
	+ [2.2: Aligning the reads to the reference](#22-aligning-the-reads-to-the-reference)
	+ [2.3: Converting SAM to BAM](#23-converting-SAM-to-BAM)
	+ [2.4: Basic alignment statistics](#24-basic-alignment-statistics)
 	+ [2.5: Coverage plot](25-coverage-plot)	 
* [3: Alignment on your own](#3-alignment-on-your-own)
* [4: Extra data](#4-extra-data)
* [5: Assembly Visualisation and Statistics Practical](#5-assembly-visualisation-and-statistics-practical)
	+ [5.1: Setup](#51-setup)
	+ [5.2: Summary Statistics with weeSAM](#52-summary-statistics-with-weeSAM)
	+ [5.3: Coverage plot on your own](#53-coverage-plot-on-your-own)
	+ [5.4: Visualisation with Tablet](#54-visualisation-with-tablet)
 
 
# 0: Overview

**YOU DO NOT NEED TO ENTER THE COMMANDS IN THIS OVERVIEW SECTION!**

In this practical, we will be aligning paired end reads to a reference sequence. Commands that you need to enter into the Terminal window (i.e. the command line) are presented in a box and with a different font, like this:

```
ls
```

Sometimes a command is long and doesn’t fit on a single line on the screen (and screen sizes vary), but it should still be entered as one single line on the computer terminal. 

```
bwa mem -t 4 my_reference_file.fasta my_read_file_1.fastq my_read_file_2.fastq > my_output_file.sam
```

A few Linux tips to remember:

1.	Use the **Tab button** to automatically complete filenames – especially long ones
2.	Use the **Up Arrow** to scroll through your previous commands, it enables you to easily re-run or re-use/change/correct old commands
3.	**Case matters**, the following file names are all different:

```
Myfile.txt
MyFile.txt
MYFILE.txt
myfile.txt
my file.txt
my_file.txt
```

Also watch out for number 1s being confused with lowercase letter L’s, and capital O’s being confused with zeroes

```
l = lower case letter L
1 = number one
O = capital letter O
0 = zero
```

# 1: Setup

**Make sure you are logged into the alpha2 server with MobaXterm.**

In this session, we will be working with two sets of Illumina paired end reads which were simulated from a SARS-CoV-2 genome; these simulated reads were created using ART (Huang et al., 2012: [10.1093/bioinformatics/btr708](10.1093/bioinformatics/btr708)). The goal now is to align these reads to a reference genome sequence, with an ultimate goal of creating a consensus sequence for mutation anlysis.

To start off, you will need to copy the data we need for the practical to your home directory. First change directory (cd) to your home directory

```
cd
```

Then copy (cp) the data folder (-r for recursive - we want the folder and all it's contents) to your current directory (which will be your home directory after the above command was entered):

```
cp -r /home4/VBG_data/Richard .
```

Then change directory to the Sim1 (Simulated sample 1) data folder

```
cd ~/Richard/Sim1/
```

Next, list the contents of the directory so you can see the files we will be working with:

```
ls
```

You should see the FASTQ paired-end read files:

**S1\_R1.fq**  
**S1\_R2.fq**

And also a FASTA reference sequence files:

**sars2_ref.fasta**  

## 1.1: Basic read statistics

We will first use a tool called prinseq to count the number of reads in each file. As these are paired end reads, there should be one read from each read pair in each file – and hence the same number of reads in each file. We will also use prinseq to output statistics on the read lengths, but prinseq itself can do much much more.

```
prinseq-lite.pl -stats_info -stats_len -fastq S1_R1.fq -fastq2 S2_R2.fq
```

***Command breakdown:***

1.	**prinseq-lite.pl** is the name of the program
2.	**-stats\_info** tells prinseq to output basic stats on the reads (number of reads and bases)
3.	**-stats\_len** tells prinseq to output basic stats on read lengths (min, max, mean etc)
4.	**-fastq S1\_R1.fq** the name of the 1st FASTQ file
5.	**-fastq2 S2\_R2.fq** the name of the 2nd FASTQ file in the pair

### Common Issue
* A common issue here is not entering the prinseq command on one line in the terminal - you should only use the enter key at the end of the command to execute it.
* Another common issue is typos - check the command carefully if you get an error - it is likely you have mispelled a file or argument

***
### Questions
**Question 1** – How many reads and bases are in the read files 1 and 2?

**Question 2** – What is the average (mean) length of the reads? 
***

The statistics are split into those for the first FASTQ file of the read pair (e.g. stats\_info, stats\_len, etc) and those for the second FASTQ file of the read pair (e.g. stats\_info2, stats\_len2, etc), and should look a bit like this (but with different numbers!):

```
stats_info	bases		48000000
stats_info	reads		320000
stats_info2	bases		48000000
stats_info2	reads		320000
stats_len	max		150
stats_len	mean		150.00
stats_len	median		150
stats_len	min		150
stats_len	mode		150
stats_len	modeval		320000
stats_len	range		1
stats_len	stddev		0.00
stats_len2	max		150
stats_len2	mean		150.00
stats_len2	median		150
stats_len2	min		150
stats_len2	mode		150
stats_len2	modeval		320000
stats_len2	range		1
stats_len2	stddev		0.00 
```

Paired read files should always have the same number of lines/reads (the ordering of the reads in each file is also critical), so if your two paired files have a different number of reads, something has gone wrong (e.g. filtering/trimming went wrong and corrupted the output, or maybe files from different samples are being used). 
 
# 2: Read Alignment

There are many tools available to align reads onto a reference sequence: bwa, bowtie2, minimap2, bbMap, to name but a few.

We will be using [BWA](http://bio-bwa.sourceforge.net) to align our paired end reads to a reference sequence and output a [SAM (Sequence Alignment Map)](https://samtools.github.io/hts-specs/SAMv1.pdf) file. The SAM file contains the result of each read’s alignment to the given reference sequence. 

## 2.1: Indexing the reference sequence

First, we need to create a BWA index of the reference sequence. Tools such as BWA need to index the sequence first to create a fast lookup (or index) of short sequence seeds within the reference sequence. This enables the tools to rapidly align millions of reads:

```
bwa index sars2_ref.fasta
```

If you list (ls) the contents of the directory, you should see the BWA index files, they will all have the prefix sars2\_ref.fasta, and will have extensions such as **.amb**, **.ann**, **.bwt**, **.pac**, and **.sa**.

```
ls
```

## 2.2: Aligning the reads to the reference

Next, we want to align our reads to the reference sequence using the BWA mem algorithm:

```
bwa mem -t 4 sars2_ref.fasta S1_R1.fq S1_R2.fq > S1.sam
```

***Command breakdown:***

1. **bwa** = the name of the program we are executing
2. **mem** = the BWA algorithm to use (recommended for illumina reads > 70nt)
3. **-t 4** = use 4 computer threads
4. **sars2\_ref.fasta** = the name (and location) of the reference genome to align to
5. **S1\_R1.fq** = the name of read file 1
6. **S2\_R2.fq** = the name of read file 2
7. **>** = direct the output into a file
8. **S1.sam** = the name of the output SAM file to create 

Overall, this command will create an output file called S1.sam in the current directory, which contains the results (in SAM format) of aligning all our reads to the reference sequence sars2\_ref.fasta.

When bwa has finished (and your prompt comes back), check that the SAM file has been created.

```
ls
```

There should now be a file called **S1.sam** in the directory.

### Common issue
A common mistake is not waiting for your previous command to finish, and entering the next command into the terminal before the prompt has returned. You need to wait until the **username@alpha2** command prompt returns before entering the next command - the bwa alignment can sometimes take a few minutes.

## 2.3: Converting SAM to BAM

Typically, a SAM file contains a single line for each read in the data set, and this line stores the alignment result of each read (reference name, alignment location, CIGAR string, the read sequence itself, quality, etc).

SAM files are in a text format (which you can open and view if you like: head S1.sam), but can take up a lot of disk storage space. It is good practice to convert your SAM files to BAM (Binary Alignment Map) files, which are compressed binary versions of the same data, and can be sorted and indexed easily to make searches faster. We will use [samtools](https://samtools.github.io) to convert our SAM to BAM, and sort and index the BAM file:

```
samtools sort S1.sam -o S1.bam
```

```
samtools index S1.bam
```

***Command breakdown:***

1.	The first command tells samtools to **sort** the SAM file, and to also output (**-o**)the sorted data in BAM format to a file called **S1.bam**
3.	We then use samtools to **index** the BAM file S1.bam (indexing [which relies on sorted data] enables faster searches downstream).


There should now be two new files in the directory called: 

**S1.bam** (the BAM file)  
**S1.bam.bai** (the BAM index file) 

Now let’s list (ls) the contents of the directory to check we have our new files, and also check out their sizes:

```
ls -lh
```

***Command breakdown:***
* **-l** tells the list (**ls**) command to give the output in a long list format, whilst the **h** tells it to provide file sizes in a human readable format, this is the 5th column, which will have the size of each file in a format such as 2.5M (M for megabytes) or 9.5G (G for gigabytes).

***

### Questions
**Question 3** – How big is the SAM file compared to the BAM file?

***

**NB:** If your SAM file is 0B (i.e. 0 bytes, empty) then something went wrong with the bwa alignment step, so restart from there. If you SAM file is fine (i.e. >0), but your BAM file is 0B (i.e. empty), then something went wrong with your SAM to BAM conversion so re-do that step. 

We don’t need our original SAM file anymore (as we have the BAM file now) so we remove (rm) the SAM file S1.sam:

```
rm S1.sam
```

## 2.4: Basic alignment statistics

One common thing to check is how many reads have been aligned (or mapped) to the reference, and how many are not aligned (or unmapped). Samtools can report this for us easily, utilising the aligner SAM flags you learnt about in the previous session.

**Reminder:** the 2nd column in the SAM file contains the flag for the read alignment. If the flag includes the number 4 flag in its makeup then the read is unmapped, if it doesn’t include the number 4 in it's makeup then it is mapped.

### Number of unmapped reads
```
samtools view -c -f4 S1.bam
```

***Command breakdown***

1.	**samtools view** = to view the file S1.bam
2.	**–c** = count the read alignments
3.	**–f4** = only include read alignments that do have the unmapped flag 4

### Number of mapped read alignments:
```
samtools view -c -F4 S1.bam
```

***Command breakdown***

1.	**samtools view** = to view the file S1.bam
2.	**–c** = count the read alignments
3.	**–F4** = skip read alignments that contain the unmapped Flag 4 

***
### Questions

**Question 4** – how many reads are mapped to the sars2_ref.fasta genome?

**Question 5** – how many reads are unmapped?
***

Technically, the above command gives the number of mapped read **alignments** not reads. A read could be mapped equally well to multiple positions (one will be called the primary alignment, and others secondary alignments [sam flag 256]), or a read could be split into two parts (e.g. spliced) with one part being the primary alignment and the others supplementary [sam flag 2048]

So to get the true number of mapped reads you need to count only the alignments that do not have flags 4 (unmapped), 256 (not primary), and 2048 (supplementary) = 4 + 256 + 2048 = 2308

### Number of mapped reads

```
samtools view -c -F4 -F256 -F2048 S1.bam
```

or summing up the F flag values together:

```
samtools view -c -F2308 S1.bam
```

For small RNA viruses, secondary and supplementary alignments tend to be rare, but it is important to know the distinction between mapped **reads** and mapped read **alignments**.

## 2.5: Coverage plots

We previously used samtools to count the number of mapped and unmapped reads (using samtools view -c commands), now let’s explore the read mapping in more detail by creating a coverage plot using a tool called weeSAM: (https://github.com/centre-for-virus-research/weeSAM)[https://github.com/centre-for-virus-research/weeSAM]

weeSAM analyses a SAM or BAM file, generates a graphical coverage plot, and reports a range of summary statistics such as:

* **Ref_Name**: The identifier of the reference.
* **Ref_Len**: The length in bases of each reference.
* **Mapped\_Reads**: Number of reads mapped to each reference.
* **Breadth**: The number of sites in the genome covered by reads.
* **%\_Covered**: The percent of sites in the genome which have coverage.
* **Min\_Depth**: Minimum read depth observed.
* **Max\_Depth**: Max read depth observed.
* **Avg\_Depth**: Mean read depth observed.
* **Std\_Dev**: Standard deviation of the mean (Avg_Depth).
* **Above\_0.2_Depth**: Percentage of sites which have greater than 0.2 * Avg_Depth.
* **Above\_1_Depth**: Percentage of sites which are above Avg_Depth.
* **Above\_1.8_Depth**: Percentage of sites which have greater than 1.8 * Avg_Depth.
* **Variation\_Coefficient**: The mean of Std_Dev of the mean.

The Average Depth (Avg_Depth) is perhaps the most important field, along with Breadth which will tell you how much of the genome is covered by aligned reads. But the fields such as Std\_Dev and Above_0.2_Depth can give an indication of the variability in the coverage across the genome.

Let’s run weeSAM on our sample:

```
weeSAM --bam S1.bam --html S1
```

An explanation of this command is:

1.	**weeSAM**: the name of the program we are using
2.	**--bam**: flag to signify input bam file
3.	**S1.bam**: the name of our bam file to analyse
4.	**--html**: flag to signify output html file
5.	**S1**: the name prefix to use for the output

If you list the contents of the directory you should see that a folder called **S1\_html\_results** has been created:

```
ls
```

Inside this folder is a HTML file that we can view in a web browser (like Firefox or Chrome), the HTML file has the summary statistics and coverage plot so lets take a look and open the html file.

```
firefox S1_html_results/S1.html
```

***RJO CHECK - will this launch via MobaXterm or should they download?***

You should see something like this:

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/1b_weesam_summ.png)

***RJO UPDATE FIGURE***

***
### Questions
**Question 9** – what is the average depth of coverage across the SARS-CoV-2 reference genome?
***

Now let’s view the coverage plot by clicking on the hyperlink (blue and underlined) in the Ref_Name column, you should see a coverage plot similar to this:

![](https://github.com/WCSCourses/GCV23/blob/main/modules/ReferenceAlignment/1b_weesam.png)

***RJO UPDATE FIGURE***

The x-axis represents the genome position, whilst the y-axis represents the Depth of Coverage at each genome position. 

**NB:** The reference sequence filename is sars2_ref.fasta, but the actual name of the sequence itself is MN908947.fasta, you can open up the file yourself to check this if you want (head –n1 sars2_ref.fasta).

**Close the weeSAM and Firefox windows before proceeding!**

***RJO CHECK - not if not launcged via MobaXterm***

### Common issue
A common issue here is due to the fact that we have launched firefox from the terminal (wihtout running it background - see advanced linux commands). In order to get our command prompt back (the username@alpha2) we need to close the firefox window down, the prompt should then return.


# 3: Alignment on your own

You now need to use bwa to align the reads for the Sim2 samples to the sars2_ref.fasta reference sequence.

You need to work out the commands yourself based on the previous commands for the Sim1 sample. 

Here is a reminder of the commands you used for Sim1 (S1) which you will need to adapt. 

**NB:** Essentially, you will want to the input FASTQ filenames and the output files from S1 to S2


```
bwa mem -t 4 sars2_ref.fasta S1_R1.fq S1_R2.fq > S1.sam
```
```
samtools sort S1.sam -o S1.bam
```
```
samtools index S1.bam
```
```
rm S1.sam
```
```
samtools view -c -f4 S1.bam
```
```
samtools view -c -F2308 S1.bam
```
```
weeSAM --bam S1.bam --html S1
```

***
### Questions

**Question 6** – how many reads are mapped to the sars2_ref.fasta genome for sample Sim2/S2?

**Question 7** – how many reads are unmapped?
***

# 4: Extra Data

If you are looking for something extra to do, there are additional data sets located in the folder:

### ~/Richard/Ebola/

You will find a set of (gzipped) FASTQ paired end read files, and a reference FASTA sequence to align them to.

The reads are from a patient from the ebola epidemic in West Africa 2014 {Gire et al, 2014} [https://www.ncbi.nlm.nih.gov/pubmed/25214632](https://www.ncbi.nlm.nih.gov/pubmed/25214632)

The reference ebola sequence is from a 2007 outbreak in Democratic Republic of Congo. 

Try aligning the reads to the reference yourself.

### ~/Richard/Noisey/

This is a real HCV sample, but the read quality is quite poor making it quite noisey. Again, two HCV ref sequences are supplied (HCV_1a and HCV_1B). Align the paired end reads to each reference and determine what subtype the sample is by comparing mapping and coverage statistics.

### ~/Richard/Mystery/

This is a mystery sample, combine all the given references sequences into one file using the “cat” command, align the reads to that combined reference and then determine what the virus in the sample is.
 
# 5: Assembly Visualisation with Tablet

[Tablet](https://ics.hutton.ac.uk/tablet/) is a tool for the visualisation of next generation sequence assemblies and alignments. It goes beyond simple coverage plots, and allows you to scroll across the genome, zoom into errors of interests, highlight mutations to the reference, and investigate the assembly.

Tablet requires three files:

1.	A bam file, e.g. 1a.bam
2.	A bam index file, e.g. 1a.bam.bai
3.	A reference sequence file: e.g. 1a\_hcv\_ref.fasta

**Tablet demonstration**

## 6: Consensus and variant calling

In this practical, rather than using our previous simulated HCV samples, we will use some real SARS-CoV-2 samples as people are more likely to be familiar with this virus (and it's ORFs and mutations). These samples have already been cleaned, aligned and trimmed to the Wuhan-Hu-1 reference sequence and have indexed BAM files that are ready for consensus and variant calling. 

### 6.1: Setup and data

First, lets move into the SARS2 data directory:

```
cd ~/Richard/SARS2
```

If you list the contents of this directory you will see
In this session, we will be working on some more Illumina paired end read data. The FASTQ data was downloaded from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser) (ENA), and there are 4 samples in total (the samples are not related to one another), with R1 and R2 FASTQ files for each:

* ERR9105817 - ARTIC primer version 4.1
* ERR9731990 - ARTIC primer version 4.1
* ERR9761275 - ARTIC primer version 4.1
* ERR9788433 - ARTIC primer version 4.1

The primer scheme to use is:

```
~/artic-ncov2019/primer_schemes/nCoV-2019/V4.1/SARS-CoV-2.scheme.bed
```

Your task is to work as a group in the breakout rooms to analyse these samples. Initial read QC (with trim_galore) is not required (but you could add it if you wanted).  You should:

* align the reads to the Wuhan-Hu-1 reference sequence
* Report the number of mapped reads
* Trim the ARTIC primers
* Call a consensus sequence
* Use Pangolin to assign a lineage
* Use SPEAR to call the mutations

This is a flexible session, and a chance to collate all the steps that you have learnt onto a single sample(s).

As a group you could:

* Analyse a sample each and collate the results. As there are only 4 samples (and groups will likely be larger than 4) - multiple people could analyse a single sample and check you get the same results
* Write a bash script to process the sample automatically. Remember all the steps to analyse a sample are the same, it is just the input/output names that are changing. Completed example bash scripts will be uploaded here after the session.

The data in this folder is from a run on an Illumina MiSeq machine. The name of the folder implies it was run on the 3rd July 2020 (200703), the machine ID is M01569, the run ID is 0148_000000000-J53HN, and this was called Batch70 locally within the [Medical Research Council-University of Glasgow Centre for Virus Research](https://www.gla.ac.uk/research/az/cvr/) (CVR) as part of a Covid-19 Genomics UK Consortium ([COG-UK](https://www.cogconsortium.uk)) sequencing run. The samples were sequenced using Version 1 (V1) of the ARTIC [nCoV-2019](https://github.com/artic-network/primer-schemes/tree/master/nCoV-2019) amplicon primers.

There are four samples in this run called:

* CVR2058
* CVR2078
* CVR2092
* CVR2101

If you list the contents of the directory you should see paired end reads (R1.fastq and R2.fastq) for each of the four samples:

```
ls
```

These samples have already been trimmed/filtered using [trim_galore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/), so we do not need to QC them. For information, this is command that was used for each sample:

```
#DO NOT ENTER THE BELOW COMMAND- IT IS JUST FOR INFORMATION 
trim_galore -q 20 --length 50 --paired SAMPLE_R1.fastq SAMPLE_R2.fastq
```
### 1.2.2: Illumina consensus generation walkthrough

We will be using the [Wuhan-Hu-1](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) SARS-CoV-2 isolate as the reference sequence, which has the GenBank accession number [MN908947](https://www.ncbi.nlm.nih.gov/nuccore/MN908947); this is also the sequence used as the basis for the SARS-CoV-2 RefSeq sequence [NC_045512](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512). We will use [bwa](https://github.com/lh3/bwa) to align our reads to the reference.

First we need to index the reference sequence that we will be aligning our reads to. Indexing enables the aligner to quickly look up the best places to start aligning a read by using small sequences call 'seeds' from within the read and looking up if and where those seeds occur in the reference genome using the index.

```
bwa index ~/SARS-CoV-2/MN908947.fasta 
```

This should of created a range of bwa index files (MN908947.fasta.amb/.ann/.bwt/.pac/.sa files), so list (```ls```) the contents of the directory to check:

```
ls ~/SARS-CoV-2/
```

You only need to index a genome once, if you are aligning many samples to the same genome sequence (as we will be on this course) you do not need to re-run the reference index step; but don't confuse this step with the BAM indexing step which does need to be done for each sample.

Now we will align the reads from sample CVR2058 to the reference sequence using bwa:


```
bwa mem -t4 ~/SARS-CoV-2/MN908947.fasta CVR2058_R1.fastq CVR2058_R2.fastq > CVR2058.sam

```
Breaking this command down:

* **bwa** = the name of the program
* **mem** = the name of the bwa algorithm to use (it is recommended for reads > 70bp)
* **-t4** = use 4 computational threads
* **~/SARS-CoV-2/MN908947.fasta** = the path to the reference file (and index)
* **CVR2058_R1.fastq** = FASTQ read pair 1
* **CVR2058_R2.fastq** = FASTQ read pair 2
* **> CVR2058.sam** = direct the output into the file CVR2058.sam rather than to the command line

We now have a SAM file, which we immediately convert into a BAM file as they are compressed and faster to work with downstream. We can sort the file at the same time as converting it:


```
samtools sort -@4 CVR2058.sam -o CVR2058.bam
```
Breaking this command down:

* **samtools** = the name of the program
* **sort** = the name of the function within samtools to use
* **-@4** = use 4 threads
* **CVR2058.sam** = the input file
* **-o CVR2058.bam** =  the output file

We no longer need the SAM file so can delete (```rm```) it:

```
rm CVR2058.sam 

```
Now we need to index the BAM file, this makes downstream analyses faster and many downstream tools will only work if the index file has been created:

```
samtools index CVR2058.bam
```

If we list the contents of the directory we should see the index file with a .bai extension has been created:

```
ls
```


As we have used ARTIC amplicons, each read will typically start and end with a primer sequence. The primer sequence does not come from the sample's viral genome (and the primer may actually be slightly different to the virus) so all primer sequences should be removed from the read ends so they don't interfere with consensus and variant calling. 

To do this we use a tool called [ivar](https://andersen-lab.github.io/ivar/html/manualpage.html) which requires a [BED](https://software.broadinstitute.org/software/igv/BED) file containing the primer coordinates on the reference genome:

```
ivar trim -i CVR2058.bam -b ~/artic-ncov2019/primer_schemes/nCoV-2019/V1/nCoV-2019.bed -p CVR2058_trim.bam
```
**NB:** Use **TAB Completion** to help enter the -b primer bed file!

Breaking this command down:

* **ivar** = the name of the program
* **trim** = the name of the function within ivar to use (trimming primers)
* **-i CVR2058.bam** = the name of the input BAM fie (it is in the current directory)
* **-b ~/artic-ncov2019/primer_schemes/nCoV-2019/V1/nCoV-2019.bed** = the path to primer BED file
* **-p CVR2058_trim.bam** = the name of the output file to create (which will be a primer trimmed version of the input BAM file)

For those who are interested ivar works by soft clipping the primer sequences from the read alignments in the BAM file (rather than actually trimming the reads sequences) based on the alignment coordinates. Soft clipping is a form of trimming that is embedded within the CIGAR string of a read alignment. The CIGAR string is the 6th field of the SAM/BAM file, if you were to examine the BAM file manually you should see lots of 'S' characters in the CIGAR field: see the [CIGAR specification](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information.

We now need to sort and then index this BAM file. Sort:

```
samtools sort -@4 CVR2058_trim.bam -o CVR2058_trim_sort.bam 
```

Rename the file back to CVR2058\_trim.bam as we don't need the unsorted BAM file:

```
mv CVR2058_trim_sort.bam CVR2058_trim.bam
```

Index the BAM:

```
samtools index CVR2058_trim.bam
```

We now have a sorted and indexed BAM file (CVR2058\_trim.bam) that contains our sample's paired end reads (CVR2058\_R1.fastq CVR2058\_R2.fastq) aligned to the Wuhan-Hu-1 (MN908947) reference genome, with the amplicon primer sequences clipped off. So we are now ready to call a consensus sequence:

```
samtools mpileup -aa -A -d 0 -Q 0 CVR2058_trim.bam | ivar consensus -p CVR0258 -t 0.4
```

Breaking this command down, there are two parts:

1. samtools [mpileup](http://www.htslib.org/doc/samtools-mpileup.html) which essentially outputs the base and indel counts for each genome position
	* **-aa** = output data for absolutely all positions (even zero coverage ones)
	* **-A** = count orphan reads (reads whose pair did not map)
	* **-d 0** = override the maximum depth (default is 8000 which is typically too low for viruses)
	* **-Q 0** = minimum base quality, 0 essentially means all the data
2. ivar [consensus](https://andersen-lab.github.io/ivar/html/manualpage.html) - this calls the consensus - the output of the samtools mpileup command is piped '|' directly into ivar
	* -p CVR2058 = prefix with which to name the output file
	* -t 0.4 = the minimum frequency threshold that a base must match to be used in calling the consensus base at a position. In this case, an ambiguity code will be used if more than one base is > 40% (0.4). See [ivar manual]

By default, ivar consensus uses a minimum depth (-m) of 10 and a minimum base quality (-q) of 20 to call the consensus; these defaults can be changed by using the appropriate arguments. If a genome position has a depth less than the minimum, an 'N' base will be used in the consensus sequence by default.

ivar will output some basic statistics to the screen such as:

```
#DO NOT ENTER THIS - IT IS THE IVAR OUTPUT YOU SHOULD SEE
Reference length: 29903
Positions with 0 depth: 121
Positions with depth below 10: 121
```
and you should see our consensus sequence (CVR0258.fa) in the directory:

```
ls
```

which you can view via the command line (again, we will be doing variants in later sessions):

```
more CVR0258.fa 
```

### 1.2.3: Exercise: Generating Illumina consensus sequences yourself

There are three other samples in the Illumina data directory:

* CVR2078
* CVR2092
* CVR2101

You should now choose atleast one sample to create a consensus sequence for yourself by running through the above steps, but adapting them for the next sample (you simply need to change the input read names, and the output file names from CVR2058 to your next sample name). A reminder that the commands used were:

```
bwa mem -t4 ~/SARS-CoV-2/MN908947.fasta CVR2058_R1.fastq CVR2058_R2.fastq > CVR2058.sam
```

```
samtools sort -@4 CVR2058.sam -o CVR2058.bam
```

```
rm CVR2058.sam 
```

```
samtools index CVR2058.bam
```

```
ivar trim -i CVR2058.bam -b ~/artic-ncov2019/primer_schemes/nCoV-2019/V1/nCoV-2019.bed -p CVR2058_trim.bam
```

```
samtools sort -@4 CVR2058_trim.bam -o CVR2058_trim_sort.bam 
```

```
mv CVR2058_trim_sort.bam CVR2058_trim.bam
```

```
samtools index CVR2058_trim.bam
```

```
samtools mpileup -aa -A -d 0 -Q 0 CVR2058_trim.bam | ivar consensus -p CVR0258 -t 0.4
```

**QUESTION** - what is number of mapped reads in each of the samples you have looked at? Hint:

```
samtools view -c -F2308 input.bam
```

Overall, you should again see that we are simply running the same set of commands over and over again for different samples but just changing the input and output names. This is where the power of simple bash scripting and bioinformatics pipelines come into play, as relatively simple scripts can be written to automate this process.

### 1.2.4: Non-amplicon Illumina samples

If the sample was shotgun rather than amplicons, simply omit the ivar trim (and subsequent BAM manipultion) steps:


```
bwa mem -t4 ~/SARS-CoV-2/MN908947.fasta CVR2058_R1.fastq CVR2058_R2.fastq > CVR2058.sam
```

```
samtools sort -@4 CVR2058.sam -o CVR2058.bam
```

```
rm CVR2058.sam 
```

```
samtools index CVR2058.bam
```

```
samtools mpileup -aa -A -d 0 -Q 0 CVR2058.bam | ivar consensus -p CVR0258 -t 0.4
```



## 4: SARS-CoV-2 Group Practical

In this session, we will be working on some more Illumina paired end read data. The FASTQ data was downloaded from the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser) (ENA), and there are 4 samples in total (the samples are not related to one another), with R1 and R2 FASTQ files for each:

* ERR9105817 - ARTIC primer version 4.1
* ERR9731990 - ARTIC primer version 4.1
* ERR9761275 - ARTIC primer version 4.1
* ERR9788433 - ARTIC primer version 4.1

The primer scheme to use is:

```
~/artic-ncov2019/primer_schemes/nCoV-2019/V4.1/SARS-CoV-2.scheme.bed
```

Your task is to work as a group in the breakout rooms to analyse these samples. Initial read QC (with trim_galore) is not required (but you could add it if you wanted).  You should:

* align the reads to the Wuhan-Hu-1 reference sequence
* Report the number of mapped reads
* Trim the ARTIC primers
* Call a consensus sequence
* Use Pangolin to assign a lineage
* Use SPEAR to call the mutations

This is a flexible session, and a chance to collate all the steps that you have learnt onto a single sample(s).

As a group you could:

* Analyse a sample each and collate the results. As there are only 4 samples (and groups will likely be larger than 4) - multiple people could analyse a single sample and check you get the same results
* Write a bash script to process the sample automatically. Remember all the steps to analyse a sample are the same, it is just the input/output names that are changing. Completed example bash scripts will be uploaded here after the session.


## 5: Warnings

I would consider this VM a good place to learn BUT not necessarily a good place to conduct 'real' analyses. The reason being is that many of the SARS-CoV-2 tools and datasets are updated very frequently which means many will be out of date on the VM already (many of the tools were installed a few months ago). Tools such as Pangolin and SPEAR do however have good update functions.


## Glossary

SAM Fields

SAM Flags

FASTQ Scores

Extract unmapped reads

Remove unmapped reads




