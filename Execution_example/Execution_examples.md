## FAQ

### Question 1: What input arguments do I need to run duplexCaller  ?

duplexCaller takes as input the following three arguments:

```
	a. a BAM file de-duplicated using MBCs (indexed)

	b. a bed-like file with the positions of interest (only chromosome TAB location)

	c. a reference genome in fasta format (indexed)
```

To make the example more interactive let us assume that we have the following BAM:

```
Example_chr17.bam
```

and few positions to screen similar to the ones in the provided file:
```
example_positions.txt
```

We cannot provide the reference and the BAM due to storage restrictions, but you can find reference genome  elsewhere 

Please contact Dimitrios to obtain the Example_chr17.bam or other toy samples.


### Question 2: How do I run duplexCaller ?

Type the following command:

```
python duplexCaller.py bamFile=Example_chr17.bam positionFile=example_positions.txt referenceGenome=myRef.fa
```

### Question 3: What is the output ?

If you execute the previous example you will see a folder named Example_chr17 

and inside this folder you will find a subfolder named RESULTS and also another subfolder named INTERM_FILES used to store intermediate files.

The intermediate files are used by the generateBam tool for visualization purposes see below.

In this current version we store many intermediate files (such as the VariantDuplexDict.txt files) to better inspect the results.

The duplexCaller results are written in a TAB-limited file named 

```
Example_chr17_VariantReport.txt
```

This report looks like:

```
chr17	37884037	C	4158	0	2077	2081	0	0	0	0,0	REF	1141,375	0,0	0,0	0,0
chr17	37873777	C	2548	0	2545	0	3	0	0	0,0	REF	0,0	1,0	0,0	0,0
chr17	37866005	C	1805	0	0	0	1805	0	0	0,0	REF	0,0	1063,331	0,0	0,0
.....

```

For each position we write the genomic coordinates (col 1,2), the reference allele (col 3), the total reads (col 4) and then the number of reads supporting A,C,G,T as well as cases with DEL and INS (col 5-10). 

The last 5 columns (named for example DistFrag_*,DUP_*), report the number of distinct fragments and duplexes per alternative allele including INS and DEL (comma delimited values)

This information is used further to accept or reject variants based on other criteria (e.g., minimum number of duplexes==2) 


### Question 4: How do I run  generateBAM  ?

First you need to execute duplexCaller and produce the results. 

Then you type something like the following:

```
python generateBam.py bamFile=Example_chr17.bam positionFile=example_positions.txt
```
The program generates a bam file named (from the Example_chr17_variants.bam in our example) in the RESULTS subfolder that contains the reads supporting the variants of interest.


### Question 5: I have problems to run the analysis, what should I do ?

You try to re-produce the example, you try to resolve the dependenices (pysam, samtools, samtools in the PATH) and....you try a bit more in general...

You can always contact Dimitrios and ask for help =)




