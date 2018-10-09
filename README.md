# duplexCaller

Circulating free DNA sequencing (cfDNA-Seq) can portray cancer genome landscapes, but highly sensitive and specific technologies are necessary to accurately detect mutations with often low variant frequencies. We developed a customizable hybrid-capture cfDNA-Seq technology using off-the-shelf molecular barcodes and a novel duplex DNA molecule identification tool for enhanced error correction.

Here we provide the duplexCaller software for variant filtering and error correction.

We also provide the program generateBam that generates a BAM file that contains only the reads
tha support the variants of interest as processed by duplexCaller. Helpful utility for debugging and visualization.

## Publication

Title: Ultra-Sensitive Mutation Detection and Genome-Wide DNA Copy Number Reconstruction by Error-Corrected Circulating Tumor DNA Sequencing

Journal: Clinical Chemistry

Published: August 2018

DOI: 10.1373/clinchem.2018.289629

## Dependencies and System Requirements

duplexCaller tool is written in Python 2.The Python compiler used for development is the Python Python 2.7.10 (default, Oct 23 2015, 19:19:21).The program has been developed in a Mac OS computer with El Capitan version 10.11.5. 

The current version works for Mac OS and Fedora systems. 

Users of CentoOS and Ubuntu please contact Dimitrios Kleftogiannis (see below) for more information. The program does not work for Windows systems. 

The program depends on pysam libraries downloaded from:
```
http://pysam.readthedocs.io/en/latest/index.html
```
The program also depends on SAMtools, so please make sure that SAMtools is installed and configured.

You might need to add SAMtools in your path, so after you install SAMtools you might need a command like:
```
PATH=$PATH:/your/path/to/Samtools
```

Please note that the project is UNDER DEVELOPMENT, so changes might appear in next versions.

The current version (13-Sep-2018), handles adaptor read-through cases and discard reads with mapping quality 0.

Please open Execution_examples to see how to use the program. 

#### Releases

17-Sep-2018: Unified version that fixes the bug between Mac and Linux systems.

## Contact

Dr Dimitrios Kleftogiannis

Comments and bug reports are welcome, email to dimitrios DOT kleftogiannis AT icr DOT ac DOT uk 

I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
You are free to modify, extend or distribute this code, as long as our copyright notice is included whole and unchanged. 

## Licence

Copyright 2017 -- Dimitrios Kleftogiannis -- The Institute of Cancer Research (ICR)

Centre for Evolution and Cancer -- Translational Oncogenomics Team & Bioinformatics Team
       			
duplexCaller is freely available software and licensed under the Educational Community License, Version 2.0 (the "License") 

You may not use this program except in compliance with the License. You may obtain a copy of the License at: https://opensource.org/licenses/ECL-2.0




