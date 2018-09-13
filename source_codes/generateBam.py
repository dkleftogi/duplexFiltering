#!/usr/local/bin/python

'''
					Ultra-sensitive mutation detection and genome-wide DNA copy 
            number reconstruction by error corrected circulating tumor DNA sequencing

BEGIN COPYRIGHT NOTICE

     generateBam code -- (c) 2017 Dimitrios Kleftogiannis -- ICR -- www.icr.ac.uk

     Copyright 2017 Dimitrios Kleftogiannis Licensed under the
     Educational Community License, Version 2.0 (the "License"); you may
     not use this file except in compliance with the License. 

     You may obtain a copy of the License at

     https://opensource.org/licenses/ECL-2.0

     Unless required by applicable law or agreed to in writing,
     software distributed under the License is distributed on an "AS IS"
     BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
     or implied. See the License for the specific language governing
     permissions and limitations under the License.

     Published reports of research using this code (or a modified version) should cite the 
     relevant article of this tool:

    "Mansukhani, Barber et. al Ultra-sensitive mutation detection and genome-wide DNA copy 
    number reconstruction by error corrected circulating tumor DNA sequencing"
     
     Comments and bug reports are welcome.
       
     Email to dimitrios.kleftogiannis@icr.ac.uk 

     I would also appreciate hearing about how you used this code, improvements that you have made to it.
 
     You are free to modify, extend or distribute this code, as long as this copyright notice is included whole and unchanged. 

END COPYRIGHT NOTICE
 
UTILITY
  This program takes as input the output results from duplexCaller.py and generates a BAM file
  with all reads supporting the variants of interests. Since it is a slow process this implementation
  is separate from the duplexCaller code that runs fast for ~50 positions at time.


INPUT ARGUMENTS
    
    1. bam file                     : a bam file from the tumour or plasma of interest
    
    2. bed-like file                : a bed-like file with the positions of interest (chrom TAB position)
                              

DEPENDENCIES
    
    This is a Python program, thus you need the Python 2 compiler to be installed in your computer.
    
    The program has been developed and tested in a Mac OS computer with El Capitan version 10.11.5
    
    The Python compiler used for development is the Python Python 2.7.10 (default, Oct 23 2015, 19:19:21)
    
    The program works for Unix-like systems but has not been tested for Windows operating systems. 
    
    The program has not been tested in Cygwing-like systems running under Windows operating systems.   

    The program depends on pysam libraries downloaded from http://pysam.readthedocs.io/en/latest/index.html

    The program also depends on samtools, so please make sure that SAMtools is installed and configured properly in your system

    You might need to add samtools in your path so after you intall SAMtools you might need a command like: 

    PATH=$PATH:/your/path/to/Samtools


RUNNING
	
	An execution example is as follows:

    python generateBam.py bamFile=Example_chr17.bam positionFile=example_positions.txt

    Please remember that you need to run the duplexCaller first. 

    This program reads the output of duplexCaller

'''

#modules we need, I might have some extra that didnt use.
#Remember that this program is under-developement so you may find block of codes used for testing.
import sys
import os
import pysam
import re
from collections import defaultdict
from itertools import groupby
import datetime
import time
import threading

exitFlag = 0

#prints information about program's execution
def printUsage():
    print('To run this program please type the following:')
    print('\tpython generateBam.py bamFile=file.bam positionFile=file.txt\n')
    print('Where:\n') 
    print('\tfile.bam is a bam file.\n')
    print('\tfile.txt is a bed-like file with your input genomic positions. Please write chromosome and start only\n')
    print('Execution example:\n') 
    print('\tpython generateBam.py bamFile=Example_chr17.bam positionFile=example_positions.txt\n')
    print('\tThe input bam need to be indexed (produce files .bai with samtools) \n')
    print('\tThis program requires execution of duplexCaller with the same input parameters.\n')
    print('Please give the arguments in the indicated order similar to the provided example!\n') 

class myThread (threading.Thread):
    def __init__(self, threadID, name, chrom,bamFilePrefix, bamFile,SAM_FILES):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.chrom = chrom
        self.fileName = bamFilePrefix
        self.inputBam = bamFile
        self.myDir = SAM_FILES
    def run(self):
        myStr="\t\tStarting " + self.name +" to process "+ self.chrom
        print('%s\n'%myStr) 
        processChrom(self.name,self.chrom,self.fileName, self.inputBam,self.myDir)

#this is the actual command
def processChrom(threadName, chrom, fileName, inputBam, myDir):
   if exitFlag:
      threadName.exit()
   #print "%s: %s" % (threadName, time.ctime(time.time()))
   command='samtools view '+inputBam+' '+chrom+' > '+myDir+'/'+fileName+'_'+chrom+'.sam' 
   os.system(command)
   #print(command)

#save the genomic positions of interest; remember this is not the panel design
def storePositionsFile(bedFile):
    #save the positions
    aDict={}
    #check if bed file exists
    if os.path.exists(bedFile):
        #the file exists
        #open the file and read it line by line
        fileIN=open(bedFile,'r')
        #read the first line and store it to the dictionary
        for eachLine in fileIN:
            line = eachLine.rstrip('\n')
            tmp=line.split("\t")
            chrom=tmp[0]
            startPos=tmp[1]
            key=chrom+'_'+startPos
            aDict[key]=int(startPos)
            #print repr(key)
        fileIN.close()
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] ERROR from function: storePositionsFile. The bed file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()
    return aDict

#save the genomic positions of interest; remember this is not the panel design
def storeChromosomes(bedFile):
    #save the positions
    aList=[]
    #check if bed file exists
    if os.path.exists(bedFile):
        #the file exists
        #open the file and read it line by line
        fileIN=open(bedFile,'r')
        #read the first line and store it to the dictionary
        for eachLine in fileIN:
            line = eachLine.rstrip('\n')
            tmp=line.split("\t")
            chrom=tmp[0]
            startPos=tmp[1]
            key=chrom
            aList.append(key)
            #print repr(key)
        fileIN.close()
    else:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] ERROR from function: storeChromosomes. The bed file does not exist!\n'%(st))
        print('************************************************************************************************************************************\n')
        sys.exit()
    #remove duplicates from the list
    return list(set(aList))


def generateMergedSAM(chromList,bamFilePrefix,SAM_FILES,RESULTS):

    outFileName=SAM_FILES+'/merged_variants.sam'
    outFile=open(outFileName,'w')
    a=0
    #process each chromosome separately
    for chrom in chromList:
        #first store all QNAME found
        print('\tProcessing: %s\n'%chrom)
        currentQNAMEfile=SAM_FILES+'/'+bamFilePrefix+'_'+chrom+'_QNAME.sam'
        #open it, read it and store the QNAMES
        fileIN=open(currentQNAMEfile,'r')
        qnameDict=defaultdict(list)
        for eachLine in fileIN:
            line = eachLine.rstrip('\n')
            tmp=line.split("\t")
            QNAME=tmp[0]
            qnameDict[QNAME].append(eachLine)
        #at this point we have all qnames so we just need to parse the big sam...
        fileIN.close()

        currentSAMfile=SAM_FILES+'/'+bamFilePrefix+'_'+chrom+'.sam'
        fileIN=open(currentSAMfile,'r')
        for eachLine in fileIN:
            line = eachLine.rstrip('\n')
            tmp=line.split("\t")
            prompt_QNAME=tmp[0]
            recordsFound=qnameDict[prompt_QNAME]
            if len(recordsFound)==0:
                a=a+1
            else:
                outFile.write('%s\n'%eachLine)
        fileIN.close()
        del qnameDict
    outFile.close()

#erase the big sam files to save some space...
def cleanFiles(SAM_FILES):
    command='rm '+SAM_FILES+'/*.sam'
    os.system(command)

#main function of the program
def myMain():
    #check the number of input arguments
    if len(sys.argv)!=3:
        print('************************************************************************************************************************************\n')
        print('\t\t\t\t\tYour input arguments are not correct!\n')
        print('\t\t\t\tCEC Bioinformatics & CEC Translational Oncogenomics Team\n')
        print('\t\t\tCopyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n')
        printUsage()
    else:
        print('************************************************************************************************************************************\n')
        print('\t\tgenerateBam: Generate a BAM with all reads supporting the variants. Helpful for visualizing of the results.\n')
        print('\t\t\t\t\tCEC Bioinformatics & CEC Translational Oncogenomics Team\n')
        print('\t\tCopyright 2017 ICR -- Dimitrios Kleftogiannis -- dimitrios.kleftogiannis@icr.ac.uk\n')
        #parse the first input arguments 
        #here if the user does not write the correct argument name it gets an error and the program stops
        bamFile=sys.argv[1].split('bamFile=')
        bamFile=bamFile[1]
        #parse the second argument
        bedFile=sys.argv[2].split('positionFile=')
        bedFile=bedFile[1]
        #print the arguments given by user; is good for 'self' debugging
        print('Execution started with the following parameters:\n')
        print('1. bamFile:                  \t\t\t\t%s' % bamFile)
        print('2. bedFile:                  \t\t\t\t%s' % bedFile)
        #save the bam file prefix for further naming of files
        #save the bam file prefix for further naming of files
        bamFilePrefix=bamFile.split('/')
        bamFilePrefix=bamFilePrefix[-1]
        bamFilePrefix=bamFilePrefix[:-4]
        
        #save the bam file prefix for further naming of files
        #generate folders for storing the final results and the intermediate results
        RESULTS=bamFilePrefix+'/RESULTS'
        SAM_FILES=bamFilePrefix+'/INTERM_FILES'
        command='mkdir -p '+bamFilePrefix+' && mkdir -p '+RESULTS+' && mkdir -p '+SAM_FILES
        os.system(command)

        #store the positions for testing
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n[%s] Function storePositionsFile: store the input genomic positions'%(st))
        posList=storePositionsFile(bedFile)
        chromList=storeChromosomes(bedFile)

        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('\n[%s] Process reads supporting variants per chromosome: Releasing threads...'%(st))
        #release the threads
        threadLock = threading.Lock()
        threads = []
        i=1
        for chrom in chromList:
        # Create new threads and start them
            myStr='Thread-'+str(i)
            thread = myThread(i, myStr,chrom,bamFilePrefix,bamFile,SAM_FILES)
            thread.start()
            threads.append(thread)
            i=i+1
        #find your threads
        for t in threads:
            t.join()
        #and collect them
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] Process reads supporting variants per chromosome: All threads collected...'%(st))   

        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] Generating one SAM file containing only the variants.'%(st))
        generateMergedSAM(chromList,bamFilePrefix,SAM_FILES,RESULTS)


        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] Generating a bam file containing only the variants.'%(st))
        #simply merge the results and produce the bam file with the variants
        command='samtools view -H '+bamFile+' > '+SAM_FILES+'/header.sam && sort '+SAM_FILES+'/merged_variants.sam | uniq >'+SAM_FILES+'/_merged_variants_1.sam && cat '+SAM_FILES+'/header.sam '+SAM_FILES+'/_merged_variants_1.sam | samtools sort > '+SAM_FILES+'/_merged_variants_2.sam && samtools view -b '+SAM_FILES+'/_merged_variants_2.sam > '+RESULTS+'/'+bamFilePrefix+'_variants.bam && samtools index '+RESULTS+'/'+bamFilePrefix+'_variants.bam'
        print(command)
        os.system(command)

       #clear the intermediate files
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        print('[%s] Clear all intermediate files produced'%(st))
        #cleanFiles(SAM_FILES)

        print('************************************************************************************************************************************\n')

#this is where we start
if __name__=='__main__':
    myMain()




