# BarCounter
BarCounter is a tool designed to process Next-Generation Sequencing data and count the number of Unique Molecular Identifiers (UMIs) for each specified Antibody Derived Tag (ADT) for each cell barcode specified in a user provided barcode whitelist.  

This tool was designed for use with 10X Genomics v3 3' Gene Expression sequencing data produced via Cell Hashing or CITE-seq protocols. Tag counts for each whitelist cell barcode with observed counts will be written to a comma separated values (.csv) file with one barcode per line.  

Column 1 will contain the barcode sequence, column 2 will contain the total counts from all tags, and all subsequent columns will contain the counts for each tag in the order specified in the taglist.  

A log file with GMT timestamps will be created with all user-displayed messages.  

Outputs will be written to the user specified directory. If the output directory does not exist at the time of the program running, BarCounter will create it.  

Barcounter can be compiled using GCC version 6.3.0 or newer:  
```
gcc Bar_Count.c barcodes.c tags.c umis.c -lz -o barcounter
```

### Definitions:
- *barcode whitelist*: A file with either a .txt or .gz extension that lists all valid cell barcodes with one barcode per line.  
- *taglist*: A comma separated values file (.csv) containing all ADT sequences and names. One tag is listed per line in the format SEQUENCE,Name.  
    - **ex. GTCAACTCTTTAGCG,HT1**

### Arguments:
- `-w`: barcode whitelist  
- `-t`: taglist  
- `-1`: read1 fastq, comma separated list of files (ex. -1 sample1_S1_L001_R1_001.fastq.gz,sample1_S1_L002_R1_001.fastq.gz)  
- `-2`: read2 fastq, comma separated list of files (ex. -2 sample1_S1_L001_R2_001.fastq.gz,sample1_S1_L002_R2_001.fastq.gz)  
- `-o`: output directory  
- `-h`: (optional) This displays a help message with the proper usage. Inclusion of -h will immediately exit the program.  

### Assumptions:
The cell barcode is expected to be 16bp long and begin at the firt base in read1.  
Tag sequences are expected to 15bp long and begin at the first base in read2.  
All tag names are required to be unique.  
All tag sequences are required to have a minimum hamming distance of three from all other tags.  
UMIs are expected to be 12bp long and begin at base 17 in read1.  
Sequence data (read1 and read2) is expected to be in Ilumina standard gzipped fastq format.  
Fastq files are expected to follow Illumina standard naming convention (ex. sample1_S1_L001_R1_001.fastq.gz).  
Fastq file names are underscore delimited: the first field is the sample name, the fourth field is the read number.  
Both sample name and read number must be present in fastq file names for successful completion of BarCounter.  
All input fastq file names must contain the same sample name.  
All input read1 fastq file names must contain "R1", all input read2 fastq file names must contain "R2".  

### Requirements:
Required RAM increases with the number of whitelist barcodes, tags, and UMIs. However, the increase in memory usage is smaller as the size of the inputs increases. For a dataset containing ~40M reads and 30K cells 4 - 5 Gb of memory is usually sufficient.  

BarCounter is single threaded, a single CPU is sufficient.  

### Licensing
All code was written by Elliott Swanson of the Allen Institute for Immunology (elliott.swanson@alleninstitute.org).  

BarCounter is free for academic use as detailed in LICENSE.  

Please refer to the following publication for citations: https://doi.org/10.1101/2020.09.04.283887
