[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.org/ahmedmagds/WhatsGNU.svg?branch=master)](https://travis-ci.org/ahmedmagds/WhatsGNU)
# WhatsGNU
What's Gene Novelty Unit: A Tool For Identifying Proteomic Novelty.
## Introduction
WhatsGNU utilizes the natural variation in public databases to rank protein sequences based on the number of observed exact protein matches (the GNU score) in all known genomes of a certain species & can quickly create whole protein reports.<br/>
WhatsGNU compresses proteins database based on exact match to much fewer number of proteins that differ by at least one amino acid. WhatsGNU will save a copy of the compressed database in two formats; database.txt and database.pickle for faster subsequent uses.<br/>
## Available Databases
Five precompressed databases (.pickle) are available to download and use:

### Ortholog Mode:
1. [_Mycobacterium tuberculosis_](https://www.dropbox.com/sh/8nqowtd4fcf7dgs/AAAdXiqcxTsEqfIAyNE9TWwRa?dl=0) Version: 07/09/2019 (compressed 26,794,006 proteins in 6563 genomes to 434,725 protein variants). 
2. [_Pseudomonas aeruginosa_](https://www.dropbox.com/sh/r0wvoig3alsz7xg/AABPoNu6FdN7zG2PP9BFezQYa?dl=0) Version: 07/06/2019 (compressed 14,475,742 proteins in 4713 genomes to 1,288,892 protein variants).
3. [_Staphylococcus aureus_](https://www.dropbox.com/sh/p292mia4oc99hx6/AACPuv7uoYUkZ1WCBDX0XPSVa?dl=0) Version: 06/14/2019 (compressed 27,213,667 proteins in 10350 genomes to 571,848 protein variants).<br/>

### Big Data basic Mode:
1. [_Salmonella enterica_] Enterobase Version: 08/29/2019 (compressed 975,262,506 proteins in 216,642 genomes to 5,056,335 protein variants).
2. [_Staphylococcus aureus_](https://www.dropbox.com/s/bcs922768tjrwwg/Sau_Staphopia_basic_43914.pickle?dl=0) Staphopia Version: 06/27/2019 (compressed 115,178,200 proteins in 43,914 genomes to 2,228,761 protein variants).

**The five databases are available to download by visiting the link or using the wget command as follows:**

**_S. aureus_ Ortholog**
```
$wget -O Sau.zip https://www.dropbox.com/sh/p292mia4oc99hx6/AACPuv7uoYUkZ1WCBDX0XPSVa?dl=0
$unzip Sau.zip -d WhatsGNU_Sau_Ortholog
```

**_Mycobacterium tuberculosis_ Ortholog**
```
$wget -O TB.zip https://www.dropbox.com/sh/8nqowtd4fcf7dgs/AAAdXiqcxTsEqfIAyNE9TWwRa?dl=0
$unzip TB.zip -d WhatsGNU_TB_Ortholog
```
**_Pseudomonas aeruginosa_ Ortholog**
```
$wget -O Pa.zip https://www.dropbox.com/sh/r0wvoig3alsz7xg/AABPoNu6FdN7zG2PP9BFezQYa?dl=0
$unzip Pa.zip -d WhatsGNU_Pa_Ortholog
```
**_S. aureus_ Staphopia**
```
$wget -O Sau_Staphopia_basic_43914.pickle https://www.dropbox.com/s/bcs922768tjrwwg/Sau_Staphopia_basic_43914.pickle?dl=0
```
**_S. enterica_ Enterobase**
```
$wget -O Senterica_Enterobase_basic_216642.pickle
```
The three Ortholog databases contain all available genomes for the species from GenBank as per version day. To know the genomes included in each database, download [List of Genomes included](https://drive.google.com/file/d/1zJoxYznrsUjrXs5lwSU_lg0o6r8KfSWq/view?usp=sharing). The databases for these 3 Ortholog databases will be updated 3 times per year to include new sequenced genomes.

## WhatsGNU toolbox
1. ### WhatsGNU_get_GenBank_genomes.py
This script downloads genomic fna files or protein faa files from GenBank. 

2. ### WhatsGNU_database_customizer.py
This script customizes the protein faa files from GenBank, RefSeq, Prokka and RAST by adding a strain name to the start of each protein. This script can also customize the strain names for gff file to be used in Roary for pangenome analysis, if the Ortholog mode is going to be used in WhatsGNU. 

3. ### WhatsGNU_main.py
In basic mode, this script ranks protein sequences based on the number of observed exact protein matches (the GNU score) in all known genomes of a particular species. It generates a report for all the proteins in your query in seconds using exact match compression technique. In ortholog mode, the script will additionally link the different alleles of an ortholog group using the clustered proteins output file from Roary or similar pangenome analysis tools. In this mode, WhatsGNU will calculate Ortholog Variant Rarity Index (OVRI) (scale 0-1). This metric is calculated as the number of alleles in an orthologous group that have a GNU score less than or equal to the GNU score of any given allele divided by the sum of GNU scores in the orthologous group. This index represents how unusual a given GNU score is within an ortholog group by measuring how many other protein alleles in the ortholog group have that GNU score or lower. For instance, an allele of GNU=8 in an ortholog group that has 6 alleles with this distribution of GNU scores [300,20,15,8,2,1] will get an OVRI of (8+2+1)/346= 0.03. On the other hand, the allele with GNU=300 will get an OVRI of (300+20+15+8+2+1)/346= 1. An allele with an OVRI of 1 is relatively common regardless of the magnitude of the GNU score, while an allele with OVRI of 0.03 is relatively rare. This index helps distinguish between ortholog groups with high levels of diversity and ortholog groups that are highly conserved.   
 
4. ### WhatsGNU_plotter.py
This script plots:
* Heatmap of GNU scores of orthologous genes in different isolates.
* Metadata distribution bar plot of genes. 
* Histogram of the GNU scores of all proteins in a genome.
* Volcano plot showing genes with a lower average GNU score in one group (case) compared to the other (control). The x-axis is the delta average GNU score (Average_GNU_score_case – Average_GNU_score_control) in the ortholog group. Lower average GNU score in cases will have a negative value on the x-axis (red dots) while lower average GNU score in the control group will have positive value on the x-axis (green dots). The y-axis could be drawn as a -log10(P value) from Mann–Whitney-Wilcoxon test. In this case, lower average GNU score in one group (upper left for case or upper right for control) would be of interest as shown by a significant P value (-log10( P value) > 1.3). The y-axis can also be the average OVRI in the case group for negative values on the x-axis or average OVRI in the control group for positive values on the x-axis.

## Dependencies
* [Python3.x](https://www.python.org/downloads/)
* Blastp (optional for WhatsGNU_main.py and required for WhatsGNU_plotter.py)
* NumPy  (required for WhatsGNU_plotter.py)
* SciPy   (required for WhatsGNU_plotter.py)
* Matplotlib  (required for WhatsGNU_plotter.py)
## Installation
WhatsGNU is a command-line application written in Python3. Simply download and use!
```
$git clone https://github.com/ahmedmagds/WhatsGNU
$cd WhatsGNU/bin
$chmod +x *.py
$pwd 
#pwd will give you a path/to/folder/having/WhatsGNU which you will use in next command
$export PATH=$PATH:/path/to/folder/having/WhatsGNU/bin
```
If you need it permanently, you can add this last line to your .bashrc or .bash_profile. 
## Test
* Type WhatsGNU_main.py -h and it should output help screen.
* Type WhatsGNU_main.py -v and you should see an output like WhatsGNU_main.py 1.0.
## Input
1. database (precompressed (.pickle or .txt) or raw (.faa)).
2. Query protein FASTA file (.faa) or folder of query files.

Optional for _S. aureus_:
Download the CSV file of [CC/ST frequencies](https://drive.google.com/file/d/1PaxWdKAyHOO_pAM0-Knx-6G5HYKoXcKU/view?usp=sharing) in the _S. aureus_ database.
## Usage for WhatsGNU_main.py
### Use precompressed databases
```
$WhatsGNU_main.py -d Sau_Ortholog_10350.pickle -dm ortholog query.faa
or
$WhatsGNU_main.py -d Senterica_Enterobase_basic_216642.pickle -dm basic query.faa
```
You can also use a folder of multiple .faa query files as input (e.g. folder_faa/ has all .faa files to be processed)
```
$WhatsGNU.py -d TB_Ortholog_6563.pickle -dm ortholog folder_faa/
```
### Use precompressed databases with more features
You can assign output folder name using -o instead of default (WhatsGNU_results_timestamp)
```
$WhatsGNU_main.py -d Sau_Staphopia_basic_43914.pickle -dm basic -o output_results_folder query.faa
```
Create a file of each protein with all associated ids from the database (Note: large file (~ 1 Gb for 3000 proteins))
```
$WhatsGNU_main.py -d Pa_Ortholog_4713.pickle -dm ortholog -i -o output_results_folder query.faa
```
Create a file of top 10 genomes with hits
```
$WhatsGNU_main.py -d Sau_Ortholog_10350.pickle -dm ortholog -t query.faa
```
Check how many hits you get from a particular genome in the database (**It has to be used with -t**). The names of the different strains in the databases and their corresponding Genbank strain name and GCA number are available from [List of Genomes included](https://drive.google.com/file/d/1zJoxYznrsUjrXs5lwSU_lg0o6r8KfSWq/view?usp=sharing)
```
$WhatsGNU_main.py -d Sau_Ortholog_10350.pickle -dm ortholog -t -s FDAARGOS_31_GCA_001019015.2_CC8_ query.faa
```
Get Metadata (CC/ST) composition of your hits in the report (**Only for _S. aureus_ and you will need to download [CC/ST frequencies](https://drive.google.com/file/d/1PaxWdKAyHOO_pAM0-Knx-6G5HYKoXcKU/view?usp=sharing) to be used with -e**)
```
$WhatsGNU_main.py -d Sau_Ortholog_10350.pickle -dm ortholog -e metadata_frequencies.csv query.faa
```
Get a fasta (.faa) file of all proteins with GNU score of zero.
```
$WhatsGNU_main.py -d Sau_Ortholog_10350.pickle -dm ortholog -f query.faa
```
**The following options work with -dm ortholog:**

Run blastp on the proteins with GNU score of zero and modify the report with ortholog information.
```
$WhatsGNU_main.py -d Pa_Ortholog_4713.pickle -dm ortholog -b query.faa
```
Get the output report of blastp run (works with -b).
```
$WhatsGNU_main.py -d Pa_Ortholog_4713.pickle -dm ortholog -b -op query.faa
```
Select a blastp percent identity and coverage cutoff values [Default 80], range(0,100).
```
$WhatsGNU_main.py -d Sau_Ortholog_10350.pickle -dm ortholog -b –w 90 –c 50 query.faa
```
Select an OVRI cutoff value [Default 0.045], range (0-1).
```
$WhatsGNU_main.py -d TB_Ortholog_6563.pickle -dm ortholog -ri 0.09 query.faa
```
### Use all features together
```
$WhatsGNU_main.py -d Sau_Ortholog_10350.pickle -dm ortholog -o output_results_folder -i -t -s strain_name -e metadata_frequencies.csv -f -b -op –w 95 –c 40 -ri 0.09 query.faa
```
### Command line options
```
$WhatsGNU_main.py -h
```
```
usage: WhatsGNU_main.py [-h] [-m MKDATABASE | -d DATABASE] [-a] [-j]
                        [-r [ROARY_CLUSTERED_PROTEINS]] [-dm {ortholog,basic}]
                        [-ri [RARITY_INDEX]] [-o OUTPUT_FOLDER] [--force]
                        [-p PREFIX] [-t] [-s STRAINHITS] [-e METADATA] [-i]
                        [-f] [-b] [-op] [-w [PERCENT_IDENTITY]]
                        [-c [PERCENT_COVERAGE]] [-q] [-v]
                        query_faa

WhatsGNU v1.0 utilizes the natural variation in public databases to rank
protein sequences based on the number of observed exact protein matches
(the GNU score) in all known genomes of a particular species. It generates a
report for all the proteins in your query in seconds.

positional arguments:
  query_faa             Query protein FASTA file/s to analyze (.faa)

optional arguments:
  -h, --help            show this help message and exit
  -m MKDATABASE, --mkdatabase MKDATABASE
                        you have to provide path to faa file or a folder of
                        multiple faa files for compression
  -d DATABASE, --database DATABASE
                        you have to provide path to your compressed database
  -a, --pickle          Save database in pickle format [Default only txt file]
  -j, --sql             Save database in SQL format for large Databases
                        [Default only txt file]
  -r [ROARY_CLUSTERED_PROTEINS], --roary_clustered_proteins [ROARY_CLUSTERED_PROTEINS]
                        clustered_proteins output file from roary to be used
                        with -m
  -dm {ortholog,basic}, --database_mode {ortholog,basic}
                        select a mode from 'ortholog' or 'basic' to be used
                        with -d
  -ri [RARITY_INDEX], --rarity_index [RARITY_INDEX]
                        select an ortholog variant rarity index (OVRI) cutoff
                        value in range (0-1)[0.045] for ortholog mode
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        Database output prefix to be created for results
                        (default: timestamped WhatsGNU_results in the current
                        directory)
  --force               Force overwriting existing results folder assigned
                        with -o (default: off)
  -p PREFIX, --prefix PREFIX
                        Prefix for output compressed database (default:
                        WhatsGNU_compressed_database)
  -t, --topgenomes      create a file of top 10 genomes with hits
  -s STRAINHITS, --strainhits STRAINHITS
                        check how many hits you get from a particular
                        strain,it has to be used with -t
  -e METADATA, --metadata METADATA
                        get the metadata composition of your hits, use the
                        metadata_frequency.csv file produced by the WhatsGNU
                        customizer script
  -i, --ids_hits        create a file of each protein with locus_tags (ids) of
                        all hits from the database, large file (~ 1 Gb for
                        3000 pts)
  -f, --faa_GNU_0       get a fasta (.faa) file of all proteins with GNU score
                        of zero
  -b, --blastp          run blastp on the proteins with GNU score of zero and
                        modify the report with ortholog_info, blastp has to be
                        installed
  -op, --output_blastp  get the output report of blastp run, it has to be used
                        with -b
  -w [PERCENT_IDENTITY], --percent_identity [PERCENT_IDENTITY]
                        select a blastp percent identity cutoff value [80],
                        range(0,100)
  -c [PERCENT_COVERAGE], --percent_coverage [PERCENT_COVERAGE]
                        select a blastp percent coverage cutoff value [80],
                        range(0,100)
  -q, --quiet           No screen output [default OFF]
  -v, --version         print version and exit
```
## Output
### Always
**query_WhatsGNU_report_v1.txt** (tab-separated output file)

protein | length | function | sequence | GNU score
------- | ------ | -------- | -------- | ---------
protein_id | 261 | hydrolase | MKVQIYQLP | 1918

Note: If -c option is used for _S. aureus_, CC/ST percentages' columns will be added to the report.

**WhatsGNU_date_time.log** (Log file, e.g. WhatsGNU_v1_20190209_183406.log)

### Optional
Option | File | Description
------ | ---- | -----------
-b | query_WhatsGNU_hits_v1.txt | each protein with all hits_ids from the database,large file (~ 1 Gb for S. aureus)
-t | query_WhatsGNU_topgenomes_v1.txt | top 10 genomes with hits to your query

### Use your own database
#### First time use with unprocessed database (-m with one concatenated (.faa) file of all genomes of a species)
```
$WhatsGNU.py -m database.faa query.faa
```
#### Subsequent uses (-d with the compressed database)
```
$WhatsGNU.py -d database.pickle query.faa
or
$WhatsGNU.py -d database.txt query.faa
```
## Instructions for creating a database
### Simple (works for basic report)
1. Download proteomes of a bacterial species (.faa) in a Directory from GenBank FTP site (ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/)
2. cd Directory
3. gunzip *.faa.gz
4. cat *.faa > database.faa
### Advanced 
1. Download proteomes of a bacterial species (.faa) in a Directory from GenBank FTP site (ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/)
2. cd Directory
3. gunzip *.faa.gz

When Proteome files are downloaded from GenBank and unzipped, they have the protein ids and sequences as following:
```
>protein_1
MSDMF
>protein_2
MRTYZ
```
The protein ids usually have little or no information about the strain so once all proteins mixed in the database, strains' information are lost. Adding strain name to the start of each protein would help to later count the top genomes. Also it will help identify protein origin. Different genomes sometimes have the same strain name so to overcome this problem we recommend using strain_GCA#. This information can be downloaded as an excel sheet from NCBI and then in excel concatenate the two columns of strain name and GCA#. Then add it to each protein in the .faa file to look as following:
```
>strain_GCA_12345.1_protein_1
MSDMF
>strain_GCA_12345.1_protein_2
MRTYZ
```
4. cat *.faa > database.faa

At this point the database is ready to be used in WhatsGNU for the first time with -m as previously explained.
## Requests for creating a database
Requests to process a database for a specific species are welcomed and will be considered
## Bugs
Please submit via the GitHub issues page: https://github.com/ahmedmagds/WhatsGNU/issues
## Software Licence
GPLv3: https://github.com/ahmedmagds/WhatsGNU/blob/master/LICENSE
## Author
Ahmed M. Moustafa: [ahmedmagds](https://github.com/ahmedmagds)<br/>
Twitter: [Ahmed_Microbes](https://twitter.com/Ahmed_Microbes)

