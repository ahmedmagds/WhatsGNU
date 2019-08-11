[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://travis-ci.org/ahmedmagds/WhatsGNU.svg?branch=master)](https://travis-ci.org/ahmedmagds/WhatsGNU)
# WhatsGNU (Under update as of 08/11/2019, check for updated version in next few days)
What's Gene Novelty Unit: A Tool For Identifying Proteomic Novelty.
## Introduction
WhatsGNU utilizes the natural variation in public databases to rank protein sequences based on the number of observed exact protein matches (the GNU score) in all known genomes of a certain species & can quickly create whole protein reports.<br/>
WhatsGNU compresses proteins database based on exact match to much fewer number of proteins that differ by at least one amino acid. WhatsGNU will save a copy of the compressed database in two formats; database.txt and database.pickle for faster subsequent uses.<br/>
## Available Databases
Three precompressed databases (.pickle) are available to download and use:
1. [_Mycobacterium tuberculosis_](https://drive.google.com/drive/folders/1U2S6OUVJ6o3Q8dhilj2A97Kj4SHH56gT?usp=sharing) Version: 01/30/2019 (compressed 18,230,371 proteins in 4497 genomes to 443,237). 
2. [_Pseudomonas aeruginosa_](https://drive.google.com/drive/folders/1bZtgzMQWvRnrZ33aq6RAECOZKYSCKylA?usp=sharing) Version: 01/27/2019 (compressed 14,475,742 proteins in 2329 genomes to 872,836).
3. [_Staphylococcus aureus_](https://drive.google.com/drive/folders/1cusXLqOEa2K3XhnCstuuWWVcPhI9deth?usp=sharing) Version: 01/21/2019 (compressed 22,738,456 proteins in 8524 genomes to 565,843).<br/>

These three databases contain all available annotated genomes for the species from GenBank as per version day. To know the genomes included in each database, download [List of Genomes included](https://drive.google.com/file/d/1zJoxYznrsUjrXs5lwSU_lg0o6r8KfSWq/view?usp=sharing). The databases for these 3 organisms will be updated 3 times per year to include new sequenced genomes.

## Dependencies
[Python3.x](https://www.python.org/downloads/)<br/>
## Installation
WhatsGNU is a command-line application written in Python3, with no additional dependencies beyond the standard Python3.x package. Simply download and use!
```
$git clone https://github.com/ahmedmagds/WhatsGNU
$cd WhatsGNU
$chmod +x WhatsGNU.py
$pwd 
#pwd will give you a path/to/folder/having/WhatsGNU which you will use in next command
$export PATH=$PATH:/path/to/folder/having/WhatsGNU
```
If you need it permanently, you can add this last line to your .bashrc or .bash_profile. 
## Test
* Type WhatsGNU.py -h and it should output help screen.
* Type WhatsGNU.py -v and you should see an output like WhatsGNU.py 1.0.
## Input
1. database (precompressed (.pickle or .txt) or raw (.faa)).
2. Query protein FASTA file (.faa) or folder of query files.

Optional for _S. aureus_:
Download the CSV file of [CC/ST frequencies](https://drive.google.com/file/d/1PaxWdKAyHOO_pAM0-Knx-6G5HYKoXcKU/view?usp=sharing) in the _S. aureus_ database.
## Usage
### Use precompressed databases with basic function
```
$WhatsGNU.py -d TB_013019_database.pickle query.faa
or
$WhatsGNU.py -d Pa_012719_database.pickle query.faa
or
$WhatsGNU.py -d Sau_012119_database.pickle query.faa
```
You can also use a folder of multiple .faa query files as input (e.g. folder_faa has all .faa files to be processed)
```
$WhatsGNU.py -d Sau_012119_database.pickle folder_faa/
```
### Use precompressed databases with more features
You can assign output folder name using -o instead of default (WhatsGNU_results_v1_timestamp)
```
$WhatsGNU.py -d Sau_012119_database.pickle -o output_results_folder query.faa
```
Create a file of each protein with all associated ids from the database (Note: large file (~ 1 Gb for 3000 proteins))
```
$WhatsGNU.py -d Sau_012119_database.pickle -b -o output_results_folder query.faa
```
Create a file of top 10 genomes with hits
```
$WhatsGNU.py -d Sau_012119_database.pickle -t query.faa
```
Check how many hits you get from a particular genome in the database (**It has to be used with -t**). The names of the different strains in the databases and their corresponding Genbank strain name and GCA number are available from [List of Genomes included](https://drive.google.com/file/d/1zJoxYznrsUjrXs5lwSU_lg0o6r8KfSWq/view?usp=sharing)
```
$WhatsGNU.py -d Sau_012119_database.pickle -t -s USA300_FPR3757_CC8_GCA_000013465.1 query.faa
```
Get MLST CC/ST composition of your hits in the report (**Only for _S. aureus_ and you will need to download [CC/ST frequencies](https://drive.google.com/file/d/1PaxWdKAyHOO_pAM0-Knx-6G5HYKoXcKU/view?usp=sharing) to be used with -c**)
```
$WhatsGNU.py -d Sau_012119_database.pickle -c Saureus_CC_ST_names_frequencies_012119.csv query.faa
```
### Use all features together
```
$WhatsGNU.py -d Sau_012119_database.pickle -o results_WhatsGNU -b -t -s USA300_FPR3757_CC8_GCA_000013465.1 -c Saureus_CC_ST_names_frequencies_012119.csv query.faa
```
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
### Command line options
```
$WhatsGNU.py -h
```
```
usage: WhatsGNU.py [-h] [-m MKDATABASE | -d DATABASE] [-o OUTPUT_FOLDER] [-t]
                   [-s STRAINHITS] [-b] [-c [CCST_TYPING]] [-v]
                   query_faa

WhatsGNU v1.0 utilizes the natural variation in public databases to rank
protein sequences based on the number of observed exact protein matches
(the GNU score) in all known genomes of certain species & can quickly create
whole protein reports

positional arguments:
  query_faa             Query protein FASTA file to analyze (.faa)

optional arguments:
  -h, --help            show this help message and exit
  -m MKDATABASE, --mkdatabase MKDATABASE
                        you have to provide path to faa file format to create
                        compressed database in txt and pickle formats
  -d DATABASE, --database DATABASE
                        you have to provide path to your processed database
  -o OUTPUT_FOLDER, --output_folder OUTPUT_FOLDER
                        give name for output folder to be created for results
                        (default: timestamped WhatsGNU_results_v1 in the
                        current directory)
  -t, --topgenomes      create a file of top 10 genomes with hits
  -s STRAINHITS, --strainhits STRAINHITS
                        check how many hits you get from a particular
                        strain,it has to be used with -t
  -b, --hits            create a file of each protein with all hits from the
                        database,large file (~ 1 Gb for 3000 pts)
  -c [CCST_TYPING], --CCST_typing [CCST_TYPING]
                        get the CC/ST composition of your hits (Note: Works
                        only for S.aureus)
  -v, --version         print version and exit
```
## Output
### Always
**query_WhatsGNU_report_v1.txt** (tab-separated output file)

protein | length | function | sequence | GNU score
------- | ------ | -------- | -------- | ---------
protein_id | 261 | hydrolase | MKVQIYQLP | 1918

Note: If -c option is used for _S. aureus_, CC/ST percentages' columns will be added to the report.

**WhatsGNU_v1_date_time.log** (Log file, e.g. WhatsGNU_v1_20190209_183406.log)

### Optional
Option | File | Description
------ | ---- | -----------
-b | query_WhatsGNU_hits_v1.txt | each protein with all hits_ids from the database,large file (~ 1 Gb for S. aureus)
-t | query_WhatsGNU_topgenomes_v1.txt | top 10 genomes with hits to your query
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

