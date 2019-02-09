[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
# WhatsGNU
What's Gene Novelty Unit: A Tool For Identifying Proteomic Novelty.
## Introduction
WhatsGNU utilizes the natural variation in public databases to rank protein sequences based on the number of observed exact protein matches (the GNU score) in all known genomes of a certain species & can quickly create whole protein reports.<br/>
WhatsGNU compresses proteins database based on exact match to much fewer number of proteins that differ by at least one amino acid. WhatsGNU will save a copy of the compressed database in two formats; database.txt and database.pickle for faster subsequent uses.<br/>

Three precompressed databases (.pickle) are available to download and use:
1. [_Mycobacterium tuberculosis_](https://drive.google.com/drive/folders/1U2S6OUVJ6o3Q8dhilj2A97Kj4SHH56gT?usp=sharing) Version: 01/30/2019 (compressed 18,230,371 proteins in 4497 genomes to 443,237). 
2. [_Pseudomonas aeruginosa_](https://drive.google.com/drive/folders/1bZtgzMQWvRnrZ33aq6RAECOZKYSCKylA?usp=sharing) Version: 01/27/2019 (compressed 14,475,742 proteins in 2329 genomes to 872,836).
3. [_Staphylococcus aureus_](https://drive.google.com/drive/folders/1cusXLqOEa2K3XhnCstuuWWVcPhI9deth?usp=sharing) Version: 01/21/2019 (compressed 22,738,456 proteins in 8524 genomes to 565,843).<br/>

These three databases contain all available annotated genomes for the species from GenBank as per version day. The databases for these 3 organisms will be updated 3 times per year to include new sequenced genomes.

## Dependencies
[Python3.x](https://www.python.org/downloads/)<br/>
## Installation
WhatsGNU is a command-line application written in Python3, with no additional dependencies beyond the standard Python3.x package. Simply download and use!
```
$git clone https://github.com/ahmedmagds/WhatsGNU
$cd WhatsGNU
$chmod +x WhatsGNU.py
```
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
Check how many hits you get from a particular genome in the database (**It has to be used with -t**)
```
$WhatsGNU.py -d Sau_012119_database.pickle -t -s USA300_FPR3757_CC8_GCA_000013465.1 query.faa
```
Get MLST CC/ST composition of your hits in the report (**Only for _S. aureus_ and you will need to download [CC/ST database frequencies](https://drive.google.com/file/d/1PaxWdKAyHOO_pAM0-Knx-6G5HYKoXcKU/view?usp=sharing) to be used with -c**)
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
                   query.faa

WhatsGNU v1.0 utilizes the natural variation in public databases to rank
protein sequences based on the number of observed exact protein matches
(the GNU score) in all known genomes of certain species & can quickly create
whole protein reports

positional arguments:
  query.faa             Query protein FASTA file to analyze (.faa)

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
  -t, --tophits         create a file of top 10 genomes with hits
  -s STRAINHITS, --strainhits STRAINHITS
                        check how many hits you get from a particular
                        strain,it has to be used with -t
  -b, --hits            create a file of each protein with all hits from the
                        database,large file (~ 1 Gb for 3000 pts)
  -c [CCST_TYPING], --CCST_typing [CCST_TYPING]
                        get the CC/ST composition of your hits (Note: Works
                        only for Saureus)
  -v, --version         print version and exit
```
## Output
File | Description
------------ | -------------
query_WhatsGNU_report_v1.txt | tab-separated output file
query_WhatsGNU_hits_v1.txt | tab-separated output file
## Instructions for creating a database
put somehting here
## Requests for creating a database
Requests to process a database for a specific species are welcomed and will be considered
## Bugs
Please submit via the GitHub issues page: https://github.com/ahmedmagds/WhatsGNU/issues
## Software Licence
GPLv3: https://github.com/ahmedmagds/WhatsGNU/blob/master/LICENSE
## Author
Ahmed M. Moustafa: [ahmedmagds](https://github.com/ahmedmagds)<br/>
Twitter: [Ahmed_Microbes](https://twitter.com/Ahmed_Microbes)

