# WhatsGNU
WhatsGNU: a tool for identifying proteomic novelty
# Usage
```
$WhatsGNU_v1.0.py -h
```
```
usage: WhatsGNU.py [-h] [-m MKDATABASE | -d DATABASE] [-o OUTPUT_FOLDER] [-t]
                   [-s STRAINHITS] [-b] [-c [CCST_TYPING]] [-v]
                   fastafile_faa

WhatsGNU v1.0 utilizes the naturalvariation in public databases to rank
protein sequences based on the number of observed exact proteinmatches
(the GNU score) in all known genomes of certain species & can quickly create
whole protein reports

positional arguments:
  fastafile_faa         protein FASTA file to analyse (.faa)

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
                        strain,you have to use it with -t
  -b, --hits            create a file of each protein with all hits from the
                        database,large file (~ 1 Gb for 3000 pts)
  -c [CCST_TYPING], --CCST_typing [CCST_TYPING]
                        get the CC/ST composition of your hits (Note: Works
                        only for Saureus)
  -v, --version         print version and exit
  ```
