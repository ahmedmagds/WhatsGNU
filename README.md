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
2. [_Pseudomonas aeruginosa_](https://www.dropbox.com/sh/r0wvoig3alsz7xg/AABPoNu6FdN7zG2PP9BFezQYa?dl=0) Version: 07/06/2019 (compressed 14,475,742 proteins in 4712 genomes to 1,288,892 protein variants).
3. [_Staphylococcus aureus_](https://www.dropbox.com/sh/p292mia4oc99hx6/AACPuv7uoYUkZ1WCBDX0XPSVa?dl=0) Version: 06/14/2019 (compressed 27,213,667 proteins in 10350 genomes to 571,848 protein variants).<br/>

### Big Data basic Mode:
1. [_Salmonella enterica_](https://www.dropbox.com/s/gbjengikpynxo12/Senterica_Enterobase_basic_216642.pickle?dl=0) Enterobase Version: 08/29/2019 (compressed 975,262,506 proteins in 216,642 genomes to 5,056,335 protein variants).
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
The three Ortholog databases contain all available genomes for the species from GenBank as per version day. To know the genomes included in each database, download [List of Genomes included](https://www.dropbox.com/s/w7z6htvot8167ep/List_of_genomes_included_092019.xlsx?dl=0). The databases for these 3 Ortholog databases will be updated 3 times per year to include new sequenced genomes.

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
* Metadata distribution bar plot of proteins. 
* Histogram of the GNU scores of all proteins in a genome.
* Volcano plot showing proteins with a lower average GNU score in one group (case) compared to the other (control). The x-axis is the delta average GNU score (Average_GNU_score_case – Average_GNU_score_control) in the ortholog group. Lower average GNU score in cases will have a negative value on the x-axis (red dots) while lower average GNU score in the control group will have positive value on the x-axis (green dots). The y-axis could be drawn as a -log10(P value) from Mann–Whitney-Wilcoxon test. In this case, lower average GNU score in one group (upper left for case or upper right for control) would be of interest as shown by a significant P value (-log10( P value) > 1.3). The y-axis can also be the average OVRI in the case group for negative values on the x-axis or average OVRI in the control group for positive values on the x-axis.

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

## Usage for WhatsGNU_main.py
### Input
1. database (precompressed (.pickle or .txt) or raw (.faa)).
2. Query protein FASTA file (.faa) or folder of query files.

Optional for _S. aureus_:
The CSV file of Metadata (CC/ST) frequencies for the _S. aureus_ database.
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
Check how many hits you get from a particular genome in the database (**It has to be used with -t**). The names of the different strains in the databases and their corresponding Genbank strain name and GCA number are available from [List of Genomes included](https://www.dropbox.com/s/w7z6htvot8167ep/List_of_genomes_included_092019.xlsx?dl=0)
```
$WhatsGNU_main.py -d Sau_Ortholog_10350.pickle -dm ortholog -t -s FDAARGOS_31_GCA_001019015.2_CC8_ query.faa
```
Get Metadata (CC/ST) composition of your hits in the report (**Only for _S. aureus_ and you will need to use the metadata_frequencies.csv file (available to download with the database) with -e**)
```
$WhatsGNU_main.py -d Sau_Ortholog_10350.pickle -dm ortholog -e metadata_frequencies.csv query.faa
```
Get a fasta (.faa) file of all proteins with GNU score of zero.
```
$WhatsGNU_main.py -d Sau_Ortholog_10350.pickle -dm ortholog -f query.faa
```
### The following options work with -dm ortholog

Run blastp on the proteins with GNU score of zero and modify the report with ortholog information.
```
$WhatsGNU_main.py -d Pa_Ortholog_4713.pickle -dm ortholog -b query.faa
```
**Note:** If -b is used, WhatsGNU will search for compressed_db_orthologs.faa and compressed_db_orthologs_info.txt in the same path for the compressed database as they are needed for the blastp run.

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
### Output
#### Always with -m or -d
**query_WhatsGNU_report_v1.txt** (tab-separated output file)
##### Basic Mode

protein | GNU score | length | function | sequence |
------- | --------- | ------ | -------- | -------- |
strain_x_protein_1 | 2 | 3 | argG | MVM | 
##### Ortholog Mode (in addition to the previous five columns)

ortholog_group |	ortho_gp_total_sequences_number |	ortho_gp_total_variants_number	| minimum_GNU	| maximum_GNU	| average_GNU	| OVRI	| OVRI interpretation |
-------------- | ------------------------------- | ------------------------------ | ----------- | ----------- | ----------- | ---- | ------------------- |
argG | 100 | 5 | 2 | 50 | 38 | 0.02 | rare |

Explanation for the columns in the report:
For instance, if strain_x_protein_1 (sequence: MVM) belongs to argG orthologous group which has 5 protein variants (MMMM,MVVM, MVM, MVV and VVM) with GNU scores [50,35,10,3,2]:
* Column 1: protein name
* Column 2: GNU score (number of exact matches in the database)
* Column 3: protein sequence length
* Column 4: function from the database
* Column 5: protein sequence
* Column 6: name of the orthologous group
* Column 7: total number of sequences (sum of GNU scores) in the orthologous group
* Column 8: Number of protein variants (alleles) in the orthologous group
* Column 5: minimum GNU score in the orthologous group
* Column 6: maximum GNU score in the orthologous group
* Column 7: average GNU score in the orthologous group
* Column 8: Ortholog Variant Rarity Index (OVRI) (scale is 0-1) which is<br/>
(GNU score of the allele + lower GNU scores)/(Sum of GNU scores in the ortholog group)<br/>
For example for the variant that has GNU=2 it will be 2/100 = 0.02<br/>
while for variant that has GNU=10 it will be (10+3+2)/100 = 0.15<br/>
and finally the variant that has GNU=50 it will be (50+35+10+3+2)/100 = 1
* Column 9: A rare/frequent tag to the protein based on its OVRI. The default cutoff value which is arbitrary is 0.045 so anything below this value is rare and above is frequent.      

Note: If -e option is used for _S. aureus_, CC/ST percentages' columns will be added to the report.

**WhatsGNU_date_time.log** (Log file, e.g. WhatsGNU_v1_20190209_183406.log)

#### Always with -m
* compressed_db.txt (if -a, compressed_db.pickle will be created)
* compressed_db_orthologs.faa (if “-r clustered_proteins” is used with -m)
* compressed_db_orthologs_info.txt (if “-r clustered_proteins” is used with -m)

#### Optional
Option | File | Description
------ | ---- | -----------
-i | query_WhatsGNU_hits.txt | each protein with all hits_ids from the database,large file (~ 1 Gb for S. aureus)
-t | query_WhatsGNU_topgenomes.txt | top 10 genomes with hits to your query
-f | query_WhatsGNU_zeros.faa | file of all proteins with GNU score of zero
-op | query_WhatsGNU_zeros_blast_report.txt | output report of blastp run

## Usage for WhatsGNU_plotter.py
### Input
A folder of query_WhatsGNU_report.txt files.<br/>

### Heatmap
Plot a heatmap of GNU scores for these proteins in proteins.faa using this strains’ order. Assign a title using -t. Font size and figure size (w,h) are given by -f and -fs, respectively. Annotate the heatmap cells with OVRI rare tag using -r option.
```
WhatsGNU_plotter.py -hp ortholog -q proteins.faa -r -d strains_order.txt -t title -r -f 14 -fs 14 10 prefix_name WhatsGNU_reports_folder/
```

### Metadata percentage distribution
Plot a metadata percentage bar plot for the GNU scores of the proteins in proteins.faa for each WhatsGNU report.
```
WhatsGNU_plotter.py -mb basic -q proteins.faa prefix_name WhatsGNU_reports_folder/
```
### Histogram
Plot a blue histogram of the GNU scores for each WhatsGNU report using 100 bins and get a text file showing novel and conserved proteins with -p option to assign cutoffs.
```
WhatsGNU_plotter.py -x -e blue -b 100 -p 50 5000 prefix_name WhatsGNU_reports_folder/
```
### Volcano plot
Plot two scatterplots that shows either statistical significance (P value) or average OVRI versus magnitude of change (Delta_average_GNU_Score). The case/control tag is provided in isolates_case_control_tag.csv. The option -c 100 is a percentage of isolates a protein must be in to be included. A summary statistics file is also created.
```
WhatsGNU_plotter.py -st isolates_case_control_tag.csv -c 100 prefix_name WhatsGNU_reports_folder/
```
### All features together
```
WhatsGNU_plotter.py -hp ortholog -q proteins.faa -d strains_order.txt -t title -r -f 16 -fs 14 10 -mb ortholog -x -e blue -b 100 -st isolates_case_control_tag.csv -c 100 prefix_name WhatsGNU_reports_folder/
```
### Command line options
```
$WhatsGNU_plotter.py -h
usage: WhatsGNU_plotter.py [-h] [-hp {ortholog,basic}] [-l LIST_GENES]
                           [-q FASTA] [-op] [-d STRAINS_ORDER] [-r]
                           [-rc RARITY_COLOR] [-fs FIGURE_SIZE FIGURE_SIZE]
                           [-hc HEATMAP_COLOR] [-mc MASKED_COLOR]
                           [-f FONT_SIZE] [-t TITLE] [-mb {ortholog,basic}]
                           [-w] [-s SELECT_METADATA] [-x] [-e HISTOGRAM_COLOR]
                           [-b HISTOGRAM_BINS]
                           [-p NOVEL_CONSERVED NOVEL_CONSERVED]
                           [-st STRAINS_TAG_VOLCANO] [-c CUTOFF_VOLCANO]
                           [-cc CASE_CONTROL_NAME CASE_CONTROL_NAME]
                           prefix_name directory_path

WhatsGNU_plotter script for WhatsGNU v1.0.

positional arguments:
  prefix_name           prefix name for the the output folder and
                        heatmap/volcano output files
  directory_path        path to directory of WhatsGNU reports

optional arguments:
  -h, --help            show this help message and exit
  -hp {ortholog,basic}, --heatmap {ortholog,basic}
                        heatmap of GNU scores for orthologous genes in
                        multiple isolates
  -l LIST_GENES, --list_genes LIST_GENES
                        a txt file of ortholog group names from one of the
                        WhatsGNU reports for heatmap
  -q FASTA, --fasta FASTA
                        a FASTA file of sequences for the proteins of interest
                        for heatmap or metadata barplot
  -op, --output_blastp  get the output report of blastp run, it has to be used
                        with -q
  -d STRAINS_ORDER, --strains_order STRAINS_ORDER
                        list of strains order for heatmap
  -r, --rarity          Annotate heatmap cells with OVRI(default: off)
  -rc RARITY_COLOR, --rarity_color RARITY_COLOR
                        OVRI data text color in the heatmap
  -fs FIGURE_SIZE FIGURE_SIZE, --figure_size FIGURE_SIZE FIGURE_SIZE
                        heatmap width and height in inches w,h, respectively
  -hc HEATMAP_COLOR, --heatmap_color HEATMAP_COLOR
                        heatmap color
  -mc MASKED_COLOR, --masked_color MASKED_COLOR
                        missing data color in heatmap
  -f FONT_SIZE, --font_size FONT_SIZE
                        heatmap font size
  -t TITLE, --title TITLE
                        title for the heatmap [Default:WhatsGNU heatmap]
  -mb {ortholog,basic}, --metadata_barplot {ortholog,basic}
                        Metadata percentage distribution for proteins in a
                        FASTA file
  -w, --all_metadata    all metadata
  -s SELECT_METADATA, --select_metadata SELECT_METADATA
                        select some metadata
  -x, --histogram       histogram of GNU scores
  -e HISTOGRAM_COLOR, --histogram_color HISTOGRAM_COLOR
                        histogram color
  -b HISTOGRAM_BINS, --histogram_bins HISTOGRAM_BINS
                        number of bins for the histograms [10]
  -p NOVEL_CONSERVED NOVEL_CONSERVED, --novel_conserved NOVEL_CONSERVED NOVEL_CONSERVED
                        upper and lower GNU score limits for novel and
                        conserved proteins novel_GNU_upper_limit,
                        conserved_GNU_lower_limit, respectively [Default 10,
                        100]
  -st STRAINS_TAG_VOLCANO, --strains_tag_volcano STRAINS_TAG_VOLCANO
                        a csv file of the strains of the two groups to be
                        compared with (case/control) tag
  -c CUTOFF_VOLCANO, --cutoff_volcano CUTOFF_VOLCANO
                        a percentage of isolates a protein must be in [Default:
                        100]
  -cc CASE_CONTROL_NAME CASE_CONTROL_NAME, --case_control_name CASE_CONTROL_NAME CASE_CONTROL_NAME
                        case and control groups' names [Default: case control]
```
### Output
A heatmap, metadata percentage distribution bar plot, histogram and two volcano plots and summary statistics files.

## Instructions for creating a database
### Simple
1. Download proteomes of a species (.faa) in a Directory from GenBank
```
$WhatsGNU_get_GenBank_genomes.py -f GCAs.txt Species_faa
```
2. Modify the faa files to have the strains' names
```
$WhatsGNU_database_customizer.py -c -g Species_modified Species_faa/
```
3. Run WhatsGNU_main.py in basic mode
```
$WhatsGNU_main.py -m Species_modified_concatenated.faa query.faa
```
### Advanced (e.g. _S. aureus_)
1. Download genomes of a species (.fna) in a Directory from GenBank  
```
$WhatsGNU_get_GenBank_genomes.py -c GCAs.txt Sau_fna
$gunzip Sau_fna/*
```
2. Annotate the genomes using [Prokka](https://github.com/tseemann/prokka) 
**An example command for _S. aureus_ is given, change it or use any other options from Prokka**
```
$for i in `cat file_names.list`;do prokka --kingdom Bacteria --outdir prokka_$i --gcode 11 --genus Staphylococcus --species aureus --strain $i --prefix $i --locustag $i Species_fna/$i*.fna; done
$find ./ -name '*.faa' -exec cp -prv '{}' '/Sau_faa/' ';'
$find ./ -name '*.gff' -exec cp -prv '{}' '/Sau_gff/' ';'
```
3. Modify the faa and gff files to have the strains' names
```
$WhatsGNU_database_customizer.py -c -p -l strain_name_list.csv Sau_modified_faa Sau_faa/
$WhatsGNU_database_customizer.py -i -s -l strain_name_list.csv -g Sau_modified_gff Sau_gff/
```
The strain_name_list.csv is a comma-separated list of 3+ columns: file_name, old locustag, new locustag and optionally metadata. If metadata are provided, the script will concatenate the new locustag with metadata using ‘_’ as a separator. The new locustag in this case will be:  new_locustag_metadata_. In case of GenBank, RefSeq and RAST, use NA for the old locustag column in the list.csv file. 

4. Run [Roary](https://sanger-pathogens.github.io/Roary/) for pangenome analysis
**An example command for Roary is given, change it or use any other options from Roary**
```
$roary Sau_modified_gff/*.gff
```
5.Run WhatsGNU_main.py in Ortholog mode using clustered_proteins output file from Roary
```
$WhatsGNU_main.py -m Sau_modified_concatenated.faa -r clustered_proteins query.faa
```
## Command line options for WhatsGNU_get_GenBank_genomes.py
```
$WhatsGNU_get_GenBank_genomes.py -h
usage: WhatsGNU_get_GenBank_assemblies.py [-h] [-f] [-c] [-r]
                                          list output_folder

Get GenBank assemblies (faa or/and fna) for WhatsGNU v1.0

positional arguments:
  list           a list.txt file of GenBank accession numbers (GCA#.#)
  output_folder  give name for output folder to be created

optional arguments:
  -h, --help     show this help message and exit
  -f, --faa      protein faa file from GenBank
  -c, --contigs  genomic fna file from GenBank
  -r, --remove   remove assembly_summary_genbank.txt after done
```
## Command line options for WhatsGNU_database_customizer.py
```
$WhatsGNU_database_customizer.py -h
usage: WhatsGNU_database_customizer.py [-h] [-g | -p | -r | -s] [-z]
                                       [-l LIST_CSV] [-i] [-c]
                                       prefix_name directory_path

Database_customizer script for WhatsGNU v1.0.

positional arguments:
  prefix_name           prefix name for the output folder and the one
                        concatenated modified file
  directory_path        path to directory of faa, RAST txt or gff files

optional arguments:
  -h, --help            show this help message and exit
  -g, --GenBank_RefSeq  faa files from GenBank or RefSeq
  -p, --prokka          faa files from Prokka
  -r, --RAST            spreadsheet tab-separated text files from RAST
  -s, --gff_file        gff file from prokka, needed if planning to run Roary
  -z, --gzipped         compressed file (.gz)
  -l LIST_CSV, --list_csv LIST_CSV
                        a file.csv of 3+ columns: file_name, old locustag, new
                        locustag and optionally metadata
  -i, --individual_files
                        individual modified files
  -c, --concatenated_file
                        one concatenated modified file of all input files
```
## Requests for creating a database
Requests to process a database for a specific species are welcomed and will be considered
## Bugs
Please submit via the GitHub issues page: https://github.com/ahmedmagds/WhatsGNU/issues
## Software Licence
GPLv3: https://github.com/ahmedmagds/WhatsGNU/blob/master/LICENSE
## Citations
* Please cite Prokka 'Seemann 2014, Bioinformatics;30(14):2068-9' if you use WhatsGNU1.0.
* Please also cite Roary 'Page et al. 2015, Bioinformatics;31(22):3691-3693' if you use WhatsGNU1.0.
* Please also cite BLAST+ 'Camacho et al. 2009, BMC Bioinformatics;10:421' if you use WhatsGNU1.0.
* Please cite Staphopia 'Petit RA III and Read TD 2018, PeerJ;6:e5261' if you use Staphopia _S. aureus_ Database.
* Please cite Enterobase 'Alikhan NF et al. 2018, PLoS Genetics;14(4):e1007261' if you use Enterobase _S. enterica_ Database.
## Author
Ahmed M. Moustafa: [ahmedmagds](https://github.com/ahmedmagds)<br/>
Twitter: [Ahmed_Microbes](https://twitter.com/Ahmed_Microbes)

