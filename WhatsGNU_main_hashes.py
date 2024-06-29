#!/usr/bin/env python3
# PROGRAM: WhatsGNU is a Python3 program that ranks protein sequences in a genome
# faa file generated from annotation programs based on the number of observed
# exact protein matches in a public or private database.

# Copyright (C) 2019 Ahmed M. Moustafa

#########################################################################################
# LICENSE
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################################

# DATE CREATED: January 25, 2019

# AUTHOR: Ahmed M Moustafa

# CONTACT1: moustafaam@email.chop.edu
# CONTACT2: ahmedmagdy2009@hotmail.com

# AFFILIATION: Pediatric Infectious Disease Division, Children’s Hospital of Philadelphia,
# Abramson Pediatric Research Center, University of Pennsylvania, Philadelphia,
# Pennsylvania, 19104, USA

# CITATION1: Ahmed M Moustafa and Paul J Planet
# WhatsGNU: a tool for identifying proteomic novelty
# Genome Biology(2020)21:58, doi: https://doi.org/10.1186/s13059-020-01965-w.

import os
import sys
from collections import defaultdict
from collections import Counter
import time
import pickle
import argparse
import logging
import tempfile
import subprocess
from shutil import rmtree as rmt
import mmh3

START_TIME = time.time()

PARSER = argparse.ArgumentParser(
    prog="WhatsGNU_main.py",
    description="WhatsGNU v1.4 utilizes the natural\
 variation in public databases to rank protein sequences based on the number of observed exact protein\
 matches (the GNU score) in all known genomes of a particular species. It generates a report for all the\
 proteins in your query in seconds.",
)
GROUP = PARSER.add_mutually_exclusive_group()
GROUP.add_argument(
    "-d",
    "--database",
    type=str,
    help="you have to provide path to your compressed database",
)
PARSER.add_argument(
    "-o",
    "--output_folder",
    type=str,
    help="Database output prefix to be created for \
results (default: timestamped WhatsGNU_results in the current directory)",
)
PARSER.add_argument(
    "--force",
    help="Force overwriting existing results folder assigned with -o (default: off)",
    action="store_true",
)
PARSER.add_argument(
    "-p",
    "--prefix",
    type=str,
    help="Prefix for output compressed database \
(default: WhatsGNU_compressed_database)",
)
PARSER.add_argument(
    "-t",
    "--topgenomes",
    help="create a file of top N genomes with most number of exact matches to query [Default top 10 genomes]",
    action="store_true",
)
PARSER.add_argument(
    "-csv",
    help="csv file of hashed inputs",
    type=str,
)
PARSER.add_argument(
    "-tn",
    "--topgenomes_count",
    type=int,
    help="select number of closest top genomes to show [Default top 10 genomes]",
)
PARSER.add_argument(
    "-s",
    "--strainhits",
    type=str,
    help="check how many hits you get from a particular strain,\
it has to be used with -t",
)
PARSER.add_argument(
    "-i",
    "--ids_hits",
    help="create a file of each protein with locus_tags (ids) of all hits from \
the database, large file (~ 1 Gb for 3000 pts)",
    action="store_true",
)
PARSER.add_argument(
        "--accession_names",
        help="to be used with --ids_hits. If this option is selected, writes the id_hits file with the accession names.",
        action="store_true"
)
PARSER.add_argument(
        "--hash_values",
        help="to be used with --ids_hits. Default option. This options writes the id_hits file with the hashed values.",
        action="store_true"
)
PARSER.add_argument(
    "-q",
    "--quiet",
    help="No screen output [default OFF]",
    action="store_true",
)
PARSER.add_argument(
    "-v",
    "--version",
    help="print version and exit",
    action="version",
    version="%(prog)s 1.4",
)
PARSER.add_argument(
    "query_faa", type=str, help="Query protein FASTA file/s to analyze (.faa)"
)
if len(sys.argv) == 1:
    PARSER.print_help()
    sys.exit(0)

ARGS = PARSER.parse_args()
if bool(vars(ARGS)["topgenomes"]) and not bool(vars(ARGS)["csv"]):
    PARSER.exit(status=0, message="Error: You have to use -csv with -t\n")
if bool(vars(ARGS)["ids_hits"]) and not bool(vars(ARGS)["csv"]):
    PARSER.exit(status=0, message="Error: You have to use -csv with -i\n")
if bool(vars(ARGS)["accession_names"]) and not bool(vars(ARGS)["ids_hits"]):
    PARSER.exit(status=0, message="Error: You have to use --accession_names with -i\n")
if bool(vars(ARGS)["hash_values"]) and not bool(vars(ARGS)["ids_hits"]):
    PARSER.exit(status=0, message="Error: You have to use --hash_values with -i\n")
OS_SEPARATOR = os.sep

if ARGS.topgenomes:
    if ARGS.topgenomes_count:
        top_count = ARGS.topgenomes_count
    else:
        top_count = 10
#####variables######
SEQUENCES_DICT = {}
TIMESTR = time.strftime("%Y%m%d_%H%M%S")
#####create results folder######
if ARGS.output_folder:
    try:
        os.mkdir(ARGS.output_folder)
        RESULTS_FOLDER = ARGS.output_folder + OS_SEPARATOR
    except:
        if ARGS.force:
            rmt(ARGS.output_folder)
            os.mkdir(ARGS.output_folder)
            RESULTS_FOLDER = ARGS.output_folder + OS_SEPARATOR
        else:
            PARSER.exit(
                status=0,
                message="Folder exists, Please change --output_folder or use --force\n")
else:
    os.mkdir("WhatsGNU_results_{}".format(TIMESTR))
    RESULTS_FOLDER = "WhatsGNU_results_{}{}".format(TIMESTR,OS_SEPARATOR)

###############Logging##################
LOG_FILE = 'WhatsGNU_'+ TIMESTR
if ARGS.quiet:
    LOG_LIST = [
        logging.FileHandler("{0}{1}.log".format(RESULTS_FOLDER, LOG_FILE))
    ]
else:
    LOG_LIST = [
        logging.FileHandler("{0}{1}.log".format(RESULTS_FOLDER, LOG_FILE)),
        logging.StreamHandler()
    ]
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s",
    handlers=LOG_LIST)
if ARGS.force:
    logging.info("overwrote results folder({})".format(RESULTS_FOLDER))
else:
    logging.info("created results folder({})".format(RESULTS_FOLDER))
#####faa files to be processed######
QUERY = ARGS.query_faa
QUERY_LIST = []
try:
    for file in os.listdir(QUERY):
        if file.endswith(".faa"):
            QUERY_LIST.append(QUERY + file)
    logging.info(
        "You provided folder of {} faa files as queries".format(len(QUERY_LIST))
    )
    if len(QUERY_LIST) == 0:
        logging.error("The directory did not have any query faa files")
        PARSER.exit(
            status=0,
            message="The directory did not have any query faa files\n",
        )
except:
    if QUERY.endswith(".faa"):
        QUERY_LIST.append(QUERY)
        logging.info("You provided one faa file as a query")
    else:
        logging.error(
            "You did not provide single faa file or path to directory with multiple faa files"
        )
        PARSER.exit(
            status=0,
            message="You did not provide single faa file or path to directory with multiple faa files\n",
        )
#####Run database for first time#####
PROTEIN_COUNTER = 0
FILE_COUNTER = 0
#####load database file#########
if ARGS.database:
    if ARGS.database.endswith(".pickle"):
        try:
            PICKLE_IN = open(ARGS.database, "rb")
            SEQUENCES_DICT = pickle.load(PICKLE_IN)
            if bool(SEQUENCES_DICT):
                logging.info(
                    "opened previously created pickle_dict in --- {:.3f} seconds ---".format(
                        time.time() - START_TIME
                    )
                )
        except:
            logging.error(
                "The compressed database pickle_dict file you provided is empty or corrupted"
            )
            PARSER.exit(
                status=0,
                message="The compressed database pickle_dict file you provided is empty or corrupted\n",
            )
    else:
        logging.error(
            "No proper compressed database (file.pickle, file.txt or file.db) was provided using -d"
        )
        PARSER.exit(
            status=0,
            message="No proper compressed database (file.pickle, file.txt or file.db) was provided using -d\n",
        )
#############################
master_hash_values = {}
if ARGS.csv:
    for line in open(ARGS.csv):
        line = line.rstrip()
        line_list = line.split(',')
        hash_val = line_list[1]
        genome_name = line_list[0]
        master_hash_values[hash_val] = genome_name
#print(master_hash_values)
#########whatsgnu############
def compute_hash(sequence):
    """Compute a hash for a given amino acid sequence."""
    #return hashlib.sha256(sequence.encode()).hexdigest()
    hash = mmh3.hash64(sequence, 42)[0]
    if hash < 0: hash += 2**64
    return hash
##############################
REPORT_LIST = ["protein", "GNU score", "length", "function", "sequence",
"ortholog_group", "ortho_gp_total_sequences_number", "ortho_gp_total_variants_number",
"minimum_GNU", "maximum_GNU", "average_GNU", "OVRI", "OVRI interpretation"]
for QUERYFILE in QUERY_LIST:
    prtn_id_dict = {}
    QUERYFILE_OBJECT = open(QUERYFILE, "r")
    line_check = QUERYFILE_OBJECT.readline()
    if not line_check.startswith(">"):
        logging.error("Not a FASTA file: {}".format(QUERYFILE))
        PARSER.exit(status=0, message="Not a FASTA file\n")
    QUERYFILE_OBJECT.seek(0)

    file_hits = (
        RESULTS_FOLDER
        + (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
        + "_WhatsGNU_hits.txt"
    )
    file_report = (
        RESULTS_FOLDER
        + (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
        + "_WhatsGNU_report.txt"
    )
    output_file_report = open(file_report, "w")
    logging.info(
        "opened report file for {}".format(
            (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
        )
    )
    output_file_report.write("{}\n".format("\t".join(REPORT_LIST[:5])))

    if ARGS.ids_hits:  # output big (roughly 1GB) hits file
        output_file_hits = open(file_hits, "w")
        output_file_hits.write("{}\t{}\n".format("protein_query", "hits"))
        logging.info(
            "opened hits file for {} as per your request of -i".format(
                (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
            )
        )

    processed_counter = 0
    strain_name_list = []
    all_hits_list = []
    sequence_info = ""
    sequence_string = ""

    for line in QUERYFILE_OBJECT:
        line = line.rstrip()
        if line.startswith(">"):
            if len(sequence_string) > 0:
                #print(sequence_string)
                processed_counter += 1
                try:
                    if ARGS.database:
                        if ARGS.database.endswith(".db"): #SQL option
                            seq_tuple = (sequence_string,)
                            c.execute('select list_of_ids from WhatsGNU where sequence = ?', seq_tuple)
                            ids_string = c.fetchone()[0]
                            ids_list = ids_string.split('.')
                        else:
                            hashed_sequence = compute_hash(sequence_string)
                            ids_string = SEQUENCES_DICT[hashed_sequence]
                    else:
                        hashed_sequence = compute_hash(sequence_string)
                        ids_string = SEQUENCES_DICT[hashed_sequence]
                    function,ids_string2 = ids_string.split('./',1)

                    prtn_id = sequence_info.split()[0]
                    prtn_id_dict[prtn_id] = hashed_sequence
                    query_sequence_details = [str(len(sequence_string)),
                                              function, sequence_string]
                    if ARGS.topgenomes:
                        for locus_id2 in ids_string2.split('.')[0:-1]:
                            locus_id2_name = master_hash_values[locus_id2]
                            strain_name_list.append(locus_id2_name)
                    ids_list = []
                    if ARGS.ids_hits: #ids_hits is for redcarpet
                        ids_tabbed_string = ids_string2.rsplit('.',1)[0]
                        if ARGS.accession_names:
                            split_ids_tabbed_string = ids_tabbed_string.split('.')
                            for id in split_ids_tabbed_string:
                                id_name = master_hash_values[id]
                                ids_list.append(id_name)
                            names_tabbed_string = '\t'.join(ids_list) #replaced . with \t
                    protein_name_GNU = [sequence_info, ids_string.rsplit('.',1)[-1]]
                        #GNU_score = ids_list[-1]
                    #GNU -8, ortho_gp_name -7, ortho_no_seq -6, no_seq_variants -5
                    #min_GNU -4, max_GNU -3, avg_GNU -2, variant_freq -1
                    if ARGS.ids_hits:
                        if ARGS.accession_names:
                            output_file_hits.write(
                            "{}\t{}\n".format(hashed_sequence, names_tabbed_string)
                        )
                        else:
                            output_file_hits.write(
                            "{}\t{}\n".format(hashed_sequence, ids_tabbed_string)
                        )

                    output_file_report.write("{}\t{}\n".format(
                        "\t".join(protein_name_GNU),
                        "\t".join(query_sequence_details)))
                except: #If the gnu score is 0 (don't see the sequence in the dictionary)
                    query_sequence_details = [sequence_info, '0', str(len(sequence_string)),
                                              'NA', sequence_string]

                    prtn_id = sequence_info.split()[0]
                    prtn_id_dict[prtn_id] = hashed_sequence

                    if ARGS.ids_hits:
                        output_file_hits.write("{}\tno_hits\n".format(hashed_sequence)) #sequence_info
                    output_file_report.write("{}\n".format(
                        "\t".join(query_sequence_details)))
                if ARGS.database:
                    if ARGS.database.endswith(".db"):
                        logging.info(
                            "processed protein {} of {} in {:.3F}".format(
                                processed_counter,
                                (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0],
                                time.time() - START_TIME,
                            )
                        )
                sequence_string = ""
            sequence_info = line.lstrip(">")
        else:
            sequence_string += line
    processed_counter += 1
    try:
        if ARGS.database:
            if ARGS.database.endswith(".db"):
                seq_tuple = (sequence_string,)
                c.execute('select list_of_ids from WhatsGNU where sequence = ?', seq_tuple)
                ids_string = c.fetchone()[0]
                ids_list = ids_string.split('.')
            else:
                hashed_sequence = compute_hash(sequence_string)
                ids_string = SEQUENCES_DICT[hashed_sequence]
                #print(ids_string)
        else:
            hashed_sequence = compute_hash(sequence_string)
            ids_string = SEQUENCES_DICT[hashed_sequence]
        function,ids_string2 = ids_string.split('./',1)

        prtn_id = sequence_info.split()[0]
        prtn_id_dict[prtn_id] = hashed_sequence
        query_sequence_details = [str(len(sequence_string)),
                                  function, sequence_string]
        if ARGS.topgenomes:
            for locus_id2 in ids_string2.split('.')[0:-1]: #basic
                locus_id2_name = master_hash_values[locus_id2]
                strain_name_list.append(locus_id2_name)
        if ARGS.ids_hits:
            ids_list = []
            if ARGS.accession_names:
                split_ids_tabbed_string = ids_tabbed_string.split('.')
                for id in split_ids_tabbed_string:
                    id_name = master_hash_values[id]
                    ids_list.append(id_name)
                names_tabbed_string = '\t'.join(ids_list)
            else:
                ids_tabbed_string = ids_string2.rsplit('.',1)[0]
        protein_name_GNU = [sequence_info, ids_string.rsplit('.',1)[-1]]
            #GNU_score = ids_list[-1]
        #GNU -8, ortho_gp_name -7, ortho_no_seq -6, no_seq_variants -5
        #min_GNU -4, max_GNU -3, avg_GNU -2, variant_freq -1
        if ARGS.ids_hits:
            if ARGS.accession_names:
                output_file_hits.write(
                "{}\t{}\n".format(hashed_sequence, names_tabbed_string) #sequence_info
            )
            else:
                output_file_hits.write(
                "{}\t{}\n".format(hashed_sequence, ids_tabbed_string) #sequence_info
            )
        output_file_report.write("{}\t{}\n".format(
            "\t".join(protein_name_GNU),
            "\t".join(query_sequence_details)))
    except:
        query_sequence_details = [sequence_info, '0', str(len(sequence_string)),
                                  'NA', sequence_string]

        prtn_id = sequence_info.split()[0]
        prtn_id_dict[prtn_id] = hashed_sequence

        if ARGS.ids_hits:
            output_file_hits.write("{}\tno_hits\n".format(hashed_sequence))
        output_file_report.write("{}\n".format(
            "\t".join(query_sequence_details)))
    if ARGS.database:
        if ARGS.database.endswith(".db"):
            logging.info(
                "processed protein {} of {} in {:.3F}".format(
                    processed_counter,
                    (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0],
                    time.time() - START_TIME,
                )
            )
    logging.info(
        "processed {} proteins of {} in {:.3F}".format(
            processed_counter,
            (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0],
            time.time() - START_TIME,
        )
    )
    if ARGS.topgenomes:  # get top 10 genomes with hits # Have to relate back to csv file here?
        file_topgenomes = (
            RESULTS_FOLDER
            + (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
            + "_WhatsGNU_topgenomes.txt"
        )
        C = Counter(strain_name_list)
        most_occur = C.most_common(top_count)
        with open(file_topgenomes, "w") as output_file_topgenomes:
            output_file_topgenomes.write(
                "WhatsGNU found {} total hits from the database in --- {:.3f} seconds ---\n".format(
                    len(strain_name_list), time.time() - START_TIME
                )
            )
            output_file_topgenomes.write(
                "\n".join("{}\t{}".format(x[0], x[1]) for x in most_occur)
            )
            if ARGS.strainhits:
                output_file_topgenomes.write("\nNo of hits from strain {} is {}".format(
                    ARGS.strainhits, C[ARGS.strainhits]))
        logging.info("Found top {} genomes with hits".format(len(most_occur)))
    if ARGS.ids_hits:
        output_file_hits.close()
        fn = RESULTS_FOLDER + (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0] + "_prtn_id_hashes.csv"
        op_fn = open(fn, 'w')
        for k,v in prtn_id_dict.items():
            op_fn.write(f'{k},{v}\n')
    output_file_report.close()
logging.info("Done in --- {:.3f} seconds ---".format(time.time() - START_TIME))
logging.info("""Thanks for using WhatsGNU1.1, I hope you found it useful.
Please cite WhatsGNU 'Moustafa AM and Planet PJ 2020, Genome Biology;21:58'.
Please cite Prokka 'Seemann 2014, Bioinformatics;30(14):2068-9' if you use WhatsGNU.
Please also cite Roary 'Page et al. 2015, Bioinformatics;31(22):3691-3693' if you use WhatsGNU.
Please also cite BLAST+ 'Camacho et al. 2009, BMC Bioinformatics;10:421' if you use WhatsGNU.
Please cite Staphopia 'Petit RA III and Read TD 2018, PeerJ;6:e5261' if you use Staphopia S. aureus Database.
Please cite Enterobase 'Alikhan NF et al. 2018, PLoS Genetics;14(4):e1007261' if you use Enterobase S. enterica Database.
The manual is extensive and available to read at https://github.com/ahmedmagds/WhatsGNU
If you have problems, please file at https://github.com/ahmedmagds/WhatsGNU/issues""")
