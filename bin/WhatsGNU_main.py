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

# CITATION1: WhatsGNU, https://github.com/ahmedmagds/WhatsGNU
# CITATION2: Ahmed M Moustafa and Paul J Planet (2019)
# WhatsGNU: a tool for identifying proteomic novelty
# TBD TBD(TBD):TBD-TBD, doi: TBD.

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
import sqlite3
from shutil import rmtree as rmt

START_TIME = time.time()

PARSER = argparse.ArgumentParser(
    prog="WhatsGNU_main.py",
    description="WhatsGNU v1.0 utilizes the natural\
 variation in public databases to rank protein sequences based on the number of observed exact protein\
 matches (the GNU score) in all known genomes of a particular species. It generates a report for all the\
 proteins in your query in seconds.",
)
GROUP = PARSER.add_mutually_exclusive_group()
GROUP.add_argument(
    "-m",
    "--mkdatabase",
    type=str,
    help="you have to provide path to faa file or \
a folder of multiple faa files for compression",
)
GROUP.add_argument(
    "-d",
    "--database",
    type=str,
    help="you have to provide path to your compressed database",
)
PARSER.add_argument(
    "-a",
    "--pickle",
    help="Save database in pickle format [Default only txt file]",
    action="store_true",
)
PARSER.add_argument(
    "-j",
    "--sql",
    help="Save database in SQL format for large Databases [Default only txt file]",
    action="store_true",
)
PARSER.add_argument(
    "-r",
    "--roary_clustered_proteins",
    type=str,
    nargs="?",
    const="clustered_proteins",
    help="clustered_proteins output file from roary to be used with -m",
)
PARSER.add_argument(
    "-dm",
    "--database_mode",
    type=str,
    choices=['ortholog', 'basic'],
    help="select a mode from 'ortholog' or 'basic' to be used with -d",
)
PARSER.add_argument(
    "-ri",
    "--rarity_index",
    nargs="?",
    type=float,
    help="select an ortholog variant rarity index (OVRI) cutoff value in range (0-1)[0.045] for ortholog mode",
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
    help="create a file of top 10 genomes with hits",
    action="store_true",
)
PARSER.add_argument(
    "-s",
    "--strainhits",
    type=str,
    help="check how many hits you get from a particular strain,\
it has to be used with -t",
)
PARSER.add_argument(
    "-e",
    "--metadata",
    type=str,
    help="get the metadata composition of your hits, use the metadata_frequency.csv file produced by the WhatsGNU customizer script",
)
PARSER.add_argument(
    "-i",
    "--ids_hits",
    help="create a file of each protein with locus_tags (ids) of all hits from \
the database, large file (~ 1 Gb for 3000 pts)",
    action="store_true",
)
PARSER.add_argument(
    "-f",
    "--faa_GNU_0",
    help="get a fasta (.faa) file of all proteins with GNU score of zero",
    action="store_true",
)
PARSER.add_argument(
    "-b",
    "--blastp",
    help="run blastp on the proteins with GNU score of zero and modify the \
report with ortholog_info, blastp has to be installed",
    action="store_true",
)
PARSER.add_argument(
    "-op",
    "--output_blastp",
    help="get the output report of blastp run, it has to be used with -b",
    action="store_true",
)
PARSER.add_argument(
    "-w",
    "--percent_identity",
    type=float,
    nargs="?",
    help="select a blastp percent identity cutoff value [80], range(0,100)",
)
PARSER.add_argument(
    "-c",
    "--percent_coverage",
    type=float,
    nargs="?",
    help="select a blastp percent coverage cutoff value [80], range(0,100)",
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
    version="%(prog)s 1.0",
)
PARSER.add_argument(
    "query_faa", type=str, help="Query protein FASTA file/s to analyze (.faa)"
)
if len(sys.argv) == 1:
    PARSER.print_help()
    sys.exit(0)
ARGS = PARSER.parse_args()
if bool(vars(ARGS)["strainhits"]) and not bool(vars(ARGS)["topgenomes"]):
    PARSER.exit(status=0, message="Error: You have to use -s with -t\n")
if bool(vars(ARGS)["database"]) and bool(vars(ARGS)["roary_clustered_proteins"]):
    PARSER.exit(status=0, message="Error: You cannot use -r with -d\n")
if bool(vars(ARGS)["database"]) and not bool(vars(ARGS)["database_mode"]):
    PARSER.exit(status=0, message="Error: You have to specify the database mode using -dm with -d\n")
if bool(vars(ARGS)["mkdatabase"]) and bool(vars(ARGS)["database_mode"]):
    PARSER.exit(status=0, message="Error: You cannot use -dm with -m\n")
if bool(vars(ARGS)["database"]) and bool(vars(ARGS)["prefix"]):
    PARSER.exit(status=0, message="Error: You cannot use -p with -d\n")
if bool(vars(ARGS)["output_blastp"]) and bool(vars(ARGS)["database_mode"] == 'basic'):
    PARSER.exit(status=0, message="Error: You have to use -p with -dm ortholog\n")
if bool(vars(ARGS)["sql"]) and bool(vars(ARGS)["roary_clustered_proteins"]):
    PARSER.exit(status=0, message="Error: You cannot use -j with -r, -j is limited for basic mode\n")
if bool(vars(ARGS)["blastp"]) and bool(vars(ARGS)["database_mode"] == 'basic'):
    PARSER.exit(status=0, message="Error: You have to use -b with -dm ortholog\n")
if bool(vars(ARGS)["blastp"]) and bool(vars(ARGS)["mkdatabase"]) and not bool(vars(ARGS)["roary_clustered_proteins"]):
    PARSER.exit(status=0, message="Error: You have to use -b with -r (or ortholog mode using -dm)\n")
if bool(vars(ARGS)["output_blastp"]) and not bool(vars(ARGS)["blastp"]):
    PARSER.exit(status=0, message="Error: You have to use -b with -op\n")
if bool(vars(ARGS)["percent_identity"]) and not bool(vars(ARGS)["blastp"]):
    PARSER.exit(status=0, message="Error: You have to use -b with -w\n")
if bool(vars(ARGS)["percent_coverage"]) and not bool(vars(ARGS)["blastp"]):
    PARSER.exit(status=0, message="Error: You have to use -b with -c\n")
OS_SEPARATOR = os.sep
if ARGS.rarity_index:
    RARITY_INDEX_CUTOFF = ARGS.rarity_index
else:
    RARITY_INDEX_CUTOFF = 0.045
if RARITY_INDEX_CUTOFF < 0.0 or RARITY_INDEX_CUTOFF > 1.0:
    PARSER.exit(status=0, message="Error: rarity_index cutoff values should be in range [0-1]\n")
##########blastp check##############
if ARGS.blastp:
    try:
        GETVERSION = subprocess.Popen("blastp -version", shell=True, stdout=subprocess.PIPE).stdout
        VERSION = GETVERSION.read()
        print("Found blastp (version:{})".format(VERSION.decode().splitlines()[1]))
    except:
        PARSER.exit(status=0, message="Error: blastp cannot be found\n")
    if ARGS.percent_identity:
        PERCENT_IDENTITY_CUTOFF = ARGS.percent_identity
    else:
        PERCENT_IDENTITY_CUTOFF = 80.0
    if ARGS.percent_coverage:
        PERCENT_COOVERAGE_CUTOFF = ARGS.percent_coverage
    else:
        PERCENT_COOVERAGE_CUTOFF = 80.0
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
#####database faa files to be compressed######
if ARGS.mkdatabase:
    DATABASE_FILES = ARGS.mkdatabase
    DATABASE_FILES_LIST = []
    try:
        for file in os.listdir(DATABASE_FILES):
            if file.endswith(".faa"):
                DATABASE_FILES_LIST.append(DATABASE_FILES + file)
        logging.info(
            "You provided folder of {} Database faa files to be processed".format(len(DATABASE_FILES_LIST))
        )
        if len(DATABASE_FILES_LIST) == 0:
            logging.error("The Database directory did not have any faa files")
            PARSER.exit(
                status=0,
                message="The Database directory did not have any faa files\n",
            )
    except:
        if DATABASE_FILES.endswith(".faa"):
            DATABASE_FILES_LIST.append(DATABASE_FILES)
            logging.info("You provided one Database faa file to be compressed")
        else:
            logging.error(
                "You did not provide single Database faa file or path to directory with multiple faa files"
            )
            PARSER.exit(
                status=0,
                message="You did not provide single Database faa file or path to directory with multiple faa files\n",
            )
    if ARGS.prefix:
        DB_PREFIX = ARGS.prefix
    else:
        DB_PREFIX = 'WhatsGNU_compressed_database'

#####Run database for first time#####
PROTEIN_COUNTER = 0
FILE_COUNTER = 0
if ARGS.mkdatabase:
    SEQUENCES_DICT_DB = defaultdict(list)
    DB_SEQUENCE_STRING = ""
    DB_SEQUENCE_INFO = ""
    for DATABASE_FILE in DATABASE_FILES_LIST:
        if DATABASE_FILE.endswith(".faa"):
            DATABASEFILE_OBJECT = open(DATABASE_FILE, "r")
            DB_LINE_CHECK = DATABASEFILE_OBJECT.readline()
            if not DB_LINE_CHECK.startswith(">"):
                logging.error("{} is not in a FASTA format".format(DATABASE_FILE))
                PARSER.exit(status=0, message="Database is not in a FASTA format\n")
            DATABASEFILE_OBJECT.seek(0)
            try:
                for line in DATABASEFILE_OBJECT:
                    line = line.rstrip()
                    if line.startswith(">"):
                        PROTEIN_COUNTER += 1
                        if len(DB_SEQUENCE_STRING) > 0:
                            SEQUENCES_DICT_DB[DB_SEQUENCE_STRING].append(DB_SEQUENCE_INFO)
                            DB_SEQUENCE_STRING = ""
                        DB_SEQUENCE_INFO = line.lstrip(">")
                    else:
                        DB_SEQUENCE_STRING += line
                SEQUENCES_DICT_DB[DB_SEQUENCE_STRING].append(DB_SEQUENCE_INFO)
                DATABASEFILE_OBJECT.close()
                FILE_COUNTER += 1
            except:
                logging.error(
                    "Cannot process the Database faa file {} provided to a compressed database".format(DATABASE_FILE)
                )
                PARSER.exit(
                    status=0,
                    message="Cannot process the Database faa file provided to a compressed database\n",
                )
            logging.info("processed file {} in --- {:.3f} seconds ---".format(FILE_COUNTER, time.time() - START_TIME))
    if bool(SEQUENCES_DICT_DB):
        logging.info(
            "processed {} proteins in {} file/s to a compressed database of {} variants in --- {:.3f} seconds ---".format(
                PROTEIN_COUNTER, len(DATABASE_FILES_LIST), len(SEQUENCES_DICT_DB), time.time() - START_TIME
            )
        )
#########Ortholog from Roary##########
    if ARGS.roary_clustered_proteins:
        ORTHOLOG_CSV_LIST = []
        ORTHOLOG_NAMES = []
        ORTHOLOG_NAMES_USED = []
        ORTHOLOG_DICT = defaultdict(list)
        ORTHOLOG_FASTA = defaultdict(list)
        ORTHOLOG_DICT_SORTED = {}
        ORTHOLOG_DICT_VARIANTS = {}
        ORTHOLOG_FASTA_LIST = []
        ROARY_ORTHO_DICT = {}
        with open(ARGS.roary_clustered_proteins, "r") as file_handle:
            for line in file_handle:
                line = line.rstrip()
                ortholog_gp_name, ids = line.split(": ")
                listids = ids.split("\t")
                for locus_id in listids:
                    ROARY_ORTHO_DICT[locus_id] = ortholog_gp_name
                number_of_seqs = str(len(listids))
                listids.insert(0, number_of_seqs)
                listids.insert(0, ortholog_gp_name)
                ORTHOLOG_CSV_LIST.append(listids)
                ORTHOLOG_NAMES.append(ortholog_gp_name)
            logging.info(
                "processed roary clustered_proteins file in --- {:.3f} seconds ---".format(
                    time.time() - START_TIME
                )
            )
        ORTHO_COUNTER = 0
        ORPHAN_COUNTER = 0
        START_TIME2 = time.time()
        for ptn_seq in SEQUENCES_DICT_DB:
            ORTHO_COUNTER += 1
            ptn_seq_ids_list = SEQUENCES_DICT_DB[ptn_seq]
            ptn_seq_locustag_1st = ptn_seq_ids_list[0].split()[0] #whitespace split
            ptn_seq_GNU_score = len(ptn_seq_ids_list)
            ptn_seq_ids_list.append(str(ptn_seq_GNU_score)) #add GNU score
            try:
                ortho_group_name = ROARY_ORTHO_DICT[ptn_seq_locustag_1st]
                ORTHOLOG_NAMES_USED.append(ortho_group_name)
                ORTHOLOG_DICT[ortho_group_name].append(ptn_seq_GNU_score) #ortholog alleles GNU scores
                ORTHOLOG_FASTA[ortho_group_name].append(ptn_seq) #ortholog protein FASTA sequences
                SEQUENCES_DICT_DB[ptn_seq].append(ortho_group_name) #add ortholog group name
            except:
                SEQUENCES_DICT_DB[ptn_seq].append('Orphan') #if locustag deleted by roary
                ORPHAN_COUNTER += 1
        logging.info('{} proteins are Orphans with no information in the clustered_proteins file'.format(ORPHAN_COUNTER))
        ORTHOLOG_FILE_NAME = (RESULTS_FOLDER
            + DB_PREFIX
            + "_orthologs"
            + ".faa"
        )
        OUTPUT_FILE_ORTHO = open(ORTHOLOG_FILE_NAME, "w")
        for record in ORTHOLOG_FASTA:
            ortholog_sequence = max(ORTHOLOG_FASTA[record], key=len) #longest sequence
            ORTHOLOG_FASTA_LIST.append(ortholog_sequence)
            OUTPUT_FILE_ORTHO.write(
                ">{}\n{}\n".format(record, ortholog_sequence))
        #FASTA of a representative sequence from each ortho_gp needed for blastp
        OUTPUT_FILE_ORTHO.close()
        for record in ORTHOLOG_DICT:
            ortho_all = []
            ortholog_list_sorted = sorted(ORTHOLOG_DICT[record], reverse=True)
            ORTHOLOG_DICT_VARIANTS[record] = ortholog_list_sorted
            for i in ortholog_list_sorted:
                ortho_multiplied = []
                ortho_multiplied = [i] *i
                ortho_all += ortho_multiplied
            ORTHOLOG_DICT_SORTED[record] = ortho_all
        ORTHOLOG_COUNTER_2 = 0
        START_TIME3 = time.time()
        for ptn_seq in SEQUENCES_DICT_DB:
            ORTHOLOG_COUNTER_2 += 1
            ptn_seq_ids_list = SEQUENCES_DICT_DB[ptn_seq]
            try:
                record_GNU_score = int(ptn_seq_ids_list[-2]) #GNU score -8
                ortholog_group_name = ptn_seq_ids_list[-1] #ortholog_gp_name -7
                ortholog_list_GNU = ORTHOLOG_DICT_SORTED[ortholog_group_name] #otholog GNU score for each member
                ortholog_list_variants = ORTHOLOG_DICT_VARIANTS[ortholog_group_name]
                ptn_seq_ids_list.append(str(len(ortholog_list_GNU))) #total no of sequences in the group #-6
                ptn_seq_ids_list.append(str(len(ortholog_list_variants))) #number of variants in a gp #-5
                ptn_seq_ids_list.append(str(ortholog_list_GNU[-1])) #minimum GNU #-4
                ptn_seq_ids_list.append(str(ortholog_list_GNU[0])) #maximum GNU #-3
                ortholog_GNU_avg = sum(ortholog_list_GNU)/sum(ortholog_list_variants)
                ptn_seq_ids_list.append(str("{:.1f}".format(ortholog_GNU_avg))) #avg GNU #-2
                ptn_seq_ids_list.append(str("{:.4f}".format(
                    sum(ortholog_list_variants[ortholog_list_variants.index(
                    record_GNU_score):])/sum(ortholog_list_variants)))) #variant_freq #-1
                #GNU -8, ortho_gp_name -7, ortho_no_seq -6, no_seq_variants -5
                #min_GNU -4, max_GNU -3, avg_GNU -2, variant_freq -1
            except:
                logging.warning(
                    "No ortholog group information was found for {}".format(
                        ptn_seq_ids_list[0]
                    )
                )
                ptn_seq_ids_list.extend(['NA']*6) #otholog group name
        ORTHOLOG_INFO_FILE_NAME = (RESULTS_FOLDER
            + DB_PREFIX
            + "_orthologs_info"
            + ".txt"
        )
        OUTPUT_FILE_ORTHO_INFO = open(ORTHOLOG_INFO_FILE_NAME, "w")
        for sequence in ORTHOLOG_FASTA_LIST:
            sequence_record_list = SEQUENCES_DICT_DB[sequence]
            ortholog_name = sequence_record_list[-7]
            ortho_members_count = sequence_record_list[-6]
            seq_variants_count = sequence_record_list[-5]
            minimum_GNU = sequence_record_list[-4]
            maximum_GNU = sequence_record_list[-3]
            mean_GNU = sequence_record_list[-2]
            OUTPUT_FILE_ORTHO_INFO.write(
                "{}\t{}\t{}\t{}\t{}\t{}\n".format(ortholog_name,
                ortho_members_count, seq_variants_count, minimum_GNU,
                maximum_GNU, mean_GNU
                ) #ortholog groups information needed for blastp
            )
        OUTPUT_FILE_ORTHO_INFO.close()
        if bool(SEQUENCES_DICT_DB):
            logging.info(
                "Added {} ortholog groups information to the compressed database in --- {:.3f} seconds ---".format(
                    len(set(ORTHOLOG_NAMES_USED)), time.time() - START_TIME
                )
            )
    else: #basic mode
        for ptn_seq in SEQUENCES_DICT_DB:
            ptn_seq_ids_list = SEQUENCES_DICT_DB[ptn_seq]
            ptn_seq_GNU_score = len(ptn_seq_ids_list)
            ptn_seq_ids_list.append(str(ptn_seq_GNU_score)) #GNU score
    logging.info(
        "Added GNU score to the compressed database in --- {:.3f} seconds ---".format(
            time.time() - START_TIME
        )
    )
#####Save database for first time as txt file#####
    try:
        TXT_FILE_NAME = (RESULTS_FOLDER
            + DB_PREFIX
            + ".txt"
        )
        OUTPUT_FILE_DB = open(TXT_FILE_NAME, "w")
        for record in SEQUENCES_DICT_DB:
            OUTPUT_FILE_DB.write(
                "{}\t{}\n".format(record, "._/".join(SEQUENCES_DICT_DB[record]))
            )
        OUTPUT_FILE_DB.close()
        logging.info(
            "saved database ({}) of {} proteins as txt file in --- {:.3f} seconds ---".format(
                TXT_FILE_NAME, len(SEQUENCES_DICT_DB), time.time() - START_TIME
            )
        )
    except:
        logging.critical(
            "cannot save compressed database as txt file, this time will be ok as the compressed dictionary will be used"
        )
#####Save database for first time as pickle file#####
    if ARGS.pickle:
        try:
            PICKLE_FILE_NAME = (RESULTS_FOLDER
                + DB_PREFIX
                + ".pickle"
            )
            PICKLE_OUT = open(PICKLE_FILE_NAME, "wb")
            pickle.dump(SEQUENCES_DICT_DB, PICKLE_OUT)
            PICKLE_OUT.close()
            logging.info(
                "saved pickle_dict ({}) in --- {:.3f} seconds ---".format(
                    PICKLE_FILE_NAME, time.time() - START_TIME
                )
            )
        except:
            logging.warning("cannot save pickle file")
#####Save database for first time as SQLite3 file#####
    if ARGS.sql:
        try:
            SQL_FILE_NAME = (RESULTS_FOLDER
                + DB_PREFIX
                + ".db"
            )
            CONN = sqlite3.connect(SQL_FILE_NAME)
            c = CONN.cursor()
            c.execute("create table WhatsGNU (sequence text, list_of_ids text)")
            for record in SEQUENCES_DICT_DB:
                c.execute("insert into WhatsGNU values (?,?)", (record, '._/'.join(SEQUENCES_DICT_DB[record])))
            CONN.commit()
            c.close()
            CONN.close()
            logging.info(
                "saved SQL DB ({}) in --- {:.3f} seconds ---".format(
                    SQL_FILE_NAME, time.time() - START_TIME
                )
            )
        except:
            logging.warning("cannot save SQL DB file")

#####load database file#########
if ARGS.mkdatabase:
    SEQUENCES_DICT = SEQUENCES_DICT_DB
    logging.info(
        "As you just created a compressed database using -m, it will be used this time, next time provide the database using -d"
    )
elif ARGS.database:
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
            if ARGS.blastp:
                ORTHOLOG_FILE_NAME = (
                    ARGS.database.split(".pickle")[0]
                    + "_orthologs"
                    + ".faa"
                )
                ORTHOLOG_INFO_FILE_NAME = (
                    ARGS.database.split(".pickle")[0]
                    + "_orthologs_info"
                    + ".txt"
                )
        except:
            logging.error(
                "The compressed database pickle_dict file you provided is empty or corrupted"
            )
            PARSER.exit(
                status=0,
                message="The compressed database pickle_dict file you provided is empty or corrupted\n",
            )
    elif ARGS.database.endswith(".txt"):
        try:
            SEQUENCES_DICT = {}
            TXT_DB_FILEOBJECT = open(ARGS.database, "r")
            logging.info(
                "opened previously created compressed txt database in --- {:.3f} seconds ---".format(
                    time.time() - START_TIME
                )
            )
            for line in TXT_DB_FILEOBJECT:
                line = line.rstrip()
                seq, ids = line.split("\t")
                listids = ids.split("._/")
                SEQUENCES_DICT[seq] = listids
            if bool(SEQUENCES_DICT):
                logging.info(
                    "processed compressed txt database to dictionary in --- {:.3f} seconds ---".format(
                        time.time() - START_TIME
                    )
                )
            if ARGS.blastp:
                ORTHOLOG_FILE_NAME = (
                    ARGS.database.split(".txt")[0]
                    + "_orthologs"
                    + ".faa"
                )
                ORTHOLOG_INFO_FILE_NAME = (
                    ARGS.database.split(".txt")[0]
                    + "_orthologs_info"
                    + ".txt"
                )
        except:
            logging.error(
                "The compressed database txt file you provided is empty or corrupted"
            )
            PARSER.exit(
                status=0,
                message="The compressed database txt file you provided is empty or corrupted\n",
            )
    elif ARGS.database.endswith(".db"):
        try:
            CONN = sqlite3.connect(ARGS.database)
            c = CONN.cursor()
            logging.info(
                "Connected to SQL database in --- {:.3f} seconds ---".format(
                    time.time() - START_TIME
                )
            )
        except:
            logging.error(
                "Could not connect to SQL database file you provided"
            )
            PARSER.exit(
                status=0,
                message="Could not connect to SQL database file you provided\n",
            )
    else:
        logging.error(
            "No proper compressed database (file.pickle, file.txt or file.db) was provided using -d"
        )
        PARSER.exit(
            status=0,
            message="No proper compressed database (file.pickle, file.txt or file.db) was provided using -d\n",
        )
else:
    logging.error(
        "Neither you created new database using -m (file.faa) nor proper database (file.pickle, file.txt or file.db) was provided using -d"
    )
    PARSER.exit(
        status=0,
        message="Neither you created new database using -m (file.faa) nor proper database (file.pickle or file.txt) was provided using -d\n",
    )

########Ortho info#########
if ARGS.blastp:
    try:
        ORTHO_DICT = {}
        with open(ORTHOLOG_INFO_FILE_NAME, "r") as filehandle:
            for line in filehandle:
                line = line.rstrip()
                ORTHO_DICT[line.split("\t")[0]] = line.split("\t")[1:]
        logging.info("Processed Ortholog info file ({})".format(ORTHOLOG_INFO_FILE_NAME))
    except:
        logging.error(
            "Could not process Ortholog info file"
        )
        PARSER.exit(
            status=0,
            message="Could not process Ortholog info file\n",
        )
    if os.path.exists(ORTHOLOG_FILE_NAME):
        logging.info("Found Ortholog faa file ({})".format(ORTHOLOG_FILE_NAME))
    else:
        logging.error(
            "Could not find Ortholog faa file"
        )
        sys.exit(0)
########metadata database info#########
if ARGS.metadata:
    METADATA_NAMES_LIST = []
    METADATA_FREQUENCIES_LIST = []
    METADATA_NAMES_FREQ_LIST_FILE = open(ARGS.metadata, "r")
    for line in METADATA_NAMES_FREQ_LIST_FILE:
        line = line.rstrip()
        metadata, metadata_freq = line.split(",")
        METADATA_NAMES_LIST.append(metadata)
        METADATA_FREQUENCIES_LIST.append(int(metadata_freq))
    logging.info("processed metadata list")
#########whatsgnu###########
REPORT_LIST = ["protein", "GNU score", "length", "function", "sequence",
"ortholog_group", "ortho_gp_total_sequences_number", "ortho_gp_total_variants_number",
"minimum_GNU", "maximum_GNU", "average_GNU", "OVRI", "OVRI interpretation"]
for QUERYFILE in QUERY_LIST:
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
    if ARGS.blastp:
        GNU_report_tmp = tempfile.NamedTemporaryFile(mode='w+')
        file_report = (
            RESULTS_FOLDER
            + (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
            + "_WhatsGNU_report.txt"
        )
        output_file_report = open(file_report, "w")
        if ARGS.output_blastp:
            blast_report = (
                RESULTS_FOLDER
                + (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
                + "_WhatsGNU_zeros_blast_report.txt"
            )
            logging.info(
                "opened blast output file for {}".format(
                    (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
                )
            )
        else:
            blast_report_tmp = tempfile.NamedTemporaryFile(mode='w+')
            logging.info(
                "opened temporary file for blast results for {}".format(
                    (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
                )
            )
        if ARGS.faa_GNU_0:
            ZEROS_file_hits = (
                RESULTS_FOLDER
                + (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
                + "_WhatsGNU_zeros.faa"
            )
            OUTPUT_FILE_Zeros = open(ZEROS_file_hits, "w")
            logging.info(
                "opened fasta file for proteins with GNU_score of zero for {}".format(
                    (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
                )
            )
        else:
            OUTPUT_FILE_Zeros = tempfile.NamedTemporaryFile(mode='w+')
            logging.info(
                "opened temporary fasta file for proteins with GNU_score of zero for {}".format(
                    (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
                )
            )
    else:
        file_report = (
            RESULTS_FOLDER
            + (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
            + "_WhatsGNU_report.txt"
        )
        output_file_report = open(file_report, "w")
        if ARGS.faa_GNU_0:
            ZEROS_file_hits = (
                RESULTS_FOLDER
                + (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
                + "_WhatsGNU_zeros.faa"
            )
            OUTPUT_FILE_Zeros = open(ZEROS_file_hits, "w")
            logging.info(
                "opened fasta file for proteins with GNU_score of zero for {}".format(
                    (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
                )
            )
    logging.info(
        "opened report file for {}".format(
            (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
        )
    )
    if (ARGS.roary_clustered_proteins) or (ARGS.database_mode == 'ortholog'):
        if ARGS.metadata:
            if ARGS.blastp:#ortho normal blastp
                GNU_report_tmp.write("{}\t{}\n".format(
                    "\t".join(REPORT_LIST),
                    "\t".join(METADATA_NAMES_LIST)))
                output_file_report.write("{}\t{}\n".format(
                    "\t".join(REPORT_LIST),
                    "\t".join(METADATA_NAMES_LIST)))
            else:#ortho normal normal
                output_file_report.write("{}\t{}\n".format(
                    "\t".join(REPORT_LIST),
                    "\t".join(METADATA_NAMES_LIST)))
        else:#ortho normal
            if ARGS.blastp:#ortho normal blastp
                GNU_report_tmp.write("{}\n".format("\t".join(REPORT_LIST)))
                output_file_report.write("{}\n".format("\t".join(REPORT_LIST)))
            else:#ortho normal normal
                output_file_report.write("{}\n".format("\t".join(REPORT_LIST)))
    else:#ARGS.database_mode == 'basic' or just -m without -roary
        if ARGS.metadata:#basic metadata
            output_file_report.write("{}\t{}\n".format(
                "\t".join(REPORT_LIST[:5]),
                "\t".join(METADATA_NAMES_LIST)))
        else:#basic normal
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
    metadata_count_list = []
    metadata_percentage_list = []
    strain_name_list = []
    all_hits_list = []
    if ARGS.metadata:
        metadata_zeros_list = ["0.00"] * (len(METADATA_NAMES_LIST))
    ortho_NA_list = ["NA"] * 6
    sequence_info = ""
    sequence_string = ""
    for line in QUERYFILE_OBJECT:
        line = line.rstrip()
        if line.startswith(">"):
            if len(sequence_string) > 0:
                processed_counter += 1
                try:
                    metadata_count_list = []
                    metadata_percentage_list = []
                    if ARGS.database:
                        if ARGS.database.endswith(".db"):
                            seq_tuple = (sequence_string,)
                            c.execute('select list_of_ids from WhatsGNU where sequence = ?', seq_tuple)
                            ids_string = c.fetchone()[0]
                            ids_list = ids_string.split('._/')
                        else:
                            ids_list = SEQUENCES_DICT[sequence_string]
                    else:
                        ids_list = SEQUENCES_DICT[sequence_string]
                    try:
                        function = ids_list[0].split(' ', 1)[1]
                    except:
                        function = ids_list[0]
                    query_sequence_details = [str(len(sequence_string)),
                                              function, sequence_string]
                    if (ARGS.roary_clustered_proteins) or (ARGS.database_mode == 'ortholog'):
                        if ARGS.topgenomes:
                            for locus_id2 in ids_list[0:-8]: #ortho
                                strain_name_list.append(locus_id2.split('|')[0])
                        if ARGS.ids_hits or ARGS.metadata:
                            ids_tabbed_string = "\t".join(ids_list[0:-8])
                        rarity_value = float(ids_list[-1])
                        if rarity_value <= RARITY_INDEX_CUTOFF:
                            variant_frequency = 'rare'
                        else:
                            variant_frequency = 'frequent'
                        protein_name_GNU = [sequence_info, ids_list[-8]]
                        report_values = [ids_list[-7], ids_list[-6],
                                ids_list[-5], ids_list[-4], ids_list[-3],
                                ids_list[-2], str(rarity_value), variant_frequency]
                    else:
                        if ARGS.topgenomes:
                            for locus_id2 in ids_list[0:-1]: #basic
                                strain_name_list.append(locus_id2.split('|')[0])
                        if ARGS.ids_hits or ARGS.metadata:
                            ids_tabbed_string = "\t".join(ids_list[0:-1])
                        protein_name_GNU = [sequence_info, ids_list[-1]]
                        #GNU_score = ids_list[-1]
                    #GNU -8, ortho_gp_name -7, ortho_no_seq -6, no_seq_variants -5
                    #min_GNU -4, max_GNU -3, avg_GNU -2, variant_freq -1
                    if ARGS.ids_hits:
                        output_file_hits.write(
                            "{}\t{}\n".format(sequence_info, ids_tabbed_string)
                        )
                    if ARGS.metadata:
                        for metadatum in METADATA_NAMES_LIST:
                            metadata_count_list.append(
                                ids_tabbed_string.count("_" + metadatum + "_")
                            )
                        metadata_percentage_list = [
                            str("{:.2f}".format(metadatum_count * 100 / metadatum_freq))
                            for metadatum_count, metadatum_freq in zip(
                                metadata_count_list, METADATA_FREQUENCIES_LIST
                            )
                        ]
                        ##############write to report#########
                    if (ARGS.roary_clustered_proteins) or (ARGS.database_mode == 'ortholog'):#ortho
                        if ARGS.metadata:#ortho metadata
                            if ARGS.blastp:#ortho metadata blastp
                                GNU_report_tmp.write("{}\t{}\t{}\t{}\n".format(
                                    "\t".join(protein_name_GNU),
                                    "\t".join(query_sequence_details),
                                    "\t".join(report_values),
                                    "\t".join(metadata_percentage_list)))
                            else:#ortho metadata nobp
                                output_file_report.write("{}\t{}\t{}\t{}\n".format(
                                    "\t".join(protein_name_GNU),
                                    "\t".join(query_sequence_details),
                                    "\t".join(report_values),
                                    "\t".join(metadata_percentage_list)))
                        else:#ortho nometa
                            if ARGS.blastp:#ortho nometa blastp
                                GNU_report_tmp.write("{}\t{}\t{}\n".format(
                                    "\t".join(protein_name_GNU),
                                    "\t".join(query_sequence_details),
                                    "\t".join(report_values)))
                            else:#ortho nometa nobp
                                output_file_report.write("{}\t{}\t{}\n".format(
                                    "\t".join(protein_name_GNU),
                                    "\t".join(query_sequence_details),
                                    "\t".join(report_values)))
                    else:#ARGS.database_mode == 'basic' or just -m without -roary
                        if ARGS.metadata:#basic metadata
                            output_file_report.write("{}\t{}\t{}\n".format(
                                "\t".join(protein_name_GNU),
                                "\t".join(query_sequence_details),
                                "\t".join(metadata_percentage_list)))
                        else:#basic nometa
                            output_file_report.write("{}\t{}\n".format(
                                "\t".join(protein_name_GNU),
                                "\t".join(query_sequence_details)))
                except:
                    query_sequence_details = [sequence_info, '0', str(len(sequence_string)),
                                              'NA', sequence_string]
                    if ARGS.ids_hits:
                        output_file_hits.write("{}\tno_hits\n".format(sequence_info))
                    if (ARGS.roary_clustered_proteins) or (ARGS.database_mode == 'ortholog'):#ortho
                        if ARGS.metadata:#ortho metadata
                            if ARGS.blastp:#ortho metadata blastp
                                GNU_report_tmp.write("{}\t{}\t0\trare\t{}\n".format(
                                    "\t".join(query_sequence_details),
                                    "\t".join(ortho_NA_list),
                                    "\t".join(metadata_zeros_list)))
                                OUTPUT_FILE_Zeros.write(
                                    ">{}\n{}\n".format(sequence_info, sequence_string))
                            else:#ortho metadata nobp
                                output_file_report.write("{}\t{}\t0\trare\t{}\n".format(
                                    "\t".join(query_sequence_details),
                                    "\t".join(ortho_NA_list),
                                    "\t".join(metadata_zeros_list)))
                                if ARGS.faa_GNU_0:
                                    OUTPUT_FILE_Zeros.write(
                                        ">{}\n{}\n".format(sequence_info, sequence_string))
                        else:#ortho nometa
                            if ARGS.blastp:#ortho nometa blastp
                                GNU_report_tmp.write("{}\t{}\t0\trare\n".format(
                                    "\t".join(query_sequence_details),
                                    "\t".join(ortho_NA_list)))
                                OUTPUT_FILE_Zeros.write(
                                    ">{}\n{}\n".format(sequence_info, sequence_string))
                            else:#ortho nometa nobp
                                output_file_report.write("{}\t{}\t0\trare\n".format(
                                    "\t".join(query_sequence_details),
                                    "\t".join(ortho_NA_list)))
                                if ARGS.faa_GNU_0:
                                    OUTPUT_FILE_Zeros.write(
                                        ">{}\n{}\n".format(sequence_info, sequence_string))
                    else:#ARGS.database_mode == 'basic' or just -m without -roary
                        if ARGS.faa_GNU_0:
                            OUTPUT_FILE_Zeros.write(
                                ">{}\n{}\n".format(sequence_info, sequence_string))
                        if ARGS.metadata:#basic metadata
                            output_file_report.write("{}\t{}\n".format(
                                "\t".join(query_sequence_details),
                                "\t".join(metadata_zeros_list)))
                        else:#basic nometa
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
        metadata_count_list = []
        metadata_percentage_list = []
        if ARGS.database:
            if ARGS.database.endswith(".db"):
                seq_tuple = (sequence_string,)
                c.execute('select list_of_ids from WhatsGNU where sequence = ?', seq_tuple)
                ids_string = c.fetchone()[0]
                ids_list = ids_string.split('._/')
            else:
                ids_list = SEQUENCES_DICT[sequence_string]
        else:
            ids_list = SEQUENCES_DICT[sequence_string]
        try:
            function = ids_list[0].split(' ', 1)[1]
        except:
            function = ids_list[0]
        query_sequence_details = [str(len(sequence_string)),
                                  function, sequence_string]
        if (ARGS.roary_clustered_proteins) or (ARGS.database_mode == 'ortholog'):
            if ARGS.topgenomes:
                for locus_id2 in ids_list[0:-8]: #ortho
                    strain_name_list.append(locus_id2.split('|')[0])
            if ARGS.ids_hits or ARGS.metadata:
                ids_tabbed_string = "\t".join(ids_list[0:-8])
            rarity_value = float(ids_list[-1])
            if rarity_value <= RARITY_INDEX_CUTOFF:
                variant_frequency = 'rare'
            else:
                variant_frequency = 'frequent'
            protein_name_GNU = [sequence_info, ids_list[-8]]
            report_values = [ids_list[-7], ids_list[-6],
                    ids_list[-5], ids_list[-4], ids_list[-3],
                    ids_list[-2], str(rarity_value), variant_frequency]
        else:
            if ARGS.topgenomes:
                for locus_id2 in ids_list[0:-1]: #basic
                    strain_name_list.append(locus_id2.split('|')[0])
            if ARGS.ids_hits or ARGS.metadata:
                ids_tabbed_string = "\t".join(ids_list[0:-1])
            protein_name_GNU = [sequence_info, ids_list[-1]]
            #GNU_score = ids_list[-1]
        #GNU -8, ortho_gp_name -7, ortho_no_seq -6, no_seq_variants -5
        #min_GNU -4, max_GNU -3, avg_GNU -2, variant_freq -1
        if ARGS.ids_hits:
            output_file_hits.write(
                "{}\t{}\n".format(sequence_info, ids_tabbed_string)
            )
        if ARGS.metadata:
            for metadatum in METADATA_NAMES_LIST:
                metadata_count_list.append(
                    ids_tabbed_string.count("_" + metadatum + "_")
                )
            metadata_percentage_list = [
                str("{:.2f}".format(metadatum_count * 100 / metadatum_freq))
                for metadatum_count, metadatum_freq in zip(
                    metadata_count_list, METADATA_FREQUENCIES_LIST
                )
            ]
            ##############write to report#########
        if (ARGS.roary_clustered_proteins) or (ARGS.database_mode == 'ortholog'):#ortho
            if ARGS.metadata:#ortho metadata
                if ARGS.blastp:#ortho metadata blastp
                    GNU_report_tmp.write("{}\t{}\t{}\t{}\n".format(
                        "\t".join(protein_name_GNU),
                        "\t".join(query_sequence_details),
                        "\t".join(report_values),
                        "\t".join(metadata_percentage_list)))
                else:#ortho metadata nobp
                    output_file_report.write("{}\t{}\t{}\t{}\n".format(
                        "\t".join(protein_name_GNU),
                        "\t".join(query_sequence_details),
                        "\t".join(report_values),
                        "\t".join(metadata_percentage_list)))
            else:#ortho nometa
                if ARGS.blastp:#ortho nometa blastp
                    GNU_report_tmp.write("{}\t{}\t{}\n".format(
                        "\t".join(protein_name_GNU),
                        "\t".join(query_sequence_details),
                        "\t".join(report_values)))
                else:#ortho nometa nobp
                    output_file_report.write("{}\t{}\t{}\n".format(
                        "\t".join(protein_name_GNU),
                        "\t".join(query_sequence_details),
                        "\t".join(report_values)))
        else:#ARGS.database_mode == 'basic' or just -m without -roary
            if ARGS.metadata:#basic metadata
                output_file_report.write("{}\t{}\t{}\n".format(
                    "\t".join(protein_name_GNU),
                    "\t".join(query_sequence_details),
                    "\t".join(metadata_percentage_list)))
            else:#basic nometa
                output_file_report.write("{}\t{}\n".format(
                    "\t".join(protein_name_GNU),
                    "\t".join(query_sequence_details)))
    except:
        query_sequence_details = [sequence_info, '0', str(len(sequence_string)),
                                  'NA', sequence_string]
        if ARGS.ids_hits:
            output_file_hits.write("{}\tno_hits\n".format(sequence_info))
        if (ARGS.roary_clustered_proteins) or (ARGS.database_mode == 'ortholog'):#ortho
            if ARGS.metadata:#ortho metadata
                if ARGS.blastp:#ortho metadata blastp
                    GNU_report_tmp.write("{}\t{}\t0\trare\t{}\n".format(
                        "\t".join(query_sequence_details),
                        "\t".join(ortho_NA_list),
                        "\t".join(metadata_zeros_list)))
                    OUTPUT_FILE_Zeros.write(
                        ">{}\n{}\n".format(sequence_info, sequence_string))
                else:#ortho metadata nobp
                    output_file_report.write("{}\t{}\t0\trare\t{}\n".format(
                        "\t".join(query_sequence_details),
                        "\t".join(ortho_NA_list),
                        "\t".join(metadata_zeros_list)))
                    if ARGS.faa_GNU_0:
                        OUTPUT_FILE_Zeros.write(
                            ">{}\n{}\n".format(sequence_info, sequence_string))
            else:#ortho nometa
                if ARGS.blastp:#ortho nometa blastp
                    GNU_report_tmp.write("{}\t{}\t0\trare\n".format(
                        "\t".join(query_sequence_details),
                        "\t".join(ortho_NA_list)))
                    OUTPUT_FILE_Zeros.write(
                        ">{}\n{}\n".format(sequence_info, sequence_string))
                else:#ortho nometa nobp
                    output_file_report.write("{}\t{}\t0\trare\n".format(
                        "\t".join(query_sequence_details),
                        "\t".join(ortho_NA_list)))
                    if ARGS.faa_GNU_0:
                        OUTPUT_FILE_Zeros.write(
                            ">{}\n{}\n".format(sequence_info, sequence_string))
        else:#ARGS.database_mode == 'basic' or just -m without -roary
            if ARGS.faa_GNU_0:
                OUTPUT_FILE_Zeros.write(
                    ">{}\n{}\n".format(sequence_info, sequence_string))
            if ARGS.metadata:#basic metadata
                output_file_report.write("{}\t{}\n".format(
                    "\t".join(query_sequence_details),
                    "\t".join(metadata_zeros_list)))
            else:#basic nometa
                output_file_report.write("{}\n".format(
                    "\t".join(query_sequence_details)))
    logging.info(
        "processed {} proteins of {} in {:.3F}".format(
            processed_counter,
            (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0],
            time.time() - START_TIME,
        )
    )
    #blast step here
    if ARGS.blastp:
        OUTPUT_FILE_Zeros.seek(0)
        logging.info("Running blast for proteins with GNU score of Zero for {}".format((
            QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]))
        if ARGS.output_blastp:
            if ARGS.faa_GNU_0:
                blast_results = os.system("blastp -query {} -subject {} -evalue 1e-06 -max_target_seqs 5 -max_hsps 1 -outfmt '6 qseqid sacc evalue qcovs pident' -out {}".format(ZEROS_file_hits, ORTHOLOG_FILE_NAME, blast_report))
                logging.info("Saving the blast output as {}".format(blast_report.rsplit(OS_SEPARATOR, 1)[-1]))
            else:
                blast_results = os.system("blastp -query {} -subject {} -evalue 1e-06 -max_target_seqs 5 -max_hsps 1 -outfmt '6 qseqid sacc evalue qcovs pident' -out {}".format(OUTPUT_FILE_Zeros.name, ORTHOLOG_FILE_NAME, blast_report))
                logging.info("Saving the blast output as {}".format(blast_report.rsplit(OS_SEPARATOR, 1)[-1]))
        else:
            if ARGS.faa_GNU_0:
                blast_results = os.system("blastp -query {} -subject {} -evalue 1e-06 -max_target_seqs 5 -max_hsps 1 -outfmt '6 qseqid sacc evalue qcovs pident' -out {}".format(ZEROS_file_hits, ORTHOLOG_FILE_NAME, blast_report_tmp.name))
                logging.info("Saving the blast output as a temporary file")
            else:
                blast_results = os.system("blastp -query {} -subject {} -evalue 1e-06 -max_target_seqs 5 -max_hsps 1 -outfmt '6 qseqid sacc evalue qcovs pident' -out {}".format(OUTPUT_FILE_Zeros.name, ORTHOLOG_FILE_NAME, blast_report_tmp.name))
                logging.info("Saving the blast output as a temporary file")
        if blast_results == 0:
            logging.info(
                "blast ran successfully & now editing the GNU report with the blast results"
            )
            #####Ortho info file parsed before opening query file to ORTHO_DICT#####
            #####Blast results Parser#####
            blast_dict = {}
            if ARGS.output_blastp:
                blast_report_file = open(blast_report, "r")
                for line in blast_report_file:
                    line = line.rstrip()
                    line_info = line.split("\t")
                    if float(line_info[-1]) > PERCENT_IDENTITY_CUTOFF:
                        if float(line_info[-2]) > PERCENT_COOVERAGE_CUTOFF:
                            if  line_info[0] not in blast_dict:
                                blast_dict[line_info[0]] = [line_info[1]] + ORTHO_DICT[line_info[1]]
                blast_report_file.close()
            else:
                blast_report_tmp.seek(0)
                for line in blast_report_tmp:
                    line = line.rstrip()
                    line_info = line.split("\t")
                    if float(line_info[-1]) > PERCENT_IDENTITY_CUTOFF:
                        if float(line_info[-2]) > PERCENT_COOVERAGE_CUTOFF:
                            if  line_info[0] not in blast_dict:
                                blast_dict[line_info[0]] = [line_info[1]] + ORTHO_DICT[line_info[1]]
                blast_report_tmp.close()
            #########GNU report parser##########
            GNU_report_dict = {}
            GNU_report_tmp.seek(0)
            GNU_report_tmp.readline()
            for line in GNU_report_tmp:
                line = line.rstrip()
                protein_name = line.split("\t")[0]
                GNU_info = line.split("\t")[1:]
                protein_id, protein_info = protein_name.split(" ", 1)
                GNU_info.insert(0, protein_info)
                GNU_report_dict[protein_id] = GNU_info
            #########GNU report update with blast results##########
            for record in GNU_report_dict:
                try:
                    blast_list = GNU_report_dict[record]
                    blast_list[5] = blast_dict[record][0]
                    blast_list[6] = blast_dict[record][1]
                    blast_list[7] = blast_dict[record][2]
                    blast_list[8] = blast_dict[record][3]
                    blast_list[9] = blast_dict[record][4]
                    blast_list[10] = blast_dict[record][5]
                except:
                    pass
            #########Write New GNU report ##########
            for record in GNU_report_dict:
                record_joined_name = record + ' ' + GNU_report_dict[record][0]
                record_GNU_list = GNU_report_dict[record][1:]
                output_file_report.write(
                    "{}\t{}\n".format(
                        record_joined_name, "\t".join(record_GNU_list)
                                    )
                )
        else:
            logging.critical("Error: No blast results")
            GNU_report_tmp.seek(0)
            GNU_report_tmp.readline()
            for line in GNU_report_tmp:
                output_file_report.write(line)
        GNU_report_tmp.close()
    if ARGS.topgenomes:  # get top 10 genomes with hits
        file_topgenomes = (
            RESULTS_FOLDER
            + (QUERYFILE.rsplit(OS_SEPARATOR, 1)[-1]).split(".faa")[0]
            + "_WhatsGNU_topgenomes.txt"
        )
        C = Counter(strain_name_list)
        most_occur = C.most_common(10)
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
    output_file_report.close()
    if ARGS.faa_GNU_0:
        OUTPUT_FILE_Zeros.close()
    QUERYFILE_OBJECT.close()
logging.info("Done in --- {:.3f} seconds ---".format(time.time() - START_TIME))
logging.info("""Thanks for using WhatsGNU1.0, I hope you found it useful.
Please cite Prokka 'Seemann 2014, Bioinformatics;30(14):2068-9' if you use WhatsGNU1.0.
Please also cite Roary 'Page et al. 2015, Bioinformatics;31(22):3691-3693' if you use WhatsGNU1.0.
Please also cite BLAST+ 'Camacho et al. 2009, BMC Bioinformatics;10:421' if you use WhatsGNU1.0.
Please cite Staphopia 'Petit RA III and Read TD 2018, PeerJ;6:e5261' if you use Staphopia S. aureus Database.
Please cite Enterobase 'Alikhan NF et al. 2018, PLoS Genetics;14(4):e1007261' if you use Enterobase S. enterica Database.
The manual is extensive and available to read at https://github.com/ahmedmagds/WhatsGNU
If you have problems, please file at https://github.com/ahmedmagds/WhatsGNU/issues""")
