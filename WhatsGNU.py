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
from collections import defaultdict
from collections import Counter
import time
import pickle
import re
import argparse
import logging

START_TIME = time.time()

PARSER = argparse.ArgumentParser(
    prog="WhatsGNU.py",
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
    help="you have to provide path to faa file format to \
create compressed database in txt and pickle formats",
)
GROUP.add_argument(
    "-d",
    "--database",
    type=str,
    help="you have to provide path to your processed database",
)
PARSER.add_argument(
    "-o",
    "--output_folder",
    type=str,
    help="give name for output folder to be created for \
results (default: timestamped WhatsGNU_results_v1 in the current directory)",
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
    "-b",
    "--hits",
    help="create a file of each protein with all hits from the database,\
large file (~ 1 Gb for 3000 pts)",
    action="store_true",
)
PARSER.add_argument(
    "-c",
    "--CCST_typing",
    type=str,
    nargs="?",
    const="Saureus_CC_ST_names_frequencies_012119.csv",
    help="get the CC/ST composition of your hits (Note: Works only for S.aureus)",
)
PARSER.add_argument(
    "-v",
    "--version",
    help="print version and exit",
    action="version",
    version="%(prog)s 1.0",
)
PARSER.add_argument(
    "query_faa", type=str, help="Query protein FASTA file to analyze (.faa)"
)
ARGS = PARSER.parse_args()

if bool(vars(ARGS)["strainhits"]) and not bool(vars(ARGS)["topgenomes"]):
    PARSER.exit(status=0, message="Error: You have to use -s with -t\n")

#####variables######
SEQUENCES_DICT = {}
######Logging##################
TIMESTR = time.strftime("%Y%m%d_%H%M%S")
LOG_FILE = "WhatsGNU_v1_" + TIMESTR + ".log"
logging.basicConfig(
    filename=LOG_FILE,
    format="%(asctime)s:%(levelname)s:%(message)s",
    level=logging.DEBUG,
)

#####faa files to be processed######
QUERY = ARGS.query_faa
QUERY_LIST = []
try:
    for file in os.listdir(QUERY):
        if file.endswith(".faa"):
            QUERY_LIST.append(QUERY + file)
    print("You provided folder of {} faa files to be processed".format(len(QUERY_LIST)))
    logging.info(
        "You provided folder of {} faa files to be processed".format(len(QUERY_LIST))
    )
    if len(QUERY_LIST) == 0:
        logging.error("The directory did not have any faa files")
        PARSER.exit(
            status=0,
            message="The directory did not have any faa files\n",
        )
except:
    if QUERY.endswith(".faa"):
        QUERY_LIST.append(QUERY)
        print("You provided one faa file to be processed")
        logging.info("You provided one faa file to be processed")
    else:
        logging.error(
            "You did not provide single faa file or path to directory with multiple faa files"
        )
        PARSER.exit(
            status=0,
            message="You did not provide single faa file or path to directory with multiple faa files\n",
        )

#####Run database for first time#####
if ARGS.mkdatabase:
    SEQUENCES_DICT_DB = defaultdict(list)
    DB_SEQUENCE_STRING = ""
    DB_SEQUENCE_INFO = ""
    if ARGS.mkdatabase.endswith(".faa"):
        DATABASEFILE_OBJECT = open(ARGS.mkdatabase, "r")
        DB_LINE_CHECK = DATABASEFILE_OBJECT.readline()
        if not DB_LINE_CHECK.startswith(">"):
            logging.error("{} is not in a FASTA format".format(ARGS.mkdatabase))
            PARSER.exit(status=0, message="Database is not in a FASTA format\n")
        DATABASEFILE_OBJECT.seek(0)
        try:
            for line in DATABASEFILE_OBJECT:
                line = line.rstrip()
                if line.startswith(">"):
                    if len(DB_SEQUENCE_STRING) > 0:
                        SEQUENCES_DICT_DB[DB_SEQUENCE_STRING].append(DB_SEQUENCE_INFO)
                        DB_SEQUENCE_STRING = ""
                    DB_SEQUENCE_INFO = line.lstrip(">")
                else:
                    DB_SEQUENCE_STRING += line
            SEQUENCES_DICT_DB[DB_SEQUENCE_STRING].append(DB_SEQUENCE_INFO)
            if bool(SEQUENCES_DICT_DB):
                print(
                    "processed database of {} proteins to compressed dictionary in --- {:.3f} seconds ---".format(
                        len(SEQUENCES_DICT_DB), time.time() - START_TIME
                    )
                )
                logging.info(
                    "processed database of {} proteins to compressed dictionary in --- {:.3f} seconds ---".format(
                        len(SEQUENCES_DICT_DB), time.time() - START_TIME
                    )
                )
                DATABASEFILE_OBJECT.close()
        except:
            logging.error(
                "Cannot process the faa file provided to a compressed dictionary"
            )
            PARSER.exit(
                status=0,
                message="Cannot process the faa file provided to a compressed dictionary\n",
            )
        #####Save database for first time#####
        TIMESTR2 = time.strftime("%Y%m%d_%H%M")
        try:
            TXT_FILE_NAME = (
                ARGS.mkdatabase.split(".faa")[0]
                + "_compressed_database_"
                + TIMESTR2
                + ".txt"
            )
            OUTPUT_FILE_DB = open(TXT_FILE_NAME, "w")
            for record in SEQUENCES_DICT_DB:
                OUTPUT_FILE_DB.write(
                    "{}\t{}\n".format(record, "._/".join(SEQUENCES_DICT_DB[record]))
                )
            OUTPUT_FILE_DB.close()
            print(
                "saved database ({}) of {} proteins as txt file in --- {:.3f} seconds ---".format(
                    TXT_FILE_NAME, len(SEQUENCES_DICT_DB), time.time() - START_TIME
                )
            )
            logging.info(
                "saved database ({}) of {} proteins as txt file in --- {:.3f} seconds ---".format(
                    TXT_FILE_NAME, len(SEQUENCES_DICT_DB), time.time() - START_TIME
                )
            )
        except:
            print(
                "cannot save compressed database as txt file, this time will be ok as the compressed dictionary will be used"
            )
            logging.critical(
                "cannot save compressed database as txt file, this time will be ok as the compressed dictionary will be used"
            )
        #########################################
        try:
            PICKLE_FILE_NAME = (
                ARGS.mkdatabase.split(".faa")[0]
                + "_compressed_database_"
                + TIMESTR2
                + ".pickle"
            )
            PICKLE_OUT = open(PICKLE_FILE_NAME, "wb")
            pickle.dump(SEQUENCES_DICT_DB, PICKLE_OUT)
            PICKLE_OUT.close()
            print(
                "saved pickle_dict ({}) in --- {:.3f} seconds ---".format(
                    PICKLE_FILE_NAME, time.time() - START_TIME
                )
            )
            logging.info(
                "saved pickle_dict ({}) in --- {:.3f} seconds ---".format(
                    PICKLE_FILE_NAME, time.time() - START_TIME
                )
            )
        except:
            print("cannot save pickle file")
            logging.warning("cannot save pickle file")
    else:
        logging.error(
            "The file you provided to be processed as a database is not having extension .faa file"
        )
        PARSER.exit(
            status=0,
            message="The file you provided to be processed as a database\
        is not having extension .faa file\n",
        )

#####load database file#########
if ARGS.mkdatabase:
    SEQUENCES_DICT = SEQUENCES_DICT_DB
    print(
        "As you just created a compressed database using -m, it will be used this time, next time provide the database using -d"
    )
    logging.info(
        "As you just created a compressed database using -m, it will be used this time, next time provide the database using -d"
    )
elif ARGS.database:
    if ARGS.database.endswith(".pickle"):
        try:
            PICKLE_IN = open(ARGS.database, "rb")
            SEQUENCES_DICT = pickle.load(PICKLE_IN)
            if bool(SEQUENCES_DICT):
                print(
                    "opened previously created pickle_dict in --- {:.3f} seconds ---".format(
                        time.time() - START_TIME
                    )
                )
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
    elif ARGS.database.endswith(".txt"):
        try:
            SEQUENCES_DICT = {}
            TXT_DB_FILEOBJECT = open(ARGS.database, "r")
            print(
                "opened previously created compressed txt database in --- {:.3f} seconds ---".format(
                    time.time() - START_TIME
                )
            )
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
                print(
                    "processed compressed txt database to dictionary in --- {:.3f} seconds ---".format(
                        time.time() - START_TIME
                    )
                )
                logging.info(
                    "processed compressed txt database to dictionary in --- {:.3f} seconds ---".format(
                        time.time() - START_TIME
                    )
                )
        except:
            logging.error(
                "The compressed database txt file you provided is empty or corrupted"
            )
            PARSER.exit(
                status=0,
                message="The compressed database txt file you provided is empty or corrupted\n",
            )
    else:
        logging.error(
            "No proper compressed database (file.pickle or file.txt) was provided using -d"
        )
        PARSER.exit(
            status=0,
            message="No proper compressed database (file.pickle or file.txt) was provided using -d\n",
        )
else:
    logging.error(
        "Neither you created new database using -m (file.faa) nor proper database (file.pickle or file.txt) was provided using -d"
    )
    PARSER.exit(
        status=0,
        message="Neither you created new database using -m (file.faa) nor proper database (file.pickle or file.txt) was provided using -d\n",
    )

#####create results folder######
if ARGS.output_folder:
    try:
        os.mkdir(ARGS.output_folder)
        RESULTS_FOLDER = ARGS.output_folder + "/"
        print("created folder ({}) as requested".format(RESULTS_FOLDER))
        logging.info(
            "created folder ({}) as requested".format(RESULTS_FOLDER)
        )
    except:
        os.mkdir("WhatsGNU_results_v1_" + TIMESTR)
        RESULTS_FOLDER = "./WhatsGNU_results_v1_{}/".format(TIMESTR)
        print(
            "the folder you specified exists so created default results folder ({})".format(
                RESULTS_FOLDER
            )
        )
        logging.warning(
            "the folder you specified exists so created default results folder ({})".format(
                RESULTS_FOLDER
            )
        )
else:
    os.mkdir("WhatsGNU_results_v1_" + TIMESTR)
    RESULTS_FOLDER = "./WhatsGNU_results_v1_{}/".format(TIMESTR)
    print("created default results folder({})".format(RESULTS_FOLDER))
    logging.info("created default results folder({})".format(RESULTS_FOLDER))
#########whatsgnu###########
for QUERYFILE in QUERY_LIST:
    QUERYFILE_OBJECT = open(QUERYFILE, "r")
    line_check = QUERYFILE_OBJECT.readline()
    if not line_check.startswith(">"):
        logging.error("Not a FASTA file: {}".format(QUERYFILE))
        PARSER.exit(status=0, message="Not a FASTA file\n")
    QUERYFILE_OBJECT.seek(0)
    file_hits = (
        RESULTS_FOLDER
        + (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0]
        + "_WhatsGNU_hits_v1.txt"
    )
    file_report = (
        RESULTS_FOLDER
        + (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0]
        + "_WhatsGNU_report_v1.txt"
    )
    output_file_report = open(file_report, "w")
    print(
        "opened report file for {}".format(
            (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0]
        )
    )
    logging.info(
        "opened report file for {}".format(
            (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0]
        )
    )
    CC_ST_names_list = []
    CC_ST_frequencies_list = []
    if ARGS.CCST_typing:  # get CC/ST composition of hits
        CC_ST_names_freq_list_file = open(ARGS.CCST_typing, "r")
        for line in CC_ST_names_freq_list_file:
            line = line.rstrip()
            CC_ST, CC_ST_freq = line.split(",")
            CC_ST_names_list.append(CC_ST)
            CC_ST_frequencies_list.append(int(CC_ST_freq))
        output_file_report.write(
            "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                "protein",
                "length",
                "function",
                "sequence",
                "GNU score",
                "\t".join(CC_ST_names_list),
            )
        )
        print("processed CC/ST list")
        logging.info("processed CC/ST list")
    else:
        output_file_report.write(
            "{}\t{}\t{}\t{}\t{}\n".format(
                "protein", "length", "function", "sequence", "GNU score"
            )
        )
    if ARGS.hits:  # output big (roughly 1GB) hits file
        output_file_hits = open(file_hits, "w")
        output_file_hits.write("{}\t{}\n".format("protein_query", "hits"))
        print(
            "opened hits file for {} as per your request of -b".format(
                (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0]
            )
        )
        logging.info(
            "opened hits file for {} as per your request of -b".format(
                (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0]
            )
        )
    processed_counter = 0
    CC_ST_count_list = []
    CC_ST_percentage_list = []
    all_hits_list = []
    zeros_list = ["0.00"] * 57
    sequence_info = ""
    sequence_string = ""
    if ARGS.topgenomes:  # get top 10 genomes with hits
        compiled_pattern = re.compile(r"\w+_GCA_\d+.\d_")
    for line in QUERYFILE_OBJECT:
        line = line.rstrip()
        if line.startswith(">"):
            if len(sequence_string) > 0:
                processed_counter += 1
                try:
                    found_strain_name = []
                    CC_ST_count_list = []
                    CC_ST_percentage_list = []
                    ids_list = SEQUENCES_DICT[sequence_string]
                    ids_tabbed_string = "\t".join(ids_list)
                    try:
                        function = re.search(r".+_GCA_\d+.\d+_(.+)", ids_list[0])
                        function_cleaned = function.group(1)
                    except:
                        function_cleaned = ids_list[0]
                    if ARGS.topgenomes:
                        all_hits_list.extend(
                            re.findall(compiled_pattern, ids_tabbed_string)
                        )
                    if ARGS.hits:
                        output_file_hits.write(
                            "{}\t{}\n".format(sequence_info, ids_tabbed_string)
                        )
                    if ARGS.CCST_typing:
                        for CC_ST in CC_ST_names_list:
                            CC_ST_count_list.append(
                                ids_tabbed_string.count("_" + CC_ST + "_")
                            )
                        CC_ST_percentage_list = [
                            str("{:.2f}".format(CC_ST_count * 100 / CC_ST_freq))
                            for CC_ST_count, CC_ST_freq in zip(
                                CC_ST_count_list, CC_ST_frequencies_list
                            )
                        ]
                        output_file_report.write(
                            "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                sequence_info,
                                len(sequence_string),
                                function_cleaned,
                                sequence_string,
                                len(ids_list),
                                "\t".join(CC_ST_percentage_list),
                            )
                        )
                    else:
                        output_file_report.write(
                            "{}\t{}\t{}\t{}\t{}\n".format(
                                sequence_info,
                                len(sequence_string),
                                function_cleaned,
                                sequence_string,
                                len(ids_list),
                            )
                        )
                    print(
                        "processed protein {} of {} in {:.3F}".format(
                            processed_counter,
                            (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0],
                            time.time() - START_TIME,
                        )
                    )
                    logging.info(
                        "processed protein {} of {} in {:.3F}".format(
                            processed_counter,
                            (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0],
                            time.time() - START_TIME,
                        )
                    )
                except:
                    if ARGS.hits:
                        output_file_hits.write("{}\tno_hits\n".format(sequence_info))
                    if ARGS.CCST_typing:
                        output_file_report.write(
                            "{}\t{}\tNA\t{}\t{}\n".format(
                                sequence_info,
                                str(len(sequence_string)),
                                sequence_string,
                                "\t".join(zeros_list),
                            )
                        )
                    else:
                        output_file_report.write(
                            "{}\t{}\tNA\t{}\t0\n".format(
                                sequence_info,
                                str(len(sequence_string)),
                                sequence_string,
                            )
                        )
                    print(
                        "processed protein {} of {} in {:.3F}".format(
                            processed_counter,
                            (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0],
                            time.time() - START_TIME,
                        )
                    )
                    logging.info(
                        "processed protein {} of {} in {:.3F}".format(
                            processed_counter,
                            (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0],
                            time.time() - START_TIME,
                        )
                    )
                sequence_string = ""
            sequence_info = line.lstrip(">")
        else:
            sequence_string += line
    processed_counter += 1
    try:
        found_strain_name = []
        CC_ST_count_list = []
        CC_ST_percentage_list = []
        ids_list = SEQUENCES_DICT[sequence_string]
        ids_tabbed_string = "\t".join(ids_list)
        try:
            function = re.search(
                r".+_GCA_\d+.\d+_(.+)", ids_list[0]
            )  # for programs other than prokka (like RAST)
            function_cleaned = function.group(1)
        except:
            function_cleaned = ids_list[0]
        if ARGS.topgenomes:
            all_hits_list.extend(re.findall(compiled_pattern, ids_tabbed_string))
        if ARGS.hits:
            output_file_hits.write("{}\t{}\n".format(sequence_info, ids_tabbed_string))
        if ARGS.CCST_typing:
            for CC_ST in CC_ST_names_list:
                CC_ST_count_list.append(ids_tabbed_string.count("_" + CC_ST + "_"))
            CC_ST_percentage_list = [
                str("{:.2f}".format(CC_ST_count * 100 / CC_ST_freq))
                for CC_ST_count, CC_ST_freq in zip(
                    CC_ST_count_list, CC_ST_frequencies_list
                )
            ]
            output_file_report.write(
                "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    sequence_info,
                    len(sequence_string),
                    function_cleaned,
                    sequence_string,
                    len(ids_list),
                    "\t".join(CC_ST_percentage_list),
                )
            )
        else:
            output_file_report.write(
                "{}\t{}\t{}\t{}\t{}\n".format(
                    sequence_info,
                    len(sequence_string),
                    function_cleaned,
                    sequence_string,
                    len(ids_list),
                )
            )
        print(
            "processed protein {} of {} in {:.3F}".format(
                processed_counter,
                (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0],
                time.time() - START_TIME,
            )
        )
        logging.info(
            "processed protein {} of {} in {:.3F}".format(
                processed_counter,
                (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0],
                time.time() - START_TIME,
            )
        )
    except:
        if ARGS.hits:
            output_file_hits.write("{}\tno_hits\n".format(sequence_info))
        if ARGS.CCST_typing:
            output_file_report.write(
                "{}\t{}\tNA\t{}\t{}\n".format(
                    sequence_info,
                    str(len(sequence_string)),
                    sequence_string,
                    "\t".join(zeros_list),
                )
            )
        else:
            output_file_report.write(
                "{}\t{}\tNA\t{}\t0\n".format(
                    sequence_info, str(len(sequence_string)), sequence_string
                )
            )
        print(
            "processed protein {} of {} in {:.3F}".format(
                processed_counter,
                (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0],
                time.time() - START_TIME,
            )
        )
        logging.info(
            "processed protein {} of {} in {:.3F}".format(
                processed_counter,
                (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0],
                time.time() - START_TIME,
            )
        )
    if ARGS.topgenomes:  # get top 10 genomes with hits
        file_topgenomes = (
            RESULTS_FOLDER
            + (QUERYFILE.rsplit("/", 1)[-1]).split(".faa")[0]
            + "_WhatsGNU_topgenomes_v1.txt"
        )
        C = Counter(all_hits_list)
        if ARGS.strainhits:
            print(
                "No of hits from strain {} is {}".format(
                    ARGS.strainhits, C[ARGS.strainhits + "_"]
                )
            )
            logging.info(
                "No of hits from strain {} is {}".format(
                    ARGS.strainhits, C[ARGS.strainhits + "_"]
                )
            )
        most_occur = C.most_common(10)
        with open(file_topgenomes, "w") as output_file_topgenomes:
            output_file_topgenomes.write(
                "WhatsGNU found {} total hits from the database in --- {:.3f} seconds ---\n".format(
                    len(all_hits_list), time.time() - START_TIME
                )
            )
            output_file_topgenomes.write(
                "\n".join("{}\t{}".format(x[0].rstrip("_"), x[1]) for x in most_occur)
            )
        print("Found top 10 genomes with hits")
        logging.info("Found top 10 genomes with hits")
    if ARGS.hits:
        output_file_hits.close()
    output_file_report.close()
    QUERYFILE_OBJECT.close()
print("Done in --- {:.3f} seconds ---".format(time.time() - START_TIME))
logging.info("Done in --- {:.3f} seconds ---".format(time.time() - START_TIME))
