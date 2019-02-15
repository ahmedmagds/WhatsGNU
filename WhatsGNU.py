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

import sys
import os
from collections import defaultdict
from collections import Counter
import time
import pickle
import re
import argparse
import logging

start_time = time.time()

parser = argparse.ArgumentParser(prog='WhatsGNU.py', description="WhatsGNU v1.0 utilizes the natural\
 variation in public databases to rank protein sequences based on the number of observed exact protein\
 matches (the GNU score) in all known genomes of a particular species. It generates a report for all the\
 proteins in your query in seconds")
group = parser.add_mutually_exclusive_group()
group.add_argument("-m", "--mkdatabase", type=str, help="you have to provide path to faa file format to \
create compressed database in txt and pickle formats")
group.add_argument("-d", "--database", type=str, help="you have to provide path to your processed database")
parser.add_argument("-o", "--output_folder", type=str, help="give name for output folder to be created for \
results (default: timestamped WhatsGNU_results_v1 in the current directory)")
parser.add_argument("-t", "--topgenomes", help="create a file of top 10 genomes with hits", action="store_true")
parser.add_argument("-s", "--strainhits", type=str, help="check how many hits you get from a particular strain,\
it has to be used with -t")
parser.add_argument("-b", "--hits", help="create a file of each protein with all hits from the database,\
large file (~ 1 Gb for 3000 pts)", action="store_true")
parser.add_argument("-c", "--CCST_typing", type=str, nargs='?', const="Saureus_CC_ST_names_frequencies_012119.csv",
help="get the CC/ST composition of your hits (Note: Works only for S.aureus)")
parser.add_argument("-v", "--version", help="print version and exit", action="version", version='%(prog)s 1.0')
parser.add_argument("query_faa", type=str, help="Query protein FASTA file to analyze (.faa)")
args = parser.parse_args()

if bool(vars(args)['strainhits']) and not bool(vars(args)['topgenomes']):
    parser.exit(status=0, message='Error: You have to use -s with -t\n')

#####variables######
_format = "fasta"
sequences_dict = {}
sequences_dict_d = defaultdict(list)
db_sequence_string = ''
db_sequence_info = ''
######Logging##################
timestr = time.strftime("%Y%m%d_%H%M%S")
log_file = 'WhatsGNU_v1_'+timestr+'.log'
logging.basicConfig(filename=log_file,format='%(asctime)s:%(levelname)s:%(message)s',level=logging.DEBUG)

#####faa files to be processed######
query = args.query_faa
query_list = []
try:
    for file in os.listdir(query):
        if file.endswith(".faa"):
            query_list.append(query+file)
    print("You provided folder of {} faa files to be processed".format(len(query_list)))
    logging.info("You provided folder of {} faa files to be processed".format(len(query_list)))
except:
    if query.endswith(".faa"):
        query_list.append(query)
        print("You provided one faa file to be processed")
        logging.info("You provided one faa file to be processed")
    else:
        logging.error("You did not provide single faa file or path to directory with multiple faa files")
        parser.exit(status=0, message='You did not provide single faa file or path to directory with \
        multiple faa files\n')

#####Run database for first time#####
if args.mkdatabase:
    if args.mkdatabase.endswith(".faa"):
        try:
            databasefile_object = open(args.mkdatabase,'r')
            for line in databasefile_object:
                line = line.rstrip()
                if line.startswith('>'):
                    if len(db_sequence_string) > 0:
                        sequences_dict_d[db_sequence_string].append(db_sequence_info)
                        db_sequence_string = ''
                    db_sequence_info = line.lstrip('>')
                else:
                    db_sequence_string += line
            sequences_dict_d[db_sequence_string].append(db_sequence_info)
            if bool(sequences_dict_d):
                print("processed database of {} proteins to compressed dictionary in --- {:.3f} seconds ---".format(
                len(sequences_dict_d),time.time() - start_time))
                logging.info("processed database of {} proteins to compressed dictionary in --- {:.3f} seconds ---".format(
                len(sequences_dict_d),time.time() - start_time))
                databasefile_object.close()
        except:
            logging.error("Cannot process the faa file provided to a compressed dictionary")
            parser.exit(status=0, message='Cannot process the faa file provided to a compressed dictionary\n')
        #####Save database for first time#####
        timestr2 = time.strftime("%Y%m%d_%H%M")
        try:
            txt_file_name= args.mkdatabase.split(".faa")[0]+"_compressed_database_"+timestr2+".txt"
            output_file_3 = open(txt_file_name,'w')
            for record in sequences_dict_d:
                output_file_3.write("{}\t{}\n".format(record, '._/'.join(sequences_dict_d[record])))
            output_file_3.close()
            print("saved database ({}) of {} proteins as txt file in --- {:.3f} seconds ---".format(
            txt_file_name,len(sequences_dict_d),time.time() - start_time))
            logging.info("saved database ({}) of {} proteins as txt file in --- {:.3f} seconds ---".format(
            txt_file_name,len(sequences_dict_d),time.time() - start_time))
        except:
            print("cannot save compressed database as txt file,\ 
            this time will be ok as the compressed dictionary will be used")
            logging.critical("cannot save compressed database as txt file,\
            this time will be ok as the compressed dictionary will be used")
        #########################################
        try:
            pickle_file_name = args.mkdatabase.split(".faa")[0]+"_compressed_database_"+timestr2+".pickle"
            pickle_out = open(pickle_file_name,"wb")
            pickle.dump(sequences_dict_d, pickle_out)
            pickle_out.close()
            print("saved pickle_dict ({}) in --- {:.3f} seconds ---".format(pickle_file_name,time.time() - start_time))
            logging.info("saved pickle_dict ({}) in --- {:.3f} seconds ---".format(
            pickle_file_name,time.time() - start_time))
        except:
            print("cannot save pickle file")
            logging.warning("cannot save pickle file")
    else:
        logging.error("The file you provided to be processed as a database is not having extension .faa file")
        parser.exit(status=0, message='The file you provided to be processed as a database\
        is not having extension .faa file\n')

#####load database file#########
if args.mkdatabase:
    sequences_dict = sequences_dict_d
    print("""As you just created a compressed database using -m, it will be used this time, next time provide the database
    using -d""")
    logging.info("""As you just created a compressed database using -m, it will be used this time, next time provide the
    database using -d""")
elif args.database:
    if args.database.endswith(".pickle"):
        try:
            pickle_in = open(args.database,"rb")
            sequences_dict = pickle.load(pickle_in)
            if bool(sequences_dict):
                print("opened previously created pickle_dict in --- {:.3f} seconds ---".format(time.time() - start_time))
                logging.info("opened previously created pickle_dict in --- {:.3f} seconds ---".format(
                time.time() - start_time))
        except:
            logging.error("The compressed database pickle_dict file you provided is empty or corrupted")
            parser.exit(status=0, message='The compressed database pickle_dict file you provided is empty or corrupted\n')
    elif args.database.endswith(".txt"):
        try:
            sequences_dict = {}
            database_fo = open(args.database, 'r')
            print('opened previously created compressed txt database in --- {:.3f} seconds ---'.format(
            time.time() - start_time))
            logging.info('opened previously created compressed txt database in --- {:.3f} seconds ---'.format(
            time.time() - start_time))
            for line in database_fo:
                line = line.rstrip()
                seq,ids = line.split('\t')
                listids = ids.split('._/')
                sequences_dict[seq] = listids
            if bool(sequences_dict):
                print("processed compressed txt database to dictionary in --- {:.3f} seconds ---".format(
                time.time() - start_time))
                logging.info("processed compressed txt database to dictionary in --- {:.3f} seconds ---".format(
                time.time() - start_time))
        except:
            logging.error("The compressed database txt file you provided is empty or corrupted")
            parser.exit(status=0, message='The compressed database txt file you provided is empty or corrupted\n')
    else:
        logging.error("No proper compressed database (file.pickle or file.txt) was provided using -d")
        parser.exit(status=0, message='No proper compressed database (file.pickle or file.txt) was provided using -d\n')
else:
    logging.error("""Neither you created new database using -m (file.faa) nor proper database (file.pickle or file.txt)
    was provided using -d""")
    parser.exit(status=0, message="""Neither you created new database using -m (file.faa) nor proper database (file.pickle
    or file.txt) was provided using -d\n""")

#####create results folder######
if args.output_folder:
    try:
        os.mkdir(args.output_folder)
        results_folder = "./"+args.output_folder+"/"
        print("created folder ({}) in the current directory".format(results_folder))
        logging.info("created folder ({}) in the current directory".format(results_folder))
    except:
        os.mkdir("WhatsGNU_results_v1_"+timestr)
        results_folder = "./WhatsGNU_results_v1_{}/".format(timestr)
        print("the folder you specified exists so created default results folder ({})".format(results_folder))
        logging.warning("the folder you specified exists so created default results folder ({})".format(results_folder))
else:
    os.mkdir("WhatsGNU_results_v1_"+timestr)
    results_folder = "./WhatsGNU_results_v1_{}/".format(timestr)
    print("created default results folder({})".format(results_folder))
    logging.info("created default results folder({})".format(results_folder))
#########whatsgnu###########
for queryfile in query_list:
    queryfile_object = open(queryfile,'r')
    line_check = queryfile_object.readline()
    if not line_check.startswith(">"):
        raise TypeError("Not a FASTA file: ", queryfile)
    queryfile_object.seek(0)
    file_hits = results_folder + (queryfile.rsplit('/', 1)[-1]).split(".faa")[0] + "_WhatsGNU_hits_v1.txt"
    file_report = results_folder + (queryfile.rsplit('/', 1)[-1]).split(".faa")[0] + "_WhatsGNU_report_v1.txt"
    output_file_2 = open(file_report,'w')
    print("opened report file for {}".format((queryfile.rsplit('/', 1)[-1]).split(".faa")[0]))
    logging.info("opened report file for {}".format((queryfile.rsplit('/', 1)[-1]).split(".faa")[0]))
    CC_ST_names_list = []
    CC_ST_frequencies_list = []
    if args.CCST_typing: #get CC/ST composition of hits
        CC_ST_names_freq_list_file = open(args.CCST_typing, "r")
        for line in CC_ST_names_freq_list_file:
            line = line.rstrip()
            CC_ST,CC_ST_freq = line.split(',')
            CC_ST_names_list.append(CC_ST)
            CC_ST_frequencies_list.append(int(CC_ST_freq))
        output_file_2.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
        'protein', 'length', 'function','sequence','GNU score', '\t'.join(CC_ST_names_list)))
        print("processed CC/ST list")
        logging.info("processed CC/ST list")
    output_file_2.write("{}\t{}\t{}\t{}\t{}\n".format('protein', 'length', 'function','sequence','GNU score'))
    if args.hits: #output big (roughly 1GB) hits file
        output_file = open(file_hits, 'w')
        output_file.write("{}\t{}\n".format("protein_query", "hits"))
        print("opened hits file for {} as per your request of -b".format((queryfile.rsplit('/', 1)[-1]).split(".faa")[0]))
        logging.info("opened hits file for {} as per your request of -b".format(
        (queryfile.rsplit('/', 1)[-1]).split(".faa")[0]))
    processed_counter = 0
    CC_ST_count_list = []
    CC_ST_percentage_list = []
    all_hits_list = []
    zeros_list = ['0.00'] * 57
    sequence_info =''
    sequence_string = ''
    if args.topgenomes: #get top 10 genomes with hits
        compiled_pattern = re.compile(r"\w+_GCA_\d+.\d_")
    for line in queryfile_object:
        line = line.rstrip()
        if line.startswith('>'):
            if len(sequence_string) > 0:
                processed_counter += 1
                try:
                    found_strain_name = []
                    CC_ST_count_list = []
                    CC_ST_percentage_list = []
                    ids_list = sequences_dict[sequence_string]
                    ids_tabbed_string = '\t'.join(ids_list)
                    try:
                        function = re.search(r".+_GCA_\d+.\d+_(.+)", ids_list[0])
                        function_cleaned = function.group(1)
                    except:
                        function_cleaned = ids_list[0]
                    if args.topgenomes:
                        all_hits_list.extend(re.findall(compiled_pattern,ids_tabbed_string))
                    if args.hits:
                        output_file.write("{}\t{}\n".format(sequence_info, ids_tabbed_string))
                    if args.CCST_typing:
                        for CC_ST in CC_ST_names_list:
                            CC_ST_count_list.append(ids_tabbed_string.count('_'+CC_ST+'_'))
                        CC_ST_percentage_list = [str('{:.2f}'.format(
                        CC_ST_count*100/CC_ST_freq)) for CC_ST_count,CC_ST_freq in zip(
                        CC_ST_count_list,CC_ST_frequencies_list)]
                        output_file_2.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        sequence_info,len(sequence_string),function_cleaned,sequence_string,len(ids_list),
                        '\t'.join(CC_ST_percentage_list)))
                    else:
                        output_file_2.write("{}\t{}\t{}\t{}\t{}\n".format(
                        sequence_info,len(sequence_string),function_cleaned,sequence_string,len(ids_list)))
                    print("processed protein {} of {} in {:.3F}".format(
                    processed_counter,(queryfile.rsplit('/', 1)[-1]).split(".faa")[0],time.time() - start_time))
                    logging.info("processed protein {} of {} in {:.3F}".format(
                    processed_counter,(queryfile.rsplit('/', 1)[-1]).split(".faa")[0],time.time() - start_time))
                except:
                    if args.hits:
                        output_file.write("{}\tno_hits\n".format(sequence_info))
                    if args.CCST_typing:
                        output_file_2.write("{}\t{}\tNA\t{}\t{}\n".format(
                        sequence_info,str(len(sequence_string)),sequence_string,'\t'.join(zeros_list)))
                    else:
                        output_file_2.write("{}\t{}\tNA\t{}\t0\n".format(
                        sequence_info,str(len(sequence_string)),sequence_string))
                    print("processed protein {} of {} in {:.3F}".format(
                    processed_counter,(queryfile.rsplit('/', 1)[-1]).split(".faa")[0], time.time() - start_time))
                    logging.info("processed protein {} of {} in {:.3F}".format(
                    processed_counter,(queryfile.rsplit('/', 1)[-1]).split(".faa")[0], time.time() - start_time))
                sequence_string = ''
            sequence_info = line.lstrip('>')
        else:
            sequence_string += line
    processed_counter += 1
    try:
        found_strain_name = []
        CC_ST_count_list = []
        CC_ST_percentage_list = []
        ids_list = sequences_dict[sequence_string]
        ids_tabbed_string = '\t'.join(ids_list)
        try:
            function = re.search(r".+_GCA_\d+.\d+_(.+)", ids_list[0]) #for programs other than prokka (like RAST)
            function_cleaned = function.group(1)
        except:
            function_cleaned = ids_list[0]
        if args.topgenomes:
            all_hits_list.extend(re.findall(compiled_pattern,ids_tabbed_string))
        if args.hits:
            output_file.write("{}\t{}\n".format(sequence_info, ids_tabbed_string))
        if args.CCST_typing:
            for CC_ST in CC_ST_names_list:
                CC_ST_count_list.append(ids_tabbed_string.count('_'+CC_ST+'_'))
            CC_ST_percentage_list = [str('{:.2f}'.format(
            CC_ST_count*100/CC_ST_freq)) for CC_ST_count,CC_ST_freq in zip(CC_ST_count_list,CC_ST_frequencies_list)]
            output_file_2.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
            sequence_info,len(sequence_string),function_cleaned,sequence_string,len(ids_list),
            '\t'.join(CC_ST_percentage_list)))
        else:
            output_file_2.write("{}\t{}\t{}\t{}\t{}\n".format(
            sequence_info,len(sequence_string),function_cleaned,sequence_string,len(ids_list)))
        print("processed protein {} of {} in {:.3F}".format(
        processed_counter,(queryfile.rsplit('/', 1)[-1]).split(".faa")[0],time.time() - start_time))
        logging.info("processed protein {} of {} in {:.3F}".format(
        processed_counter,(queryfile.rsplit('/', 1)[-1]).split(".faa")[0],time.time() - start_time))
    except:
        if args.hits:
            output_file.write("{}\tno_hits\n".format(sequence_info))
        if args.CCST_typing:
            output_file_2.write("{}\t{}\tNA\t{}\t{}\n".format(
            sequence_info,str(len(sequence_string)),sequence_string,'\t'.join(zeros_list)))
        else:
            output_file_2.write("{}\t{}\tNA\t{}\t0\n".format(sequence_info,str(len(sequence_string)),sequence_string))
        print("processed protein {} of {} in {:.3F}".format(
        processed_counter,(queryfile.rsplit('/', 1)[-1]).split(".faa")[0], time.time() - start_time))
        logging.info("processed protein {} of {} in {:.3F}".format(
        processed_counter,(queryfile.rsplit('/', 1)[-1]).split(".faa")[0], time.time() - start_time))
    if args.topgenomes: #get top 10 genomes with hits
        file_topgenomes = results_folder + (queryfile.rsplit('/', 1)[-1]).split(".faa")[0] + "_WhatsGNU_topgenomes_v1.txt"
        C = Counter(all_hits_list)
        if args.strainhits:
            print("No of hits from strain {} is {}".format(args.strainhits,C[args.strainhits+'_']))
            logging.info("No of hits from strain {} is {}".format(args.strainhits,C[args.strainhits+'_']))
        most_occur = C.most_common(10)
        with open(file_topgenomes, 'w') as output_file_4:
            output_file_4.write("WhatsGNU found {} total hits from the database in --- {:.3f} seconds ---\n".format(
            len(all_hits_list),time.time() - start_time))
            output_file_4.write('\n'.join('{}\t{}'.format(x[0].rstrip('_'),x[1]) for x in most_occur))
        print("Found top 10 genomes with hits")
        logging.info("Found top 10 genomes with hits")
    if args.hits:
        output_file.close()
    output_file_2.close()
    queryfile_object.close()
print("Done in --- {:.3f} seconds ---".format(time.time() - start_time))
logging.info("Done in --- {:.3f} seconds ---".format(time.time() - start_time))
