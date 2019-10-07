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

# DATE CREATED: June 27, 2019

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
import shutil
import urllib.request as request
from contextlib import closing
import argparse

PARSER = argparse.ArgumentParser(
    prog="WhatsGNU_get_GenBank_assemblies.py",
    description="Get GenBank assemblies (faa or/and fna) for WhatsGNU v1.0",
)
PARSER.add_argument(
    "-f",
    "--faa",
    help="protein faa file from GenBank",
    action="store_true",
)
PARSER.add_argument(
    "-c",
    "--contigs",
    help="genomic fna file from GenBank",
    action="store_true",
)
PARSER.add_argument(
    "-r",
    "--remove",
    help="remove assembly_summary_genbank.txt after done",
    action="store_true",
)
PARSER.add_argument(
    "list",
    type=str,
    help="a list.txt file of GenBank accession numbers (GCA#.#)",
)
PARSER.add_argument(
    "output_folder",
    type=str,
    help="give name for output folder to be created",
)
if len(sys.argv) == 1:
    PARSER.print_help()
    sys.exit(0)
ARGS = PARSER.parse_args()
OS_SEPARATOR = os.sep
#####create results folder######
try:
    os.mkdir(ARGS.output_folder)
    RESULTS_FOLDER = ARGS.output_folder + OS_SEPARATOR
    print("created folder ({}) as requested".format(RESULTS_FOLDER))
except:
    PARSER.exit(
        status=0,
        message="the folder you specified exists, please specify another folder name\n",)
if ARGS.faa and ARGS.contigs:
    CONTIGS_FOLDER = RESULTS_FOLDER + 'contigs'
    PROTEINS_FOLDER = RESULTS_FOLDER + 'proteins'
    os.mkdir(CONTIGS_FOLDER)
    os.mkdir(PROTEINS_FOLDER)
###########NCBI FTP List#############
TXT_FILE_NAME = (RESULTS_FOLDER + "assembly_summary_genbank.txt")
TXT_FILE_NAME_OBJECT = open(TXT_FILE_NAME, "wb")
with closing(request.urlopen(
    'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt')) as file_link:
    shutil.copyfileobj(file_link, TXT_FILE_NAME_OBJECT)
TXT_FILE_NAME_OBJECT.close()
print("Downloaded the assembly_summary_genbank.txt")
TXT_FILE_NAME_OBJECT = open(TXT_FILE_NAME, "r", encoding='utf-8')
GCA_DICT = {}
for line in TXT_FILE_NAME_OBJECT:
    if line.startswith("#"):
        continue
    GCA_DICT[line.split("\t")[0]] = line.split("\t")[19] #folder ftp url
TXT_FILE_NAME_OBJECT.close()
##############GCA list#################
LIST_FILE_OBJECT = open(ARGS.list, 'r')
GCA_LIST = []
for line in LIST_FILE_OBJECT:
    line = line.rstrip()
    GCA_LIST.append(line)

########################################
for GCA in GCA_LIST:
    try:
        GCA_ftp_url = GCA_DICT[GCA] #url from assembly_summary_genbank.txt
        GCA_url, GCA_acc = GCA_ftp_url.rsplit("/", 1) #splitting url to get accession
        GCA_file_url = GCA_url + '/'+  GCA_acc +'/'+ GCA_acc #joining them again
        print('Working on {}'.format(GCA))
        if ARGS.faa and not ARGS.contigs:
            faa_file = RESULTS_FOLDER + GCA_acc + '.faa.gz'
            faa_url = GCA_file_url + '_protein.faa.gz'
            faa_file_object = open(faa_file, 'wb')
            with closing(request.urlopen(faa_url)) as file_link:
                shutil.copyfileobj(file_link, faa_file_object)
            faa_file_object.close()
        if ARGS.contigs and not ARGS.faa:
            contigs_file = RESULTS_FOLDER + GCA_acc + '.fna.gz'
            contigs_url = GCA_file_url + '_genomic.fna.gz'
            contigs_file_object = open(contigs_file, 'wb')
            with closing(request.urlopen(contigs_url)) as file_link:
                shutil.copyfileobj(file_link, contigs_file_object)
            contigs_file_object.close()
        if ARGS.faa and ARGS.contigs:
            faa_file = PROTEINS_FOLDER + OS_SEPARATOR + GCA_acc + '.faa.gz'
            faa_url = GCA_file_url + '_protein.faa.gz'
            faa_file_object = open(faa_file, 'wb')
            with closing(request.urlopen(faa_url)) as file_link:
                shutil.copyfileobj(file_link, faa_file_object)
            faa_file_object.close()
            contigs_file = CONTIGS_FOLDER + OS_SEPARATOR + GCA_acc + '.fna.gz'
            contigs_url = GCA_file_url + '_genomic.fna.gz'
            contigs_file_object = open(contigs_file, 'wb')
            with closing(request.urlopen(contigs_url)) as file_link:
                shutil.copyfileobj(file_link, contigs_file_object)
            contigs_file_object.close()
    except:
        print("Could not get the ftp url for {} from assembly_summary_genbank.txt".format(GCA))
if ARGS.remove:
    os.remove(TXT_FILE_NAME)
