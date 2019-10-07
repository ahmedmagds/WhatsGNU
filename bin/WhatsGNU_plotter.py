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

# DATE CREATED: May 30, 2019

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
import argparse
import tempfile
import subprocess
from collections import defaultdict
from math import ceil as cl
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats as st


PARSER = argparse.ArgumentParser(
    prog="WhatsGNU_plotter.py", description="WhatsGNU_plotter script for WhatsGNU v1.0.",)
PARSER.add_argument("-hp", "--heatmap", type=str,
choices=['ortholog', 'basic'], help="heatmap of GNU scores for orthologous genes in multiple isolates",)
PARSER.add_argument("-l", "--list_genes", type=str, help="a txt file of ortholog group names from one of the WhatsGNU reports for heatmap",)
PARSER.add_argument("-q", "--fasta", type=str, help="a FASTA file of sequences for the proteins of interest for heatmap or metadata barplot",)
PARSER.add_argument("-op", "--output_blastp", help="get the output report of blastp run, it has to be used with -q", action="store_true",)
PARSER.add_argument("-d", "--strains_order", type=str, help="list of strains order for heatmap",)
PARSER.add_argument("-r", "--rarity", help="Annotate heatmap cells with OVRI(default: off)", action="store_true",)
PARSER.add_argument("-rc", "--rarity_color", type=str, help="OVRI data text color in the heatmap")
PARSER.add_argument("-fs", "--figure_size", nargs=2, type=float, help="heatmap width and height in inches w,h, respectively",)
PARSER.add_argument("-hc", "--heatmap_color", type=str, help="heatmap color")
PARSER.add_argument("-mc", "--masked_color", type=str, help="missing data color in heatmap")
PARSER.add_argument("-f", "--font_size", type=int, help="heatmap font size")
PARSER.add_argument("-t", "--title", type=str, help="title for the heatmap [Default:WhatsGNU heatmap]")
PARSER.add_argument("-mb", "--metadata_barplot", type=str, choices=['ortholog', 'basic'], help="Metadata percentage distribution for proteins in a FASTA file",)
PARSER.add_argument("-w", "--all_metadata", help="all metadata", action="store_true",)
PARSER.add_argument("-s", "--select_metadata", type=str, help="select some metadata",)
PARSER.add_argument("-x", "--histogram", help="histogram of GNU scores", action="store_true",)
PARSER.add_argument("-e", "--histogram_color", type=str, help="histogram color")
PARSER.add_argument("-b", "--histogram_bins", type=int, help="number of bins for the histograms [10]")
PARSER.add_argument("-p", "--novel_conserved", nargs=2, type=int, help="upper and lower GNU score limits for novel and conserved proteins novel_GNU_upper_limit, conserved_GNU_lower_limit, respectively [Default 10, 100]",)
PARSER.add_argument("-st", "--strains_tag_volcano", type=str, help="a csv file of the strains of the two groups to be compared with (case/control) tag",)
PARSER.add_argument("-c", "--cutoff_volcano", type=int, help="a percentage of isolates a protein must be in [Default: 100]",)
PARSER.add_argument("-cc", "--case_control_name", nargs=2, type=str, help="case and control groups' names [Default: case control]",)
PARSER.add_argument("prefix_name", type=str, help="prefix name for the the output folder and heatmap/volcano output files")
PARSER.add_argument("directory_path", type=str, help="path to directory of WhatsGNU reports")
if len(sys.argv) == 1:
    PARSER.print_help()
    sys.exit(0)
ARGS = PARSER.parse_args()
if bool(vars(ARGS)["heatmap"] == 'ortholog') and not (bool(vars(ARGS)["fasta"]) or bool(vars(ARGS)["list_genes"])):
    PARSER.exit(status=0, message="Error: You have to use -hp ortholog with -q or -l\n")
if bool(vars(ARGS)["heatmap"] == 'basic') and not bool(vars(ARGS)["fasta"]):
    PARSER.exit(status=0, message="Error: You have to use -hp basic with -q\n")
if bool(vars(ARGS)["output_blastp"]) and not bool(vars(ARGS)["fasta"]):
    PARSER.exit(status=0, message="Error: You have to use -op with -q\n")
if bool(vars(ARGS)["strains_order"]) and not bool(vars(ARGS)["heatmap"]):
    PARSER.exit(status=0, message="Error: You have to use -d with -hp\n")
if bool(vars(ARGS)["rarity"]) and not bool(vars(ARGS)["heatmap"] == 'ortholog'):
    PARSER.exit(status=0, message="Error: You have to use -r with -hp ortholog\n")
if bool(vars(ARGS)["rarity_color"]) and not bool(vars(ARGS)["rarity"]):
    PARSER.exit(status=0, message="Error: You have to use -rc with -r\n")
if bool(vars(ARGS)["figure_size"]) and not bool(vars(ARGS)["heatmap"]):
    PARSER.exit(status=0, message="Error: You have to use -fs with -hp\n")
if bool(vars(ARGS)["heatmap_color"]) and not bool(vars(ARGS)["heatmap"]):
    PARSER.exit(status=0, message="Error: You have to use -hc with -hp\n")
if bool(vars(ARGS)["masked_color"]) and not bool(vars(ARGS)["heatmap"]):
    PARSER.exit(status=0, message="Error: You have to use -mc with -hp\n")
if bool(vars(ARGS)["font_size"]) and not bool(vars(ARGS)["heatmap"]):
    PARSER.exit(status=0, message="Error: You have to use -f with -hp\n")
if bool(vars(ARGS)["title"]) and not bool(vars(ARGS)["heatmap"]):
    PARSER.exit(status=0, message="Error: You have to use -t with -hp\n")
if bool(vars(ARGS)["metadata_barplot"]) and not bool(vars(ARGS)["fasta"]):
    PARSER.exit(status=0, message="Error: You have to use -mb with -q\n")
if bool(vars(ARGS)["all_metadata"]) and not bool(vars(ARGS)["metadata_barplot"]):
    PARSER.exit(status=0, message="Error: You have to use -w with -mb\n")
if bool(vars(ARGS)["select_metadata"]) and not bool(vars(ARGS)["metadata_barplot"]):
    PARSER.exit(status=0, message="Error: You have to use -s with -mb\n")
if bool(vars(ARGS)["histogram_color"]) and not bool(vars(ARGS)["histogram"]):
    PARSER.exit(status=0, message="Error: You have to use -e with -x\n")
if bool(vars(ARGS)["histogram_bins"]) and not bool(vars(ARGS)["histogram"]):
    PARSER.exit(status=0, message="Error: You have to use -b with -x\n")
if bool(vars(ARGS)["novel_conserved"]) and not bool(vars(ARGS)["histogram"]):
    PARSER.exit(status=0, message="Error: You have to use -p with -x\n")
if bool(vars(ARGS)["cutoff_volcano"]) and not bool(vars(ARGS)["strains_tag_volcano"]):
    PARSER.exit(status=0, message="Error: You have to use -c with -st\n")
if bool(vars(ARGS)["case_control_name"]) and not bool(vars(ARGS)["strains_tag_volcano"]):
    PARSER.exit(status=0, message="Error: You have to use -cc with -st\n")
OS_SEPARATOR = os.sep
################################
if ARGS.fasta:
    try:
        GETVERSION = subprocess.Popen("blastp -version", shell=True, stdout=subprocess.PIPE).stdout
        VERSION = GETVERSION.read()
        print("Found blastp (version:{})".format(VERSION.decode().splitlines()[1]))
    except:
        PARSER.exit(status=0, message="Error: blastp cannot be found\n")
################################
if ARGS.list_genes:
    GENES_LIST = []
    GENES_FILE_OBJECT = open(ARGS.list_genes, 'r')
    for line in GENES_FILE_OBJECT:
        line = line.rstrip()
        GENES_LIST.append(line)
###################################
if ARGS.novel_conserved:
    NOVEL_LIMIT = ARGS.novel_conserved[0]
    CONSERVED_LIMIT = ARGS.novel_conserved[1]
else:
    NOVEL_LIMIT = 10
    CONSERVED_LIMIT = 100
#########process files in a directory to a list##################
DIRECTORY_PATH = ARGS.directory_path
FILES_LIST = []
STRAINS_LIST = []
for file in os.listdir(DIRECTORY_PATH):
    if file.endswith("WhatsGNU_report.txt"):
        FILES_LIST.append(DIRECTORY_PATH + file)
####################OUTPUT_FOLDER_name##################
FOLDER_PATH = os.path.abspath(DIRECTORY_PATH).rsplit(OS_SEPARATOR, 1)[0] + OS_SEPARATOR
OUTPUT_FOLDER = FOLDER_PATH + ARGS.prefix_name + OS_SEPARATOR
try:
    os.mkdir(OUTPUT_FOLDER)
except:
    PARSER.exit(status=0, message="Error: Folder exists, choose another prefix_name\n")
##########################HEATMAP#################################
##########strains_order##########
if ARGS.heatmap:
    if ARGS.strains_order:
        STRAINS_LIST = []
        STRAINS_FILE_OBJECT = open(ARGS.strains_order, 'r')
        for line in STRAINS_FILE_OBJECT:
            line = line.rstrip()
            STRAINS_LIST.append(line)
    else:
        for file_name in FILES_LIST:
            strain_name = (file_name.rsplit(OS_SEPARATOR, 1)[-1]).split("_WhatsGNU_report.txt")[0]
            STRAINS_LIST.append(strain_name)
    ###########fasta#############
    if ARGS.fasta:
        GENES_DICT = {}
        HISTOGRAM_DICT = {}
        NOVEL_DICT = {}
        CONSERVED_DICT = {}
        GENES_LIST = []
        FASTA_FILE = open(ARGS.fasta, 'r')
        for line in FASTA_FILE:
            line = line.rstrip()
            if line.startswith(">"):
                GENES_LIST.append(line.lstrip('>'))
                GENES_DICT[line.lstrip('>')] = {}
        FASTA_FILE.close()
        for file_name in FILES_LIST:
            novel = []
            conserved = []
            strain_name = (file_name.rsplit(OS_SEPARATOR, 1)[-1]).split("_WhatsGNU_report.txt")[0]
            file_object = open(file_name, 'r')
            if ARGS.heatmap == 'ortholog':
                if ((file_object.readline().rstrip()).split('\t'))[5] != 'ortholog_group':
                    PARSER.exit(status=0, message="Error: This report '{}' was created in WhatsGNU basic mode\n".format(strain_name))
            else:
                file_object.readline()
            if ARGS.output_blastp:
                blast_report = (OUTPUT_FOLDER +
                    (file_name.rsplit(OS_SEPARATOR, 1)[-1]).split("_WhatsGNU_report.txt")[0]
                    + "_WhatsGNU_plotter_blast_report.txt"
                )
            else:
                blast_report_file = tempfile.NamedTemporaryFile(mode='w+')
            OUTPUT_FILE_SEQ = tempfile.NamedTemporaryFile(mode='w+')
            GNU_scores = []
            for line in file_object:
                line = line.rstrip()
                GNU_report_list = line.split('\t')
                if ARGS.histogram:
                    GNU_scores.append(int(GNU_report_list[1]))
                    if int(GNU_report_list[1]) < NOVEL_LIMIT:
                        novel.append(GNU_report_list[0])
                    elif int(GNU_report_list[1]) > CONSERVED_LIMIT:
                        conserved.append(GNU_report_list[0])
                if ARGS.heatmap == 'ortholog':
                    record_id = (GNU_report_list[0].split(' ')[0] + '__' + GNU_report_list[1]
                            + '__' + GNU_report_list[5] + '__' + GNU_report_list[-2]
                            + '__' + GNU_report_list[-1])
                else:
                    record_id = (GNU_report_list[0].split(' ')[0] + '__' + GNU_report_list[1])
                record_seq = GNU_report_list[4]
                OUTPUT_FILE_SEQ.write('>'+record_id+'\n'+record_seq+'\n')
            if ARGS.histogram:
                HISTOGRAM_DICT[strain_name] = GNU_scores
                NOVEL_DICT[strain_name] = novel
                CONSERVED_DICT[strain_name] = conserved
            OUTPUT_FILE_SEQ.seek(0)
            if ARGS.output_blastp:
                blast_results = os.system("blastp -query {} -subject {} -evalue 1e-06 -max_target_seqs 5 -outfmt '6 qseqid sacc evalue qcovs pident' -out {}".format(ARGS.fasta, OUTPUT_FILE_SEQ.name, blast_report))
                blast_report_file = open(blast_report, "r")
            else:
                blast_results = os.system("blastp -query {} -subject {} -evalue 1e-06 -max_target_seqs 5 -outfmt '6 qseqid sacc evalue qcovs pident' -out {}".format(ARGS.fasta, OUTPUT_FILE_SEQ.name, blast_report_file.name))
                blast_report_file.seek(0)
            for line in blast_report_file:
                line = line.rstrip()
                blast_report_list = line.split('\t')
                if float(blast_report_list[-1]) > 80.0:
                    subject_record_list = blast_report_list[1].split('__')
                    try:#capturing paralogs or broken proteins with highest GNU score
                        if int(subject_record_list[1]) > int(GENES_DICT[blast_report_list[0]][strain_name][1]):
                            GENES_DICT[blast_report_list[0]][strain_name] = subject_record_list
                    except:
                        GENES_DICT[blast_report_list[0]][strain_name] = subject_record_list
        OUTPUT_FILE_NAME = OUTPUT_FOLDER + ARGS.prefix_name + '_WhatsGNU_heatmap' + '.txt'
        OUTPUT_FILE_OBJECT = open(OUTPUT_FILE_NAME, 'w')
        OUTPUT_FILE_OBJECT.write("\t")
        OUTPUT_FILE_OBJECT.write('\t'.join(STRAINS_LIST))
        FINAL_HEATMAP_LIST = []
        FINAL_HEATMAP_SIGNIFICANCE = []
        for i in GENES_LIST:
            heatmap_GNU_list = []
            heatmap_Significance = []
            OUTPUT_FILE_OBJECT.write('\n')
            OUTPUT_FILE_OBJECT.write(i)
            for strain_name in STRAINS_LIST:
                OUTPUT_FILE_OBJECT.write("\t")
                try:
                    GNU_score = GENES_DICT[i][strain_name][1]
                    OUTPUT_FILE_OBJECT.write(GNU_score)
                    heatmap_GNU_list.append(int(GNU_score))
                    if ARGS.heatmap == 'ortholog':
                        Significance = GENES_DICT[i][strain_name][-1][0]
                        if Significance == 'r':
                            heatmap_Significance.append(Significance)
                        else:
                            heatmap_Significance.append('')
                except:
                    OUTPUT_FILE_OBJECT.write('NA')
                    heatmap_GNU_list.append(np.nan)
                    if ARGS.heatmap == 'ortholog':
                        heatmap_Significance.append('')
            FINAL_HEATMAP_LIST.append(heatmap_GNU_list)
            if ARGS.heatmap == 'ortholog':
                FINAL_HEATMAP_SIGNIFICANCE.append(heatmap_Significance)
    ###########GENES_LIST#####################
    if ARGS.list_genes:
        HISTOGRAM_DICT = {}
        for file_name in FILES_LIST:
            strain_name = (file_name.rsplit(OS_SEPARATOR, 1)[-1]).split("_WhatsGNU_report.txt")[0]
            file_object = open(file_name, 'r')
            if ((file_object.readline().rstrip()).split('\t'))[5] != 'ortholog_group':
                PARSER.exit(status=0, message="Error: This report '{}' was created in WhatsGNU basic mode\nUse -l for reports that were created with the ortholog mode\n".format(strain_name))
            GNU_scores = []
            for line in file_object:
                line = line.rstrip()
                GNU_report_list = line.split('\t')
                if ARGS.histogram:
                    GNU_scores.append(int(GNU_report_list[1]))
                    if int(GNU_report_list[1]) < NOVEL_LIMIT:
                        novel.append(GNU_report_list[0])
                    elif int(GNU_report_list[1]) > CONSERVED_LIMIT:
                        conserved.append(GNU_report_list[0])
            if ARGS.histogram:
                HISTOGRAM_DICT[strain_name] = GNU_scores
                NOVEL_DICT[strain_name] = novel
                CONSERVED_DICT[strain_name] = conserved
            file_object.close()
        #################Extract genes from GNU report#######################
        GENES_DICT = {}
        for i in GENES_LIST:
            GENES_DICT[i] = {}
            for file_name in FILES_LIST:
                strain_name = (file_name.rsplit(OS_SEPARATOR, 1)[-1]).split("_WhatsGNU_report.txt")[0]
                file_object = open(file_name, 'r')
                file_object.readline()
                for line in file_object:
                    line = line.rstrip()
                    GNU_report_list = line.split('\t')
                    GNU_significance_list = []
                    if i == GNU_report_list[5]:
                        GNU_significance_list.append(GNU_report_list[1])
                        GNU_significance_list.append(GNU_report_list[-1])
                        GENES_DICT[i][strain_name] = GNU_significance_list # dictionary of dictionary
                        break
                file_object.close()
        ###############heatmap_matplotlib######################
        OUTPUT_FILE_NAME2 = OUTPUT_FOLDER + ARGS.prefix_name + '_WhatsGNU_heatmap' + '.txt'
        OUTPUT_FILE_OBJECT2 = open(OUTPUT_FILE_NAME2, 'w')
        OUTPUT_FILE_OBJECT2.write("\t")
        OUTPUT_FILE_OBJECT2.write('\t'.join(STRAINS_LIST))
        FINAL_HEATMAP_LIST = []
        FINAL_HEATMAP_SIGNIFICANCE = []
        for i in GENES_LIST:
            heatmap_GNU_list = []
            heatmap_Significance = []
            OUTPUT_FILE_OBJECT2.write('\n')
            OUTPUT_FILE_OBJECT2.write(i)
            for strain_name in STRAINS_LIST:
                OUTPUT_FILE_OBJECT2.write("\t")
                try:
                    GNU_score = GENES_DICT[i][strain_name][0]
                    Significance = GENES_DICT[i][strain_name][-1][0]
                    OUTPUT_FILE_OBJECT2.write(GNU_score)
                    heatmap_GNU_list.append(int(GNU_score))
                    if Significance == 'r':
                        heatmap_Significance.append(Significance)
                    else:
                        heatmap_Significance.append('')
                except:
                    OUTPUT_FILE_OBJECT2.write('NA')
                    heatmap_GNU_list.append(np.nan)
                    heatmap_Significance.append('')
            FINAL_HEATMAP_LIST.append(heatmap_GNU_list)
            FINAL_HEATMAP_SIGNIFICANCE.append(heatmap_Significance)
    ###########convert lists to array###############
    FINAL_HEATMAP_ARRAY = np.array(FINAL_HEATMAP_LIST)
    FINAL_HEATMAP_ARRAY = np.ma.masked_invalid(FINAL_HEATMAP_ARRAY)
    if ARGS.heatmap == 'ortholog':
        FINAL_HEATMAP_SIGNIFICANCE_ARRAY = np.array(FINAL_HEATMAP_SIGNIFICANCE)
    ####################plot heatmap########################
    fig, ax = plt.subplots()
    #heatmap color
    if ARGS.heatmap_color:
        HEATMAP_COLOR = ARGS.heatmap_color
    else:
        HEATMAP_COLOR = 'Reds'
    cmap = plt.get_cmap(HEATMAP_COLOR)
    im = ax.imshow(FINAL_HEATMAP_ARRAY, cmap)
    if ARGS.masked_color:
        MASK_COLOR = ARGS.masked_color
    else:
        MASK_COLOR = '#D3D3D3'
    cmap.set_bad(color=MASK_COLOR)
    #font size
    if ARGS.font_size:
        FONT_SIZE = ARGS.font_size
    else:
        FONT_SIZE = 8
    # Show all ticks...
    ax.set_xticks(np.arange(len(STRAINS_LIST)))
    ax.set_yticks(np.arange(len(GENES_LIST)))
    # Label ticks with the respective list entries
    ax.set_xticklabels(STRAINS_LIST, fontsize=FONT_SIZE)
    ax.set_yticklabels(GENES_LIST, fontsize=FONT_SIZE)
    plt.ylabel('Proteins', fontsize=FONT_SIZE)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    # Loop over data dimensions and create text annotations.
    if ARGS.heatmap == 'ortholog':
        if ARGS.rarity:
            if ARGS.rarity_color:
                RARITY_TEXT_COLOR = ARGS.rarity_color
            else:
                RARITY_TEXT_COLOR = "k"
            for i in range(len(GENES_LIST)):
                for j in range(len(STRAINS_LIST)):
                    text = ax.text(j, i, FINAL_HEATMAP_SIGNIFICANCE_ARRAY[i, j],
                                   ha="center", va="center", color=RARITY_TEXT_COLOR, fontsize=FONT_SIZE)
    #Heatmap Title
    if ARGS.title:
        ax.set_title(ARGS.title)
    else:
        ax.set_title("WhatsGNU heatmap")
    #colorbar legend
    cbar = ax.figure.colorbar(im, ax=ax, pad=0.02)
    cbar.ax.set_ylabel("GNU score", rotation=-90, va="bottom", fontsize=FONT_SIZE)
    fig.tight_layout()
    #figure size
    if ARGS.figure_size:
        FIGURE_SIZE = tuple(ARGS.figure_size)
    else:
        FIGURE_SIZE = tuple([10, 5])
    fig.set_size_inches(FIGURE_SIZE, forward=True) # width and height
    #adjust for lengthy gene list
    if len(GENES_LIST) > 10:
        plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
    # save heatmap file
    HEATMAP_FILE_NAME = OUTPUT_FOLDER + ARGS.prefix_name + '_WhatsGNU_heatmap' + '.pdf'
    plt.savefig(HEATMAP_FILE_NAME, dpi=300, bbox_inches='tight')
##########################metadata##############################
if ARGS.metadata_barplot:
    GNU_REPORT_DICT = {}
    if ARGS.fasta:#parse the fasta sequence file and get the ptn names
        GENES_DICT = {}
        HISTOGRAM_DICT = {}
        GENES_LIST = []
        FASTA_FILE = open(ARGS.fasta, 'r')
        for line in FASTA_FILE:
            line = line.rstrip()
            if line.startswith(">"):
                GENES_LIST.append(line.lstrip('>'))
                GENES_DICT[line.lstrip('>')] = {}
        FASTA_FILE.close()
        if ARGS.select_metadata:#get user selected metadata and capture length
            HEADER_LINE = []
            SELECT_FILE = open(ARGS.select_metadata, 'r')
            for line in SELECT_FILE:
                line = line.rstrip()
                HEADER_LINE.append(line)
            METADATA_LEN = len(HEADER_LINE)#length of metadata
            SELECT_FILE.close()
        for file_name in FILES_LIST:
            metadata_blast_dict = {}
            strain_name = (file_name.rsplit(OS_SEPARATOR, 1)[-1]).split("_WhatsGNU_report.txt")[0]
            file_object = open(file_name, 'r')
            if ARGS.select_metadata:#get the index of each metadata
                metadata_index = []
                for i in HEADER_LINE:
                    file_object.seek(0)
                    metadata_index.append(((file_object.readline().rstrip()).split('\t')).index(i))
            else:
                if ARGS.metadata_barplot == 'ortholog':#length of metadata
                    METADATA_LEN = len((file_object.readline().rstrip()).split('\t')[13:])
                else:
                    METADATA_LEN = len((file_object.readline().rstrip()).split('\t')[5:])
            file_object.seek(0)
            if ARGS.all_metadata:
                if ARGS.metadata_barplot == 'ortholog':#capture metadata
                    HEADER_LINE = (file_object.readline().rstrip()).split('\t')[13:]
                else:
                    HEADER_LINE = (file_object.readline().rstrip()).split('\t')[5:]
            elif ARGS.select_metadata:
                HEADER_LINE = []
                SELECT_FILE = open(ARGS.select_metadata, 'r')
                for line in SELECT_FILE:
                    line = line.rstrip()
                    HEADER_LINE.append(line)
                METADATA_LEN = len(HEADER_LINE)
                SELECT_FILE.close()
            else:
                if ARGS.metadata_barplot == 'ortholog':
                    HEADER_LINE = (file_object.readline().rstrip()).split('\t')[13:]
                else:
                    HEADER_LINE = (file_object.readline().rstrip()).split('\t')[5:]
                if METADATA_LEN > 10:#keep 9 only
                    HEADER_LINE = HEADER_LINE[:9]
                    HEADER_LINE.append('Others')
            if ARGS.output_blastp:
                blast_report = (OUTPUT_FOLDER +
                    (file_name.rsplit(OS_SEPARATOR, 1)[-1]).split("_WhatsGNU_report.txt")[0]
                    + "_WhatsGNU_barplot_blast_report.txt"
                )
            else:
                blast_report_file = tempfile.NamedTemporaryFile(mode='w+')
            OUTPUT_FILE_SEQ = tempfile.NamedTemporaryFile(mode='w+')
            for line in file_object:
                line = line.rstrip()
                GNU_report_list = line.split('\t')
                record_id = GNU_report_list[0].split(' ')[0]
                record_seq = GNU_report_list[4]
                if ARGS.metadata_barplot == 'ortholog':
                    GNU_REPORT_DICT[record_id] = '__'.join(GNU_report_list[13:])
                else:
                    GNU_REPORT_DICT[record_id] = '__'.join(GNU_report_list[5:])
                OUTPUT_FILE_SEQ.write('>'+record_id+'\n'+record_seq+'\n')
            OUTPUT_FILE_SEQ.seek(0)
            if ARGS.output_blastp:
                blast_results = os.system("blastp -query {} -subject {} -evalue 1e-06 -max_target_seqs 5 -outfmt '6 qseqid sacc evalue qcovs pident' -out {}".format(ARGS.fasta, OUTPUT_FILE_SEQ.name, blast_report))
                blast_report_file = open(blast_report, "r")
            else:
                blast_results = os.system("blastp -query {} -subject {} -evalue 1e-06 -max_target_seqs 5 -outfmt '6 qseqid sacc evalue qcovs pident' -out {}".format(ARGS.fasta, OUTPUT_FILE_SEQ.name, blast_report_file.name))
                blast_report_file.seek(0)
            for line in blast_report_file:
                line = line.rstrip()
                blast_report_list = line.split('\t')
                if  blast_report_list[0] not in metadata_blast_dict:
                    subject_metadata = []
                    for i in (GNU_REPORT_DICT[blast_report_list[1]]).split('__'):
                        subject_metadata.append(float(i))
                    metadata_blast_dict[blast_report_list[0]] = subject_metadata
            for record in metadata_blast_dict:
                if ARGS.all_metadata:
                    CC_included = metadata_blast_dict[record]
                    if len(CC_included) > 10:
                        FONT_SIZE2 = 8
                    else:
                        FONT_SIZE2 = 10
                elif ARGS.select_metadata:
                    CC_included = []
                    CC_all = metadata_blast_dict[record]
                    if ARGS.metadata_barplot == 'ortholog':
                        for i in metadata_index:
                            CC_included.append(CC_all[(i - 13)])
                    else:
                        for i in metadata_index:
                            CC_included.append(CC_all[(i - 5)])
                    if len(CC_included) > 10:
                        FONT_SIZE2 = 8
                    else:
                        FONT_SIZE2 = 10
                else:
                    if METADATA_LEN > 10:
                        CC_percent = metadata_blast_dict[record]
                        other_metadata = CC_percent[9:]
                        total_others = 0
                        for i in other_metadata:
                            total_others += i
                        if total_others > 100:
                            total_others = 100
                        CC_included = CC_percent[0:9]
                        CC_included.append(total_others)
                        FONT_SIZE2 = 8
                    else:
                        CC_included = metadata_blast_dict[record]
                        FONT_SIZE2 = 10
                for metadatum_name, metadatum_freq in zip(HEADER_LINE, CC_included):
                    plt.bar(metadatum_name, metadatum_freq, width=0.8, align='center')
                plt.xlabel('Metadata', fontsize=12, fontname="sans-serif", fontweight="bold")
                plt.ylabel('Percentage', fontsize=12, fontname="sans-serif", fontweight="bold")
                plt.xticks(HEADER_LINE, fontsize=FONT_SIZE2, rotation=90, fontname="sans-serif")
                plt.yticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100], fontsize=FONT_SIZE2, fontname="sans-serif")
                plt.title(record, fontname="sans-serif", fontweight="bold")
                figure_name = OUTPUT_FOLDER + strain_name + 'WhatsGNU_barplot_' +record+ '.png'
                plt.savefig(figure_name, format='png', bbox_inches='tight', dpi=300)
                plt.close()
###############histogram######################
if ARGS.histogram:
    if not ARGS.heatmap:
        HISTOGRAM_DICT = {}
        NOVEL_DICT = {}
        CONSERVED_DICT = {}
        for file_name in FILES_LIST:
            strain_name = (file_name.rsplit(OS_SEPARATOR, 1)[-1]).split("_WhatsGNU_report.txt")[0]
            STRAINS_LIST.append(strain_name)
            file_object = open(file_name, 'r')
            file_object.readline()
            GNU_scores = []
            novel = []
            conserved = []
            for line in file_object:
                line = line.rstrip()
                GNU_report_list = line.split('\t')
                GNU_scores.append(int(GNU_report_list[1]))
                if int(GNU_report_list[1]) < NOVEL_LIMIT:
                    novel.append(GNU_report_list[0])
                elif int(GNU_report_list[1]) > CONSERVED_LIMIT:
                    conserved.append(GNU_report_list[0])
            HISTOGRAM_DICT[strain_name] = GNU_scores
            NOVEL_DICT[strain_name] = novel
            CONSERVED_DICT[strain_name] = conserved
            file_object.close()
    if ARGS.histogram_color:
        HISTO_COLOR = ARGS.histogram_color
    else:
        HISTO_COLOR = 'b'
    if ARGS.histogram_bins:
        HISTO_BINS = ARGS.histogram_bins
    else:
        HISTO_BINS = 10
    for strain_name in STRAINS_LIST:
        GNU_scores_list = HISTOGRAM_DICT[strain_name]
        NOVEL_GENES_LIST = NOVEL_DICT[strain_name]
        CONSERVED_GENES_LIST = CONSERVED_DICT[strain_name]
        OUTPUT_FILE_NAME3 = OUTPUT_FOLDER + strain_name + '_WhatsGNU_histogram.txt'
        OUTPUT_FILE_OBJECT3 = open(OUTPUT_FILE_NAME3, 'w')
        OUTPUT_FILE_OBJECT3.write('GNU<{}'.format(NOVEL_LIMIT))
        for i in NOVEL_GENES_LIST:
            OUTPUT_FILE_OBJECT3.write('\t')
            OUTPUT_FILE_OBJECT3.write(i)
        OUTPUT_FILE_OBJECT3.write('\nGNU>{}'.format(CONSERVED_LIMIT))
        for i in CONSERVED_GENES_LIST:
            OUTPUT_FILE_OBJECT3.write('\t')
            OUTPUT_FILE_OBJECT3.write(i)
        OUTPUT_FILE_OBJECT3.close()
        fig2, ax2 = plt.subplots()
        y, x, _ = plt.hist(GNU_scores_list, facecolor=HISTO_COLOR, bins=HISTO_BINS)#it accepts hex
        max_GNU = x.max()
        LINE_LOCATION = -0.12*y.max()
        WORD_LOCATION = (-0.12*y.max()) - (0.5*(0.12*y.max()))
        ax2.annotate('', xy=(0, LINE_LOCATION), xytext=(max_GNU, LINE_LOCATION), arrowprops=dict(arrowstyle='<->'), annotation_clip=False)
        ax2.annotate('Novelty', xy=(0, WORD_LOCATION), xytext=(0, WORD_LOCATION), annotation_clip=False)
        ax2.annotate('Ultraconserved', xy=(max_GNU, WORD_LOCATION), xytext=(max_GNU, WORD_LOCATION), horizontalalignment='right', annotation_clip=False)
        plt.margins(x=0)
        plt.xlabel('GNU score')
        plt.ylabel('Number of Proteins')
        plt.title('Histogram of GNU scores of {}'.format(strain_name))
        histogram_file_name = OUTPUT_FOLDER + strain_name +'_WhatsGNU_histogram.pdf'
        plt.savefig(histogram_file_name, dpi=300, bbox_inches='tight')
####################################Volcano###########################################
if ARGS.strains_tag_volcano:
    DIRECTORY_PATH = ARGS.directory_path
    CASE_COUNTER = 0
    CONTROL_COUNTER = 0
    if ARGS.strains_tag_volcano:
        CD_DICT = {}
        FILES_LIST = []
        CD_FILE_OBJECT = open(ARGS.strains_tag_volcano, 'r')
        for line in CD_FILE_OBJECT:
            line = line.rstrip()
            strain, CD_tag = line.split(',')
            if CD_tag == 'case':
                CASE_COUNTER += 1
            elif CD_tag == 'control':
                CONTROL_COUNTER += 1
            else:
                PARSER.exit(status=0, message="You are using tag other than case and control (case-sensitive)\n")
            FILES_LIST.append(DIRECTORY_PATH + strain)
            CD_DICT[strain] = CD_tag

    ###########################
    UNTARGETED = defaultdict(list)
    ORTHO_STRAIN_FUNCTION = {}
    for file_name in FILES_LIST:
        logp_list = []
        strain_name = file_name.rsplit(OS_SEPARATOR, 1)[-1]
        file_object = open(file_name, 'r')
        if ((file_object.readline().rstrip()).split('\t'))[5] != 'ortholog_group':
            PARSER.exit(status=0, message="Error: This report '{}' was created in WhatsGNU basic mode\nVolcano plot works for reports that were created with the ortholog mode\n".format(strain_name))
        ortho_strain_GNU = defaultdict(list)
        ortho_strain_RI = defaultdict(list)
        CD_tag = CD_DICT[strain_name]
        for line in file_object:
            line = line.rstrip()
            GNU_report_list = line.split('\t')
            if GNU_report_list[5] != 'NA':#to exclude unannotated GNU=0
                ortho_strain_GNU[GNU_report_list[5]].append(int(GNU_report_list[1])) #GNU score
                ortho_strain_RI[GNU_report_list[5]].append(float(GNU_report_list[-2])) #variant frequency
                ORTHO_STRAIN_FUNCTION[GNU_report_list[5]] = GNU_report_list[3] #function
        for record in ortho_strain_GNU:
            record_details = (CD_tag + '__' + str(max(ortho_strain_GNU[record])) + '__' + str(max(ortho_strain_RI[record])))
            UNTARGETED[record].append(record_details)
        file_object.close()
    #################################################
    OUTPUT_FILE_NAME4 = OUTPUT_FOLDER + ARGS.prefix_name + '_WhatsGNU_volcano' + '.txt'
    OUTPUT_FILE_OBJECT4 = open(OUTPUT_FILE_NAME4, 'w')
    OUTPUT_FILE_NAME5 = OUTPUT_FOLDER + ARGS.prefix_name + '_ptns_present_in_one_gp' + '.txt'
    OUTPUT_FILE_OBJECT5 = open(OUTPUT_FILE_NAME5, 'w')
    if ARGS.cutoff_volcano:
        CUTOFF_VALUE = ARGS.cutoff_volcano
        CASE_CUTOFF = int((CUTOFF_VALUE * CASE_COUNTER)/100)
        CONTROL_CUTOFF = int((CUTOFF_VALUE * CONTROL_COUNTER)/100)
    else:
        CUTOFF_VALUE = 100
        CASE_CUTOFF = int((CUTOFF_VALUE * CASE_COUNTER)/100)
        CONTROL_CUTOFF = int((CUTOFF_VALUE * CONTROL_COUNTER)/100)
    if ARGS.case_control_name:
        CASE_NAME = 'Lower GNU in {}'.format(ARGS.case_control_name[0])
        CONTROL_NAME = 'Lower GNU in {}'.format(ARGS.case_control_name[1])
        OUTPUT_FILE_OBJECT4.write('ortho_group\tfunction\tDELTA_AVG_GNU\tAVG_OVRI\tMann-Whitney U statistic\tMann–Whitney–Wilcoxon pvalue\t-log10 pvalue\t{}\t{}\n'.format(ARGS.case_control_name[0], ARGS.case_control_name[1]))
        OUTPUT_FILE_OBJECT5.write('ortho_group\tfunction\tAVG_GNU\tAVG_OVRI\t{}\t{}\n'.format(ARGS.case_control_name[0], ARGS.case_control_name[1]))
    else:
        CASE_NAME = 'Lower GNU in case'
        CONTROL_NAME = 'Lower GNU in control'
        OUTPUT_FILE_OBJECT4.write('ortho_group\tfunction\tDELTA_AVG_GNU\tAVG_OVRI\tMann-Whitney U statistic\tMann–Whitney–Wilcoxon pvalue\t-log10 pvalue\tcase\tcontrol\n')
        OUTPUT_FILE_OBJECT5.write('ortho_group\tfunction\tDELTA_AVG_GNU\tAVG_OVRI\tcase\tcontrol\n')
    VOLCANO = {}
    DELTA_AVG_GNU_LIST = []
    for record in UNTARGETED:
        all_list = UNTARGETED[record]
        function = ORTHO_STRAIN_FUNCTION[record]
        CD_list = []
        OVRI_list = []
        OVRI_case = []
        OVRI_control = []
        GNU_case = []
        GNU_control = []
        for i in all_list:
            CD, GNU, OVRI = i.split('__')
            CD_list.append(CD)
            if CD == 'case':
                GNU_case.append(int(GNU))
                OVRI_case.append(float(OVRI))
            elif CD == 'control':
                GNU_control.append(int(GNU))
                OVRI_control.append(float(OVRI))
            else:
                PARSER.exit(status=0, message="You are using tag other than case and control (case-sensitive)\n")
        if CD_list.count('case') >= CASE_CUTOFF or CD_list.count('control') >= CONTROL_CUTOFF:
            if len(GNU_case) == 0:
                DELTA_AVG_GNU = (sum(GNU_control)/len(GNU_control))
                AVG_RI = sum(OVRI_control)/len(OVRI_control)
                OUTPUT_FILE_OBJECT5.write(record+'\t'+function+'\t'+str(DELTA_AVG_GNU)+'\t'+str(AVG_RI)+'\t'+str(CD_list.count('case'))+'\t'+str(CD_list.count('control'))+'\n')
            elif len(GNU_control) == 0:
                DELTA_AVG_GNU = (sum(GNU_case)/len(GNU_case))
                AVG_RI = sum(OVRI_case)/len(OVRI_case)
                OUTPUT_FILE_OBJECT5.write(record+'\t'+function+'\t'+str(DELTA_AVG_GNU)+'\t'+str(AVG_RI)+'\t'+str(CD_list.count('case'))+'\t'+str(CD_list.count('control'))+'\n')
            else:
                DELTA_AVG_GNU = (sum(GNU_case)/len(GNU_case)) - (sum(GNU_control)/len(GNU_control))
                try:
                    wtstatistics, pvalue = st.mannwhitneyu(GNU_case, GNU_control, alternative='two-sided')
                    logp = (-1) * (np.log10(pvalue))
                except:
                    pvalue = 1.0
                    logp = 0.0
                    wtstatistics = CASE_COUNTER * CONTROL_COUNTER
                if DELTA_AVG_GNU < 0:
                    AVG_RI = sum(OVRI_case)/len(OVRI_case)
                elif DELTA_AVG_GNU > 0:
                    AVG_RI = sum(OVRI_control)/len(OVRI_control)
                else:
                    AVG_RI = (sum(OVRI_case)+sum(OVRI_control))/(len(OVRI_control)+len(OVRI_case))
                DELTA_AVG_GNU_LIST.append(DELTA_AVG_GNU)
                logp_list.append(logp)
                VOLCANO[record] = (DELTA_AVG_GNU, AVG_RI, pvalue, logp)
                OUTPUT_FILE_OBJECT4.write(record+'\t'+function+'\t'+str(DELTA_AVG_GNU)+'\t'+str(AVG_RI)+'\t'+str(wtstatistics)+'\t'+str(pvalue)+'\t'+str(logp)+'\t'+str(CD_list.count('case'))+'\t'+str(CD_list.count('control'))+'\n')
    BORDER_LINE = (max(DELTA_AVG_GNU_LIST)-min(DELTA_AVG_GNU_LIST))/4
    OUTPUT_FILE_OBJECT4.close()
    OUTPUT_FILE_OBJECT5.close()
    #############DELTA_AVG_GNU and AVG_OVRI##############
    VOLCANO_OVRI_FILE_NAME = OUTPUT_FOLDER + ARGS.prefix_name + '_WhatsGNU_volcano_OVRI' + '.png'
    fig, ax = plt.subplots()
    LINE_LOCATION = -0.22
    WORD_LOCATION = -0.27
    for record in VOLCANO:
        DELTA_AVG_GNU, AVG_RI, pvalue, logp = VOLCANO[record]
        if DELTA_AVG_GNU > BORDER_LINE: #capture dots above -log10 pvalue and borderline
            dot_color = 'green'
        elif DELTA_AVG_GNU < (-1*BORDER_LINE):
            dot_color = 'red'
        else:
            dot_color = 'grey'
        plt.scatter(DELTA_AVG_GNU, AVG_RI, c=dot_color)
    plt.xlabel('DELTA_AVG_GNU', fontsize=12, fontname="sans-serif", fontweight="bold")
    plt.ylabel('AVG_OVRI', fontsize=12, fontname="sans-serif", fontweight="bold")
    plt.xticks(fontsize=12, fontname="sans-serif")
    plt.yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], fontsize=12, fontname="sans-serif")
    ax.annotate('', xy=(min(DELTA_AVG_GNU_LIST), LINE_LOCATION), xytext=(max(DELTA_AVG_GNU_LIST), LINE_LOCATION), arrowprops=dict(arrowstyle='<->'), annotation_clip=False)
    ax.annotate(CASE_NAME, xy=(min(DELTA_AVG_GNU_LIST), WORD_LOCATION), xytext=(min(DELTA_AVG_GNU_LIST), WORD_LOCATION), annotation_clip=False)
    ax.annotate(CONTROL_NAME, xy=(max(DELTA_AVG_GNU_LIST), WORD_LOCATION), xytext=(max(DELTA_AVG_GNU_LIST), WORD_LOCATION), horizontalalignment='right', annotation_clip=False)
    plt.savefig(VOLCANO_OVRI_FILE_NAME, format='png', bbox_inches='tight', dpi=300)
    plt.close()
    ###############Mann–Whitney–Wilcoxon#################
    VOLCANO_WMU_FILE_NAME = OUTPUT_FOLDER + ARGS.prefix_name + '_WhatsGNU_volcano_WMu' + '.png'
    YTICKS_LIST = []
    for i in range(cl(max(logp_list))+1):
        YTICKS_LIST.append(float(i))
    fig, ax = plt.subplots()
    LINE_LOCATION = -((max(YTICKS_LIST))/4)
    WORD_LOCATION = LINE_LOCATION - 0.3
    for record in VOLCANO:
        DELTA_AVG_GNU, AVG_RI, pvalue, logp = VOLCANO[record]
        if logp > 1.3:
            if DELTA_AVG_GNU > BORDER_LINE:
                dot_color = 'green'
            elif DELTA_AVG_GNU < (-1*BORDER_LINE):
                dot_color = 'red'
            else:
                dot_color = 'grey'
        else:
            dot_color = 'grey'
        plt.scatter(DELTA_AVG_GNU, logp, c=dot_color)
    plt.xlabel('DELTA_AVG_GNU', fontsize=12, fontname="sans-serif", fontweight="bold")
    plt.ylabel('-log10 pvalue', fontsize=12, fontname="sans-serif", fontweight="bold")
    plt.xticks(fontsize=12, fontname="sans-serif")
    plt.yticks(np.arange(len(YTICKS_LIST)), YTICKS_LIST, fontsize=12, fontname="sans-serif")
    x2, y2 = [min(DELTA_AVG_GNU_LIST), max(DELTA_AVG_GNU_LIST)], [1.3, 1.3]
    plt.plot(x2, y2, ls='--')
    ax.annotate('', xy=(min(DELTA_AVG_GNU_LIST), LINE_LOCATION), xytext=(max(DELTA_AVG_GNU_LIST), LINE_LOCATION), arrowprops=dict(arrowstyle='<->'), annotation_clip=False)
    ax.annotate(CASE_NAME, xy=(min(DELTA_AVG_GNU_LIST), WORD_LOCATION), xytext=(min(DELTA_AVG_GNU_LIST), WORD_LOCATION), annotation_clip=False)
    ax.annotate(CONTROL_NAME, xy=(max(DELTA_AVG_GNU_LIST), WORD_LOCATION), xytext=(max(DELTA_AVG_GNU_LIST), WORD_LOCATION), horizontalalignment='right', annotation_clip=False)
    plt.savefig(VOLCANO_WMU_FILE_NAME, format='png', bbox_inches='tight', dpi=300)
    plt.close()
