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

# DATE CREATED: January 29, 2022

# AUTHOR: Ahmed M Moustafa

# CONTACT1: moustafaam@email.chop.edu
# CONTACT2: ahmedmagdy2009@hotmail.com

# AFFILIATION: Pediatric Infectious Disease Division, Children’s Hospital of Philadelphia,
# Abramson Pediatric Research Center, University of Pennsylvania, Philadelphia,
# Pennsylvania, 19104, USA

# CITATION1: Ahmed M Moustafa and Paul J Planet
# WhatsGNU: a tool for identifying proteomic novelty
# Genome Biology(2020)21:58, doi: https://doi.org/10.1186/s13059-020-01965-w.

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
    prog="WhatsGNU_allele_scoring.py", description="WhatsGNU_allele_scoring script for WhatsGNU v1.5.",)
PARSER.add_argument("traits", type=str, help="traits file")
PARSER.add_argument("directory_path", type=str, help="path to directory of WhatsGNU reports")
if len(sys.argv) == 1:
    PARSER.print_help()
    sys.exit(0)
ARGS = PARSER.parse_args()
OS_SEPARATOR = os.sep
#########process files in a directory to a list##################
DIRECTORY_PATH = ARGS.directory_path.rstrip(OS_SEPARATOR) + OS_SEPARATOR
FILES_LIST = []
STRAINS_LIST = []
for file in os.listdir(DIRECTORY_PATH):
    if file.endswith("WhatsGNU_report.txt"):
        FILES_LIST.append(DIRECTORY_PATH + file)
############traits file##############

# reads in positive or negative values (1's and 0's for the 'trait' (in this case,
# being in the Chilean clade or not))

traits_obj = open(ARGS.traits,'r')
traits_header = traits_obj.readline().rstrip().split(',')[1:]
trait_dict = {i:{'count':[],'+ves':[],'-ves':[]} for i in traits_header}
strain_count = 0
for line in traits_obj:
    strain_count +=1
    strain_traits = line.rstrip().split(',')
    strain_name = strain_traits[0]
    traits_hits = strain_traits[1:] #[1,0,1,1]
    for t,z in zip(traits_header,traits_hits):
        trait_dict[t][strain_name] = z
        trait_dict[t]['count'].append(z)
        if z == "1":
            trait_dict[t]['+ves'].append(strain_name)
        elif z == "0":
            trait_dict[t]['-ves'].append(strain_name)
        else:
            sys.exit('please use 1 or 0 in traits file')

#count cases and controls for each  trait
#print(trait_dict)
for t in trait_dict:
    t_lst = trait_dict[t]['count']
    c_1 = t_lst.count('1')
    c_0 = t_lst.count('0')
    trait_dict[t]['count'] = [c_1, c_0]
print(trait_dict)
#########alleles dict###############
# alleles_dict seems to be structured such that:
# the 'outermost' key is the AA sequence for an allele (so aaa in the example)
# then the values for those AA sequence keys are their own dictionaries where key is either
# 'genes' or 'isolates'

# 'genes' seems to store gene names which are pulled from the first column and
# seem to be the prokka-assigned CDS name from the database

# 'isolates' seems to store the name of the file something was found in


{'aaa':{'genes':['gene1','gene2'],'isolates':['isolate1','isolate2']},
'bbb':{'genes':['gene1','gene2'],'isolates':['isolate1','isolate2']}}
alleles_dict = {}
"""protein	GNU score	length	function	sequence
ortholog_group	ortho_gp_total_sequences_number	ortho_gp_total_variants_number
minimum_GNU	maximum_GNU	average_GNU	OVRI	OVRI interpretation"""

for f in FILES_LIST:
    strain_name = (f.rsplit(OS_SEPARATOR, 1)[-1]).split("_WhatsGNU_report.txt")[0]
    file_obj = open(f, 'r')
    file_obj.readline()
    for line in file_obj:
        # For a given item in the whatsGNU report, it makes an item like 'HCFJDM_06285 Plasmid recombination enzyme type 2__2745__420__group_15538'
        allele_lst = line.rstrip().split('\t')
        ptn_name = allele_lst[0]
        #gnu score, length, function, ortho group
        ptn_gnu = allele_lst[1]+'__'+allele_lst[2]+'__'+(allele_lst[3]).replace(',','_')+'__'+allele_lst[5]
        ptn_seq = allele_lst[4]
        # checks if the protein sequence is already in the alleles dictionary for this genome
        # if it is, add this protein name (which includes the unique prokka tag) to the list of genes associated
        # with that sequence

        # if it's not, make a new entry with that sequence as the highest level key
        # so in the case of group_999 and C2128, this would make the genes list length 2
        # So that should be fine, but then the question is what it does with this list?

        # I think it should also check here if strain_name is in alleles_dict[ptn_seq]
        # and do something accordingly if it is

        # says whether we are adding a strain to a ptn_seq entry for the first time
        # default is false
        # but that changes if either the protein seq itself has to be added (so doesn't have any strains yet)
        # or if not(strain_name in alleles_dict[ptn_seq]['isolates']) is met.

        newstrain=False

        if ptn_seq in alleles_dict:
            alleles_dict[ptn_seq]['genes'].append(ptn_name)
            if not(strain_name in alleles_dict[ptn_seq]['isolates']):
                newstrain=True
                alleles_dict[ptn_seq]['isolates'].append(strain_name)
        else:
            newstrain=True
            alleles_dict[ptn_seq] = {'genes':[ptn_name],'isolates':[strain_name],'gnu':ptn_gnu}
            # here it also iterates through the traits (so like InClade, OutOfClade in my case)
            # and then it adds another item to the alleles_dict[ptn_seq] entry that's an empty list.
            for trait in traits_header:
                alleles_dict[ptn_seq][trait] = []
        # then, it adds to a vector for that protein sequence and that trait
        # that's the value for the trait for that strain
        # which means that right now that's getting added twice for C2128

    # adds traindict[trait][strain_name] (value for that strain of that trait)
        # will it do this twice for a duplicate protein?

        # At the time of this for loop, alleles_dict[ptn_seq][trait] for a given trait
        # is either an empty list []
        # or it's a list of 1s and 0s of length(number of instances of ptn_seq that have been seen so far)
        # really though it should only be appended if the strain is new to the allele_dict[]

        # Should only append if the strain is not already represented for that ptn_seq
        if newstrain:
            for trait in traits_header:
                alleles_dict[ptn_seq][trait].append((trait_dict[trait][strain_name]))

################################
FOLDER_PATH = os.path.abspath(DIRECTORY_PATH).rsplit(OS_SEPARATOR, 1)[0] + OS_SEPARATOR
fisher_values = {}
repeat_count2 = 0
for trait in trait_dict:#(tpap,tpan,tnap,tnan)
    trait_op = open('{}{}_allele_scoring_AEC.csv'.format(FOLDER_PATH,trait),'w')
    trait_op.write(','.join(['allele','GNU score','length', 'function', 'ortholog_group',
    'trait+_allele+','trait+_allele-','trait-_allele+','trait-_allele-',
    'sensitivity','specificity','p_value','odds_ratio','Bonferroni_Pv',
    'trait+_allele+_isolates','trait+_allele-_isolates','trait-_allele+_isolates',
    'trait-_allele-_isolates','Hits\n']))
    trait_results = {}
    trait_fisher_values = {}
    repeat_count = 0

    # tp = Total number of positives (read in from the traits file)
    # tn = Total number of positives (read in from the traits file)
    tp, tn = trait_dict[trait]['count']

    # '+ves' should be a full list of isolates for which trait is + (according to traits file)
    # '-ves' should be a full list of isolates for which trait is - (according to traits file)
    # positives = list of unique strain names that are positive for the trait
    # negatives = list of unique strain names that are negative for the trait
    positives = set(trait_dict[trait]['+ves'])
    negatives = set(trait_dict[trait]['-ves'])

    # Iterates through each allele (protein seq)
    for allele in alleles_dict:
        # ptn_names is the full list of unique protein
        ptn_names = alleles_dict[allele]['genes']

        # so the total # of allele-positive isolates are correct because that's using set() on the list of strains
        ptn_isolates = set(alleles_dict[allele]['isolates'])

        # and then ptn_1s_0s is the length of the list of proteins matching each sequence
        # (including duplicates per genome)
        # Not anymore
        ptn_1s_0s = alleles_dict[allele][trait]

        # so, tpap here is overcounting if two proteins from one genome are in alleles_dict[allele]['genes']
        # and therefore alleles_dict[allele]['trait']
        tpap = ptn_1s_0s.count('1')
        tnap = ptn_1s_0s.count('0')

        # this is calculated correctly but would be underestimated if there were duplicates
        # how could it be negative?

        # it doesnt seem possible for the trait-positive, allele-positive genomes to be greater
        # than the total number of trait-positive genomes
        tpan = tp - tpap

        # list of isolates in ptn_isolates that are positive for the allele
        # (should be correct)
        tpap_isolates = ' '.join(list(ptn_isolates - negatives))#t+a+

        # list of isolates in 'positives' that are negative for the allele
        # (not in the list of isolates from alelle_dict)
        tpan_isolates = ' '.join(list(positives - ptn_isolates))#t+a-

        # list of isolates in ptn_isolates that are negative for the trait
        tnap_isolates = ' '.join(list(ptn_isolates - positives))#t-a+

        # "true negatives": negative trait and negative allele
        tnan_isolates = ' '.join(list(negatives - ptn_isolates))#t-a-
        # ^^^ these should all be correct despite the duplicate
        # this should also be correct
        isolates_results = [tpap_isolates,tpan_isolates,tnap_isolates,tnan_isolates, (' '.join(ptn_names).replace(',','_'))]

        # this should also be correct
        alleles_dict[allele]['isolates_results'] = isolates_results

        # if number of trait-positive, allele-positives > number of trait positives??
        # not possible
        # shouldn't need this if statement after the bug is fixed

        #################
        if tpan < 0:
            tpan = 0
        #################

        # trait negatives - trait negatives with positive alleles (should always be >=0)
        tnan = tn - tnap

        #################
        if tnan < 0:
            tnan = 0
        try:
            sen = (float(tpap) / tp * 100)
        except:
            sen = 0.0
        try:
            spc = (float(tnan) / tn * 100)
        except:
            spc = 0.0
        cont_table = [[tpap,tpan],
                     [tnap,tnan]]
        cont_tuple = (tpap,tpan,tnap,tnan)
        if cont_tuple in fisher_values:
            odds_ratio = fisher_values[cont_tuple][0]
            p_value = fisher_values[cont_tuple][1]
            repeat_count2 += 1
            if cont_tuple in trait_fisher_values:
                repeat_count += 1
            else:
                trait_fisher_values[cont_tuple] = fisher_t
        else:
            try:
                fisher_t = st.fisher_exact(cont_table)
            except:
                print(cont_table)
            p_value = fisher_t[1]
            odds_ratio = fisher_t[0]
            fisher_values[cont_tuple] = fisher_t
            trait_fisher_values[cont_tuple] = fisher_t
        trait_results[allele] = [str(i) for i in [tpap,tpan,tnap,tnan,sen,spc,p_value,odds_ratio]]
    #print('trait repeat_count',repeat_count)
    #print('total repeat_count',repeat_count2)
    #print('alleles',len(alleles_dict))
    number_of_tests = len(alleles_dict) - repeat_count
    #print(number_of_tests)
    for allele2 in trait_results:
        trait_results[allele2].append(str(min(
            float(trait_results[allele2][6]) * number_of_tests , 1.0)))
        trait_op.write(allele2+','+','.join((alleles_dict[allele2]['gnu']).split('__'))+
        ','+','.join(trait_results[allele2])+','+','.join(alleles_dict[allele2]['isolates_results'])+'\n')
    trait_op.close()
    #print(list(trait_fisher_values.keys()))
#print(len(list(fisher_values.keys())))
################################
def Setup_results(genedic, traitsdic, collapse):
    """
    This is the method that actually does all the counting of genes,
    calculation of p-values and post-test adjustment of p-values.
    The majority of running time is currently spent doing Fisher's exact
    test. The running time is somewhat shortened by storing the p-values
    of known 2x2 cell distributions, so that identical Fisher's tests
    will not have to be run.
    """
    # Need to create one results dictionary for each trait

    all_traits = {}
    gene_trait_combinations = {}
    fisher_calculated_values = {}
    for trait in traitsdic:
        all_traits[trait] = {}
        # Also need a list of all p-values to perform stepwise FDR corr
        p_value_list = []
        log.info("Gene-wise counting and Fisher's exact tests for trait: "
              "%s" % str(trait))

        # We also need a number of tests variable for each genedic.
        # (The number of tests may not be the same, depending on how
        # many genes are in 0 or all strains.)
        number_of_tests = len(genedic)
        # Need initial number of tests for status bar
        initial_number_of_tests = number_of_tests

        # Additionally, we need a dictionary to, for each gene
        # hold the gene-trait status (e.g "AB" or "ab" of each strain)
        gene_trait_combinations[trait] = {}
        distributions = {}

        i = 0  # Progress
        for gene in genedic:
            # Status:
            sys.stdout.write("\r{:.2%}".format(float(i) /
                                               initial_number_of_tests))
            sys.stdout.flush()
            i += 1  # Progress

            stats = Perform_statistics(traitsdic[trait], genedic[gene])
            stat_table = stats["statistics"]

            num_pos = stat_table["tpgp"] + stat_table["tpgn"]
            num_neg = stat_table["tngp"] + stat_table["tngn"]

            if (stat_table["tpgp"] + stat_table["tngp"]) == 0:
                # No included isolates have the gene
                # Remove 1 from the number of tests
                number_of_tests -= 1
                # proceed to the next gene
                continue

            if (stat_table["tpgn"] + stat_table["tngn"]) == 0:
                # All included isolates have the gene
                number_of_tests -= 1
                continue

            # IF gene_trait_combination exists already, ie the current
            # gene is perfectly correlated with another gene, collapse
            # the two into a single unit

            # Should consider adding a softer correlation required than
            # 100%, i.e. genes that are highly co-distributed should
            #  perhaps be merged earlier?
            if stats["hash"] in distributions and collapse:
                # Find out which gene is correlated
                corr_gene = distributions[stats["hash"]]
                # Collapse the two genes
                newgenename = corr_gene + "--" + gene
                distributions[stats["hash"]] = newgenename
                # Remove 1 from the number of tests
                number_of_tests -= 1
                mergedgenes = True

                del(gene_trait_combinations[trait][corr_gene])
                gene_trait_combinations[trait][newgenename] = \
                                                     stats["gene_trait"]
            else:
                distributions[stats["hash"]] = gene
                mergedgenes = False
                gene_trait_combinations[trait][gene] = \
                                                     stats["gene_trait"]

            obs_table = [[stat_table["tpgp"],
                          stat_table["tpgn"]],
                         [stat_table["tngp"],
                          stat_table["tngn"]]]
            obs_tuple = (stat_table["tpgp"],
                         stat_table["tpgn"],
                         stat_table["tngp"],
                         stat_table["tngn"])
            if obs_tuple in fisher_calculated_values:
                odds_ratio = fisher_calculated_values[obs_tuple][0]
                p_value = fisher_calculated_values[obs_tuple][1]
            else:
                fisher = ss.fisher_exact(obs_table)
                p_value = fisher[1]
                odds_ratio = fisher[0]
                fisher_calculated_values[obs_tuple] = fisher

            if not mergedgenes:
                all_traits[trait][gene] = {
                "NUGN": genedic[gene]["Non-unique Gene name"],
                "Annotation": genedic[gene]["Annotation"],
                "tpgp": stat_table["tpgp"],
                "tngp": stat_table["tngp"],
                "tpgn": stat_table["tpgn"],
                "tngn": stat_table["tngn"],
                "sens": ((float(stat_table["tpgp"]) / num_pos * 100) if
                    num_pos > 0 else 0.0),
                "spes": ((float(stat_table["tngn"]) / num_neg * 100) if
                    num_neg > 0 else 0.0),
                "OR": odds_ratio,
                "p_v": p_value}
                p_value_list.append((gene, p_value))
            else:
                temp = all_traits[trait][corr_gene]
                del(all_traits[trait][corr_gene])
                all_traits[trait][newgenename] = {
                "NUGN": temp["NUGN"] + "--" +
                        genedic[gene]["Non-unique Gene name"],
                "Annotation": temp["Annotation"] + "--" +
                              genedic[gene]["Annotation"],
                "tpgp": stat_table["tpgp"],
                "tngp": stat_table["tngp"],
                "tpgn": stat_table["tpgn"],
                "tngn": stat_table["tngn"],
                "sens": ((float(stat_table["tpgp"]) / num_pos * 100) if
                    num_pos > 0 else 0.0),
                "spes": ((float(stat_table["tngn"]) / num_neg * 100) if
                    num_neg > 0 else 0.0),
                "OR": odds_ratio,
                "p_v": p_value}
                p_value_list.append((newgenename, p_value))


        sys.stdout.write("\r100.00%")
        sys.stdout.flush()
        sys.stdout.write("\n")
        log.info("Adding p-values adjusted for testing multiple hypotheses")

        # Now calculate Benjamini-Hochberg and Bonferroni p-values
        # Sorted list of tuples: (gene, p-value)
        # s_p_v = abbreviation for sorted_p_values
        s_p_v = sorted(p_value_list, key=lambda x: x[1])

        # Find out which p-values are ties
        # Note: Changed from step-down to step-up
        # (from least significant to most significant)
        tie = [s_p_v[i-1][1] == s_p_v[i][1] for i in
            xrange(1, len(s_p_v))]
        # Initialize dics of corrected p values
        # bh_c_p_v = abbreviation for bh_corrected_p_values
        bh_c_p_v = {}
        # The least significant gene is entered into the dic
        bh_c_p_v[s_p_v[len(s_p_v)-1][0]] = last_bh = s_p_v[len(s_p_v)-1][1]

        for (ind, (gene, p)) in reversed(list(enumerate(s_p_v[:-1]))):
            bh_c_p_v[gene] = min([last_bh, p*number_of_tests/(ind+1.0)]\
                                ) if not tie[ind] else last_bh
            last_bh = bh_c_p_v[gene]

        # Now add values to dictionaries:
        for gene in all_traits[trait]:
            all_traits[trait][gene]["B_p"] = min(
                all_traits[trait][gene]["p_v"] * number_of_tests , 1.0)
            all_traits[trait][gene]["BH_p"] = min(bh_c_p_v[gene], 1.0)

    return {"Results": all_traits,
            "Gene_trait_combinations": gene_trait_combinations}

def Perform_statistics(traits, genes):
    """
    Helper function for Setup_results.
    The method that adds the presence/absence status for each
    trait-status combination.
    """
    r = {"tpgp": 0, "tpgn": 0, "tngp": 0, "tngn": 0}
    # tpgn = trait positive, gene negative, etc
    gene_trait = {}
    distribution_hash = ""
    for t in traits:  # For each strain
        # If missing data, continue without including isolate
        # Note, this should no longer be invoked, since missing-value
        # isolates are not in the traits dic to begin with
        if traits[t] in ["NA","-","."]:
            gene_trait[t] = "NA"
            distribution_hash += str(genes[t])
            continue

        try:
            if int(traits[t]) == 1 and genes[t] == 1:
                r["tpgp"] += 1
                distribution_hash += "1"
                gene_trait[t] = "AB"
            elif int(traits[t]) == 1 and genes[t] == 0:
                r["tpgn"] += 1
                distribution_hash += "0"
                gene_trait[t] = "aB"
            elif int(traits[t]) == 0 and genes[t] == 1:
                r["tngp"] += 1
                distribution_hash += "1"
                gene_trait[t] = "Ab"
            elif int(traits[t]) == 0 and genes[t] == 0:
                r["tngn"] += 1
                distribution_hash += "0"
                gene_trait[t] = "ab"
            else:
                sys.exit("There was a problem with comparing your "
                "traits and gene presence/absence files. Make sure you "
                "have formatted the traits file to specification and "
                "only use 1s and 0s, as well as NA, - or . for missing "
                "data. Also make sure the Roary file contains empty "
                "cells for non-present genes and non-empty text cells "
                "for present genes.")
        except KeyError:
            sys.stdout.write("\n")
            log.critical("CRITICAL: Could not find %s in the genes "
            "file." % str(t) )
            #log.warning("CRITICAL: Could not find %s in the genes file" % str(t))
            sys.exit("Make sure strains are named the same in your "
            "traits file as in your gene presence/absence file")
    return {"statistics": r, "hash": int(distribution_hash, 2),
            "gene_trait": gene_trait}
