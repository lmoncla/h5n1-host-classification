import sys, subprocess, glob, os, shutil, re, importlib, Bio
from subprocess import call
from Bio import SeqIO
import random
from random import shuffle
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('--input_alignment', type=str, help='alignment file from augur align')
parser.add_argument('--alignment_name', type=str, help='name of input alignment')
parser.add_argument('--output_file', type=str)
parser.add_argument('--template_xml', type=str, help='xml template')
parser.add_argument('--template_alignment_name')

args = parser.parse_args()
input_alignment = args.input_alignment
alignment_name = args.alignment_name
output_file = args.output_file
template_xml = args.template_xml
template_alignment_name = args.template_alignment_name


def read_in_alignment(input_alignment): 
    alignment_dict = {}
    for seq in SeqIO.parse(input_alignment, "fasta"): 
        alignment_dict[seq.id] = str(seq.seq)
    return(alignment_dict)



def generate_alignment_block(alignment, alignment_name): 
    seqs_list = []    
    for sequence in alignment:
        dat1 = "        <sequence id="+"\"seq_"+sequence+"\" taxon=\""+sequence+"\" totalcount=\"4\" value=\""+alignment[sequence]+"\"/>"
        seqs_list.append(dat1)
    
    new_line2 = "\n".join(seqs_list)
    new_line1 = "id=\""+alignment_name + "\"\nspec=\"Alignment\"\nname=\"alignment\">\n"
    alignment_block = new_line1+ new_line2 + "\n"
        
    return(alignment_block)




def generate_trait_block(alignment, alignment_name, trait_location, trait_block_identifier):
                                        
    traits_list = []
    for strain in alignment:
        traits_list.append(strain + "=" + strain.split("|")[trait_location])
    combined_traits = ",".join(traits_list)
    traits_string = alignment_name + "\" spec=\"mascot.util.InitializedTraitSet\" traitname=\"type\" value=\""
    new_data = traits_string + combined_traits
    trait_block = trait_block_identifier + new_data + "\">\n"
    
    return(trait_block)




def generate_date_block(alignment, alignment_name, decimal_date_location, long_date_location, date_block_identifier): 

    intro_string = date_block_identifier + alignment_name + "\" spec=\"beast.base.evolution.tree.TraitSet\" traitname=\"date\" value=\""
                    
    dates_list = []
    for strain in alignment:
        decimal_date = strain.split("|")[decimal_date_location]
        long_date = strain.split("|")[long_date_location]
        
        # test for ambiguous dates that have only year; set those to halfway through year 
        if "-XX-XX" in long_date:
            year = float(long_date.split("-")[0])
            date = long_date.split("-")[0] + ".5"
        else: 
            date = decimal_date 
            
        dates_list.append(strain + "=" + date)
    
    combined_dates = ",".join(dates_list)
    date_block = intro_string + combined_dates + "\">\n"

    return(date_block)


"""the following blocks are for estimating tip dates"""
def determine_most_recent_tip(alignment, decimal_date_location):
    all_dates = []
    for strain in alignment: 
        decimal_date = strain.split("|")[decimal_date_location]
        decimal_date = float(decimal_date)
        all_dates.append(decimal_date)
    
    maxdate = round(max(all_dates), 3)
        
    return(maxdate)




"""the prior on the tip dates should be set in years, not in years from present"""

def determine_date_bounds(most_recent_tip_date, tip_year):
    
    # for a tip with a year-only tip, the maximum furthest back in time it could be is January 1st of that year
    #upper_bound = str(round((most_recent_tip_date - tip_year), 4))
    upper_bound = str(tip_year + 0.9959)
    
    # the most recent it could be is December 31st of that year; in decimal dates, that is year.9959
    #lower_bound = str(round((most_recent_tip_date - (tip_year + 0.9959)), 4))
    lower_bound = str(tip_year)
    
    return(upper_bound, lower_bound)




def generate_mrca_distribution_block(strain, alignment_name, upper_prior_bound_year, lower_prior_bound_year, count):
    
    strain_name = strain.split("|")[0]
    strain_prior_id = strain_name + ".prior"
    
    line1 = "            <distribution id=\"insert_prior_id\" spec=\"beast.base.evolution.tree.MRCAPrior\" tipsonly=\"true\" tree=\"@Tree.t:insert_tree_id\">"
    line2 = "                <taxonset id=\"insert_id\" spec=\"TaxonSet\">"
    line3 = "                    <taxon id=\"insert_taxon_id\" spec=\"Taxon\"/>"
    line4 = "                </taxonset>"
    line5 = "                <Uniform id=\"Uniform"
    line6 = "            </distribution>"
    
    line1 = line1.replace("insert_prior_id",strain_prior_id)
    line1 = line1.replace("insert_tree_id",alignment_name)
    line2 = line2.replace("insert_id",strain_name)
    line3 = line3.replace("insert_taxon_id",strain)
    line5 = line5.replace("                <Uniform id=\"Uniform", "                <Uniform id=\"Uniform."+ str(count) + "\" lower=\"lower_prior_bound_year\" name=\"distr\" upper=\"upper_prior_bound_year\"/>")
    line5 = line5.replace("lower_prior_bound_year",lower_prior_bound_year)
    line5 = line5.replace("upper_prior_bound_year",upper_prior_bound_year)

    block = "\n".join([line1,line2,line3,line4,line5,line6])
    
    return(block)




def generate_operator(strain, alignment_name):
    strain_name = strain.split("|")[0]
    
    operator_line = "        <operator id=\"tipDatesSampler"+strain_name+"\" spec=\"TipDatesRandomWalker\" taxonset=\"@"+strain_name+"\" tree=\"@Tree.t:"+alignment_name+"\" weight=\"1.0\" windowSize=\"1.0\"/>"    
    return(operator_line)



def generate_logger(strain):
    strain_name = strain.split("|")[0]
    logger_line = "            <log idref=\""+strain_name+".prior\"/>"
    return(logger_line)



def generate_estimated_tip_dates_blocks(alignment, alignment_name, most_recent_tip): 
    
    count = 10 
    tipdate_estimation_blocks = []
    tipdate_estimation_loggers = []
    tipdate_estimation_operators = []
    
    for strain in alignment:
        decimal_date = strain.split("|")[decimal_date_location]
        long_date = strain.split("|")[long_date_location]
        
        if "-XX-XX" in long_date:
            count += 1
            year = float(long_date.split("-")[0])
            
            upper_bound_date, lower_bound_date = determine_date_bounds(most_recent_tip, year)
            block = generate_mrca_distribution_block(strain, alignment_name, upper_bound_date, lower_bound_date, count)
            tipdate_estimation_blocks.append(block)

            operator = generate_operator(strain, alignment_name)
            tipdate_estimation_operators.append(operator)

            logger = generate_logger(strain)
            tipdate_estimation_loggers.append(logger)

    return(tipdate_estimation_blocks, tipdate_estimation_loggers, tipdate_estimation_operators)



# read in the current date 
from datetime import date
today = date.today()
current_date = str(today.strftime("%Y-%m-%d"))


alignment_dict = read_in_alignment(input_alignment)


# index of the trait in the string
trait_location = 12
decimal_date_location = 1
long_date_location = 2
strain_location = 0

most_recent_tip = determine_most_recent_tip(alignment_dict, decimal_date_location)


tipdate_blocks, tipdate_loggers, tipdate_operators = generate_estimated_tip_dates_blocks(alignment_dict, alignment_name, most_recent_tip)




trait_block_identifier = "                        <typeTrait id=\"typeTraitSet.t:"
date_block_identifier = "                <trait id=\"dateTrait.t:"
lines_to_skip = ("name=\"alignment","spec=\"Alignment","        <sequence id")

with open(output_file, "w") as outfile: 
    outfile.write("")

with open(template_xml, "r") as infile: 
    for line in infile: 
        if line.startswith("id=\""):
            alignment_block = generate_alignment_block(alignment_dict, alignment_name)
            with open(output_file, "a") as outfile: 
                outfile.write(alignment_block)
            
        elif line.startswith(trait_block_identifier):
            trait_block = generate_trait_block(alignment_dict, alignment_name, trait_location, trait_block_identifier)
            with open(output_file, "a") as outfile: 
                outfile.write(trait_block)

        elif line.startswith(date_block_identifier): 
            date_block = generate_date_block(alignment_dict, alignment_name, decimal_date_location, long_date_location, date_block_identifier)
            with open(output_file, "a") as outfile: 
                outfile.write(date_block)

        elif "<!--  add tip calibration here -->" in line:
            tipdate_block = "\n".join(tipdate_blocks)
            with open(output_file, "a") as outfile: 
                outfile.write(line + "\n" + tipdate_block + "\n")
            
        elif "<!-- insert tip calibration operators -->" in line:
            tipdate_operators_block = "\n".join(tipdate_operators)
            with open(output_file, "a") as outfile: 
                outfile.write(line + tipdate_operators_block + "\n")
                
        elif "<!-- insert tip calibration loggers -->" in line:
            tipdate_loggers_block = "\n".join(tipdate_loggers)
            with open(output_file, "a") as outfile: 
                outfile.write(line + tipdate_loggers_block + "\n")

        
        elif line.startswith(lines_to_skip):
            pass
        
        else: 
            line = line.replace(template_alignment_name, alignment_name)
            with open(output_file, "a") as outfile: 
                outfile.write(line)
