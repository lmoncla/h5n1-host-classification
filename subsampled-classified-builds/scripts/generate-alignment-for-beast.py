import Bio 
from Bio import SeqIO
import json
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--input_alignment', type=str, help='alignment file from augur align')
parser.add_argument('--metadata_file', type=str, help='metadata file output from clade annotations and classification')
parser.add_argument('--output_file', type=str)
parser.add_argument('--furin_annotations_file', type=str, help='json of furin annotations (yes/no)')
parser.add_argument('--furin_seqs', type=str, help='json of furin annotations')

args = parser.parse_args()
input_alignment = args.input_alignment
metadata_file = args.metadata_file
output_file = args.output_file
furin_annotations_file = args.furin_annotations_file
furin_seqs = args.furin_seqs


def read_metadata_file(metadata_file):
    output_dict = {}
    with open(metadata_file, "r") as infile: 
        for line in infile:
            
            if line.startswith("strain"):
                fields = line.split("\t")
            
            else:
                # list.index['string'] returns the index of 'string' in list
                strain = line.split("\t")[fields.index("strain")]
                host_group =  line.split("\t")[fields.index("host")].lower()
                host_species = line.split("\t")[fields.index("host_species_standardized")]
                domestic_wild = line.split("\t")[fields.index("domestic_wild")]

                accession = line.split("\t")[fields.index("isolate_id")]
                insdc_accession = line.split("\t")[fields.index("INSDC_accession")]
                clade = line.split("\t")[fields.index("h5_label_clade")]
                annotation_method = line.split("\t")[fields.index("annotation_method")]
                
                region = line.split("\t")[fields.index("region")]
                country = line.split("\t")[fields.index("country")]
                division = line.split("\t")[fields.index("division")]
                
                date = line.split("\t")[fields.index("date")]
                
                output_dict[strain] = {"accession":accession,"insdc_accession":insdc_accession,
                                       "region":region, "country":country,"division":division,
                                       "domestic_wild":domestic_wild,"host_group":host_group, 
                                       "host_species":host_species,"date":date, 
                                       "annotation_method":annotation_method, "clade":clade}
    return(output_dict)




"""read in the json file with the data on furin cleavage site"""

def read_in_cleavage_site_annotations(furin_cleavage_annotation_json):
    with open(furin_cleavage_annotation_json) as json_file:
        furin_annotations_dict = json.load(json_file)['nodes']
    return(furin_annotations_dict)




"""read in the json file with the data on furin cleavage site"""

def read_in_cleavage_site_sequences(furin_cleavage_seq_json):
    with open(furin_cleavage_seq_json) as json_file:
        furin_seq_dict = json.load(json_file)['nodes']
    return(furin_seq_dict)




"""I would like to try running these analyses with only GsGd lineages H5 viruses. The non-GsGd lineages are
incredibly diverse and also appear to be sampled differently. To see if this is causing issues with my analyses
I would like to remove them. This will allow me to keep a list of sequences that I want to not be included in the 
final alignment if I so desire."""

def read_in_unwanted_sequences(sequences_to_exclude_list):

    seqs_to_exclude = []
    with open(sequences_to_exclude_list, "r") as infile: 
        for line in infile:
            strain_name = line.split("|")[0]
            accession = line.split("|")[6]
            seqs_to_exclude.append(strain_name)
    
    seqs_to_exclude = list(set(seqs_to_exclude))
    return(seqs_to_exclude)




"""this came from here: https://stackoverflow.com/questions/6451655/how-to-convert-python-datetime-dates-to-decimal-float-years"""
from datetime import datetime as datetime
import time

def toYearFraction(date_string):
    if len(date_string) == 4:
        date = datetime.strptime(date_string, '%Y')
    elif len(date_string) == 7:
        date = datetime.strptime(date_string, '%Y-%m')
    else:
        date = datetime.strptime(date_string, '%Y-%m-%d')

    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction



def write_new_alignment_files(input_alignment, output_filename, metadata_dict, furin_annotation_dict, furin_seq_dict):
    
    with open(output_filename, "w") as outfile: 
        outfile.write("")
        
    for seq in SeqIO.parse(input_alignment, "fasta"):
        sequence = str(seq.seq)
        strain_name = seq.description
        
        # don't write out sequences we want to exclude, like non-GsGd lineage H5 virusess
            
        metadata = metadata_dict[strain_name]
        host_group = metadata_dict[strain_name]['host_group'].replace(" ","_")
        region = metadata_dict[strain_name]['region'].replace(" ","_")
        country = metadata_dict[strain_name]['country'].replace(" ","_").replace("?","")
        genbank_id = metadata_dict[strain_name]['insdc_accession'].replace(" ","_").replace("?","")
        accession = metadata_dict[strain_name]['accession'].replace(" ","_").replace("?","")
        host_species = metadata_dict[strain_name]['host_species'].replace(" ","_")
        domestic_wild = metadata_dict[strain_name]['domestic_wild'].replace(" ","_")
        date = metadata_dict[strain_name]['date']
        clade = metadata_dict[strain_name]['clade']
        
        decimal_date = toYearFraction(date.replace("-XX",""))

        # look up furin cleavage site info
        furin_site_annotation = furin_annotation_dict[strain_name]['furin_cleavage_motif']
        furin_site_sequence = furin_seq_dict[strain_name]['cleavage_site_sequence']

        new_header = "|".join([strain_name, str(decimal_date), date, host_group, region, country, genbank_id, accession, furin_site_annotation, furin_site_sequence, clade, host_species, domestic_wild])

        with open(output_filename, "a") as outfile: 
            outfile.write(">" + new_header + "\n" + sequence + "\n")





# read in the current date 

metadata_dict = read_metadata_file(metadata_file)
furin_annotations_dict = read_in_cleavage_site_annotations(furin_annotations_file)
furin_sequences_dict = read_in_cleavage_site_sequences(furin_seqs)


write_new_alignment_files(input_alignment, output_file, metadata_dict, furin_annotations_dict,furin_sequences_dict)





