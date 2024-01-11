#!/bin/env python

import sys
import csv
csv.field_size_limit(sys.maxsize)
import re 
import argparse

from collections import defaultdict


# Constants
REGEX_PG_HEADER = r"[A-Z0-9-_]+?\(.*?\)"
REGEX_GET_VARMODS = r"VARMOD\[.*?\](?:,|)"
STRINGS_CONSIDERED_FEATURES = ["CONFLICT", "SIGNAL", "INIT_MET", "PROPEP", "PEPTIDE", "MUTAGEN", "VARIANT", "CHAIN"]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input_ident_tsv", help="The input-MGF, from which the queries for querying protein-graphs should be generated")
    parser.add_argument("-out_summary", help="Output summary as tsv")
    parser.add_argument("-out_unique", help="Output unique psms as tsv")
    parser.add_argument("-out_shared", help="Output shares psms as tsv")
    parser.add_argument("-out_ft_unique", help="Output unqiue psms which can be only explained by features as tsv")
    parser.add_argument("-out_ft_shared", help="Output shared psms which can only be explained by features as tsv")
    parser.add_argument("-accession_header", choices=["fasta_acc", "fasta_desc"], help="Selecting the header-part, where the Protein information is located. For Normal searches, it is 'fasta_accc'. For Protein-Graph-FASTA, the information can be retrieved in 'fasta_desc'", default="fasta_acc")
    parser.add_argument("-same_protein_as_unique", help="In some cases it may happen that a Protein has a shared peptide across itself. In such cases, we can either report as a unique hit (true) or as a shared hit (false). In Protein-Graphs-FASTAs, (especially when using variable modifications, ) this parameter should be set to true (which is also the default)", default=True, type=bool)
    parser.add_argument("-remove_varmods", help="Remove variable modifications from the description of ProteinGraph-Headers. If set to true (default), variable modifications (which lead to haven 2 same Proteins annoting a single sequence) will be removed and such entries will be considered as unique. If set to false, for ProteinGraph-Headers, only non-modifiable peptides are consideres as unique", default=True, type=bool)
    return parser.parse_args()


def fits_in_list_fasta_acc(cell, same_protein_as_unique, _remove_varmods):
    ''' 
    Classify a simple match from the fasta_acc cell into unqiue or shared 
    (the last two do not apply and are always false)
    '''
    matches = cell.split(",")
    if same_protein_as_unique:
        # Count unique proteins instead (in case of same sequence in 1 proteins)
        matches = set(matches)

    # Decide wether it is unique or shared
    is_unique = len(matches) == 1
    is_shared = len(matches) > 1
    return is_unique, is_shared, False, False


def fits_in_list_fasta_desc(cell, same_protein_as_unique, remove_varmods):
    '''
    Classify each cell from fasta_desc (ProtGraph) into unique shared and its feature unique/shared
    '''
    # First remove all spaces 
    cell = cell.replace(" ", "")

    # Parse all accession and their descriptions (e.g. "PXXXXX(....),PYYYYY(....)"  )
    matches_regex = re.finditer(REGEX_PG_HEADER, cell, re.MULTILINE)
    prot_matches = [x.group() for x in matches_regex]
    processed_matches = prot_matches

    # Special Case, since FASTA-files CANNOT encode modifications. 
    if remove_varmods:  
        # We consider the VARMODs (since FIXMODS are globally present in the description)
        # from the same protein as unique. 
        # E.G.: >pg|ID_XXXX|P68871(42:60,mssclvg:0,),P68871(42:60,mssclvg:0,VARMOD[56:56,M:15.994915])
        # would be considered unique, since only a variable modification was applied.
        # However, for other modifications/applications it could be interesting to consider them as unique 

        # Remove varmods via regex and replace all ','
        processed_matches = [re.sub(REGEX_GET_VARMODS, "", m).replace(",", "") for m in prot_matches]
        processed_matches = set(processed_matches)

    if same_protein_as_unique:
        # Count unique proteins instead (in case of same sequence in 1 proteins)
        processed_matches = set([m[:m.find("(")] for m in prot_matches])

    # Check if it is shared or unique
    is_unique = len(processed_matches) == 1
    is_shared = len(processed_matches) > 1

    # Check if shared or unique about features
    # We have a KEYWORD-List which is checked agains. If one fits it is a feature psm
    feature_check = [any([y in x for y in STRINGS_CONSIDERED_FEATURES]) for x in prot_matches]
    is_feature_unique = (len(processed_matches) == 1) and feature_check[0]
    is_feature_shared = len(processed_matches) > 1 and all(feature_check)

    return is_unique, is_shared, is_feature_unique, is_feature_shared


def get_psms(l: list, num_idx: int, source_name_idx: int):
    ''' Helper Function to sort psms into source_name '''
    top_psms = [x for x in l if int(x[num_idx]) == 1]
    top_psms_d =  defaultdict(lambda: 0)
    for p in top_psms:
        top_psms_d[p[source_name_idx]] +=1

    return top_psms_d


if __name__ == "__main__":
    # Parse Arguments
    args = parse_args()

    # The lists where we want to sort the psms into
    unique_psms = []
    shared_psms = []
    unique_feature_psms = []
    shared_feature_psms = []

    # Read identifications file
    with open(args.input_ident_tsv, "r") as in_file:
        csv_in = csv.reader(in_file, delimiter = "\t")

        # Read header
        header = next(csv_in)

        # Set index for reading protein_info
        fasta_acc_desc_index = header.index(args.accession_header)
        source_name_index = header.index("source_name")  # For the summary file to count number of top psm hits
        num_index = header.index("num")  # To check only with the top hit

        # Get num of psm (we only want to consider the top ones)
        fasta_acc_desc_index = header.index(args.accession_header)

        # Set the classify psm method (ProtGraph-Header or simply use the accession)
        classify_psm = fits_in_list_fasta_acc if args.accession_header == "fasta_acc" else fits_in_list_fasta_desc

        # Iterate over each line and sort them into the unique, feature_unique, shared and only_feature_shared
        for l in csv_in:
            # Check if PSM is Top Hit
            if int(l[num_index]) != 1:
                # If not skip
                continue

            # Classify PSM line 
            unique, shared, feature_unique, feature_shared = classify_psm(l[fasta_acc_desc_index], args.same_protein_as_unique, args.remove_varmods)

            # Sort the PSM into the coresponding list
            if unique or feature_unique:
                if feature_unique:
                    unique_feature_psms.append(l)
                elif unique:
                    unique_psms.append(l)
            elif shared or feature_shared:
                if feature_shared:
                    shared_feature_psms.append(l)
                elif shared:
                    shared_psms.append(l)

    # Write output summary
    # Get all psms
    d_all_psms = get_psms(unique_psms + shared_psms + unique_feature_psms + shared_feature_psms, num_index, source_name_index)
    d_unique_psms = get_psms(unique_psms, num_index, source_name_index)
    d_shared_psms = get_psms(shared_psms, num_index, source_name_index)
    d_unique_ft_psms = get_psms(unique_feature_psms , num_index, source_name_index)
    d_shared_ft_psms = get_psms(shared_feature_psms, num_index, source_name_index)

    # Write into the summary tsv
    with open(args.out_summary, "w") as summary_out:
        csv_out = csv.writer(summary_out, delimiter="\t")
        csv_out.writerow(["source_name", "count_idents", "count_unique", "count_shared", "count_unique_with_features", "count_shared_with_only_features"])
        for key in d_all_psms.keys():
            csv_out.writerow(
                [
                    key,
                    d_all_psms[key],
                    d_unique_psms[key],
                    d_shared_psms[key],
                    d_unique_ft_psms[key],
                    d_shared_ft_psms[key],
                ]
            )

    # Write each psm list individually
    for psms, f_out in [
        (unique_psms, args.out_unique), (shared_psms, args.out_shared), 
        (unique_feature_psms, args.out_ft_unique), (shared_feature_psms, args.out_ft_shared)
    ]:
        if len(psms) != 0:
            with open(f_out, "w") as csv_out_file:
                csv_out = csv.writer(csv_out_file, delimiter="\t")
                csv_out.writerow(header)
                csv_out.writerows(psms)
