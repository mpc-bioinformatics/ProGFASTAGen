#!/bin/env python

import argparse
import sys
import csv
csv.field_size_limit(sys.maxsize)
import re
import itertools
import os

import pandas as pd


# Constants
HYDROGEN_MONO_MASS = 1.007825035
FASTA_HEADER_PATTERN = r"(.*?[pg|tr|sp|lcl|ref])\|([a-zA-Z0-9\-\_]+)\|"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-comet_txt", help="The txt-output from Comet (always needed)")
    parser.add_argument("-perc_tsv", help="The pin-output from Comet (optional, if not provided)")
    parser.add_argument("-use_n_hits", help="Use second/third/n-th hits for fdr caluclation (if set higher, the FDR will be skewed due to n!=1 hits which have a a higher score then the n=1 hits, possibly yielding less identied spectra)", default=1, type=int)
    parser.add_argument("-fasta", help="The FASTA-File which was used for identification in Comet")
    parser.add_argument("-decoy_string", help="The decoy string which was used during identification in Comet ", default="DECOY_")
    parser.add_argument("-fdr", help="The FDR (based on the qvalue) where to cut-off identification results from Comet", type=float)
    parser.add_argument("-out_all_tsv", help="Output table of all psms including decoys", default=None)
    parser.add_argument("-out_no_decoys_tsv", help="Output table with no decoys removed", default=None)
    parser.add_argument("-out_fdr_cutoff_tsv", help="Output of the cutoff-table with no decoys on the fdr", default=None)
    return parser.parse_args()


def apply_fdr(entries: list, protein_idx: int, consider_entry_based_on_num = lambda x: True, decoy_string = "DECOY_"):
    '''
    Calculate the qvalue. Append this in additional_information
    
    Here the following columns are generated: "qvalue", "is_decoy" and "fasta_acc"

    The function "consider_entry_based_on_num" takes an entry as input and can decide whether an entry should 
    be considered for FDR-Calculation (e.g. for multiple peptide hits on one spectrum, of only the best hit or other schould be considered)
    '''

    # Initialize lists, which will be returned, containing the extra information per entry
    additional_headers = ["qvalue", "is_decoy", "fasta_acc"]
    additional_info = []
        
    # Counters for targets and decoys
    num_hits = 0
    num_decs = 0

    for entry in entries:
        # Iterate over each entry

        # Get all accessions and protein_hits associated to this PSM
        matches = re.finditer(FASTA_HEADER_PATTERN, entry[protein_idx], re.MULTILINE)

        # Get for each fitting protein_hit for this psm its accession and whole line (for decoy check) and save it
        is_decoy = []
        accessions = []
        for matchNum, match in enumerate(matches, start=1):
            is_decoy.append(True if decoy_string in match.groups()[0] else False)
            accessions.append(match.groups()[1])

        # Check if this entry should influence the FDR-Calculation
        if consider_entry_based_on_num(entry):
            # If yes, update the #decoys and #targets
            if all(is_decoy):
                num_decs += 1
            else:
                num_hits += 1 

            # Append Info
            additional_info.append([
                num_decs / (num_decs + num_hits),  # qvalue: This is currently not the qvalue and will be updated later
                all(is_decoy),  # is decoy
                ",".join(accessions)  # fasta_acc
            ])
        else:
            # If not, simply ignore this value, but append information
            additional_info.append([
                additional_info[-1][0],  # qvalue: We simply copy from the last entry
                None,  # is_decoy: Denote that we ignored this entry
                ",".join(accessions)  # fasta_acc
            ])

    # In qvalue, we need to update the entries to actually have the qvalue
    # We simply iterate the list from low to high and take the lowest value of qvalue and add it into the next entries 
    qvalue = additional_info[-1][0]
    for add_in in additional_info[::-1]:
        if add_in[0] > qvalue:
            add_in[0] = qvalue
        else:
            qvalue = add_in[0]

    # Return the additonal column headers and its values (fitting to the entries order)
    return additional_headers, additional_info


def get_header_desc_from_fasta_accessions(protein_ids: set, fasta: str, decoy_string: str):
    ''' Memory-efficient solution of getting the FASTA-Header descriptions '''
    fasta_d = dict()
    # Open FASTA-file
    with open(fasta, "r") as in_fasta:
        # Iterate over each line
        for l in in_fasta:
            # Check if entry is a header and not a decoy 
            if l.startswith(">") and not l.startswith(">" + decoy_string):
                # Only get the Accession
                line_split = l.split("|", 2)
                # 
                if line_split[1] in protein_ids:
                    fasta_d[line_split[1]] = line_split[-1][:-1]
    return fasta_d


def get_exp_mass_to_charge(mass: float, charge: float):
    ''' Convert mass and charge to mass_to_charge '''
    return (mass + (HYDROGEN_MONO_MASS*charge)) / charge


def get_from_comet_txt(df, entry_id: list, query: str):
    ''' Wrapper to query the dataframe returning a single value, based on query and entry_id (--> [scan, charge, num]) '''
    # Get if available:
    res = df[
            (df["scan"] == int(entry_id[0])) & 
            (df["charge"] == int(entry_id[1])) & 
            (df["num"] == int(entry_id[2]))
        ][query]
    if len(res) == 0 and query in ("exp_neutral_mass", "retention_time_sec"):
        # If the exact entry is not found, it is also okay to look for the 
        # retention time and exp_neutral_mass in the other entries (since those are across all
        # nums equal in comet)
        # Could also be a BUG in comet
        res = df[
            (df["scan"] == int(entry_id[0])) & 
            (df["charge"] == int(entry_id[1]))
        ][query]
    
    return str(res.iloc[0])


if __name__ == "__main__":
    '''
    Opens the identification results either from Comet or Percolator (from Comet input), applys the 
    specified FDR (e.g 0.01) and append the additional columns at the end: 
    "source_file" --> The filename, from which file the PSM was retrieved (base name)
    "qvalue"  --> The calcuated qvalue (via #Decoys / #Target)
    "fasta_acc" --> The FASTA-Accession (>XX|accession|XXXXX)
    "fasta_desc" --> The FASTA-Description (>XX|XXXXX|Description)
    "used_score" --> The value of the used score (for Comet: xcorr, for Percolator: svm_score)
    "charge" --> The charge state of the PSM
    "retention_time" --> The RT of the PSM
    "exp_mass_to_charge" --> THe experimental mass to charge (percursor_mass + Hydrogen / charge)
    '''
    # Parse arguments
    args = parse_args()

    # Get the psms either from comet (or if provided from percolator)
    if args.perc_tsv is None:
        # No percolator tsv provided
        with open(args.comet_txt, "r") as in_comet:
            # Load in comet results
            next(in_comet)  # Skip the comet header

            # Load as tsv
            csv_in = csv.reader(in_comet, delimiter="\t")

            # Get headers to retrieve specific headers indices
            headers = next(csv_in)

            # Sort entries based on the xcorr
            score_idx = headers.index("xcorr")
            entries = sorted([x for x in csv_in], key=lambda x: float(x[score_idx]), reverse=True)
            
            # Apply FDR and calculate the qvalue 
            protein_idx = headers.index("protein")  # Protein-Header (contains the <XX|"accession"|XXXX (comma seperated) in FASTA)

            # Provide a filter to inlcude only num=1 hits (best hit from Comet)
            num_idx = headers.index("num")
            qval_headers, qval_entries = apply_fdr(
                entries, protein_idx=protein_idx, decoy_string=args.decoy_string,
                consider_entry_based_on_num= lambda x: True if int(x[num_idx]) < args.use_n_hits + 1 else False
            )

            ### Adding additonal columns
            # In Comet most additionaly information is already available in the txt file
            additional_headers = ["fasta_desc", "used_score", "retention_time", "exp_mass_to_charge"]  # We add these additionaly columns. "scan", "num" and "charge" is left out since, Comet, already provides such columns
            additional_entries = []
            
            # First retrieve FASTA Accession <-> Description for all found accessions
            fasta_acc_idx = qval_headers.index("fasta_acc")
            fasta_acc_desc = get_header_desc_from_fasta_accessions(
                set(itertools.chain(*[x[fasta_acc_idx].split(",") for x in qval_entries])), 
                fasta=args.fasta, decoy_string=args.decoy_string
            )

            # Fill the additional entries
            retention_time_idx = headers.index("retention_time_sec")
            charge_idx = headers.index("charge")
            exp_neutral_mass_idx = headers.index("exp_neutral_mass")
            for entry, qval_entry in zip(entries, qval_entries):
                additional_entries.append(
                    [
                        ",".join([fasta_acc_desc[x] for x in qval_entry[fasta_acc_idx].split(",")]),  # fasta_desc: As in accession, seperated by ,
                        entry[score_idx],  # used_score: To show which was used to calculate to FDR
                        entry[retention_time_idx],  # retention_time: Copied from retention_time_sec
                        get_exp_mass_to_charge(float(entry[exp_neutral_mass_idx]), int(entry[charge_idx])) # exp_mass_to_charge: Calcuated based on mass and charge
                    ]
                )

    else: 
        # Percolater tsv is provided, use this one instead!
        with open(args.perc_tsv, "r") as in_perc:

            # Load as tsv
            csv_in = csv.reader(in_perc, delimiter="\t")

            # Get headers to retrieve specific headers indices
            headers = next(csv_in)

            # Sort entries based on the svm_score
            score_idx = headers.index("svm_score")  # We use the xcorr for sorting
            entries = sorted([x for x in csv_in], key=lambda x: float(x[score_idx]), reverse=True)

            # Apply FDR and calculate the qvalue 
            protein_idx = headers.index("protein_id")  # Protein-Header (contains the <XX|"accession"|XXXX (comma seperated) in FASTA)

            # Provide a filter to inlcude only num=1 hits (best hit from Comet)
            psm_id_idx = headers.index("psm_id")
            qval_headers, qval_entries = apply_fdr(
                entries, protein_idx=protein_idx, decoy_string=args.decoy_string,
                consider_entry_based_on_num= lambda x: True if int(x[psm_id_idx].split("_")[-1]) < args.use_n_hits + 1 else False 
            )

            ### Adding additional columns
            # In Percolator we need to add more information (and supplement with info from Comet)
            additional_headers = ["scan", "num", "charge", "fasta_desc", "used_score", "retention_time", "exp_mass_to_charge"]  # We add these additionaly columns. "scan", "num" and "charge" is left out since, Comet, already provides such columns
            additional_entries = []
            
            # First retrieve FASTA Accession <-> Description for all found accessions
            fasta_acc_idx = qval_headers.index("fasta_acc")
            fasta_acc_desc = get_header_desc_from_fasta_accessions(
                set(itertools.chain(*[x[fasta_acc_idx].split(",") for x in qval_entries])), 
                fasta=args.fasta, decoy_string=args.decoy_string
            )

            # Fill the additional entries
            # We need to query the txt-file for some entries, therefor we load it in pandas
            comet_results_df = pd.read_csv(args.comet_txt, sep="\t", header=1)
            for entry, qval_entry in zip(entries, qval_entries):
                # Get scan, charge and num from the psm_id
                scan, charge, num = entry[psm_id_idx].split("_")[-3:]
                additional_entries.append(
                    [
                        scan,  # scan
                        num,  # num
                        charge,  # charge
                        ",".join([fasta_acc_desc[x] for x in qval_entry[fasta_acc_idx].split(",")]),  # fasta_desc: As in accession, seperated by ,
                        entry[score_idx],  # used_score: To show which was used to calculate to FDR
                        get_from_comet_txt(comet_results_df, [scan, charge, num], "retention_time_sec"),  # retention_time: Retrieved from the comet-results
                        get_exp_mass_to_charge(
                            float(get_from_comet_txt(comet_results_df, [scan, charge, num], "exp_neutral_mass")), 
                            int(charge)
                        )  # exp_mass_to_charge: Calcuated based on mass and charge (mass, retrieved from Comet results)
                    ]
                )

    # Add source_name
    sn_header = ["source_name", ]
    sn_entries = [".".join(args.comet_txt.split(os.sep)[-1].split(".")[:-1])] * len(entries)

    ### Write out results
    # Write all with decoys
    if args.out_all_tsv is not None: 
        with open(args.out_all_tsv, "w") as out_file:
            csv_out = csv.writer(out_file, delimiter="\t")
            csv_out.writerow(sn_header + headers + qval_headers + additional_headers)
            for a, b, c, d in zip(sn_entries, entries, qval_entries, additional_entries):
                csv_out.writerow(
                        [a] + b + c + d
                    )

    # Write all without decoys
    if args.out_no_decoys_tsv is not None: 
        with open(args.out_no_decoys_tsv, "w") as out_file:
            decoy_idx = qval_headers.index("is_decoy")
            csv_out = csv.writer(out_file, delimiter="\t")
            csv_out.writerow(sn_header + headers + qval_headers + additional_headers)
            for a, b, c, d in zip(sn_entries, entries, qval_entries, additional_entries):
                if c[decoy_idx] is not True:
                    csv_out.writerow(
                        [a] + b + c + d
                    )

    # Write all without decoys up to fdr
    if args.out_fdr_cutoff_tsv is not None: 
        with open(args.out_fdr_cutoff_tsv, "w") as out_file:
            decoy_idx = qval_headers.index("is_decoy")
            qvalue_idx = qval_headers.index("qvalue")
            csv_out = csv.writer(out_file, delimiter="\t")
            csv_out.writerow(sn_header + headers + qval_headers + additional_headers)
            for a, b, c, d in zip(sn_entries, entries, qval_entries, additional_entries):
                if c[decoy_idx] is not True and c[qvalue_idx] < args.fdr:
                    csv_out.writerow(
                        [a] + b + c + d
                    )
