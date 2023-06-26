#!/bin/env python

import sys
import csv
csv.field_size_limit(sys.maxsize)
import argparse


import xlsxwriter


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-in_psm_tsv", help="Input tsv of psms containing the columns 'source_name', 'plain_peptide', and 'fasta_desc'.")
    parser.add_argument("-out_xlsx", help="Output xlsx retruning a large heatmap, merging same found peptides per MGF, using the above mentiioned columns.")
    return parser.parse_args()


if __name__ == "__main__":
    # Parse Arguments
    args = parse_args()

    # Read identification results (containing the above mentioned columns)
    with open(args.in_psm_tsv, "r") as in_file:
        csv_in = csv.reader(in_file, delimiter="\t")

        # Read header and get the corresponding indices
        header = next(csv_in)
        source_file_idx = header.index("source_name")
        peptide_idx = header.index("plain_peptide")
        fasta_desc_idx = header.index("fasta_desc")

        # Map peptide to source_file for counting
        source_files = set()  # Tracking all source files which are added into the psm-list 
        peptide_to_source_file_dict = dict()
        for l in csv_in:
            # Go over each PSM
            source_files.add(l[source_file_idx])  # Simply add to set

            # If peptide is not tracked, initialize this entry with an empty dictionary
            if l[peptide_idx] not in peptide_to_source_file_dict:
                peptide_to_source_file_dict[l[peptide_idx]] = dict()

            # if source file is not tracked in this dicitonary of peptides (in dictionaries), initialize it as empty 
            if l[source_file_idx] not in peptide_to_source_file_dict[l[peptide_idx]]:
                peptide_to_source_file_dict[l[peptide_idx]][l[source_file_idx]] = [[], 0]

            # Map Peptide -> Source File -> (PSM, #Hits)
            peptide_to_source_file_dict[l[peptide_idx]][l[source_file_idx]][0].append(l)
            peptide_to_source_file_dict[l[peptide_idx]][l[source_file_idx]][1] += 1


        # Create rows for the final HeatMap
        source_files = sorted(list(source_files))
        summarized = [
            (
                sum([x[1] for x in source_count.values()]),  # Number of PSM matches
                len(source_count.keys()),  # Number of matches in  source files
                pep,  # Matched peptide
                next(iter(source_count.values()))[0][0][fasta_desc_idx],  # FASTA-Description fitting to peptide
                [source_count[x][1] if x in source_count else 0 for x in source_files]  # For each source file the merged peptide_count count
            ) 
            for pep, source_count in peptide_to_source_file_dict.items()
        ]
        flatted_summarized = sorted([list(x[0:4]) + x[4] for x in summarized], key=lambda x: x[0], reverse=True )

        # Set the Header for each column
        new_header = [
                "#PSMs",
                "#Source_files",
                "plain_peptide",
                "fasta_desc",
                *source_files
        ]

        # Generate the final Excel-Sheet
        workbook = xlsxwriter.Workbook(args.out_xlsx)
        worksheet1 = workbook.add_worksheet()
        worksheet1.freeze_panes(1, 0)

        # Write each row (first headers, then the data)
        for j in range(len(new_header)):
            worksheet1.write(0, j, new_header[j])
        for i in range(len(flatted_summarized)):
            for j in range(len(new_header)):
                worksheet1.write(i+1,j, flatted_summarized[i][j])

        # Heatmap Coloring
        worksheet1.conditional_format(1, 4, len(flatted_summarized)+1, len(new_header), {'type': '3_color_scale'})

        # Close the Excel-Sheet
        workbook.close()
