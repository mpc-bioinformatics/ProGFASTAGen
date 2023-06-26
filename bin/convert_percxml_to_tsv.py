#!/bin/env python

import sys
import csv
import argparse
from collections import defaultdict

from lxml import etree


csv.field_size_limit(sys.maxsize)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-pout_xml", help="Percolator XML-output file")
    parser.add_argument("-out_tsv", help="The converted output_file in tsv format")
    return parser.parse_args()


if __name__ == "__main__":
    # Parse arguments
    args = parse_args()

    psms = []  # A list of all retrieved information per psm
    # Iterate over the XML (tag-wise)
    for event, elem in etree.iterparse(args.pout_xml, events=("end",), recover=True):
        qname = etree.QName(elem)
        elem.tag = qname.localname
        namespace = qname.namespace

        if elem.tag == "psm":
            # We have a PSM. Initialize the default dict where we collect all information
            info_dict = defaultdict(lambda: list())
            info_dict["psm_id"].append(elem.attrib["{{{}}}psm_id".format(elem.nsmap["p"])])

            # Collect all the information from the subelements
            for x in elem:
                x.tag = etree.QName(x).localname

                if x.tag in ("svm_score", "q_value", "pep", "exp_mass", "calc_mass", "protein_id", "p_value" ):
                    # Get selected information from the text attribute
                    info_dict[x.tag].append(x.text)

                elif x.tag == "peptide_seq":
                    # Get peptide and remove any modifications from the sequence
                    peptide = x.attrib["seq"]
                    info_dict["peptide_seq"].append(peptide)
                    while "[" in peptide:
                        peptide = peptide[0:peptide.index("[")] + peptide[peptide.index("]")+1:]
                    info_dict["plain_peptide"].append(peptide)

            # Append to psms
            psms.append(info_dict)

    # Get header line
    header = list(psms[0].keys())

    # Get each line of a psm as a list
    psms_lists = []
    for p in psms:
        new_entry = []
        for h in header:
            new_entry.append(",".join(p[h]))
        psms_lists.append(new_entry)

    # Write to tsv
    with open(args.out_tsv, "w") as out_file:
        csv_out = csv.writer(out_file, delimiter="\t")
        csv_out.writerows([header, *psms_lists])
