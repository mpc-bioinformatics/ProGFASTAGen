#!/bin/env python

import argparse
import csv
import os
import sys


csv.field_size_limit(sys.maxsize)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-in_query_csv", help="Input-CSV of queries which should be optimized (merged together)")
    parser.add_argument("-out_query_csv", help="Output-CSV of merges queries. These can be higher then X ppm, due to overlapping. They can be queried in Protein-Graphs at once.")
    return parser.parse_args()


if __name__ == "__main__":
    # Parse arguments
    args = parse_args()

    # Open in- and output file
    with open(args.out_query_csv, "w") as out_file, open(args.in_query_csv, "r") as in_file:
        # Initialize CSV writer
        csv_out = csv.writer(out_file)
        csv_in = csv.reader(in_file)

        # Read all queries in input
        queries = [[float(y) for y in x] for x in csv_in]
        
        # Sort them ascending
        queries = sorted(queries)
        
        # Find overlaps and combine them
        num_queries = len(queries)
        index = 0
        while index < len(queries) - 1:
            if queries[index][1] > queries[index+1][0]:  # Overlap found!
                queries[index][1] = max(queries[index][1], queries[index+1][1])
                del queries[index+1]
            else:
                index += 1

        print("#Queries reduced to: {}%".format(len(queries)*100/num_queries))

        # Save optimized queries in output
        for q in queries:
            csv_out.writerow(q)
