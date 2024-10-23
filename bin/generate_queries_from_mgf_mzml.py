#!/bin/env python

import argparse

import pyopenms

# Constants
HYDROGEN_MONO_MASS = 1.007825035
WATER_MASS = 18.0105647


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-input_mgf_mzml", help="The input-MGF, from which the queries for querying protein-graphs should be generated")
    parser.add_argument("-ppm", help="The same parameter as in search engines. The tolereance (in ppm) of the MS2-precursor (default: 5 ppm)", default=5, type=float)
    parser.add_argument("-max_limit_da", help="The highest MS2-Precursor (converted to Da) which should be allowed (default: 4500 Da).", default=4500, type=float)
    parser.add_argument("-out_csv", help="The output-query file (sorted from loweset to highest query). Each line is a query here.")
    return parser.parse_args()


def parse_mgf(in_file, ppm):
    '''Parse a MGF file and return all its MS2 precursors'''
    # Save all queries in entries
    entries = []

    # Iterate linewise
    in_entry = False
    retrieved_pep_mass = False
    retrieved_charge = False
    for line in in_file:

        # Set in_entry, since we are in a MS2-Spectrum-Entry
        if line.startswith("BEGIN IONS"):
            in_entry = True
            continue

        # If in_entry, extract the charge and peptide mass
        if in_entry:
            if line.startswith("PEPMASS="):
                pepmass = float(line[len("PEPMASS="):-1].split(" ")[0])
                retrieved_pep_mass = True
                continue
            if line.startswith("CHARGE="):
                charge = int(line[len("CHARGE="):-1].replace("+", ""))
                retrieved_charge = True
                continue

        # If we go out of an entry, peptide mass and charge is definitly retrieved (since it is required in an entry)
        if line.startswith("END IONS"):
            # However, we recheck if it was actually retrieved
            if retrieved_charge and retrieved_charge:
                # If yes: we save the query:
                # Convert peptide mass to Da
                da = (float(pepmass) * float(charge)) - (HYDROGEN_MONO_MASS * float(charge))
                da = da - WATER_MASS  # Subsract H2O-Mass, due to Protein-Graphs not encoding the water mass 

                # Get lower and upper limit
                lower = da - (da / 1000000) * ppm
                upper = da + (da / 1000000) * ppm

                # Append to entries
                entries.append((lower, upper))

            # Reset state
            in_entry = False
            retrieved_pep_mass = False
            retrieved_charge = False

    return entries


def parse_mzml(mzml_file, ppm):
    ''' Read all MS2 precursors from mzML '''
    # Load MZML
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(mzml_file, exp)

    entries = []
    for spectrum in exp.getSpectra():
        if spectrum.getMSLevel() == 2:
            precs = spectrum.getPrecursors()
            if len(precs) != 1:
                raise Exception("Unexpected number of precursors in a MS2-Spectrum")
            else:
                # Get mz and Charge
                pepmass = precs[0].getMZ()
                charge = precs[0].getCharge()

                # Convert peptide mass to Da
                da = (float(pepmass) * float(charge)) - (HYDROGEN_MONO_MASS * float(charge))
                da = da - WATER_MASS  # Subsract H2O-Mass, due to Protein-Graphs not encoding the water mass 

                # Get lower and upper limit
                lower = da - (da / 1000000) * ppm
                upper = da + (da / 1000000) * ppm

                # Append to entries
                entries.append((lower, upper))

    return entries


if __name__ == "__main__":
    # Parse arguments
    args = parse_args()

    # Open Spectra file and read the precursors
    if args.input_mgf_mzml.endswith(".mgf"):
        with open(args.input_mgf_mzml, "r") as in_file:
            entries = parse_mgf(in_file, args.ppm)
    elif args.input_mgf_mzml.endswith(".mzML"):
        entries = parse_mzml(args.input_mgf_mzml, args.ppm)
    else:
        raise Exception("File extension is not supported: {}".format(args.input_mgf_mzml.split(".")[-1]))

    # Output queries (and limit by Da)
    with open(args.out_csv, "w") as out_file:
        # Write only the sorted and limited queries
        for l, u in sorted(entries):
            if l < args.max_limit_da:
                out_file.write(str(l) + "," + str(u) + "\n")
