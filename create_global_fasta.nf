#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.cgf_sp_embl_file = "proteins.txt"  // Database-file in SP-EMBL-Format. E.G.: This can be retrieved from UniProt via the txt-export
params.cgf_outdir = "$PWD/results"  // Output-Directory of the FASTA-results
params.cgf_export_data = "true"  // Boolean, if true, will export data into cgf_outdir

// Optional Parameters for Protein-Graph Generation
params.cgf_features_in_graphs = "-ft None"  // Features added to the Protein-Graph (See ProtGraph for more info)
params.cgf_peptide_limits = "--pep_miscleavages 2 --pep_min_pep_length 5 --pep_max_weight 5000" // Limits used for Exporting peptides from Protein-Graphs (see ProtGraph for more info)


// Standalone Workflow
workflow {
    // get the SP-EMBL-File
    sp_embl_file = Channel.fromPath(params.cgf_sp_embl_file)
    create_global_fasta(sp_embl_file)
}

// Importable Workflow
workflow create_global_fasta {
    take:
        // Takes a single SP-EMBL-File
        sp_embl_file
    main:
        // Create Protein-Graphs, statistics and the sqlite database
        create_sqlite_fasta_database(sp_embl_file)
        // Export sqlite into fasta
        convert_pepsqlite_to_fasta(create_sqlite_fasta_database.out[0])
    emit: 
        // Returns the corresponding Global-Peptides-FASTA
        convert_pepsqlite_to_fasta.out
}


process create_sqlite_fasta_database {
    publishDir "${params.cgf_outdir}/", mode:'copy', enabled:"${params.cgf_export_data}"

    input:
    path input_sp_embl

    output:
    path "peptides.db"
    path "proteins_statistics.csv"

    """
    PYTHONUNBUFFERED=1 protgraph -epepsqlite --pep_sqlite_database peptides.db -eo . ${params.cgf_features_in_graphs} ${params.cgf_peptide_limits} -cnp -cnpm -amw -o proteins_statistics.csv ${input_sp_embl}
    """
}

process convert_pepsqlite_to_fasta {
    publishDir "${params.cgf_outdir}/", mode:'copy', enabled:"${params.cgf_export_data}"

    input:
    path database

    output:
    path "global_peptides.fasta"

    """
     PYTHONUNBUFFERED=1 protgraph_pepsqlite_to_fasta ${database} -o global_peptides.fasta
    """
}
