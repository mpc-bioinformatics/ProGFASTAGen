#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.cmf_mgf_files = "$PWD/mgfs"  // Input-Directory of MGF-Files, which should be used to generate FASTA-files
params.cmf_sp_embl_file = "proteins.txt"  // Database-file in SP-EMBL-Format. E.G.: This can be retrieved from UniProt via the txt-export
params.cmf_outdir = "$PWD/results"  // Output-Directory of the generated FASTA-files and other intermediate files
params.cmf_export_data = "true"  // Boolean, if true, will export data into cmf_outdir

// Optional Parameters
params.cmf_one_fasta = true  // Flag to generate a fasta for the whole dataset OR a fasta per MGF-file (true --> only a single FASTA should be generated )

// For Graph-Generation
params.cmf_max_precursor_da = 5000  // Upper limit of a query in Da. MS2-Precursors higher then this value will be ommitted
params.cmf_query_ppm = 5  // The ppm (tolereance) which should be used to generate the queries. This should be the same as set in the search engine.
params.cmf_pg_additional_params = "-ft ALL -fm 'C:57.021464' -vm 'M:15.994915'"  // Additional Parameters which ProtGraph should consider. Add here the features (-ft) which should be included, the fixed and variable modifications (as specified in the search engine via "-fm", "-vm") and also include (if needed) the digestion (-d) or other things. See ProtGraph for a whole list of possible parameters
params.cmf_elbpcsr = 32 // Number of PDB-Intervals per Node in a Protein-Graph. This can be arbitrary.

// For Traversal and Limits
params.cmf_num_procs_traversal = Runtime.runtime.availableProcessors().div(2)  // Number of processes used for the graph traversal to generate a FASTA. Set this number to maximum of MaxProcess/2, due to the binary-search and to prevent DeadLocks which infrequently occur, if set higher
params.cmf_number_of_bins = 128 // Set the number of bins (they will be max_precursor/number_of_bins large)
params.cmf_timeout_for_single_query = 5 // Timeout for s single query in a bin in protein graphs (in seconds). For ~10k queries with a timeout of 5s this could take ~14 hours. This is not an upper limit and also may take loknter on how long the FASTA-Generation will need. It is also influenced by the number of applied variable modifications and applied features in a protein-graph (in ProtGraph, -ft and -vm).
params.cmf_proteins_to_limit = "__all__" // Limit the variants for the specified proteins (either "__all__" to limit all, iff they take longer then the timeout, "__none__" for none  or "PXXXXX,PXXXXX,PXXXXX" as a comma list to limit only specific ones. E.G.: In a human database, it could be interesting to only limit "P04637,P68871")
params.cmf_maximum_variant_limit = 5 // Maximum limit of variants applied on a Protein-Graph on a bin. E.G. if we found in the binary search that P53 has the following limits: -1,-1,3,1,1,1, setting this vallue would give the follwoing limits 5,5,3,1,1,1. This paramter could be used to set an upper limit of variants in a peptide. Set to -1 to allow infinite many. Set lower to reduce the size of the final FASTA-file. A limit of 5 seems reasonable.
params.cmf_use_floats = 0  // Bool wheather to use floats or integers for the masses of aminoacids (1 --> use floats, 0 --> use integers). Depending on the architeture the one or the other could be faster. Defaults to use integers.


// Standalone Workflow
workflow {
    // Get all MGF-files
    mgfs = Channel.fromPath(params.cmf_mgf_files + "/*.mgf")

    // Get SP-EMBL file
    sp_embl_file = Channel.fromPath(params.cmf_sp_embl_file)

    create_precursor_specific_fasta(mgfs, sp_embl_file)
}

// Importable Workflow
workflow create_precursor_specific_fasta {
    take:
        // Takes MGF-files
        mgf_files
        sp_embl_file
    main:
        // Convert the file to MGF
        generate_query_csvs(mgf_files)

        // Decide whether we want a single FASTA or a FASTA per file
        if (params.cmf_one_fasta) {
            // Single FASTA
            concat_query_csvs(generate_query_csvs.out.collect())
            queries = optimize_query(concat_query_csvs.out)
        } else {
            // FASTA per MGF
            queries = optimize_query(generate_query_csvs.out)
        }

        // Create Protein-Graphs using ProtGraph
        create_protein_graphs(sp_embl_file)

        // Create timout limits per Protein-Graph (and bin). 
        // This process is executed on its own without any other process, to estimate the limits without the influence of other processes
        determine_limits_using_binary_search(create_protein_graphs.out[0])


        // Combine the queries with the Protein-Graphs
        pgs_limits_and_query = create_protein_graphs.out[0].combine(determine_limits_using_binary_search.out[0]).combine(queries)

        // Create the ms2-specific FASTA.
        // This process is executed on its own without any other process, to reduce the possibility of DeadLocks
        create_precursor_specific_fasta_via_protgraphcpp(pgs_limits_and_query)

        // Merge duplicated entries into a single entry to ensure that the FASTA is "Sequence-Unique" (--> IOW: Each sequence only occurs once in the FASTA)
        compact_fasta(create_precursor_specific_fasta_via_protgraphcpp.out)
    emit:
        // Retruns each MGF, converted from a RAW-file
        compact_fasta.out
}


process generate_query_csvs {
    input:
    path input_mgf

    output:
    path "${input_mgf.baseName}.csv"

    """
    PYTHONUNBUFFERED=1 generate_queries_from_mgf.py \\
        -input_mgf ${input_mgf} \\
        -out_csv ${input_mgf.baseName}.csv \\
        -max_limit_da ${params.cmf_max_precursor_da} \\
        -ppm ${params.cmf_query_ppm}
    """
}

process concat_query_csvs {
    input:
    path input_mgf_csv

    output:
    path "all_mgfs.csv"

    """
    cat ${input_mgf_csv} > all_mgfs.csv
    """
}

process optimize_query {
    input:
    path input_csv

    output:
    path "${input_csv.baseName}_optimized.csv"

    """
    PYTHONUNBUFFERED=1 optimize_query_csv.py \\
        -in_query_csv ${input_csv} \\
        -out_query_csv ${input_csv.baseName}_optimized.csv
    """
}

process create_protein_graphs {
    publishDir "${params.cmf_outdir}/", mode:'copy', enabled:"${params.cmf_export_data}"

    input:
    path input_sp_embl

    output:
    path "database.bpcsr"
    path "proteins_statistics.csv"

    """
    PYTHONUNBUFFERED=1 protgraph \\
        -elbpcsr -eo . \\
        -elbpcsr_pdbs ${params.cmf_elbpcsr} \\
        ${params.cmf_pg_additional_params} \\
        -cnp -cnpvar -amw \\
        -o proteins_statistics.csv \\
        ${input_sp_embl}
    """
}

process determine_limits_using_binary_search {
    publishDir "${params.cmf_outdir}/", mode:'copy', enabled:"${params.cmf_export_data}"
    cpus Runtime.runtime.availableProcessors() // Tell Nextflow, that it uses all processors, to ensure that this step is not distrubed by other processes

    input:
    path database

    output:
    path "traversal_limits_cpp.csv"
    path "traversal_limits_detailed.csv"

    """
    if [ ${params.cmf_use_floats} -eq 1 ]
    then
        PYTHONUNBUFFERED=1 binary_search_on_protein_graphs.py \\
            -pgcpp_exe \$(get_cur_bin_dir.sh)/ProtGraphTraverseFloatSourceDryRun/build/protgraphtraversefloatdryrun \\
            -num_processes ${params.cmf_num_procs_traversal} \\
            -protein_graphs_bpcsr ${database} \\
            -out_detailed_statistics traversal_limits_detailed.csv \\
            -out_limits traversal_limits_cpp.csv  \\
            -max_precursor_da ${params.cmf_max_precursor_da} \\
            -bins ${params.cmf_number_of_bins} \\
            -ppm ${params.cmf_query_ppm} \\
            -timeout ${params.cmf_timeout_for_single_query} \\
            -proteins_to_limit ${params.cmf_proteins_to_limit} \\
            -max_limit ${params.cmf_maximum_variant_limit} \\
            -apply_smooting_method median
    else
        PYTHONUNBUFFERED=1 binary_search_on_protein_graphs.py \\
            -pgcpp_exe \$(get_cur_bin_dir.sh)/ProtGraphTraverseIntSourceDryRun/build/protgraphtraverseintdryrun \\
            -num_processes ${params.cmf_num_procs_traversal} \\
            -protein_graphs_bpcsr ${database} \\
            -out_detailed_statistics traversal_limits_detailed.csv \\
            -out_limits traversal_limits_cpp.csv  \\
            -max_precursor_da ${params.cmf_max_precursor_da} \\
            -bins ${params.cmf_number_of_bins} \\
            -ppm ${params.cmf_query_ppm} \\
            -timeout ${params.cmf_timeout_for_single_query} \\
            -proteins_to_limit ${params.cmf_proteins_to_limit} \\
            -max_limit ${params.cmf_maximum_variant_limit} \\
            -apply_smooting_method median
    fi
    """
}

process create_precursor_specific_fasta_via_protgraphcpp {
    cpus Runtime.runtime.availableProcessors() // Tell Nextflow, that it uses all processors, to ensure that this step is not distrubed by other processes

    input:
    tuple path(database), path(traversal_limits), path(csv_query)

    output:
    path "${csv_query.baseName}.fasta"

    """
    if [ ${params.cmf_use_floats} -eq 1 ]
    then
        \$(get_cur_bin_dir.sh)/ProtGraphTraverseFloatSourceVarLimitter/build/protgraphtraversefloatvarlimitter \\
            ${database} ${csv_query} ${params.cmf_num_procs_traversal} ${csv_query.baseName}.fasta ${traversal_limits}
        # TODO finish implementation
    else
        \$(get_cur_bin_dir.sh)/ProtGraphTraverseIntSourceVarLimitter/build/protgraphtraverseintvarlimitter \\
            ${database} ${csv_query} ${params.cmf_num_procs_traversal} ${csv_query.baseName}.fasta ${traversal_limits}
    fi
    """
}

process compact_fasta {
    publishDir "${params.cmf_outdir}/", mode:'copy', enabled:"${params.cmf_export_data}"

    input:
    path input_fasta

    output:
    path "${input_fasta.baseName}_final.fasta"

    """
    PYTHONUNBUFFERED=1 protgraph_compact_fasta ${input_fasta} -o ${input_fasta.baseName}_final.fasta
    """
}
