#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.sir_identified_files_glob = "$PWD/results/*qvalue_no_decoys_fdr_0.01.tsv"  // Input-blob of the identification tsvs, containing the followoing columns: "fasta_acc", "fasta_desc", "plain_peptide", "num", (needed to filter out only the top hits) and "source_file"
params.sir_outdir = "$PWD/results" // Output-Directory of the statistics
params.sir_export_data = "true"  // Boolean, if true, will export data into sir_outdir

// Optional Parameters
params.sir_post_fix = "" // Postfix for outputting into seperate folder (if executed multiple times)
params.sir_identification_from_protgraph = true // Flag, if the identification results are from ProtGraph-generated FASTAs. If true we use the "fasta_desc", if false we instead use the "fasta_acc".
params.sir_remove_variable_modifications = true // Flag if peptides with variable PTMs are removed (--> counted as unique). This is only considered for ProtGraph-Headers/Identification-Results (since the headers in the FASTA can distinguish between modified and not modified peptides).  
params.sir_count_same_protein_as_unique = true // Flag if proteins are counted as unique (in case of multiple sequences in the same protein). Since it can be only inferred to one protein, true should be set


// Standalone Workflow
workflow {
    // Collect all idents and prepare for workflow input
    identi_results = Channel.fromPath(params.sir_identified_files_glob)
    post_fix_idents_single = identi_results.combine(Channel.of(params.sir_post_fix))
    pf_idents_grouped = post_fix_idents_single.groupTuple(by: 1).map {tuple(it[1], it[0])}

    summarize_ident_results(
        pf_idents_grouped
    )
}

// Importable Workflow
workflow summarize_ident_results {
    take:
        // Takes identification results (as list, second tuple element) containing the followoing columns: 
        //"fasta_acc", "fasta_desc", "plain_peptide", "num" and "source_file"
        // NOTE: only add identification results of the same FDR-cutoff here (first tuple element)
        post_fix_identification_results
    main:
        // Split into seperate channels
        post_fix = post_fix_identification_results.map { it[0] }
        identification_results = post_fix_identification_results.map { it[1] }
        
        // Concatenate all identification results (already cutoff by an FDR) into a single file
        concat_to_single_identification_file(post_fix, identification_results)

        // Seperate into shared/unique/etc.. PSMs
        seperate_psms(post_fix, concat_to_single_identification_file.out)

        // Create HeatMaps of same found PSMs per MGF on unique/shared/etc..
        post_fix_single_ident = post_fix.combine(
            seperate_psms.out[1].concat(
                seperate_psms.out[2],
                seperate_psms.out[3],
                seperate_psms.out[4],
                concat_to_single_identification_file.out
            )
        )
        create_heatmaps(post_fix_single_ident)
}


process concat_to_single_identification_file {
    publishDir "${params.sir_outdir}${post_fix}", mode:'copy', enabled:"${params.sir_export_data}"

    input:
    val post_fix
    path tsv_files

    output: 
    path "all_identification_results.tsv"

    """
    { for filename in ${tsv_files}; do head -n 1 \$filename | while read line; do echo -e "\${line}"; done; break; done;
    for file in ${tsv_files}
    do
    tail -n +2 \$file | while read line; do echo -e "\${line}"; done
    done
    } > all_identification_results.tsv
    """  
}

process seperate_psms {
    publishDir "${params.sir_outdir}${post_fix}", mode:'copy', enabled:"${params.sir_export_data}"
    
    input:
    val post_fix
    path fdr_filtered_tsv 

    output:
    path("${fdr_filtered_tsv.baseName}_____all_psms_summary.tsv")
    path("${fdr_filtered_tsv.baseName}_____unique_psms.tsv")
    path("${fdr_filtered_tsv.baseName}_____shared_psms.tsv")
    path("${fdr_filtered_tsv.baseName}_____unique_ft_psms.tsv"), optional: true
    path("${fdr_filtered_tsv.baseName}_____shared_ft_psms.tsv"), optional: true

    """
    if [ ${params.sir_identification_from_protgraph} = true ]
    then
        PYTHONUNBUFFERED=1 count_and_seperate_psms.py \\
            -input_ident_tsv ${fdr_filtered_tsv} \\
            -out_summary ${fdr_filtered_tsv.baseName}_____all_psms_summary.tsv \\
            -out_unique ${fdr_filtered_tsv.baseName}_____unique_psms.tsv \\
            -out_shared ${fdr_filtered_tsv.baseName}_____shared_psms.tsv \\
            -out_ft_unique ${fdr_filtered_tsv.baseName}_____unique_ft_psms.tsv \\
            -out_ft_shared ${fdr_filtered_tsv.baseName}_____shared_ft_psms.tsv \\
            -accession_header "fasta_desc" \\
            -same_protein_as_unique ${params.sir_count_same_protein_as_unique} \\
            -remove_varmods ${params.sir_remove_variable_modifications}
    else
        PYTHONUNBUFFERED=1 count_and_seperate_psms.py \\
            -input_ident_tsv ${fdr_filtered_tsv} \\
            -out_summary ${fdr_filtered_tsv.baseName}_____all_psms_summary.tsv \\
            -out_unique ${fdr_filtered_tsv.baseName}_____unique_psms.tsv \\
            -out_shared ${fdr_filtered_tsv.baseName}_____shared_psms.tsv \\
            -out_ft_unique ${fdr_filtered_tsv.baseName}_____unique_ft_psms.tsv \\
            -out_ft_shared ${fdr_filtered_tsv.baseName}_____shared_ft_psms.tsv \\
            -accession_header "fasta_acc" \\
            -same_protein_as_unique ${params.sir_count_same_protein_as_unique} \\
            -remove_varmods ${params.sir_remove_variable_modifications}
    fi
    """
}

process create_heatmaps {
    publishDir "${params.sir_outdir}${post_fix}", mode:'copy', enabled:"${params.sir_export_data}"

    input:
    tuple val(post_fix), path(tsv_file)

    output: 
    path "${tsv_file.baseName}_____heatmap.xlsx"

    """
    PYTHONUNBUFFERED=1 export_psms_to_heatmap.py \\
        -in_psm_tsv ${tsv_file} \\
        -out_xlsx ${tsv_file.baseName}_____heatmap.xlsx
    """  
}
