#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.idc_mgf_folder = "$PWD/MGFs"  // Input-Directory of MGF-FIles, which should be used for identification
params.idc_fasta_file = "peptides.fasta"  // Database (FASTA-file) used for identification. Can be with decoys, they should then be prefixed via "DECOY_". E.g. "">sp|XXXXXX|XXXXXX" should be ">DECOY_sp|XXXXXX|XXXXXX"
params.idc_search_parameter_file = "$PWD/example_configurations/comet_config_high_res.txt"  // Search Parameters for Comet. Defaults can be found in "example_configurations"
params.idc_outdir = "$PWD/results"  // Output-Directory of the Identification-Results (mulitple output files per input file, extended with fdr_filtered)
params.idc_export_data = "true"  // Boolean, if true, will export data into idc_outdir


// Optional Parameters for Comet, Usage of Percolator and FDR
params.idc_tda = 1  // Comet-Parameter: 0 -->  Do not append database with decoys. (Please then provide a FASTA with DECOY_ prefixed to indicate decoy entries in the FASTA) | 1 --> Append decoys to database (if set to 1 it generates the reversed sequence as decoys using "DECOY_" as the prefix)
params.idc_num_parallel_threads_per_search = 4   // Number of threads used in Comet per search. From experience, a value between 4 to 16 is usually enough.
params.idc_use_percolator = 1  // 0 --> Do not use, 1 --> Use percolator for rescoring the identification results
params.idc_use_n_hits = 1  // Number of hits per spectrum to be used in FDR-Calulcation (CAUTION: can skew the results, if set higher)
params.idc_fdr = "0.01" // FDR/qvalue used for cutoff of the identification results. You can specify multiple ones using |. E.G.: "0.01|0.05"


// Standalone Workflow
workflow {
    // Get all MGF-files which should be identified
    mgfs = Channel.fromPath(params.idc_mgf_folder + "/*.mgf")

    // Get FASTA-file
    fasta_file = Channel.fromPath(params.idc_fasta_file)

    // Get Parameters-file for Comet (Search-Parameters)
    parameters_file = Channel.fromPath(params.idc_search_parameter_file)

    // Execute Workflow
    identification_via_comet(mgfs, fasta_file, parameters_file)
}

// Importable Workflow
workflow identification_via_comet {
    take:
        mgfs
        fasta_file
        parameters_file
    main:
        // Combined channel replicated the indexed fasta for each MZML to be reused
        combined_channel = fasta_file
            .combine(parameters_file)
            .combine(mgfs)
            
        // Start search
        comet_search_mgf(combined_channel)

        // Split multiple FDRs (if set) and generate a channel from it
        single_fdrs = String.valueOf(params.idc_fdr).split("\\|").collect{it as String}
        single_fdrs_channel = Channel.from(single_fdrs)

        // Remove the mzid, since we either work with the pin or txt-file in this workflow
        txt_and_pin = comet_search_mgf.out.map { mzid, txt, pin -> tuple(txt, pin) }
        
        if (params.idc_use_percolator == 1) {
            // Case: Use Percolator to rescore the identification results
            search_results_with_fdr = txt_and_pin.combine(single_fdrs_channel)
            execute_percolator(search_results_with_fdr)
            combined_txt_channel = execute_percolator.out
                .combine(fasta_file)
                .map { tuple(it[0], it[1], it[3], it[2]) }
        } else {
            // Case: Do not use Percolator and just continue with the  identification results
            combined_txt_channel = txt_and_pin
                .combine(fasta_file)
                .combine(single_fdrs)
                .map { tuple(it[0], "/NO_PERC", it[2], it[3]) }
        }

        // Cutoff idntnification table
        cutoff_identification_results(combined_txt_channel)
    emit:
        // Output fdr [0] and table (in tsv) [1]
        // The table is guaranteed to have the columns:
        // "plain_peptide" "source_file" "qvalue" "fasta_acc" "fasta_desc" "used_score" "charge" "retention_time" "exp_mass_to_charge"
        // (and additional columns provided by either percolator or comet)
        cutoff_identification_results.out[0]
}


process comet_search_mgf {
    cpus params.idc_num_parallel_threads_per_search
    publishDir "${params.idc_outdir}", mode:'copy', enabled:"${params.idc_export_data}"

    input:
    tuple path(input_fasta), path(mod_file), path(mgf_file)

    output: 
    tuple path("${mgf_file.baseName}.mzid"), path("${mgf_file.baseName}.txt"), path("${mgf_file.baseName}.pin")

    """
    # Replace specific settings regulated by this workflow (or to ensure output files)
    sed 's/^decoy_search.*/decoy_search = ${params.idc_tda} /' ${mod_file} > ${mod_file.baseName}_new.txt
    sed -i 's/^output_mzidentmlfile.*/output_mzidentmlfile = 1/' ${mod_file.baseName}_new.txt
    sed -i 's/^output_txtfile.*/output_txtfile = 1/' ${mod_file.baseName}_new.txt
    sed -i 's/^output_percolatorfile.*/output_percolatorfile = 1/' ${mod_file.baseName}_new.txt
    sed -i 's/^decoy_prefix.*/decoy_prefix = DECOY_/' ${mod_file.baseName}_new.txt
    sed -i 's/^num_threads.*/num_threads = ${params.idc_num_parallel_threads_per_search} /' ${mod_file.baseName}_new.txt

    # Execute Comet
    comet.linux.exe -P${mod_file.baseName}_new.txt -D${input_fasta} ${mgf_file}
    """  
}

process execute_percolator {
    publishDir "${params.idc_outdir}", mode:'copy', enabled:"${params.idc_export_data}"

    input:
    tuple path(txt), path(comet_pin), val(fdr)

    output: 
    tuple path(txt), path("${comet_pin.baseName}_____perc_fdr_${fdr}.tsv"), val(fdr)

    """
    # Escape Pin-File so that the generated XML by Percolator is valid
    sed \\
        -e 's/</\\&#60;/g' -e 's/>/\\&#62;/g' -e 's/\\&/\\&#38;/g' -e "s/'/\\&#39;/g" -e 's/"/\\&#34;/g' \\
        ${comet_pin} > ${comet_pin}_escaped


    # Execute Percolator
    \$(get_cur_bin_dir.sh)/percolator/usr/bin/percolator \\
        --trainFDR ${fdr} --testFDR ${fdr} -P DECOY_ \\
        --only-psms \\
        --decoy-xml-output \\
        -X ${comet_pin.baseName}.xml \\
        ${comet_pin}_escaped

    # Convert output into a tsv
     PYTHONUNBUFFERED=1 convert_percxml_to_tsv.py -pout_xml ${comet_pin.baseName}.xml -out_tsv ${comet_pin.baseName}_____perc_fdr_${fdr}.tsv
    """  
}

process cutoff_identification_results {
    publishDir "${params.idc_outdir}", mode:'copy', enabled:"${params.idc_export_data}"

    input:
    tuple path(comet_txt), path(perc_tsv), path(fasta), val(fdr)

    output:
    tuple val(fdr), path("${comet_txt.baseName}*____qvalue_no_decoys_fdr_${fdr}.tsv")
    path("${comet_txt.baseName}_____qvalue.tsv"), optional: true
    path("${perc_tsv.baseName}_____qvalue.tsv"), optional: true

    """
    if [ "${perc_tsv}" == "NO_PERC" ]; then
        # No Percolator was used
         PYTHONUNBUFFERED=1 apply_fdr_comet.py \\
            -comet_txt ${comet_txt} \\
            -use_n_hits ${params.idc_use_n_hits} \\
            -fasta ${fasta} \\
            -decoy_string "DECOY_" \\
            -fdr ${fdr} \\
            -out_all_tsv ${comet_txt.baseName}_____qvalue.tsv \\
            -out_no_decoys_tsv ${comet_txt.baseName}_____qvalue_no_decoys.tsv \\
            -out_fdr_cutoff_tsv ${comet_txt.baseName}_____qvalue_no_decoys_fdr_${fdr}.tsv
    else
        PYTHONUNBUFFERED=1 apply_fdr_comet.py \\
            -comet_txt ${comet_txt} \\
            -perc_tsv ${perc_tsv} \\
            -use_n_hits ${params.idc_use_n_hits} \\
            -fasta ${fasta} \\
            -decoy_string "DECOY_" \\
            -fdr ${fdr} \\
            -out_all_tsv ${perc_tsv.baseName}_____qvalue.tsv \\
            -out_no_decoys_tsv ${perc_tsv.baseName}_____qvalue_no_decoys.tsv \\
            -out_fdr_cutoff_tsv ${perc_tsv.baseName}_____qvalue_no_decoys_fdr_${fdr}.tsv
    fi
    """
}
