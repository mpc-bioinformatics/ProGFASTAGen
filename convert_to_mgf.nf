#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.ctm_raws = "$PWD/raws"  // Directory containing RAW-files
params.ctm_outdir = "$PWD/results"  // Output-Directory where the MGFs will be stored
params.ctm_export_data = "true"  // Boolean, if true, will export data into ctm_outdir

// Optional Parameters
params.ctm_additional_params = ""
params.ctm_num_procs_conversion = Runtime.runtime.availableProcessors()  // Number of process used to convert (CAUTION: This can be very resource intensive!)


// Standalone Workflow
workflow {
    // Convert t to MGF
    rawfiles = Channel.fromPath(params.ctm_raws + "/*.raw")
    convert_raw_via_thermorawfileparser(rawfiles)
}

// Importable Workflow
workflow convert_to_mgf {
    take:
        // Takes a list of RAW-files
        raw_files 
    main:
        // Convert the file to MGF
        convert_raw_via_thermorawfileparser(raw_files)
    emit:
        // Retruns each MGF, converted from a RAW-file
        convert_raw_via_thermorawfileparser.out
}


process convert_raw_via_thermorawfileparser {
    maxForks params.ctm_num_procs_conversion
    stageInMode "copy"

    publishDir "${params.ctm_outdir}/", mode:'copy', enabled:"${params.ctm_export_data}"

    input:
    path raw

    output:
    path "${raw.baseName}.mgf"

    """
    mono \$(get_cur_bin_dir.sh)/ThermoRawFileParser/ThermoRawFileParser.exe ${params.ctm_additional_params} --format=0 --output_file=${raw.baseName}.mgf --input=${raw} 
    """
}
