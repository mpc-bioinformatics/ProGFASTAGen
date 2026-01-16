#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters (defaults are provided as an examples)
params.ctm_raws = "$PWD/raws"  // Directory containing .raw/.d-files
params.ctm_outdir = "$PWD/results"  // Output-Directory where the MGFs will be stored


// Optional Parameters
params.ctm_additional_params_thermo = "" // Additional paramaters which may be set for conversion
params.ctm_additional_params_alphatims = "" // Additional paramaters which may be set for conversion
params.ctm_num_forks_conversion = Runtime.runtime.availableProcessors()  // The maximum number of parallel conversions done at once.
params.ctm_export_data = "true"  // Boolean, if true, will export data into ctm_outdir


// Standalone Workflow
workflow {
    // Convert to .mgf
    raw_files = Channel.fromPath(params.ctm_raws + "/*.raw")
    d_files = Channel.fromPath(params.ctm_raws + "/*.d", type: "dir")
    convert_to_mgf(raw_files, d_files)
}

// Importable Workflow
workflow convert_to_mgf {
    take:
        // Takes .raw and .d
        raw_files 
        d_files 
    main:
        // Convert the files to .mgf
        convert_raw_via_thermorawfileparser(raw_files)
        convert_d_via_alphatims(d_files)
        results = convert_raw_via_thermorawfileparser.out.concat(
            convert_d_via_alphatims.out
        )
    emit:
        // Return each input as .mgf, converted from a .raw/.d-file
        results
}


process convert_raw_via_thermorawfileparser {
    maxForks params.ctm_num_forks_conversion
    cpus 2 // Currently limited to two conversions at once (due to TRFP: https://github.com/compomics/ThermoRawFileParser/issues/23 )
    label "progfastagen_thermo_conversion"

    publishDir "${params.ctm_outdir}/", mode:'copy', enabled:"${params.ctm_export_data}"

    input:
    path raw

    output:
    path "${raw.baseName}.mgf"

    """
    thermorawfileparser ${params.ctm_additional_params_thermo} --format=0 --output_file=${raw.baseName}.mgf --input=${raw} 
    """
}


process convert_d_via_alphatims {
    maxForks params.ctm_num_forks_conversion
    cpus 2 // hardcoded limit
    label "progfastagen"

    publishDir "${params.ctm_outdir}/", mode:'copy', enabled:"${params.ctm_export_data}"

    input:
    path raw

    output:
    path "${raw.baseName}.mgf"

    """
    alphatims export mgf --threads 2 ${params.ctm_additional_params_alphatims} -o . ${raw}
    """
}
