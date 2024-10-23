#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Required Parameters
params.ctm_raws = "$PWD/raws"  // Directory containing .raw/.d-files
params.ctm_outdir = "$PWD/results"  // Output-Directory where the MGFs will be stored
params.ctm_export_data = "true"  // Boolean, if true, will export data into ctm_outdir

// Optional Parameters
params.ctm_additional_params_thermo = "" // Additional paramaters which may be set for conversion
params.ctm_additional_params_tdf2mzml = "" // Additional paramaters which may be set for conversion
params.ctm_num_forks_conversion = Runtime.runtime.availableProcessors()  // The maximum number of parallel conversions done at once.


// Standalone Workflow
workflow {
    // Convert to .mzML
    raw_files = Channel.fromPath(params.ctm_raws + "/*.raw")
    d_files = Channel.fromPath(params.ctm_raws + "/*.d", type: "dir")
    convert_to_mzml(raw_files, d_files)
}

// Importable Workflow
workflow convert_to_mzml {
    take:
        // Takes .raw and .d
        raw_files 
        d_files
    main:
        // Convert the file to .mzML
        convert_raw_via_thermorawfileparser(raw_files)
        convert_d_via_tdf2mzml(d_files)
        results = convert_raw_via_thermorawfileparser.out.concat(
            convert_d_via_tdf2mzml.out
        )
    emit:
        // Return each input as .mzML, converted from a .raw/.d-file
        results
}


process convert_raw_via_thermorawfileparser {
    maxForks params.ctm_num_forks_conversion
    container 'quay.io/biocontainers/thermorawfileparser:1.4.3--ha8f3691_0'
    cpus 2 // Currently limited to two conversions at once (due to TRFP: https://github.com/compomics/ThermoRawFileParser/issues/23 )


    publishDir "${params.ctm_outdir}/", mode:'copy', enabled:"${params.ctm_export_data}"

    input:
    path raw

    output:
    path "${raw.baseName}.mzML"

    """
    thermorawfileparser ${params.ctm_additional_params_thermo} --format=2 --output_file=${raw.baseName}.mzML --input=${raw} 
    """
}


process convert_d_via_tdf2mzml {
    maxForks params.ctm_num_forks_conversion
    container "mfreitas/tdf2mzml"
    cpus 2 // hardcoded limit

    publishDir "${params.ctm_outdir}/", mode:'copy', enabled:"${params.ctm_export_data}"

    input:
    path raw

    output:
    path "${raw.baseName}.mzML"

    """    
    export MKL_NUM_THREADS=2
    export NUMEXPR_NUM_THREADS=2
    export OMP_NUM_THREADS=2
    tdf2mzml.py -i ${raw} -o ${raw.baseName}.mzML --ms1_type centroid
    """
}


        