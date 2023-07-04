# ProGFASTAGen

The ProGFASTAGen (**Pro**tein-**G**raph-**FASTA**-**Gen**erator or **Pro**t**G**raph-**FASTA**-**Gen**erator) repository contains workflows to generate so-called ms2-specific-FASTAs (using the precursors from MGF-files) including feature-peptides, like VARIANTs or CONFLICTs if desired, or global-FASTAs (as described in [ProtGraph](https://github.com/mpc-bioinformatics/ProtGraph)). The single workflow scripts have been implemented with [Nextflow-DSL-2](https://www.nextflow.io/docs/latest/dsl2.html) and are independent to each other. Each of them can be used on their own or can be imported to other workflows for other use-cases. Further, we included three main-workflows, to show how the single workflows can be chained together. The `main_workflow_protein_fasta.nf`-workflow converts Thermo-RAW-files into MGF, searches with Comet (and Percolator) and the identification results are then further summarized. The workflows `main_workflow_global_fasta.nf` and `main_workflow_ms2_specific_fasta.nf` generate specific FASTA-files before search-engine-identification. Below are example nextflow-calls, which can be used.

Regarding the ms2-specific-FASTA-generation: The source-code of the C++ implementation for traversal can be found in `bin`. There, four implementations are present: `Float/Int`-Versions as well as `DryRun/VarLimitter`-Versions of the traversal. The `Float/Int`-Versions can be faster/slower on some specific processors and both versions can be used via a flag in the `create_ms2_specific_fasta.nf`-workflow. The `DryRun`-Version does not generate a FASTA but tests the used system (depending on a query-timeout) to determine the maximum number of variants which can be used, while not timing out. The actual FASTA-generation happens in the `VarLimitter`-Version using the generated protein-graphs at hand.

in **Prerequisites** a small description of dependencies and how to set up the host system is given. **Individual steps** describes the single workflows and how they can be called, while **Main Workflow Scripts** shows example-calls of the main workflows. In **Regenerate Results from Publication**, the calls and parameters are shown, which were used in the publication. Using the same FASTA, SP-EMBL with a similar server-setting should yield similar results as used in the publication.

## Prerequisites

### Executing on Linux

This workflow can be only executed on linux (tested on Ubuntu 22.04 and ArchLinux). Before setting up the `bin`-folder, some requiered binaries need to be present on the OS. (Focusing on Ubuntu:) The following packages need to be installed on Ubuntu (via `apt`), if not already:

```text
build-essential
wget
curl
unzip
cmake
mono-complete
python3-pip (or any environment with Python3, where pip is available)
python-is-python3 (needed for ubuntu, so that python points to python3)
```

If all packages are installed (and the python environment is set up), the setup-script needs to be executed, which downloads needed dependencies and compiles the source-code located in the `bin`-folder:

```shell
chmod +x compile_and_setup_depencies.sh  # In case this file is not executable
./compile_and_setup_depencies.sh  # Downloads dependencies, compiles the C++-implementation and sets all binaries in the bin-folder as executable
```

If the script exits without errors, the provided workflows can be executed with the command `nextflow`.

### Executing in Docker

Alternatively, docker can be used. For this, please follow the [installation guide](https://docs.docker.com/engine/install/ubuntu/) for docker. After installing docker, a local docker-container can be build with all needed dependencies for the workflows. We provide a `Dockerfile` in the `docker`-folder. To build it, execute (while beeing with a shell in the root-folder of this repository) the following:

```shell
docker build -t progfastagen:local . -f docker/Dockerfile
```

This command builds a local docker container, tagging it with `progfastagen:local`, which can be later used by nextflow. To use it with nextflow, make sure that `nextflow` is installed on the host-system. For each of the workflow example calls below, the `-with-docker progfastagen:local` then needs to be appended, to let `nextflow` know to use the local docker-container.

## Individual Steps

Each step has been implemented in such a way, that it can be executed on its own. Each subsection below, provides a brief overview and an example call of the required parameters to demonstrate how the workflow can be called. If you are interested for all the available parameters within a workflow and want modify or tune them, then please refer to the source of the workflows, where each parameter is described briefly.

### Converting RAW-files to MGF

The workflow `convert_to_mgf.nf` is a wrapper around the ThermoRawFileparser and converts RAW-files to the MGF-format. The `ctm_raws` parameter needs to be set, in order to generate the MGF-files:

```text
nextflow run convert_to_mgf.nf \
    --ctm_raws < Folder containing RAW-files > \
    --ctm_outdir < Output-Folder, where the MGFs should be stored >
```

### Generating a MS2-Specific-FASTA

The workflow `crate_ms2_specific_fasta.nf` generates a ms2-specific-FASTA-file, tailored to a set of MGF-files. Here, Protein-Graphs are generated, using an SP-EMBL-file (which can be downloaded from [UniProt](https://www.uniprot.org/) by selecting `Text` as format) and a python script prepares the queries, by extracting the MS2-precursors from the MGF-files (using a tolerance, in ppm). Using the Protein-Graphs and a `DryRun`-Version of the traversal, the maximum-variant-limits are determined for each Protein-Graph (and mass-query-range) using a binary-search. These limits are then used for the actual ms2-specific-FASTA-generation in conjunction with the extracted MS2-precursors and a compacted FASTA is returned, which is tailored to the MGF-files.

Altough of the complexity, the workflow only requires the following parameters to generate such a FASTA:

```text
nextflow run convert_to_mgf.nf \
    --cmf_mgf_files < Folder containing MGF-files > \
    --cmf_sp_embl_file < Path to a SP-EMBL-File > \
    --cmf_outdir <The Output-Folder where the traversal-limits are saved and the ms2-specific-FASTA is stored >
```

The optional parameter: `cmf_pg_additional_params` is added to ProtGraph directly, allowing every parameter, ProtGraph provides to be set there (e.g. usefull if the digestion should be changed or features/PTMs should be included/excluded, etc...), allowing arbitrary settings to generate Protein-Graphs if desired. It defaults to use all features, ProtGraph can parse.

**Note regarding PTMs/Tolerance**: The FASTA is tailored to the MS2-precursors, therefore variable and fixed modifications need to be set to the same settings as for the actual identification. This workflow defaults to carbamidomethylation (C, fixed) and oxidation (M, variable). See ProtGraph (and the workflow-parameter `cmf_pg_additional_params`) to set the PTMs accordingly in the Protein-Graphs. The same applies for the MS2-precursor-tolereance which can be set with `cmf_query_ppm` and defaults to `5ppm`.

**Note regarding Limits**: This workflows defaults to allow up to 5 seconds per query and limits peptides to contain at most 5 variants (with a maximum of 5000 Da per peptide), resulting into FASTA-files which can be 15-200GB large (depending on dataset and species). Changing these settings can drastically increase/decrease the runtime/memory usage/disk usage. We advise to change those settings slightly and to pay attention on the runtime/memory usage/disk usage if run with the newly set limits (and dataset + species) the first time.

**Note regarding identification**: If digestion is enabled (default is `Trypsin`), the resulting FASTA contains already digested entries, thus searching with a search-engine, the digestion should be set to `off/no_cut`.

### Generating a Global-FASTA

This workflow generates a so called global-FASTA, using ProtGraph, an SP-EMBL and some global limits for writing out peptides/proteins. Global-FASTAs can be generated with the `create_global_fasta.nf`-workflow. To generate a global-FASTA, only a path to a single SP-EMBL-file is required. Such a file can be downloaded from [UniProt](https://www.uniprot.org/) directly, by selecting `Text` instead of `FASTA` as the download format.

```text
nextflow run create_global_fasta.nf \
    --cgf_sp_embl_file < Path to a SP-EMBL-File > \
    --cgf_outdir < The output-folder, where the gloabl-FASTA and some Protein-Graph-statistics should be saved >
```

Per default, this workflow does not export feature-peptides and is set to only export peptides with up to 5000 Da mass and maximum of two miscleavages. It is possible to generate global-FASTAs with some specific features (like containing, `SIGNAL`, `PEPTIDE` or others) and other limits. The parameters `cgf_features_in_graphs` and `cgf_peptide_limits` can be set accordingly. These are added to ProtGraph directly, hence every parameter ProtGraph provides, can be set here (including different digestion settings).

**Note**: A dry run with ProtGraph to generate statistics how many peptide would be theoretically exported is advised prior for testing. Some Protein-Graphs with some features (e.g. P53 using variants) can contain to many peptides, which could result to very long runtimes and huge FASTAs.

**Note regarding identification**: If digestion is enabled (default is `Trypsin`), the resulting FASTA contains already digested entries, thus searching with a search-engine, the digestion should be set to `off/no_cut`.

### Identification via Coment (and Percolator)

We provide an identification workflow to showcase, that the generated FASTAs can be used with search-engines. The workflow `identification_via_comet.nf` identifies MGF-files individually, using custom search-settings for Comet (and if desired rescores the results with Percolator), applies an FDR-cutoff using the q-value (for each file) and exposes the identification results into an output-folder.

Three parameters are required, to execute the workflow:

1. The MGFs which should be identified
2. The Comet-Parameter file to set the search-settings
3. The FASTA-file which should be used for identification

Below is an example call with all required parameters (Percolator is enabled by default):

```text
nextflow run identification_via_comet.nf \
    --idc_mgf_folder < Folder containing MGF-files > \
    --idc_fasta_file < The FASTA which should be used for identification > \
    --idc_search_parameter_file < The Comet-Parameters file (Search Configuration) > \
    --idc_outdir < Output-Folder where the results of the identification files are stored >
```

Here is another example call with all required parameters (this time, turning Percolator off):

```text
nextflow run identification_via_comet.nf \
    --idc_mgf_folder < Folder containing MGF-files > \
    --idc_fasta_file < The FASTA which should be used for identification > \
    --idc_search_parameter_file < The Comet-Parameters file (Search Configuration) > \
    --idc_outdir < Output-Folder where the results of the identification files are stored > \
    --idc_use_percolator 0
```

**Note**: This identification-workflow defaults to an FDR-cutoff (q-value) of `--idc_fdr "0.01|0.05"`, reporting both FDRs. Arbitrary and multiple FDR-cutoffs can be set and should be changed to the desired value.

### Summarization of results

The `summarize_ident_results.nf`-workflow genereate convenient summarization of the identification results. Here, the identification-results are binned into 4 groups:

1. Unique PSMs (a match, which can only originate from one protein)
2. Shared PSMs (a match, which can originate from multiple proteins)
3. Unique Feature PSMs (as 1., but only containing peptides, which can be explained by a features)
4. Shared Feature PSMs (as 2., but only can be explained by features from all originating proteins)

Furthermore, heatmaps are generated to provide an overview of found peptides across all MGFs/RAW-files.

To call this method, a `glob` needs to be specified in this workflow:

```text
nextflow run summarize_ident_results.nf \
    --sir_identified_files_glob < The glob matching the desired output from the identification results >
    --sir_outdir < The output directory where the summarized results should be saved >
```

In case, the identification workflow was executed using an FDR of 0.01, you could use the following `glob`:

```text
nextflow run summarize_ident_results.nf \
    --sir_identified_files_glob "<Path_to_folder>/*qvalue_no_decoys_fdr_0.01.tsv"
    --sir_outdir < The output directory where the summarized results should be saved >
```

**Note**: This step can be used only if specific columns are present in the tables. Furthermore, it distinguishes between the identification results from a FASTA by UniProt or by ProGFASTAGen. The additional parameters control, whether to bin results in group 3 and 4, decide if variable modifications should be considered as unique, as well as if a peptide, which originates multiple times to the same protein should be considered as unique. The main-workflows set these parameters accordingly and can be used as an example.

## Main Workflow Scripts

Each individual step described above, is also imported and chained into three main-workflows:

1. `main_workflow_protein_fasta.nf` (UniProt-FASTA-search)
2. `main_workflow_global_fasta.nf` (Generation of a global-FASTA and search)
3. `main_workflow_ms2_specific_fasta.nf` (Generation of a ms2-specific-FASTA and search)

generating summarized identification results across multiple RAW-files.

In each of these workflows, it is possible to modify the parameters of the imported subworkflows, by using the imported subworkflows parameters directly (as shown in the **Individual Steps** above).

For protein-FASTA identification, only three parameters are required:

```text
nextflow run main_workflow_protein_fasta.nf \
    --main_fasta_file < The FASTA-file, to be used for identification > \
    --main_raw_files_folder < The folder containing RAW-files > \
    --main_comet_params < The parameters file for comet (for identification) > \
    --main_outdir < Output-Folder where all the results from the workflows should be saved >
```

This is also true for the other two workflows, where instead of a FASTA-file, an SP-EMBL-file needs to be provided. Such a file can be downloaded from [UniProt](https://www.uniprot.org/) directly, by selecting the format `Text`.

Here are the correpsonding calls for global-FASTA and ms2-specific-FASTA generation and identification:

```text
# global-FASTA
nextflow run main_workflow_global_fasta.nf \
    --main_sp_embl_file < The SP-EMBL-file used for Protein-Graph- and FASTA-generation > \
    --main_raw_files_folder < The folder containing RAW-files > \
    --main_comet_params< The parameters file for comet (for identification) > \
    --main_outdir < Output-Folder where all the results from the workflows should be saved >

# ms2-specific-FASTA
nextflow run main_workflow_ms2_specific_fasta.nf \
    --main_sp_embl_file < The SP-EMBL-file used for Protein-Graph- and FASTA-generation > \
    --main_raw_files_folder < The folder containing RAW-files > \
    --main_comet_params < The parameters file for comet (for identification) > \
    --main_outdir < Output-Folder where all the results from the workflows should be saved >
```

**Note**: Only defining the required parameters, uses the default parameters for every other setting. For all workflows, this would mean, that the FDR-cutoff (q-value) is set to `0.01|0.05` resulting into both FDRs considered. Furthermore, the global-FASTA and ms2-specific-FASTA workflows assume Trypsin digestion. For the global-FASTA-workflow, no features are exported by default, which may not be desired, if someone whishes to search for peptide-features (like `SIGNAL`, etc..). For the ms2-specific-FASTA-workflow, the PTMs carbamidomethylation (C, fixed) and oxidation (M, variable) are assumed, which may need to be modified.

**Note regarding example calls**: Further below you can find the calls as used in the publication. These set the most minimal parameters for a correct execution on custom datasets and can be used as an example.

## Regenerate Results from Publication

In this subsection you can find the nextflow-calls which were used to execute the 3 workflows. Executing this with the same SP-EMBL-/FASTA-file should yield the similar/same results. For generated ms2-specific-FASTAs it may happen, that these are generated with slightly different variant-limits, therefore a slightly different FASTA to search with and slightly different identification results.

The FASTA/SP-EMBL used for identification can be found [here](TODO_DL). The Comet configuration files are provided in the `example_configuration`-folder. The datasets can be retrieved from [PRIDE](https://www.ebi.ac.uk/pride/).

### PXDXXXXXX

```shell
# PXDXXXXXX ms2-specific-FASTA
nextflow run main_workflow_ms2_specific_fasta.nf \
    -with-report "PXDXXXXXX_ms2_specific_fasta/nextflow_report.html" \
    -with-timeline "PXDXXXXXX_ms2_specific_fasta/nextflow_timeline.html" \
    --main_sp_embl_file 20230619_mus_musculus_proteome.txt \
    --main_raw_files_folder PXDXXXXXX \
    --main_comet_params example_configurations/PXDXXXXXX_no_dig.txt \
    --main_outdir PXDXXXXXX_ms2_specific_fasta \
    --cmf_max_precursor_da 5000 \
    --cmf_query_ppm 5 \
    --cmf_timeout_for_single_query 5 \
    --cmf_maximum_variant_limit -1 \
    --cmf_pg_additional_params "-ft ALL -fm 'C:57.021464' -vm 'M:15.9949' -vm 'Q:-17.026549' -vm 'Q:0.984016' -vm 'N:0.984016'" \
    --idc_fdr "0.01"

# PXDXXXXXX global-FASTA
nextflow run main_workflow_global_fasta.nf \
    -with-report "PXDXXXXXX_global_fasta/nextflow_report.html" \
    -with-timeline "PXDXXXXXX_global_fasta/nextflow_timeline.html" \
    --main_sp_embl_file 20230619_mus_musculus_proteome.txt \
    --main_raw_files_folder PXDXXXXXX \
    --main_comet_params example_configurations/PXDXXXXXX_no_dig.txt \
    --main_outdir PXDXXXXXX_global_fasta \
    --cgf_features_in_graphs "-ft None" \
    --cgf_peptide_limits "--pep_miscleavages 2 --pep_min_pep_length 5 --pep_max_weight 5000" \
    --idc_fdr "0.01"

# PXDXXXXXX protein-FASTA
nextflow run main_workflow_protein_fasta.nf \
    -with-report "PXDXXXXXX_protein_fasta/nextflow_report.html" \
    -with-timeline "PXDXXXXXX_protein_fasta/nextflow_timeline.html" \
    --main_fasta_file 20230619_mus_musculus_proteome.fasta \
    --main_raw_files_folder PXDXXXXXX \
    --main_comet_params example_configurations/PXDXXXXXX_trypsin_dig.txt \
    --main_outdir PXDXXXXXX_protein_fasta \
    --idc_fdr "0.01"
```

### PXD002171

```shell
# PXD002171 ms2-specific-FASTA
nextflow run main_workflow_ms2_specific_fasta.nf \
    -with-report "PXD002171_ms2_specific_fasta/nextflow_report.html" \
    -with-timeline "PXD002171_ms2_specific_fasta/nextflow_timeline.html" \
    --main_sp_embl_file 20230619_homo_sapiens_proteome.txt \
    --main_raw_files_folder PXD002171 \
    --main_comet_params example_configurations/PXD002171_no_dig.txt \
    --main_outdir PXD002171_ms2_specific_fasta \
    --cmf_max_precursor_da 5000 \
    --cmf_query_ppm 5 \
    --cmf_timeout_for_single_query 5 \
    --cmf_maximum_variant_limit 5 \
    --cmf_pg_additional_params "-ft VARIANT -ft SIGNAL -ft INIT_MET -ft CONFLICT -ft VAR_SEQ -ft PEPTIDE -ft PROPEP -ft CHAIN -vm 'M:15.994915' -vm 'C:71.037114'" \
    --idc_fdr "0.01"

# PXD002171 global-FASTA
nextflow run main_workflow_global_fasta.nf \
    -with-report "PXD002171_global_fasta/nextflow_report.html" \
    -with-timeline "PXD002171_global_fasta/nextflow_timeline.html" \
    --main_sp_embl_file 20230619_homo_sapiens_proteome.txt \
    --main_raw_files_folder PXD002171 \
    --main_comet_params example_configurations/PXD002171_no_dig.txt \
    --main_outdir PXD002171_global_fasta \
    --cgf_features_in_graphs "-ft None" \
    --cgf_peptide_limits "--pep_miscleavages 2 --pep_min_pep_length 5 --pep_max_weight 5000" \
    --idc_fdr "0.01"

# PXD002171 protein-FASTA
nextflow run main_workflow_protein_fasta.nf \
    -with-report "PXD002171_protein_fasta/nextflow_report.html" \
    -with-timeline "PXD002171_protein_fasta/nextflow_timeline.html" \
    --main_fasta_file 20230619_homo_sapiens_proteome.fasta \
    --main_raw_files_folder PXD002171 \
    --main_comet_params example_configurations/PXD002171_trypsin_dig.txt \
    --main_outdir PXD002171_protein_fasta \
    --idc_fdr "0.01"
```

### PXD028605

```shell
# PXD028605 ms2-specific-FASTA
nextflow run main_workflow_ms2_specific_fasta.nf \
    -with-report "PXD028605_ms2_specific_fasta/nextflow_report.html" \
    -with-timeline "PXD028605_ms2_specific_fasta/nextflow_timeline.html" \
    --main_sp_embl_file 20230619_homo_sapiens_proteome.txt \
    --main_raw_files_folder PXD028605 \
    --main_comet_params example_configurations/PXD028605_no_dig.txt \
    --main_outdir PXD028605_ms2_specific_fasta \
    --cmf_max_precursor_da 5000 \
    --cmf_query_ppm 20 \
    --cmf_timeout_for_single_query 5 \
    --cmf_maximum_variant_limit 5 \
    --cmf_pg_additional_params "-ft VARIANT -ft SIGNAL -ft INIT_MET -ft CONFLICT -ft VAR_SEQ -ft PEPTIDE -ft PROPEP -ft CHAIN -fm 'C:57.021464' -vm 'M:15.9949'" \
    --idc_fdr "0.01"

# PXD028605 global-FASTA
nextflow run main_workflow_global_fasta.nf \
    -with-report "PXD028605_global_fasta/nextflow_report.html" \
    -with-timeline "PXD028605_global_fasta/nextflow_timeline.html" \
    --main_sp_embl_file 20230619_homo_sapiens_proteome.txt \
    --main_raw_files_folder PXD028605 \
    --main_comet_params example_configurations/PXD028605_no_dig.txt \
    --main_outdir PXD028605_global_fasta \
    --cgf_features_in_graphs "-ft None" \
    --cgf_peptide_limits "--pep_miscleavages 2 --pep_min_pep_length 5 --pep_max_weight 5000" \
    --idc_fdr "0.01"

# PXD028605 protein-FASTA
nextflow run main_workflow_protein_fasta.nf \
    -with-report "PXD028605_protein_fasta/nextflow_report.html" \
    -with-timeline "PXD028605_protein_fasta/nextflow_timeline.html" \
    --main_fasta_file 20230619_homo_sapiens_proteome.fasta \
    --main_raw_files_folder PXD028605 \
    --main_comet_params example_configurations/PXD028605_trypsin_dig.txt \
    --main_outdir PXD028605_protein_fasta \
    --idc_fdr "0.01"
```
