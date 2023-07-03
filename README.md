# ProGFASTAGen

The ProGFASTAGen (**Pro**tein-**G**raph-**FASTA**-**Gen**erator or **Pro**t**G**raph-**FASTA**-**Gen**erator) repository contains workflows to generate so-called ms2-specific-FASTAs (using the precursor from MGF-files) including feature-peptides, like VARIANTs or CONFLICTs if desired, or global-FASTAs (as described in [ProtGraph](https://github.com/mpc-bioinformatics/ProtGraph)). The single workflow scripts have been implemented with [Nextflow-DSL-2](https://www.nextflow.io/docs/latest/dsl2.html) and are independent to each other. Each of them can be used on their own or can be imported to other workflows for other use-cases. Further, we included three main-workflows, to show how the single workflows can be chained together. The `main_workflow_protein_fasta.nf`-workflow converts Thermo-RAW-files into MGF, searches with Comet (and Percolator) and the identification results are then further summarized. The workflows `main_workflow_global_fasta.nf` and `main_workflow_ms2_specific_fasta.nf` generate specific FASTA-files before search-engine-identification. Below are example nextflow-calls, which can be used.

Regarding, the ms2-specific-FASTA-Generation: The source-code of the C++ implementation for traversal can be found in `bin`. There, four implementations are present: `Float/Int`-Versions of the traversal as well as `DryRun/VarLimitter`-Versions of the traversal. The `Float/Int`-Versions can be faster/slower on some specific processors and both versions can be used via a flag in the workflows. The `DryRun`-Version does not generate a FASTA but tests the used system (depending on a query-timeout) to determine the maximum number of variants which can be used, while not timing out. The actual FASTA-Generation happens in the `VarLimitter`-Version using the generated protein-graphs at hand.

in **Prerequisites** a small description of dependencies and how to set up an environment is given. **Individual steps** describes the single workflows and how they can be called, while **Main Workflow Scripts** shows example-calls of the main workflows. In **Regenerate Results from Publication**, the parameters are shown, which were used in the publication and using the same FASTA, SP-EMBL with a similar server-setting should yield similar results as used in the publication.

## Prerequisites

### Executing on Linux

This workflow can be only executed on linux (tested on Ubuntu 22.04 and ArchLinux). Before setting up the `bin`-folder, some requiered binaries are need to be availabe for the OS. (Focusing on Ubuntu:) The following packages need to be installed on Ubuntu (via `apt`):

```text
build-essential
wget
curl
unzip
cmake
mono-complete
python3-pip (or any environment with Python3)
python-is-python3 (needed for ubuntu)
```

If all required binaries are available, a setup-script needs to be executed, which downloads needed dependencies and compiles the source-code located in the `bin`-folder.

```shell
chmod +x compile_and_setup_depencies.sh  # In case this file is not executable
./compile_and_setup_depencies.sh  # Downloads dependencies, compiles the C++-implementation and sets all binaries in the bin-folder as executable
```

If the above worked without errors, the provided workflows can be executed with the command `nextflow`.

### Executing in Docker

Alternatively, docker can be used. For this, please follow the [installation guide](https://docs.docker.com/engine/install/ubuntu/) for docker. After docker has been installed, a local docker-container can be build with all needed dependencies. We provide a `Dockerfile`in the `docker`-folder. To build it: execute (while beeing with a shell in the root-folder of this repository) the following:

```shell
docker build -t progfastagen:local . -f docker/Dockerfile
```

Make sure that `nextflow` installed on the host-system. For each  ofthe workflow example calls below, the `-with-docker progfastagen:local` then needs to be appended.

## Individual Steps

Each step has been implemented in such a way, that it can be executed on its own. Each subsection below, provides a brief overview of the needed parameters to execute the workflow itself. If you are interested for all the available parameters and want to tune or modify them, then please refer to the source of the workflow, whre each parameters is described briefly.

### Converting RAW-files to MGF

This workflow converts RAW-files to the MGF-format. The `ctm_raws` parameter is required, in order to generate the MGF-files

```text
nextflow run convert_to_mgf.nf \
    -- ctm_raws < Folder containing RAW-files > \
    -- ctm_outdir < Output-Folder, where the MGFs are stored >
```

### Generating a MS2-Specific-FASTA

TBD

### Generating a Global-FASTA

TBD

### Identification via Coment (and Percolator)

This workflow identifies MGF-files individually, using custom search-settings, applies an FDR-cutoff using the q-value afterwards (for each file) and exposes the identification results into an output-folder.

Example call, not using Percolator:

```text
nextflow run convert_to_mgf.nf \
    --idc_mgf_folder < Folder containing MGF-files > \
    --idc_fasta_file < The FASTA which should be used for identification > \
    --idc_search_parameter_file < The Comet-Parameters file (Search Configuration) > \
    --idc_outdir < Output-Folder where the results of the identification files are stored > \
    --idc_use_percolator 0
```

Example call, using Percolator:

```text
nextflow run convert_to_mgf.nf \
    --idc_mgf_folder < Folder containing MGF-files > \
    --idc_fasta_file < The FASTA which should be used for identification > \
    --idc_search_parameter_file < The Comet-Parameters file (Search Configuration) > \
    --idc_outdir < Output-Folder where the results of the identification files are stored >
```

**NOTE**: This workflow defaults to an fdr-cutoff (q-value) of `--idc_fdr "0.01|0.05"`, reporting both FDRs. Arbitrary and multiple FDR-cutoffs can be set and should be changed to the desired value.

### Summarization of results

TBD

## Regenerate Results from Publication

TBD
Below you can find the commands how the results of this publication were generated. The Comet-Configuration-files are also provided in the `example_configuration`-folder. The TXT- and FASTA-files are not included in this repository.

### PXDXXXXXX

```shell
# PXDXXXXXX MS2-Specific
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

# PXDXXXXXX Global digested FASTA
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

# PXDXXXXXX protein_fasta
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
# PXD002171 MS2-Specific
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

# PXD002171 Global digested FASTA
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

# PXD002171 protein_fasta
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
# PXD028605 MS2-Specific
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

# PXD028605 Global digested FASTA
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

# PXD028605 protein_fasta
nextflow run main_workflow_protein_fasta.nf \
    -with-report "PXD028605_protein_fasta/nextflow_report.html" \
    -with-timeline "PXD028605_protein_fasta/nextflow_timeline.html" \
    --main_fasta_file 20230619_homo_sapiens_proteome.fasta \
    --main_raw_files_folder PXD028605 \
    --main_comet_params example_configurations/PXD028605_trypsin_dig.txt \
    --main_outdir PXD028605_protein_fasta \
    --idc_fdr "0.01"
```
