# Calls as in the Publication

## Regenerate the Identification Results from the Publication

Below you can find the command how the results of this publication were generated. The Comet-Configuration-files are also provided in the `example_configuration`-folder. The TXT- and FASTA-files are not included in this repository.

### PXDXXXXXX

```shell
# PXDXXXXXX MS2-Specific
nextflow run main_workflow_ms2_specific_fasta.nf \
    -with-report "PXDXXXXXX_results_ms2_specific/nextflow_report.html" \
    -with-timeline "PXDXXXXXX_results_ms2_specific/nextflow_timeline.html" \
    --main_sp_embl_file 20230619_mus_musculus_proteome.txt \
    --main_raw_files_folder PXDXXXXXX \
    --main_comet_params example_configurations/PXDXXXXXX_no_dig.txt \
    --main_outdir PXDXXXXXX_results_ms2_specific \
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
    -with-report "PXD002171_results_ms2_specific/nextflow_report.html" \
    -with-timeline "PXD002171_results_ms2_specific/nextflow_timeline.html" \
    --main_sp_embl_file 20230619_homo_sapiens_proteome.txt \
    --main_raw_files_folder PXD002171 \
    --main_comet_params example_configurations/PXD002171_no_dig.txt \
    --main_outdir PXD002171_results_ms2_specific \
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
    -with-report "PXD028605_results_ms2_specific/nextflow_report.html" \
    -with-timeline "PXD028605_results_ms2_specific/nextflow_timeline.html" \
    --main_sp_embl_file 20230619_homo_sapiens_proteome.txt \
    --main_raw_files_folder PXD028605 \
    --main_comet_params example_configurations/PXD028605_no_dig.txt \
    --main_outdir PXD028605_results_ms2_specific \
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
