# PXD028605 Precursor-Specific
nextflow run main_workflow_precursor_specific_fasta.nf \
    -with-report "PXD028605_results_precursor_specific/nextflow_report.html" \
    -with-timeline "PXD028605_results_precursor_specific/nextflow_timeline.html" \
    --main_sp_embl_file 20230619_homo_sapiens_proteome.txt \
    --main_raw_files_folder PXD028605 \
    --main_comet_params example_configurations/PXD028605_no_dig.txt \
    --main_outdir PXD028605_results_precursor_specific \
    --cmf_max_precursor_da 5000 \
    --cmf_query_ppm 20 \
    --cmf_timeout_for_single_query 5 \
    --cmf_maximum_variant_limit 5 \
    --cmf_pg_additional_params "-ft VARIANT -ft SIGNAL -ft INIT_MET -ft CONFLICT -ft VAR_SEQ -ft PEPTIDE -ft PROPEP -ft CHAIN -fm 'C:57.021464' -vm 'M:15.9949'" \
    --idc_fdr "0.01"

rm -rf work/
    
# PXD028605 Global digested FASTA
nextflow run main_workflow_global_fasta.nf \
    -with-report "PXD028605_global_fasta/nextflow_report.html" \
    -with-timeline "PXD028605_global_fasta/nextflow_timeline.html" \
    --main_sp_embl_file 20230619_homo_sapiens_proteome.txt \
    --main_raw_files_folder PXD028605 \
    --main_comet_params example_configurations/PXD028605_no_dig.txt \
    --main_outdir PXD028605_global_fasta \
    --cgf_features_in_graphs "-ft None" \
    --cgf_peptide_limits "--pep_miscleavages 2 --pep_min_pep_length 5" \
    --idc_fdr "0.01"

rm -rf work/
    
# PXD028605 Protein FASTA
nextflow run main_workflow_protein_fasta.nf \
    -with-report "PXD028605_protein_fasta/nextflow_report.html" \
    -with-timeline "PXD028605_protein_fasta/nextflow_timeline.html" \
    --main_fasta_file 20230619_homo_sapiens_proteome.fasta \
    --main_raw_files_folder PXD028605 \
    --main_comet_params example_configurations/PXD028605_trypsin_dig.txt \
    --main_outdir PXD028605_protein_fasta \
    --idc_fdr "0.01"
    
rm -rf work/

# PXD002171 Precursor-Specific
nextflow run main_workflow_precursor_specific_fasta.nf \
    -with-report "PXD002171_results_precursor_specific/nextflow_report.html" \
    -with-timeline "PXD002171_results_precursor_specific/nextflow_timeline.html" \
    --main_sp_embl_file 20230619_homo_sapiens_proteome.txt \
    --main_raw_files_folder PXD002171 \
    --main_comet_params example_configurations/PXD002171_no_dig.txt \
    --main_outdir PXD002171_results_precursor_specific \
    --cmf_max_precursor_da 5000 \
    --cmf_query_ppm 5 \
    --cmf_timeout_for_single_query 5 \
    --cmf_maximum_variant_limit 5 \
    --cmf_pg_additional_params "-ft VARIANT -ft SIGNAL -ft INIT_MET -ft CONFLICT -ft VAR_SEQ -ft PEPTIDE -ft PROPEP -ft CHAIN -vm 'M:15.994915' -vm 'C:71.037114'" \
    --idc_fdr "0.01"

rm -rf work/
    
# PXD002171 Global digested FASTA
nextflow run main_workflow_global_fasta.nf \
    -with-report "PXD002171_global_fasta/nextflow_report.html" \
    -with-timeline "PXD002171_global_fasta/nextflow_timeline.html" \
    --main_sp_embl_file 20230619_homo_sapiens_proteome.txt \
    --main_raw_files_folder PXD002171 \
    --main_comet_params example_configurations/PXD002171_no_dig.txt \
    --main_outdir PXD002171_global_fasta \
    --cgf_features_in_graphs "-ft None" \
    --cgf_peptide_limits "--pep_miscleavages 2 --pep_min_pep_length 5" \
    --idc_fdr "0.01"

rm -rf work/

# PXD002171 Protein FASTA
nextflow run main_workflow_protein_fasta.nf \
    -with-report "PXD002171_protein_fasta/nextflow_report.html" \
    -with-timeline "PXD002171_protein_fasta/nextflow_timeline.html" \
    --main_fasta_file 20230619_homo_sapiens_proteome.fasta \
    --main_raw_files_folder PXD002171 \
    --main_comet_params example_configurations/PXD002171_trypsin_dig.txt \
    --main_outdir PXD002171_protein_fasta \
    --idc_fdr "0.01"

rm -rf work/


