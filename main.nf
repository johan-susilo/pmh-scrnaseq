process LOAD {
  input:
  val(row)

  output:
  // emits a pair: [ "Sample1", /path/to/raw_load_Sample1.rds ]
  tuple val(row.sample_name), path("output/*.rds")

  script:
  """
  
  Rscript ${projectDir}/scripts/01_preprocessing.R \\
    --project_dir ${projectDir} \\
    --sample_name ${row.sample_name} \\
    --ident1 ${row.ident1} \\
    --ident2 ${row.ident2} \\
    --data_path ${row.data_path} \\
    --output_dir ./output

  """
}

process QC_FILTER {

  // sort Seurat data objects into a 'data' folder
  publishDir "${params.outdir}/QC/data", pattern: "output/*.rds", mode: 'copy'
  // sort QC plots into 'plots' folder
  publishDir "${params.outdir}/QC/plots", pattern: "output/*.png", mode: 'copy'
  publishDir "${params.outdir}/QC/plots", pattern: "output/*.pdf", mode: 'copy'

  input:
  tuple val(sample_name), path(raw_rds)

  output:
  // using 'emit' allows us to grab specific file types in the workflow block!
  path "output/*_clean.rds", emit: clean_rds
  path "output/*.png", emit: png_plots, optional: true
  path "output/*.pdf", emit: pdf_plots, optional: true

  script:
  """
  Rscript ${projectDir}/scripts/02_qc_filter.R \\
    --project_dir ${projectDir} \\
    --sample_name ${sample_name} \\
    --input_rds ${raw_rds} \\
    --output_dir ./output
  """

}

process INTEGRATE {

  publishDir "${params.outdir}/Integration/data", pattern: "output/*.rds", mode: "copy"
  publishDir "${params.outdir}/Integration/plots", pattern: "output/*.pdf", mode: "copy"
  publishDir "${params.outdir}/Integration/plots", pattern: "output/*.png", mode: "copy"

  input:
  path rds_files

  output:
  path "output/TN.combined_dim30.rds", emit: integrated_rds
  path "output/*.png", emit: png_plots, optional: true
  path "output/*.pdf", emit: pdf_plots, optional: true

  script:

  """
  Rscript ${projectDir}/scripts/03_integration.R \\
    --project_dir ${projectDir} \\
    --output_dir ./output \\
    --resolution 0.4 \\
    --use_sct TRUE \\
    ${rds_files}
  """

}

process ANNOTATE {
  conda "${projectDir}/env/annotation_env.yml"

  publishDir "${params.outdir}/Annotation", mode: 'copy'

  input:
  path integrated_rds

  output:
  path "annotations/res_*/consensus/consensus_annotation.tsv", emit: consensus_tsv

  script:
  """
  Rscript ${projectDir}/scripts/04_annotation.R \\
    --project_dir ${projectDir} \\
    --rds ${integrated_rds} \\
    --output_dir ./annotations \\
    --resolution 0.4 \\
    --tissue skin
  """
}

process PLOT {
  conda "${projectDir}/env/seurat_env.yml" 
  publishDir "${params.outdir}/Final_Results", mode: 'copy'

  input:
  path integrated_rds
  path consensus_file

  output:
  path "plots/**"
  path "plots/TN.combined_ANNOTATED.rds"

  script:
  """
  Rscript ${projectDir}/scripts/05_plotting.R \\
    --project_dir ${projectDir} \\
    --rds ${integrated_rds} \\
    --consensus ${consensus_file} \\
    --output_dir ./plots \\
    --resolution 0.4
  """
}



workflow {
  sample_ch = Channel
    .fromPath(params.manifest)
    .splitCsv(header: true)
  
  loaded_ch = LOAD(sample_ch)
  filtered_ch = QC_FILTER(loaded_ch)
  all_clean_rds = filtered_ch.clean_rds.collect()

  integrated_ch = INTEGRATE(all_clean_rds)

  annotated_ch = ANNOTATE(integrated_ch.integrated_rds)
  
  PLOT(integrated_ch.integrated_rds, annotated_ch.consensus_tsv)

}