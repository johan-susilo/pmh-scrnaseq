process PREPROCESS {
  conda "${projectDir}/env/scRNA_env.yml"
  // tell Nextflow where to save the final outputs
  publishDir "${projectDir}/results/PMH_preprocessing/", mode: 'copy'

  input:
  val(row)

  output:
  path "output/*.rds"

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

workflow {
  sample_ch = Channel
    .fromPath(params.manifest)
    .splitCsv(header: true)
  
  PREPROCESS(sample_ch)
}