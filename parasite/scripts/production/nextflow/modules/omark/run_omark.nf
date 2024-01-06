process RUN_OMARK {

  label 'omamer'

  input:
  file omamer_file
  val outdir
  val db


  output:
  path("omark_output/proteins.sum"), emit: summary_file
  val outdir, emit:species_outdir
  val db

  publishDir "${params.output_dir}/${outdir}/",  mode: 'copy'

  script:
  """
  omark -f ${omamer_file} -d ${params.omamer_database} -o omark_output
  """
}