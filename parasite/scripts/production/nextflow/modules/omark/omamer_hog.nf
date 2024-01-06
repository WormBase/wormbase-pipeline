process OMAMER_HOG {

  label 'omamer'

  input:
  file protein
  val outdir
  val db
  val augustus_species
  val odb

  output:
  path "proteins.omamer", emit: omamer_file
  val outdir, emit:species_outdir
  val db

  // ourdir is Salmo_trutta (production name)
  publishDir "${params.output_dir}/${outdir}/",  mode: 'copy'

  script:
  """
  omamer search --db ${params.omamer_database} --query ${protein} --score sensitive --out proteins.omamer
  """
}