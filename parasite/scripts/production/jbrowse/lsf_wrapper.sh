
species=($(ls -d $2/species/*))

for s in "${species[@]}"
do
  s="$(basename $s)"
  bsub -q production-rh7 ". /nfs/public/rw/ensembl/perlbrew/setup_gunpowder_perl.sh; PERL5LIB=$1/lib:$PERL5LIB:/nfs/public/rw/ensembl/libs/perl/share/perl5/x86_64-linux/:/nfs/public/rw/ensembl/libs/perl/share/perl5/; perl $1/create_species_config.pl --ftp_dir $2 --jbrowse_path $3 --out_dir $4 --gff3_config_file  $5 --species $s"
done

