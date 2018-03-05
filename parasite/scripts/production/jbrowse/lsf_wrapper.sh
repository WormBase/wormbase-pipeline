
species=($(ls -d $2/species/*))

for s in "${species[@]}"
do
  s="$(basename $s)"
  bsub -q production-rh7 "module load parasite_prod_rel9; export PERL5LIB=/nfs/panda/ensemblgenomes/wormbase/software/perlmodules/:/nfs/panda/ensemblgenomes/wormbase/software/perlmodules/x86_64-linux-thread-multi-ld/:/nfs/public/rw/ensembl/libs/perl/share/perl5/:$PERL5LIB; perl $1/create_species_config.pl --ftp_dir $2 --jbrowse_path $3 --out_dir $4 --species $s"
done

