FTPDIR=/nfs/ftp/public/databases/wormbase/parasite/releases/current
FTPCHECKSUMS=${FTPDIR}/CHECKSUMS
LOCAL_FTP=/hps/nobackup/flicek/wormbase/parasite/dumps/WBPS16_postrelease/FTP
LOCALCHECKSUMS=${LOCAL_FTP}/CHECKSUMS
OUTFILE=~/compare.txt

rm $OUTFILE;
touch $OUTFILE;

while read p; do
  locmd5=$(echo $p | cut -d' ' -f1);
  locfile_suffix=$(echo $p | cut -d' ' -f2 | cut -d'/' -f2-);
  locfile=${LOCAL_FTP}/$(echo $p | cut -d' ' -f2 | cut -d'/' -f2-);
  ftpfile=${FTPDIR}/species/${locfile_suffix}
  echo "Checking ${ftpfile}...;"
  if [ -f "${ftpfile}" ]; then
    ftpmd5=$(md5sum ${ftpfile} | cut -d' ' -f1);
    if [ "$locmd5" != "$ftpmd5" ]; then
      CDIFF=$(zdiff ${locfile} ${ftpfile} | head)
      if [ "$CDIFF" != "" ]; then
      printf "${locfile}\t${ftpfile}\tdifferent\n" >> $OUTFILE;
      fi;
    fi;
    else
      printf "${locfile}\t${ftpfile}\tnot_existing\n" >> $OUTFILE;
  fi;
done < $LOCALCHECKSUMS