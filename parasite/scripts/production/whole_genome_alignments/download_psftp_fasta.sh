#!/bin/bash

# URL of the file to be downloaded
url=$1
# Location to save the file
destination=$2
# URL of the online TSV file
CHECKSUMS=$3
# URL to strip from the url
fasta_name_in_CHECKSUMS=$4
# Get the basename of the file to be downloaded
echo $url
echo $destination
echo $CHECKSUMS
echo $fasta_name_in_CHECKSUMS


# Get the expected checksum value from the online TSV file
expected_checksum=$(curl -s "${CHECKSUMS}" | grep "$fasta_name_in_CHECKSUMS" | awk '{print $1}')
expected_checksum_count=$(echo "$expected_checksum" | wc -l)
if [ "$expected_checksum_count" -ne 1 ] | [ -z "$expected_checksum" ]; then
  echo "Error: expected 1 match in the $CHECKSUMS file. Got 0 or many. Exiting."
  exit 1
fi

# Download the file
wget -O "$destination" "$url"

# Check if the download was successful
if [ $? -eq 0 ]; then
  echo "Download successful."
  else echo "Download unsuccesful."
  exit 1;
fi

# Check if the file is complete and not corrupt
checksum=$(md5sum "$destination" | awk '{print $1}')
if [ "$checksum" == "$expected_checksum" ]; then
  echo "File is complete and not corrupt."
else
  echo "checksum: $checksum"
  echo "expected_checksum: $expected_checksum"
  echo "File is corrupt."
  exit 1
fi

# Unzip if needed
if [[ $destination == *.gz ]]; then
  # Uncompress the gzipped file
  echo "Decompressing..."
  gzip -df "$destination"
  # Update the destination path
  destination=${destination%.gz}
fi
echo "Done"