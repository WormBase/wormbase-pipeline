#!/bin/bash

# Function to display help message
show_help() {
    echo "IMPORTANT: The script needs an environmental variable named APOLLO_DIR to be set. This should point to the"
    echo "apollo installation path."
    echo ""
    echo "Usage: $0 [--apname <name>] [--fasta <file/url>] [--gff3 <file/url>] [--type <types>] [--repeats]"
    echo "Options:"
    echo "  --apname         Name of the Apollo instance to create."
    echo "  --fasta          Path or URL to the FASTA file used for creating the Apollo instance."
    echo "  --gff3           Path or URL to the GFF3 file used for generating the gene models tracks."
    echo "                   The GFF3 file should follow the WBPS convention."
    echo "  --type           GFF feature types to be used for the flatfile-to-json.pl command."
    echo "                   See the documentation for --type option at: https://jbrowse.org/docs/flatfile-to-json.pl.html"
    echo "                   Example input: gene,mRNA,tRNA (default: mRNA)"
    echo "  --no-repeats     Use this flag to disable the addition of repeat tracks to the Apollo instance (Default: Enabled)"
    echo "  --rnaseq_species Enter a WBPS species name and WBPS release version. If this option is used, the script will all rnaseq Jbrowse"
    echo "                   tracks from a WBPS release for this species to your apollo instance. Input Format: <species>:<release>."
    echo "                   Example: teladorsagia_circumcincta_prjna72569:18"
    echo "  --python_path    Full path to a python executable. Recommended version 3.9.5. Defaults to: python3"
}

# Check if the APOLLO_DIR environmental variable exists
if [[ -z $APOLLO_DIR ]]; then
    echo "The script needs an environmental variable named APOLLO_DIR to be set."
    echo "This should point to the apollo installation path."
    exit 1
fi

script_dir=$(dirname "$0")

# Set deault value for --repeats and --rnaseq_species
repeats_enabled=true

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" ing
        --apname)
            APNAME="$2"
            shift 2
            ;;
        --fasta)
            INPUT_FASTA_FILE="$2"
            shift 2
            ;;
        --gff3)
            INPUT_GFF_FILE="$2"
            shift 2
            ;;
        --type)
            GFF_TYPE="$2"
            shift 2
            ;;
        --repeats)
            repeats_enabled=true
            shift
            ;;
        --no-repeats)
            repeats_enabled=false
            shift
            ;;
        --rnaseq_species)
            RNASEQ_SPECIES_RELEASE="$2"
            shift 2
            ;;            
        --python_path)
            PYTHON_PATH="$2"
            shift 2
            ;;                        
        --help)
            show_help
            exit 0
            ;;
        *)
            echo "Error: Unknown option '$1'"
            show_help
            exit 1
            ;;
    esac
done

# Check if required arguments are provided
if [[ -z $APNAME || -z $INPUT_FASTA_FILE || -z $INPUT_GFF_FILE ]]; then
    echo "Error: Missing required arguments."
    show_help
    exit 1
fi

# Set default value for --type if not provided
if [[ -z $GFF_TYPE ]]; then
    GFF_TYPE="mRNA"
fi

# Set default value for --rnaseq_species if not provided
if [[ -z $RNASEQ_SPECIES_RELEASE ]]; then
    RNASEQ_SPECIES_RELEASE=false
fi

# Set default value for --type if not provided
if [[ -z $PYTHON_PATH ]]; then
    PYTHON_PATH=python3
fi

# Use the provided arguments in your script logic
echo "Apollo Name: $APNAME"
echo "FASTA File: $INPUT_FASTA_FILE"
echo "GFF3 File: $INPUT_GFF_FILE"
echo "GFF Feature Types: $GFF_TYPE"
echo "Repeats Enabled: $repeats_enabled"
echo "RNASEQ Jbrowse input: $RNASEQ_SPECIES_RELEASE"

if [[ "${RNASEQ_SPECIES_RELEASE}" != false ]]; then
    if [[ $RNASEQ_SPECIES_RELEASE =~ ^([^:]+):([^:]+)$ ]]; then
        RNASEQ_SPECIES=${BASH_REMATCH[1]}
        RNASEQ_RELEASE=${BASH_REMATCH[2]}
        echo "WBPS Jbrowse Species: $RNASEQ_SPECIES"
        echo "WBPS Release: $RNASEQ_RELEASE"
    else
        echo "Invalid format for --rnaseq_species: $RNASEQ_SPECIES_RELEASE"
        exit 1
    fi
fi

# Apollo software directory
APBIN=${APOLLO_DIR}/bin

# Apollo instance path
APDIR=${APOLLO_DIR}/data/${APNAME}

# Final FASTA/GFF3 file paths
FASTA_FILE=${APDIR}/genome.fa
GFF_FILE=${APDIR}/annotation.gff3

if [[ -d "${APDIR}" ]]; then
    read -p "${APDIR} already exists. Do you want to overwrite it? (y/n): " answer
    if [[ "$answer" == [Yy] ]]; then
        rm -rf "${APDIR}"
        mkdir "${APDIR}"
        echo "Directory overwritten."
    else
        echo "Operation canceled. Directory not overwritten."
        exit 1
    fi
else
    mkdir -p "${APDIR}"
    echo "Directory created: ${APDIR}"
fi

# Function to download and move a FILE/URL to a specified place.
check_and_download() {
    local infile="$1"
    local outfile="$2"

    # Check if the input is a URL 
    if [[ $infile == http://* || $infile == https://* || $infile == ftp://* ]]; then
        echo "$infile is a URL. Downloading..."
        curl -s "$infile" | { [ "${infile##*.}" = "gz" ] && gzip -d || cat; } > $outfile

        # Check the exit status of the curl command
        if [ $? -ne 0 ]; then
            echo "Download failed for $infile. Exiting..."
            exit 1
        fi
    else
        # Check if the input is an existing file path
        if [[ -f $infile ]]; then
            echo "$infile is a file path. Not downloading."

            # Check if the file ends with ".gz" extension
            if [[ $infile == *.gz ]]; then
                echo "Decompressing file..."
                gzip -c -d "$infile" > "$outfile"  # Decompresses and saves to the outfile
            else
                mv "$infile" "$outfile"  # Renames the file with the specified outfile name or the same infile name
                echo "File renamed to $outfile"
            fi
        else
            echo "$infile is neither a URL nor a file path."
        fi
    fi
}

# Function to check if a specific feature type exists in the GFF file
check_feature_exists() {
    local feature_type="$1"
    local count=$(grep -c -P "\t$feature_type\t" "$GFF_FILE")
    if [ "$count" -gt 0 ]; then
        return 0  # Feature type exists
    else
        return 1  # Feature type does not exist
    fi
}

check_and_download "$INPUT_FASTA_FILE" "$FASTA_FILE"
check_and_download "$INPUT_GFF_FILE" "$GFF_FILE"

echo "Adding the reference genome..."
${APBIN}/prepare-refseqs.pl --fasta $FASTA_FILE --out ${APDIR}

if [ $? -ne 0 ]; then
    echo "Failed to add the reference genome..."
    exit 1
fi

echo "Adding the gene models..."
${APBIN}/flatfile-to-json.pl --gff $GFF_FILE --type $GFF_TYPE --trackLabel gene_models --key "Gene models" --out ${APDIR}

if [ $? -ne 0 ]; then
    echo "Failed to add the gene models..."
    exit 1
fi

echo "Index gene names..."
${APBIN}/generate-names.pl --tracks gene_models --out ${APDIR}

if [ $? -ne 0 ]; then
    echo "Failed to add the gene models..."
    exit 1
fi

if [ "$repeats_enabled" = true ]; then
    echo "Adding the repeat tracks..."
    # Check if the "low_complexity_region" feature exists
    if check_feature_exists "low_complexity_region"; then
        echo "Adding the dust.low_complexity_region repeats..."
        ${APBIN}/flatfile-to-json.pl --gff "$GFF_FILE" --trackLabel "dust.low_complexity_region" --key "Low complexity region (Dust)" --type low_complexity_region --out ${APDIR}
        if [ $? -ne 0 ]; then
            echo "Failed to add dust.low_complexity_region repeats."
            exit 1
        fi
    fi

    # Check if the "repeat_region" feature exists
    if check_feature_exists "repeat_region"; then
        echo "Adding the repeatmasker.repeat_region..."
        ${APBIN}/flatfile-to-json.pl --gff "$GFF_FILE" --trackLabel "repeatmasker.repeat_region" --key "Repetitive region (RepeatMasker)" --type repeat_region --out ${APDIR}
        if [ $? -ne 0 ]; then
            echo "Failed to add repeatmasker.repeat_region."
            exit 1
        fi
    fi

    # Check if the "tandem_repeat" feature exists
    if check_feature_exists "tandem_repeat"; then
        echo "Adding the trf.tandem_repeats..."
        ${APBIN}/flatfile-to-json.pl --gff "$GFF_FILE" --trackLabel "trf.tandem_repeat" --key "Tandem repeat (TRF)" --type tandem_repeat --out ${APDIR}
        if [ $? -ne 0 ]; then
            echo "Failed to add trf.tandem_repeats."
            exit 1
        fi
    fi
fi;

if [[ "${RNASEQ_SPECIES_RELEASE}" != false ]]; then
    echo "Adding RNASeq tracks from WBPS $RNASEQ_RELEASE Jbrowse"
    $PYTHON_PATH ${script_dir}/add_wbps_jbrowse_rnaseq_tracks_to_apollo.py -c $RNASEQ_SPECIES -a $APDIR -p $RNASEQ_RELEASE
fi;

# Exiting script successfully
exit 0