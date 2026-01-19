#
# Part of Kourami HLA typer/assembler
# (c) 2017 by  Heewook Lee, Carl Kingsford, and Carnegie Mellon University.
# See LICENSE for licensing.
#

#!/bin/bash

set -e # Exit immediately if a command exits with a non-zero status.

echo "#---------------------------------------------------------------------"
echo "# This script is adopted from run-gen-ref in bwa.kit 0.7.12 by Heng Li"
echo "# available from https://github.com/lh3/bwa/tree/master/bwakit"
echo "#"
echo "# run this to download GRCh38 reference."
echo "#"
echo "#---------------------------------------------------------------------"
echo

# --- Determine script directory and set default resources_dir ---
pushd `dirname $0` > /dev/null
SCRIPTD=`pwd`
popd > /dev/null
default_resources_dir="$SCRIPTD/../resources"
resources_dir="$default_resources_dir" # Initialize with default

# --- URL definitions ---
url38NoAltDecoy="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
url38d="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.gz"
url38="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"

# --- Function to print usage ---
function print_usage {
    echo "Usage: $0 [options] <hs38|hs38D|hs38DH|hs38NoAltDH>"
    echo ""
    echo "Options:"
    echo "  -r, --resources-dir DIR   Specify the directory for resource files."
    echo "                            Defaults to: \"$default_resources_dir\""
    echo "  -h, --help                Display this help message and exit."
    echo ""
    echo "The file containing the HLA sequences is located at <resources_dir>/hla_bwa.kit.fna.gz"
    echo ""
    echo "Analysis sets:"
    echo "  hs38         primary assembly of GRCh38 (incl. chromosomes, unplaced and unlocalized contigs) and EBV"
    echo "  hs38D        hs38 + decoy contigs"
    echo "  hs38DH       hs38 + ALT contigs + decoy contigs + HLA genes (recommended for GRCh38 mapping)"
    echo "  hs38NoAltDH  hs38 + decoy contigs + HLA alleles from [bwa.kit] (recommended for HLA typing)"
    echo ""
    echo "Note: This script downloads human reference genomes from NCBI ftp server."
    exit 1;
}

# --- Parse command line options ---
# Handling long options and short options
# Adapted from https://stackoverflow.com/a/28466267/1930608
TEMP=$(getopt -o hr: --long help,resources-dir: -n "$0" -- "$@")
if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi
eval set -- "$TEMP" # Note the quotes
unset TEMP

while true; do
  case "$1" in
    -h | --help )
      print_usage
      exit 0
      ;;
    -r | --resources-dir )
      resources_dir="$2"
      shift 2 # Past argument and value
      ;;
    -- ) # No more options
      shift
      break
      ;;
    * ) # No more options or an error
      break
      ;;
  esac
done

# --- Check for required positional argument (genome build) ---
if [ $# -ne 1 ]; then
    echo "ERROR: Genome build argument is missing or too many arguments provided."
    echo
    print_usage
fi

genome_build="$1"

# --- Create resources directory if it doesn't exist ---
# This is good practice, especially if a custom path is given.
if [ ! -d "$resources_dir" ]; then
    echo "Creating resources directory: $resources_dir"
    mkdir -p "$resources_dir"
    if [ $? -ne 0 ]; then
        echo "ERROR: Could not create resources directory: $resources_dir"
        exit 1
    fi
fi

# --- Main logic for downloading and preparing files ---
target_fa_file="$resources_dir/$genome_build.fa"
hla_file="$resources_dir/hla_bwa.kit.fna.gz" # Define once

if [ ! -e "$target_fa_file" ]; then
    echo "Target file $target_fa_file does not exist. Proceeding with download..."
    case "$genome_build" in
        "hs38DH")
            echo "Preparing hs38DH..."
            if [ ! -e "$hla_file" ]; then
                echo "ERROR: HLA file $hla_file not found in resources directory!"
                echo "Please ensure it's available before running for hs38DH or hs38NoAltDH."
                exit 1
            fi
            (wget --progress=bar:force:noscroll -O- "$url38d" | gzip -dc; gzip -dc "$hla_file") > "$target_fa_file"
            ;;
        "hs38")
            echo "Preparing hs38..."
            (wget --progress=bar:force:noscroll -O- "$url38" | gzip -dc) > "$target_fa_file"
            ;;
        "hs38D")
            echo "Preparing hs38D..."
            (wget --progress=bar:force:noscroll -O- "$url38d" | gzip -dc) > "$target_fa_file"
            ;;
        "hs38NoAltDH")
            echo "Preparing hs38NoAltDH..."
            if [ ! -e "$hla_file" ]; then
                echo "ERROR: HLA file $hla_file not found in resources directory!"
                echo "Please ensure it's available before running for hs38DH or hs38NoAltDH."
                exit 1
            fi
            (wget --progress=bar:force:noscroll -O- "$url38NoAltDecoy" | gzip -dc; gzip -dc "$hla_file") > "$target_fa_file"
            ;;
        *)
            echo "ERROR: unknown genome build: $genome_build"
            echo
            print_usage
            ;;
    esac
    echo "Successfully created $target_fa_file."
else
    echo "Target file $target_fa_file already exists. Skipping download."
fi

if [ ! -f "$target_fa_file.bwt" ]; then
    echo -e "\nPlease run 'bwa index \"$target_fa_file\"' to create the BWA index."
else
    echo "BWA index for $target_fa_file already exists."
fi

echo "Script finished."
