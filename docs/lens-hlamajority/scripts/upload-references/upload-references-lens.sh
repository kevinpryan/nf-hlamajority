#!/usr/bin/env bash
set -euo pipefail

BASE_DIR="/path-to/lens-install/references"
BUCKET="s3://bucket-dir/references/cloud"
DRYRUN="--dryrun"   # REMOVE after validation
#DRYRUN=""

while read -r path; do
    [[ -z "$path" ]] && continue

    clean="${path%/}"

    case "$clean" in

        # ----------------------------
        # Flatten into ROOT
        # ----------------------------
        dummy_file|\
        homo_sapiens/annot|\
        homo_sapiens/fasta|\
        homo_sapiens/beds|\
        homo_sapiens/cta_self|\
        homo_sapiens/erv|\
        homo_sapiens/lohhla|\
        homo_sapiens/protein|\
        homo_sapiens/vcfs|\
        viral)
            dest="$BUCKET/"
            ;;

        # ----------------------------
        # Explicit renames
        # ----------------------------
        homo_sapiens/binfotron)
            dest="$BUCKET/binfotron/"
            ;;

        homo_sapiens/neosplice/peptidome.homo_sapiens)
            dest="$BUCKET/peptidome.homo_sapiens/"
            ;;

        homo_sapiens/snaf/snaf-data)
            dest="$BUCKET/snaf-data/"
            ;;

        homo_sapiens/snpeff/GRCh38.GENCODEv37)
            dest="$BUCKET/GRCh38.GENCODEv37/"
            ;;

        homo_sapiens/starfusion/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play)
            dest="$BUCKET/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/"
            ;;

        # ----------------------------
        # Keep top-level dirs
        # ----------------------------
        antigen.garnish)
            dest="$BUCKET/antigen.garnish/"
            ;;

        mhcflurry)
            dest="$BUCKET/mhcflurry/"
            ;;

        *)
            echo "Unknown path: $clean"
            exit 1
            ;;
    esac

    echo "--------------------------------------"
    echo "SOURCE : $clean"
    echo "DEST   : $dest"
    echo "--------------------------------------"

    aws s3 sync \
        "$BASE_DIR/$clean/" \
        "$dest" \
        $DRYRUN

done < dirs-cp-contents.txt

aws s3 cp \
       dummy_file \
       $BUCKET/ \
       $DRYRUN
