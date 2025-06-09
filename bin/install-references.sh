#!/bin/bash

# to run: go to bin directory
# bash install_references.sh /path/to/dir/with/singularity/images
# last command (HLA-LA prepareGraph) can take a few hours and take up to 40G of memory

SINGULARITY_CACHE=$1
bindir=$(pwd)
cd $SINGULARITY_CACHE

# PULL SINGULARITY IMAGES
singularity pull docker://kevinr9525/cancerit-kourami:latest
singularity pull quay.io/biocontainers/bwakit:0.7.18.dev1--hdfd78af_0
singularity pull docker://kevinr9525/genome-seek_hla:v1.0.4
singularity pull https://depot.galaxyproject.org/singularity/bwa%3A0.7.18--he4a0461_1

# BUILD/DOWNLOAD KOURAMI REFERENCES
mkdir -p ${bindir}/../references; cd ${bindir}/../references
git clone https://github.com/Kingsford-Group/kourami.git
cd kourami/scripts
bash download_grch38.sh hs38NoAltDH
cd ../resources/
singularity run -B $(pwd) ${SINGULARITY_CACHE}/cancerit-kourami_latest.sif bwa index hs38NoAltDH.fa
cd ../
singularity run --cleanenv -B $(pwd) ${SINGULARITY_CACHE}/bwa%3A0.7.18--he4a0461_1 bash scripts/download_panel.sh

# BUILD/DOWNLOAD BWAKIT REFERENCES
cd ${bindir}/../references/
mkdir -p bwakit; cd bwakit
singularity run -B $(pwd) ${SINGULARITY_CACHE}/quay.io-biocontainers-bwakit-0.7.18.dev1--hdfd78af_0.img run-gen-ref hs38DH
singularity run -B $(pwd) ${SINGULARITY_CACHE}/quay.io-biocontainers-bwakit-0.7.18.dev1--hdfd78af_0.img bwa index hs38DH.fa

# BUILD/DOWNLOAD HLA-LA REFERENCES
cd ${bindir}/../references/
mkdir -p hla-la; cd hla-la
wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz
tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz
rm PRG_MHC_GRCh38_withIMGT.tar.gz
singularity run -B $(pwd) ${SINGULARITY_CACHE}/kevinr9525-genome-seek_hla-v1.0.4.img HLA-LA --action prepareGraph --PRG_graph_dir PRG_MHC_GRCh38_withIMGT
