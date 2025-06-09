SINGULARITY_CACHE=$1
bindir=$(pwd)
cd $SINGULARITY_CACHE
singularity pull docker://kevinr9525/cancerit-kourami:latest
singularity pull quay.io/biocontainers/bwakit:0.7.18.dev1--hdfd78af_0
singularity pull docker://kevinr9525/genome-seek_hla:v1.0.4

cd ${bindir}/../references
git clone https://github.com/Kingsford-Group/kourami.git
cd kourami/scripts
bash download_grch38.sh hs38NoAltDH
cd ../resources/
singularity run -B $(pwd) ${SINGULARITY_CACHE}/cancerit-kourami_latest.sif bwa index hs38NoAltDH.fa

