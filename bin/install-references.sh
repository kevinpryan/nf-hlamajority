SINGULARITY_CACHE=$1
bindir=$(pwd)
cd $SINGULARITY_CACHE

# PULL SINGULARITY IMAGES
singularity pull docker://kevinr9525/cancerit-kourami:latest
singularity pull quay.io/biocontainers/bwakit:0.7.18.dev1--hdfd78af_0
singularity pull docker://kevinr9525/genome-seek_hla:v1.0.4

# BUILD/DOWNLOAD KOURAMI REFERENCES

cd ${bindir}/../references
git clone https://github.com/Kingsford-Group/kourami.git
cd kourami/scripts
bash download_grch38.sh hs38NoAltDH
cd ../resources/
singularity run -B $(pwd) ${SINGULARITY_CACHE}/cancerit-kourami_latest.sif bwa index hs38NoAltDH.fa
cd ../
singularity run -B $(pwd) ${SINGULARITY_CACHE}/cancerit-kourami_latest.sif bash scripts/download_panel.sh

# BUILD/DOWNLOAD BWAKIT REFERENCES
cd ${bindir}/../references/
mkdir -p bwakit; cd bwakit
singularity run -B $(pwd) ${SINGULARITY_CACHE}/quay.io-biocontainers-bwakit-0.7.18.dev1--hdfd78af_0.img run-gen-ref hs38DH
singularity run -B $(pwd) ${SINGULARITY_CACHE}/quay.io-biocontainers-bwakit-0.7.18.dev1--hdfd78af_0.img bwa index hs38DH.fa

# BUILD/DOWNLOAD HLA-LA REFERENCES
cd ${bindir}/../references/
mkdir -p hla-la
wget http://www.well.ox.ac.uk/downloads/PRG_MHC_GRCh38_withIMGT.tar.gz
tar -xvzf PRG_MHC_GRCh38_withIMGT.tar.gz
singularity run -B $(pwd) ${SINGULARITY_CACHE}/kevinr9525-genome-seek_hla-v1.0.4.img HLA-LA --action prepareGraph --PRG_graph_dir PRG_MHC_GRCh38_withIMGT
