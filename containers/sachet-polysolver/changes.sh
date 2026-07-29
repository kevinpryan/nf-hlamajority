mkdir /data

sed -i 's.#!/bin/sh.#!/bin/bash.g' /home/polysolver/scripts/shell_call_hla_type
sed -i 's.\$SAMTOOLS_DIR./home/polysolver/binaries.g' /home/polysolver/scripts/shell_call_hla_type
sed -i 's.6:29941260-29945884.chr6:29941260-29945884.g' /home/polysolver/scripts/shell_call_hla_type
sed -i 's.6:31353872-31357187.chr6:31353872-31357187.g' /home/polysolver/scripts/shell_call_hla_type
sed -i 's.6:31268749-31272105.chr6:31268749-31272105.g' /home/polysolver/scripts/shell_call_hla_type

# change /home/polysolver/binaries to $NOVOALIGN_DIR in align_fork_fh.pl call
sed -i 's|\(align_fork_fh.pl.* 0 \)/home/polysolver/binaries|\1\$NOVOALIGN_DIR|' \
    /home/polysolver/scripts/shell_call_hla_type

# add NOVOALIGN_DIR variable to shell_call_hla_type
sed -i '1a NOVOALIGN_DIR=$(dirname "$(command -v novoalign)")' \
    /home/polysolver/scripts/shell_call_hla_type

sed -i 's|TMP_DIR=/home/polysolver|TMP_DIR=$outDir|' /home/polysolver/scripts/shell_call_hla_type

rm -f /home/polysolver/binaries/novoalign

