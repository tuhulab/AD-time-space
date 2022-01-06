#!/bin/bash
module load ascp/3.9.6
cd /home/projects/ku_00015/data/rna-seq
files_to_upload=$( cat /home/projects/ku_00015/people/tuhu/multiomics-ad-phd/data/geo/files_to_upload.txt )
for F in $files_to_upload
do 
    ascp -i sra-1.ssh.priv -QT -l 500m -k1 read/$F asp-sra@upload.ncbi.nlm.nih.gov:incoming
    echo $F is done
done