#!/bin/bash
cd /home/projects/ku_00015/people/tuhu/multiomics-ad-phd
module load tools
module load intel/perflibs/2020_update4 gcc/9.4.0 R/4.1.0
Rscript R/de_time.R
