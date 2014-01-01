#!/bin/bash

time nice blastn \
  -db /home/db/ncbi/nt \
  -query Pelobates_cultripes.fasta \
  -out Pelobates_cultripes-v-nt.blastn \
  -outfmt '7 std sstrand qlen slen staxids salltitles' \
  -evalue 1e-5 \
  -soft_masking true &
