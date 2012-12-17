#!/bin/sh

# raf predict ...
raf2bpseq.py $1.raf > $1.bpseq
make_coords $1.bpseq $1.coords
plot_rna $1.bpseq $1.coords --png $1.png
open $1.png
