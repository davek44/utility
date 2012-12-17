#!/bin/sh

contrafold predict $1 --bpseq $1.bpseq --parens $1.parens
make_coords $1.bpseq $1.coords
plot_rna $1.bpseq $1.coords --png $1.png
open $1.png
