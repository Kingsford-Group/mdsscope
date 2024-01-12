#!/bin/bash

K=$1
S=$2
REPEAT=$3
shift 3

# Cross product: set x input files
for f in mykkeltveit_set champarnaud_set; do
    echo ": foreach ${@} | ${f} ../sketch_histo |> %2i -a ACGT -f %1i < %f > %o |> %B_${f}.histo"
    echo ": foreach ${@} | ${f} ../sketch_histo |> %2i -a ACGT -f %1i -c < %f > %o |> %B_${f}_c.histo"
done

# Random seeds
for i in $(seq 1 $REPEAT); do
    echo ": ../create_seed |> %f %o |> seed_${i}"
done

# Syncmer
# Cross product: seed x t values
frac=$(bc <<< "scale=5; 1 / ($K - $S + 1)")
for i in $(seq 1 $REPEAT); do
    for t in 0 1 $((K/2)); do
        echo ": ../syncmer_set | seed_${i} |> %f -i %1i -s 3 -t ${t} > %o |> syncmer_set_i${i}_t${t}"
        echo ": foreach ${@} | syncmer_set_i${i}_t${t} ../sketch_histo |> %2i -a ACGT -f %1i < %f > %o |> %B_syncmer_set_i${i}_t${t}.histo"
        echo ": foreach ${@} | syncmer_set_i${i}_t${t} ../sketch_histo |> %2i -a ACGT -f %1i -c < %f > %o |> %B_syncmer_set_i${i}_t${t}_c.histo"
    done

    echo ": ../frac_set | seed_${i} |> %f -i %1i -f ${frac} > %o |> %B_i${i}"
    echo ": foreach ${@} | frac_set_i${i} ../sketch_histo |> %2i -a ACGT -f %1i < %f > %o |> %B_frac_set_i${i}.histo"
    echo ": foreach ${@} | frac_set_i${i} ../sketch_histo |> %2i -a ACGT -f %1i -c < %f > %o |> %B_frac_set_i${i}_c.histo"
done
