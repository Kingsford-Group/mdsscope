#!/bin/bash

K=$1
# S=$2
S=$((K / 2 - 1))
REPEAT=$3
shift 3

# Random seeds
for i in $(seq 1 $REPEAT); do
    echo ": ../create_seed |> %f %o |> seed_${i}"
done

sets=("mykkeltveit_set" "champarnaud_set")
for f in "${sets[@]}"; do
    for switch in s c u; do
        echo ": ${f} | ../sketch_components |> %1i -f %f -${switch} > %o |> ${f}_${switch}.scc {all_sccs}"
    done
done

frac=$(bc <<< "scale=5; 1 / ($K - $S + 1)")
for i in $(seq 1 $REPEAT); do
    for t in 0 1 $((K/2)); do
        f="syncmer_${t}_set_i${i}"
        echo ": ../syncmer_set | seed_${i} |> %f -i %1i -s ${S} -t ${t} > %o |> ${f}"
        for switch in s c u; do
            echo ": ${f} | ../sketch_components |> %1i -f %f -${switch} > %o |> ${f}_${switch}.scc {all_sccs}"
        done
        sets+=($f)
    done

    f="frac_set_i${i}"
    echo ": ../frac_set | seed_${i} |> %f -i %1i -f ${frac} > %o |> ${f}"
    sets+=($f)
    for switch in s c u; do
        echo ": ${f} | ../sketch_components |> %1i -f %f -${switch} > %o |> ${f}_${switch}.scc {all_sccs}"
    done
done

for f in "${sets[@]}"; do
    for switch in s c u; do
        echo ": foreach ${@} | ${f} ../sketch_histo |> %2i -a ACGT -f %1i -${switch} < %f > %o |> %B_${f}_${switch}.histo {all_histos}"
    done
done
