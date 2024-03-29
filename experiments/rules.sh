#!/bin/bash

K=$1
S=$((K / 2 - 1))
ILP_PYTHON=$2
REPEAT=$3
TRANSCRIPTS=$4
shift 4

# Random seeds
for i in $(seq 1 $REPEAT); do
    echo ": ../create_seed |> %f %o |> seed_${i}"
done

sets=("mykkeltveit_set" "champarnaud_set")
[ -n "$ILP_PYTHON" ] && sets+=("ilp_set")
for f in "${sets[@]}"; do
    for switch in s c u; do
        echo ": ${f} | ../sketch_components |> %1i -f %f -${switch} > %o |> ${f}_${switch}.scc {all_sccs}"
        if [ -n "$TRANCRIPTS" ]; then
            echo ": $TRANSCRIPTS | ${f} ../sketch_histo |> %2i -a ACGT -f %1i -${switch} < %f > %o |> transcripts_${f}_${switch}.histos {transcripts_histos}"
        fi
    done
done

frac=$(bc <<< "scale=5; 1 / ($K - $S + 1)")
for i in $(seq 1 $REPEAT); do
    for t in 0 1 $((K/2)); do
        f="syncmer_${t}_set_i${i}"
        echo ": ../syncmer_set | seed_${i} |> %f -i %1i -s ${S} -t ${t} > %o |> ${f}"
        for switch in s c u; do
            echo ": ${f} | ../sketch_components |> %1i -f %f -${switch} > %o |> ${f}_${switch}.scc {all_sccs}"
            if [ -n "$TRANSCRIPTS" ]; then
                echo ": $TRANSCRIPTS | ${f} ../sketch_histo seed_${i} |> %2i -a ACGT -f %1i -i %3i -${switch} < %f > %o |> transcripts_${f}_${switch}.histos {transcripts_histos}"
            fi
        done
        sets+=($f)
    done

    f="frac_set_i${i}"
    echo ": ../frac_set | seed_${i} |> %f -i %1i -f ${frac} > %o |> ${f}"
    sets+=($f)
    for switch in s c u; do
        echo ": ${f} | ../sketch_components |> %1i -f %f -${switch} > %o |> ${f}_${switch}.scc {all_sccs}"
        if [ -n "$TRANCRIPTS" ]; then
            echo ": $TRANSCRIPTS | ${f} ../sketch_histo seed_${i} |> %2i -a ACGT -f %1i -i %3i -${switch} < %f > %o |> transcripts_${f}_${switch}.histos {transcripts_histos}"
        fi
    done
done

for f in "${sets[@]}"; do
    for switch in s c u; do
        echo ": foreach ${@} | ${f} ../sketch_histo |> %2i -a ACGT -f %1i -${switch} < %f > %o |> %B_${f}_${switch}.histo {all_histos}"
    done
done
