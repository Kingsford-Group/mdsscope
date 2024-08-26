#!/usr/bin/env bash

K=$1
S=$2
REPEAT=$3
TRANSCRIPTS=$4
shift 4

# Random seeds
for i in $(seq 1 $REPEAT); do
    echo ": ../create_seed |> %f %o |> seed_${i}"
done

sets=("mykkeltveit" "champarnaud")
for set in "${sets[@]}"; do
    f="${set}_set"
    for switch in s c u; do
        echo ": foreach ${@} | ../sketch_histo |> %1i  -a ACGT --${set} -${switch} < %f > %o |> %B_${f}_${switch}.histo {all_histos}"
        if [ -n "$TRANSCRIPTS" ]; then
            echo ": $TRANSCRIPTS | ../sketch_histo |> %1i  -a ACGT --${set} -${switch} < %f > %o |> transcripts_${f}_${switch}.histos {transcripts_histos}"
        fi
    done
done

if [ -n "$S" ]; then
    frac=$(bc <<< "scale=5; 1 / ($K - $S + 1)")
    for i in $(seq 1 $REPEAT); do
        for switch in s c u; do
            for t in 0 1 $((K/2)); do
                f="syncmer_${t}_set_i${i}"
                echo ": foreach ${@} | ../sketch_histo seed_${i} |> %1i  -a ACGT --syncmer ${t} --syncmer-s ${S} -i %2i  -${switch} < %f > %o |> %B_${f}_${switch}.histo {all_histos}"
                if [ -n "$TRANSCRIPTS" ]; then
                    echo ": $TRANSCRIPTS | ../sketch_histo seed_${i} |> %1i  -a ACGT --syncmer ${t} --syncmer-s "$S" -${switch} -i %2i < %f > %o |> transcripts_${f}_${switch}.histos {transcripts_histos}"
                fi

            done

            f="frac_set_i${i}"
            echo ": foreach ${@} | ../sketch_histo seed_${i} |> %1i  -a ACGT --frac ${frac} -${switch} -i %2i < %f > %o |> %B_${f}_${switch}.histo {all_histos}"
            if [ -n "$TRANSCRIPTS" ]; then
                echo ": $TRANSCRIPTS | ../sketch_histo seed_${i} |> %1i  -a ACGT --frac ${frac} -${switch} -i %2i < %f > %o |> transcripts_${f}_${switch}.histos {transcripts_histos}"
            fi
        done
    done
fi
