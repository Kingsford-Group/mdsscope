#!/bin/bash

K=$1
S=$2
REPEAT=$3
shift 3

# Random seeds
for i in $(seq 1 $REPEAT); do
    echo ": ../create_seed |> %f %o |> seed_${i}"
done

sets=("mykkeltveit_set" "champarnaud_set")
frac=$(bc <<< "scale=5; 1 / ($K - $S + 1)")
for i in $(seq 1 $REPEAT); do
    for t in 0 1 $((K/2)); do
        f="syncmer_${t}_set_i${i}"
        echo ": ../syncmer_set | seed_${i} |> %f -i %1i -s 3 -t ${t} > %o |> ${f}"
        sets+=($f)
    done

    f="frac_set_i${i}"
    echo ": ../frac_set | seed_${i} |> %f -i %1i -f ${frac} > %o |> ${f}"
    sets+=($f)
done

for f in "${sets[@]}"; do
    cat <<EOF
: foreach ${@} | ${f} ../sketch_histo |> %2i -a ACGT -f %1i < %f > %o |> %B_${f}_s.histo {all_histos}
: foreach ${@} | ${f} ../sketch_histo |> %2i -a ACGT -f %1i -c < %f > %o |> %B_${f}_c.histo {all_histos}
: foreach ${@} | ${f} ../sketch_histo |> %2i -a ACGT -f %1i -u < %f > %o |> %B_${f}_u.histo {all_histos}
: foreach ${@} | ${f} ../sketch_components |> %2i -f %1i > %o |> %B_${f}.scc {all_sccs}
EOF
done

echo ": {all_histos} |> ./compute_stats --histos %f > %o 2>&1 |> histos"
echo ": {all_sccs} |> ./compute_stats --sccs %f > %o 2>&1 |> sccs"
