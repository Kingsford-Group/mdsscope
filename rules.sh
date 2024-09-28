#! /bin/bash

for n in "${@}"; do
    echo ": ${n}.cc |> !cxx |> {${n}_objs}"
    echo ": {${n}_objs} common.ar |> !lxxd |> ${n}"
done
