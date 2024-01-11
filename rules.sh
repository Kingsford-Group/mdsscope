#! /bin/bash

for n in "${@}"; do
    echo ": ${n}.yaggo |> !yaggo |>"
    echo ": ${n}.cc | ${n}.hpp |> !cxx |> {${n}_objs}"
    echo ": {${n}_objs} common.ar |> !lxxd |> ${n}"
done
