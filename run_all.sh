#!/bin/bash

start_th=1
end_th=10
out_file="stats.tsv"
input_file="./gen_data.txt"


# проверяет существование входного файла
if [ ! -f "$input_file" ]; then
    echo "Ошибка: входной файл $input_file не найден" >&2
    exit 1
fi

> $out_file; # почистили

for th in $(seq $start_th $end_th ); do
    echo "th=$th... "
    result=$(./n_body --th "$th" --file "$input_file")
    echo -e "$th\t$result" >> "$out_file"
done;