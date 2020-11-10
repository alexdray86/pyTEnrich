#!/bin/bash

file=$1

tail -n +4 $file | tr -s " " | sed -e 's/^[ \t]*//' | tr " " "\t" | cut -f5,6,7,9,10,11 | awk '{ if($4=="C"){ $4="-" } ; print $1, $2, $3, "_", "_", $4, $6, $5, "NA"}' | tr " " "\t" | sed 's/\bC/-/g'
