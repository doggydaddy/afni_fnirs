#!/bin/bash
input_csv=$1

sed -i 's/,lrv,/,PSD+AGG,/g' $input_csv
sed -i 's/,pcon,/,PSD-AGG,/g' $input_csv
sed -i 's/,con,/,HC,/g' $input_csv
