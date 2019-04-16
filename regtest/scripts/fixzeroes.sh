#! /bin/bash

file=$1

sed -e "s/-0.000/ 0.000/g" -e "s/-nan/ nan/g" $file > $file.$$.tmp
mv $file.$$.tmp $file

