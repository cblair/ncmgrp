#!/bin/sh

for file in `ls`; do (echo "========$file========" && cat $file) | mail well0358@vandals.uidaho.edu -s "$file contents" ; done
