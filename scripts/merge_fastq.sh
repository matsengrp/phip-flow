#!/bin/bash

cd $1
for file in *R1*; do
  [[ ! -f $file ]] && continue      # pick up only regular files

  otherfile="$2/$file"
  [[ ! -f $otherfile ]] && continue # skip if there is no matching file in folder 2
  zcat "$file" "$otherfile" > "$3/$file"
done
