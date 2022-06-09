#!/bin/bash
script_dir=`dirname "$0"`
for file in `find $script_dir/../dynamics/spam -name "*.[h|cpp]" `; do 
        #echo $file
	clang-format --style=llvm -i $file
done
