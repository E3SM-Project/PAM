#!/bin/bash
script_dir=`dirname "$0"`
for file in `find $script_dir/../dynamics/spam -type f \( -name \*.h -o -name \*.cpp \) -not -path "*/build/*"`; do
        #echo $file
	clang-format --style=llvm -i $file
done
