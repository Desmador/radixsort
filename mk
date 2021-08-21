#!/usr/bin/env tcsh

if ($# != 1) then
	echo "usage: mk <binary_name>"
	exit -1
endif

# we remove ending (just to be sure)
set binary_name=$1:r

gcc -D_GNU_SOURCE -Wall -std=c99 -O3 -o ${binary_name} ${binary_name}.c -pthread