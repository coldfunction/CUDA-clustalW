#!/bin/sh
TEST_DATA_DIR_CPU=./test_data/cpu
TEST_DATA_DIR_GPU=./test_data/gpu
if [ $# == 2 ]; then
diff $TEST_DATA_DIR_CPU/$1_$2.aln $TEST_DATA_DIR_GPU/$1_$2.aln 
else
echo "$0 [number of sequences] [length of sequences]"
echo ""
echo "Example:"
echo "    $0 100  97"
echo "    $0 100  498"
echo "    $0 100  1002"
echo "    $0 100  1523"
echo "    $0 1000 97"
echo "    $0 1000 498"
echo "    $0 1000 1002"
echo "    $0 1000 1523"
echo ""
fi

