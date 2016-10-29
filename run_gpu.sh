#!/bin/sh
TEST_DATA_DIR=./test_data/cpu
EXECUTE_FILE=./GPU/cuda_clustalw_final/src/clustalw2

if [ $# == 2 ]; then
echo ""
echo ""
echo "============================================================================"
echo "Start run $TEST_DATA_DIR/$1_$2.fa, and will output $TEST_DATA_DIR/$1_$2.aln"
$EXECUTE_FILE -INFILE=./$TEST_DATA_DIR/$1_$2.fa -ALIGN -QUIET

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
