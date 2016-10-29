#!/bin/sh
function MakeCPU(){
    echo ""
    echo ""
    echo ""
    echo "Compile CPU code ..."
    echo ""
    echo ""
    echo ""
    cd ./CPU/clustalw-2.0.11
    ./configure && make
    cd ../../
}


function MakeGPU(){
    echo ""
    echo ""
    echo ""
    echo "Compile GPU code ..."
    echo ""
    echo ""
    echo ""
    cd ./GPU/cuda_clustalw_final
    ./configure && make
    cd ../../
}

function MakeClean(){
    cd ./CPU/clustalw-2.0.11
    make clean
    cd ../../

    cd ./GPU/cuda_clustalw_final
    make clean
    cd ../../
}

if [ $# == 0 ]; then
    MakeCPU;
    MakeGPU;
else 
    if [ "$1" == "clean" ]; then
        MakeClean;
    elif [ "$1" == "GPU" ]; then
        MakeGPU;
    elif [ "$1" == "CPU" ]; then
        MakeCPU;
    fi
fi

