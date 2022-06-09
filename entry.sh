#!/bin/bash

py='/Users/vita/anaconda3/envs/dna/bin/python'
TSS='/Source/scripts/data/deal/all_TSS.bed'
DATE=`date +%Y-%m-%d_%H-%M-%S`
NTSS_LENGTH=1000

init(){
    SCRIPTPATH=$(dirname "$SCRIPT")
    cd $SCRIPTPATH
    export PYTHONPATH=$PYTHONPATH:`pwd`/utils

    TEMPDIR="$DATE"
    # TSS_LOW='/Source/scripts/data/deal/all_TSS_low.bed'
    TSS_LOW="./build/$TEMPDIR/none_tss.bed"

    echo -e "\n--${DATE}--\nTSS: ${TSS}\nTSS_LOW: ${TSS_LOW}\nNTSS_LENGTH: ${NTSS_LENGTH}\n" >> ./build/log.txt

    mkdir -p build/$TEMPDIR/features
    mkdir -p build/$TEMPDIR/figures
    
}

run(){
    echo ">>> $1"
    echo "$1" >> ./build/log.txt
    eval $1
}

init
run "$py ./cmd/pick_none_tss.py $TSS -o $TSS_LOW -l $NTSS_LENGTH"
run "$py extract_features.py -t $TSS -n $TSS_LOW -o ./build/$TEMPDIR/features"
run "$py predict.py ./build/$DATE/features/features_df_combined.csv -o ./build/$TEMPDIR/figures/"
# optional
run "$py ./cmd/plot_features.py ./build/$TEMPDIR/features/features_df_combined.csv -o ./build/$TEMPDIR/"
