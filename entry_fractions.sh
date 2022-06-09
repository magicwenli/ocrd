#!/bin/bash
set -e

py='/Users/vita/anaconda3/envs/dna/bin/python'
TSS='/Source/scripts/data/HK.all.txt.bed'
DATE=`date +%Y-%m-%d_%H-%M-%S`
NTSS_LENGTH=1000

init(){
    SCRIPTPATH=$(dirname "$SCRIPT")
    cd $SCRIPTPATH
    export PYTHONPATH=$PYTHONPATH:`pwd`/utils

    TEMPDIR="$DATE"
    # TSS_LOW='/Source/scripts/build/2022-04-23_23-22-13/none_tss.bed'
    TSS_LOW="./build/$TEMPDIR/none_tss.bed"
    fraction=(0.01 0.02 0.05 0.15 0.30 1.0)

    echo -e "\n--${DATE}--\nTSS: ${TSS}\nTSS_LOW: ${TSS_LOW}\nNTSS_LENGTH: ${NTSS_LENGTH}\n" >> ./build/log.txt

    for f in ${fraction[@]}
    do
        mkdir -p build/$TEMPDIR/features_${f}
        mkdir -p build/$TEMPDIR/figures_${f}
    done
}

run(){
    echo ">>> $1"
    echo "$1" >> ./build/log.txt
    eval $1
}

init
run "$py ./cmd/pick_none_tss.py $TSS -o $TSS_LOW -l $NTSS_LENGTH"

for f in ${fraction[@]}
do
    echo "doing $TSS with fraction $f"
    run "$py extract_features.py -t $TSS -n $TSS_LOW -f $f -o ./build/$TEMPDIR/features_${f}"
    run "$py predict.py ./build/$DATE/features_${f}/features_df_combined.csv -o ./build/$TEMPDIR/figures_${f}/ -f $f"
done
# optional
# run "$py ./cmd/plot_features.py ./build/$TEMPDIR/features_${f}/features_df_combined.csv -o ./build/$TEMPDIR/"
