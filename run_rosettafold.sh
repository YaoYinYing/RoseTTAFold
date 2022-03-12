#!/bin/bash

# make the script stop when error (non-true exit code) is occured
set -e

usage() {
        echo ""
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        echo "-i <input fasta>          Path to directory of supporting data"
        echo "Optional Parameters:"
        echo "-o <output_dir>           Path to a directory that will store the results."
        echo "-j <nproc>                Maximun number of proccessor to be used, 8 cores for default."
        echo "-m <memory>               Maximum volume memory to be used, 64 GB for default"
        echo "-t <type>                 Run type for RoseTTAFold, pyrosetta(default) or e2e "
        exit 1
}

while getopts ":i:o:j:t:m:" x; do
        case "${x}" in
        i)
                input_fasta=$OPTARG
        ;;
        o)
                output_dir=$OPTARG
        ;;
        j)
                nproc=$OPTARG
        ;;
        m)
                memory=$OPTARG
        ;;
        t)
                type=$OPTARG
        ;;
        *)
                echo What the hell?
        ;;
        esac
done

############################################################
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
# <<< conda initialize <<<
############################################################

banner="############################################################"
SCRIPT=$(dirname $0)
export PIPEDIR=$(readlink -f "$SCRIPT")

if [[ "$input_fasta" == "" ]] ; then
    usage
fi

IN=$(readlink -f "$input_fasta")
IN_DIR=$(dirname "$input_fasta")
IN_DIR=$(readlink -f "$IN_DIR")
IN_fasta=$(basename  "$input_fasta")               # input.fasta

if [[ "$output_dir" == "" ]]; then
    WDIR=$IN_DIR/${IN_fasta%.fasta}
else
    WDIR=`realpath -s $output_dir`
fi
echo Use WDIR: $WDIR

if [[ "$nproc" == "" ]] ; then
    nproc=16
fi

if [[ "$memory" == "" ]] ; then
    memory=64
fi

if [[ "$type" == "" ]] ; then
    type=pyrosetta
fi
if [[ "$type" != "pyrosetta" && "$type" != "e2e" ]] ; then
    usage
fi


CPU="$nproc"  # number of CPUs to use
MEM="$memory" # max memory (in GB)


LEN=`tail -n1 $IN | wc -m`

mkdir -p $WDIR


conda activate RoseTTAFold
############################################################
# 1. generate MSAs
############################################################
echo $banner
if [ ! -s $WDIR/t000_.msa0.a3m ]
then
    echo "Running HHblits"
    date
    cmd="$PIPEDIR/input_prep/make_msa.sh $IN $WDIR $CPU $MEM"
    echo "$cmd"
    eval "$cmd"
fi
echo $banner

############################################################
# 2. predict secondary structure for HHsearch run
############################################################
echo $banner
if [ ! -s $WDIR/t000_.ss2 ]
then
    echo "Running PSIPRED"
    cmd="$PIPEDIR/input_prep/make_ss.sh $WDIR/t000_.msa0.a3m $WDIR/t000_.ss2"
    echo "$cmd"
    eval "$cmd"
fi
echo $banner

############################################################
# 3. search for templates
############################################################
DB="$(/bin/sh $PIPEDIR/find_my_db.sh 'pdb100')"
echo $banner
if [ ! -s $WDIR/t000_.hhr ]
then
    echo "Running hhsearch"
    date
    HH="hhsearch -b 50 -B 500 -z 50 -Z 500 -mact 0.05 -cpu $CPU -maxmem $MEM -aliw 100000 -e 100 -p 5.0 -d $DB"
    cat $WDIR/t000_.ss2 $WDIR/t000_.msa0.a3m > $WDIR/t000_.msa0.ss2.a3m
    cmd="$HH -i $WDIR/t000_.msa0.ss2.a3m -o $WDIR/t000_.hhr -atab $WDIR/t000_.atab -v 0"
    echo "$cmd"
    eval "$cmd"
fi
echo $banner

WEIGHTS="$(/bin/sh $PIPEDIR/find_my_db.sh 'weight')"
if [[ "$type" == 'e2e' ]];then
    ############################################################
    # 4a. end-to-end prediction
    ############################################################
    echo $banner
    if [ ! -s $WDIR/t000_.3track.npz ]
    then
        echo "Running end-to-end prediction"
        date
        cmd="python $PIPEDIR/network/predict_e2e.py \
        -m $WEIGHTS \
        -i $WDIR/t000_.msa0.a3m \
        -o $WDIR/t000_.e2e \
        --hhr $WDIR/t000_.hhr \
        --atab $WDIR/t000_.atab \
        --db $DB"
    echo "$cmd"
    eval "$cmd"
    fi

elif [[ "$type" == "pyrosetta" ]]; then
    ############################################################
    # 4b. predict distances and orientations
    ############################################################
    echo $banner
    if [ ! -s $WDIR/t000_.3track.npz ]
    then
        echo "Predicting distance and orientations"
        cmd="python $PIPEDIR/network/predict_pyRosetta.py \
            -m $WEIGHTS \
            -i $WDIR/t000_.msa0.a3m \
            -o $WDIR/t000_.3track \
            --hhr $WDIR/t000_.hhr \
            --atab $WDIR/t000_.atab \
            --db $DB"
        echo "$cmd"
        eval "$cmd"
    fi
    echo $banner
    ############################################################
    # 5. perform modeling
    ############################################################
    echo $banner
    mkdir -p $WDIR/pdb-3track
    conda deactivate
    conda activate folding

    for m in 0 1 2
    do
        for p in 0.05 0.15 0.25 0.35 0.45
        do
            for ((i=0;i<1;i++))
            do
                if [ ! -f $WDIR/pdb-3track/model${i}_${m}_${p}.pdb ]; then
                    echo "python -u $PIPEDIR/folding/RosettaTR.py --roll -r 3 -pd $p -m $m -sg 7,3 $WDIR/t000_.3track.npz $IN $WDIR/pdb-3track/model${i}_${m}_${p}.pdb"
                fi
            done
        done
    done > $WDIR/parallel.fold.list

    N=`cat $WDIR/parallel.fold.list | wc -l`
    if [ "$N" -gt "0" ]; then
        echo "Running parallel RosettaTR.py"
        cmd="parallel -j $CPU < $WDIR/parallel.fold.list"
        echo "$cmd"
        eval "$cmd"
    fi
    echo $banner
    ############################################################
    # 6. Pick final models
    ############################################################
    echo $banner
    count=$(find $WDIR/pdb-3track -maxdepth 1 -name '*.npz' | grep -v 'features' | wc -l)
    if [ "$count" -lt "15" ]; then
        # run DeepAccNet-msa
        echo "Running DeepAccNet-msa"
        cmd="python $PIPEDIR/DAN-msa/ErrorPredictorMSA.py --roll -p $CPU $WDIR/t000_.3track.npz $WDIR/pdb-3track $WDIR/pdb-3track"
        echo "$cmd"
        eval "$cmd"
    fi

    if [ ! -s $WDIR/model/model_5.crderr.pdb ]
    then
        echo "Picking final models"
        cmd="python -u -W ignore $PIPEDIR/DAN-msa/pick_final_models.div.py \
            $WDIR/pdb-3track $WDIR/model $CPU"
        echo "$cmd"
        eval "$cmd"
        echo "Final models saved in: $2/model"
    fi

fi
echo "Done"
echo $banner