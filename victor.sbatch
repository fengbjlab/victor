#! /usr/bin/env bash

CURRENT_DIR=`pwd`

if [ $# -lt 1 ] ; then
  echo "This script submit one job to a computer cluster."
  echo "Usage: $0 [slurm_options] <your_slurm_script>"
  exit
elif [ $# -eq 1 ] ; then
  SCRIPT=$1
  OPTION=
  HAS_ARRAY=0
elif [ $# -eq 2 ] ; then
  SCRIPT=$2
  OPTION=$1
  if [[ "$1" == *--array=* ]]; then HAS_ARRAY=1; else HAS_ARRAY=0; fi
  if [ $HAS_ARRAY -eq 1 ]; then
    ARRAY_ID=${1##*--array=}
    ARRAY_ID=${ARRAY_ID%% *}
    if [ "$ARRAY_ID" == "" ]; then echo "--array= must followed by an integer"; exit 1; fi
  fi
fi

if [[ "$SCRIPT" == */* ]]; then
  D=${SCRIPT%/*}
  F=${SCRIPT##*/}
  cd $D
else
  F=$SCRIPT
fi
if [ ! -f "$F" ]; then echo "$F not exist or is not a file"; exit 1; fi
if [ ! -s "$F" ]; then echo "$F is empty"; exit 1; fi
USE_ARRAY=`grep '$SLURM_JOB_NAME\|%x' "$F" | wc -l`
SET_ARRAY=`grep '^#SBATCH --array=' "$F" | wc -l`
if [ $USE_ARRAY -gt 0 ] && [ $SET_ARRAY -eq 0 ] && [ $HAS_ARRAY -eq 0 ]; then 
  echo "!!! Error: $F uses array without setting it. Please use --array, e.g.: victor.sbatch --array=1 $F"
  exit 1
fi
if [ $HAS_ARRAY -eq 1 ]; then
  if [ -f "$F".run_$ARRAY_ID.start ]; then
    echo "--array=$ARRAY_ID for $F has been done before. Please use another number."
    exit 1
  fi
  # PREVIOUS_ID=$((ARRAY_ID-1))
  # if [ "$PREVIOUS_ID" -gt 0 ] && [ ! -f "$F".run_$PREVIOUS_ID.start ]; then
  #   echo "--array=$PREVIOUS_ID for $F has NOT been done. The array ID must be an integer starting from 1, and increase by 1 each analysis."
  #   exit 1
  # fi
fi
JOBID=`sbatch "$1" "$F" | cut -d' ' -f 4`
echo "# JOBID=$JOBID" > "$F".account
echo "# HOSTNAME=$HOSTNAME" >> "$F".account
echo "# to abort this job: cat \""$F".account\" | bash" >> "$F".account
echo "scancel ${JOBID%%.*}" >> "$F".account
cd $CURRENT_DIR
