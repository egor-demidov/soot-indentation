#!/bin/bash

FACE_NUM=0
TIP_POINT=45
BOND_TO_BREAK=27
AGGREGATE_FILE="../aggregates/size_50/aggregate_1.txt"
DUMP_DIR="../cmake-build-debug/run"
DATA_FILE="../plots/production/01a_small_aggregate.dat"
taskset 0x3 ../cmake-build-debug/01_afm_aggregate $FACE_NUM $TIP_POINT $BOND_TO_BREAK $AGGREGATE_FILE $DUMP_DIR $DATA_FILE
