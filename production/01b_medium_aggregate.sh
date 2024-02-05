#!/bin/bash

FACE_NUM=1
TIP_POINT=140
BOND_TO_BREAK=27
AGGREGATE_FILE="../aggregates/size_150/aggregate_1.txt"
DUMP_DIR="../cmake-build-debug/run2"
DATA_FILE="../plots/production/01b_medium_aggregate.dat"
taskset 0x3FFC ../cmake-build-debug/01_afm_aggregate $FACE_NUM $TIP_POINT $BOND_TO_BREAK $AGGREGATE_FILE $DUMP_DIR $DATA_FILE
