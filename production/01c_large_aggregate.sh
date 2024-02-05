#!/bin/bash

FACE_NUM=10
TIP_POINT=45
BOND_TO_BREAK=27
AGGREGATE_FILE="../aggregates/size_300/aggregate_1.txt"
DUMP_DIR="../cmake-build-debug/run3"
DATA_FILE="../plots/production/01c_large_aggregate.dat"
taskset 0x3FC0 ../cmake-build-debug/01_afm_aggregate $FACE_NUM $TIP_POINT $BOND_TO_BREAK $AGGREGATE_FILE $DUMP_DIR $DATA_FILE
