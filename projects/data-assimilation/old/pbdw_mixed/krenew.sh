#!/bin/bash
DIRNAME=/home/local/fgalarce/4d-flow-ultrasound/projects/pbdw_mixed/
exec &> $DIRNAME/krenew.out
cd $DIRNAME
python noise_pbdw_mixed.py -n=20 &> krenew.out
