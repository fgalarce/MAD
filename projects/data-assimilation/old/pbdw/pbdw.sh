#!/bin/bash
DIRNAME=/local/fgalarce/4d-flow-ultrasound/projects/pbdw/
exec &> $DIRNAME/noise_pbdw.out
cd $DIRNAME
python noise_pbdw.py -n=40 &> noise_pbdw.out
