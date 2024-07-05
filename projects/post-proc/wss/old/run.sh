#!/bin/bash
DIRNAME=/local/fgalarce/4d-flow-ultrasound/projects/wallShearStress/
exec &> $DIRNAME/krenew.out
cd $DIRNAME
python run.py &> run.out
