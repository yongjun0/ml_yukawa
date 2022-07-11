#!/bin/bash
kappa="0.5 1 1.5 2 2.5 3"
g_for_k05="11 183 61"
g_for_k10="11 14 144 217 234 72"
g_for_k15="116 22 328"
g_for_k20="158 26 31 320 449 476"
g_for_k25="328 60 890"
g_for_k30="100 1000 1174 1510 503 62"

kappa="2.1 2.2 2.3 2.4 2.5"
gamma="60 70 80 90 100"
dest="data_sample/sqw/"
mkdir -p $dest
for ik in $kappa; do
  #if [ "$ik" == "0.5"  ]; then
  # gamma=$g_for_k05
  #elif [ "$ik" == "1"  ]; then
  # gamma=$g_for_k10
  #elif [ "$ik" == "1.5"  ]; then
  # gamma=$g_for_k15
  #elif [ "$ik" == "2"  ]; then
  # gamma=$g_for_k20
  #elif [ "$ik" == "2.5"  ]; then
  # gamma=$g_for_k25
  #elif [ "$ik" == "3"  ]; then
  # gamma=$g_for_k30
  #fi
  for ig in $gamma; do
    python grnn.py ${ik} ${ig} &
  done
done


