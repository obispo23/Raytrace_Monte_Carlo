#!/bin/bash 
cd comm
g++ -std=c++11 -w -o remesh ../../../remeshing/remesh.cxx
g++ -std=c++11 -w -o raytrace ../../../source/raytrace.cxx 
g++ -std=c++11 -w -o avecalc ../../../scripts/avecalc.cxx
cd ..
for j in {1..10} #number of iterations
do
    mkdir ${1}_${j}
    cd ${1}_${j}/
    for i in {1..6} #number of cores used -1
    do
	if [ $i -lt 6 ]
	then
	    echo $i
	    ../comm/run1.sh $1 $j $i & #execute raytrace
	fi
	if [ $i -gt 5 ]
	then
	    wait
	    ../comm/applog.sh #combine files
	    ../comm/avecalc filename=logfile par=5 Nstep=50000 #calculate in, out, and yield
	    ../comm/run2.sh $j #execute remeshing
	fi
    done
    cd ../
done
wait
exit
