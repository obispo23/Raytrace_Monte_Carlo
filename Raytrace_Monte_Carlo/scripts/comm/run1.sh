#!/bin/bash
#echo $2
mkdir $3
cd $3
../../comm/raytrace filename=../../foam_$(($2-1)).stl Nsteps=50000 ParentE=$1 random=yes cutoff=8.68 bvh=8  distribution=cosine species1=Ar species2=W > logfile
cd ../

