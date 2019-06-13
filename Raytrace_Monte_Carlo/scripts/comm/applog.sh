#!/bin/bash


for i in {1..20}
do
wait
cat $i/logfile >> logfile
cat $i/in_out.dat >> in_out.dat
cat $i/outangles.dat >> outangles.dat
cat $i/penetration.dat >> penetration.dat
cat $i/angle_dis.dat >> angle_dis.dat
rm $i/penetration.dat
rm $i/outangles.dat
rm $i/angle_dis.dat
rm $i/outputs.dat
done
