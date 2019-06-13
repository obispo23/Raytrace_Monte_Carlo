#/bin/bash
for i in 300
do
    qsub -cwd -V -N foam046_$i -j y -pe dc* 5 -l  h_data=512M,h_rt=1:00:00  ./run0.sh $i
done