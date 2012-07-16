#!/bin/sh
#PBS -l nodes=1:ppn=8
#PBS -l walltime=01:00:00
#PBS -o output.txt
#PBS -j oe
#PBS -m bea
#PBS -M edz3fz@virginia.edu

cd $PBS_O_WORKDIR

NODES=`cat $PBS_NODEFILE | uniq`

echo $NODES > nodefile.txt

PORT=23335

for n in $NODES
do 
	ssh -f edz3fz@$n "cd $PBS_O_WORKDIR; mkdir -p $TMPDIR; export TMPDIR=$TMPDIR; ppserver.py -p $PORT &" &
done

python temprxGO.py -f GO_2QUG.pdb -p GO_2QUG.param -n 1000 -s 10 -k 100 -r 8 -t 300 450 -a $PORT
