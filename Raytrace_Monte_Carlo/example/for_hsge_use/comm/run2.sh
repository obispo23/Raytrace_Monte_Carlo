#!/bin/bash
i=$(($1-1))
../comm/remesh filename=../foam_$(($1-1)).stl Nsteps=50000 atom_volume=1.583E-17 flux=1E13 scale=250000 smooth=yes del=yes island=20 > remeshlog
wait
mv newmesh.stl ../foam_$1.stl



