#! /bin/bash

./main3d.gnu.TEST.MPI.ex ./inputs.3d.linear.askew-yz 

fcompare plot-3d-askew-yz/ plot-3d-askew-yz-analytic/
