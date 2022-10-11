#!/bin/bash
for k in 0.1 0.2 0.5 1.0 2.0 5.0 10.0 20.0 50.0 100
do
geometric-neb --engine psi4 --coords chain_0064.xyz --align --nebk $k psi4_64.in > $k.out
done

