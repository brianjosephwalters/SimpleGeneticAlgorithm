Programming Assignment #1 - Simple Genetic Algorithim
Brian J. Walters
2418606

Dependencies:
MPICH2 v1.4.1
GCC v4.6.3

Compiled on:
Ubuntu 12.04.3 LTS
Kernel 3.8.0-31-generic

Simple Genetic Algorithm
------------------------
Usage: mpirun [OPTIONS] ./sga [OPTIONS]
  See man mpirun for MPI usage.
Options:
  -f {file name}    output file name
  -g {value}        number of generations   0 < value
  -p {value}        size of population      0 < value
  -c {value}        crossover rate:         0 < value < 1
  -m {value}        mutation rate:          0 < value < 1
  -e {value}        elite rate:             0 < value < 1

  -d                display default values
  -h                display help

----- Default Values -----
Generations: 2000
Chromosomes: 200
Genes per Chromosome: 10
Elite rate: +0.05
Crossover Rate: +0.80
Mutation rate: +0.10
Per Generation:
  elites: 10
  crossed chromosomes: 160
  remaining chromosomes: 30
  mutations: 20 genes

Output printed to: output.txt
