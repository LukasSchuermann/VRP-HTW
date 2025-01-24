# VRP-HTW

This project contains the implementation of the branch-cut-and-price algorithm for the vehicle routing problem with heterogeneous time windows of my PhD thesis.
It is implementend in C using the MINLP solver SCIP.


## Requirements

To run, this codes needs an installation of SCIP (we used version 9.1.1). Env_Variable "SCIP_DIR" shall contain the path to the SCIP installation.
See the website for an installation guide:
https://scipopt.org/#scipoptsuite



## Usage

```markdown
$ mkdir build
$ cd build
$ cmake .. -DSCIP_DIR=/path/to/scip/installation
$ make
$ ./vrp <path-to-instance-json-file.json> 
    [-b <branching rule (options: 1 vehicle assignment, 2 arc flow, 3 vehicle arc; default: 1)>] 
    [-s <seed (options >= 0)>]
    [-sol <path-to-solution-file>]
    [-out <path-to-output-file>]
    [-src <maximum numbers of subset row cuts (options >= 0)>]
    [-nkpc <deactivate k-path cut separation>]
    [-bf <branching decision value (options in [0.0, 1.0]; 1.0 for random branching)>]
    [-rcf <reduced cost fixing properties (decay, max_depth, maxNumSRC)>]
    [-nR <deactivate root reduced cost fixing>]
    [-timeout <set time limit>]
```
