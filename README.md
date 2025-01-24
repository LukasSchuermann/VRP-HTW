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
$ ./bin/vrp <path-to-instance-json-file.json> 
    [-w <appointment window length in sec (optional; default 0)>] 
    [-o <output json file (optional; default used if no name is provided)>]
    [-s <input solution json file (optional)>]
    [-a <objective function parameters (delay, travel time, encounter probability) (default: 1, 0, 0) (mandatory) >]
    [-g <value for Gamma (optional; default 0)>]
    [-v <activate vehicle-assignment branching (optional)>]
    [-output <path to file for printing stats (optional)>]
```
