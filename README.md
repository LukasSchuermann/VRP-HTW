# VRP-HTW

This project contains the implementation of the branch-cut-and-price algorithm for the vehicle routing problem with heterogeneous time windows of my PhD thesis.
It is implementend in C using the MINLP solver SCIP.

	- In "OriginalCode" you can find all implemented techniques from the PhD thesis.
	- In "RCFC" we implemented the simple and iterative RCFC method for the VRP-HTW. This demands a manipulation of the SCIP installation, see README in the directory. Also it is based on an older version of the VRP code, e.g.~without vehicle-arc branching.
