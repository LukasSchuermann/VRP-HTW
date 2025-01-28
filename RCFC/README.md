# RCFC

This project contains the implementation of the RCFC methods for MILPs of my PhD thesis.
It is implementend in C++ using the MINLP solver SCIP.

## Requirements

For an efficient implementation of the RCFC method, we want to stop the (primal) simplex algorithm as soon as the considered basic variable exceeds a value.
As a workaround, we called the simplex algorithm with a limited number of iterations and checked in between each run.
For that, we need to add a function to SCIP that yields the primal solution of a column, even when the simplex is not terminated yet.

The following steps need to be followed to use our code:
  1. Download the source code of the SCIP Optimization Suite from https://scipopt.org/index.php#download (we used version 9.1.1)
  2. In "scip/src/lpi/" we need to adapt "lpi.h" and the corresponding LP-solver file ("lpi_spx2.cpp" in our case) accordingly:
     ```markdown
     # Add this function to lpi.h
     SCIP_EXPORT
     SCIP_Real SCIPlpiColGetNewLPvalRCFC(
        SCIP_LPI*           lpi,
        int                 colIndex
     );
     ```
     ```markdown
     # Add this declaration to lpi_spx2.cpp
     SCIP_Real SCIPlpiColGetNewLPvalRCFC(
        SCIP_LPI*           lpi,
        int                 colIndex
     ){
        return lpi->spx->getPrimalRealIndex(colIndex);
     }
     ```
  3. In "scip/src/scip/" we need to append the following code to "scip_probing.h" and "scip_probing.c":
     ```markdown
     # Append these functions to scip_probing.h (e.g. at line 547)
     SCIP_EXPORT
     int SCIPGetNlpiColsRCFC(
        SCIP*           scip
     );    
     SCIP_EXPORT
     SCIP_COL** SCIPGetlpiColsRCFC(
        SCIP*           scip
     ); 
     ```
     ```markdown
     # Append these declarations to scip_probing.c (e.g. at line 1328)
     int SCIPGetNlpiColsRCFC(
        SCIP*           scip
     ){
        return getNlpiColsRCFC(scip->lp);
     }
     SCIP_COL** SCIPGetlpiColsRCFC(
        SCIP*           scip
     ){
        return getlpiColsRCFC(scip->lp);
     }
     ```
  4. Additionally, in "scip/src/scip/" we need to append the following code to "lp.h" and "lp.c":
     ```markdown
     # Append these functions to lp.h (e.g. at line 1638)
     int getNlpiColsRCFC(
        SCIP_LP*           lp
     );
     SCIP_COL** getlpiColsRCFC(
        SCIP_LP*           lp
     );
     ```
     ```markdown
     # Append these declarations to lp.c (e.g. at line 18945)
     int getNlpiColsRCFC(
        SCIP_LP*               lp
     ){
        return lp->nlpicols;
     }
     SCIP_COL** getlpiColsRCFC(
        SCIP_LP*                lp
     ){
        return lp->lpicols;
     }
     ```
     ```markdown
     # Exchange line 12444 in lp.c
     if( lp->flushed && lp->solved )
     # by:
     if( lp->flushed && lp->solved && !(lp->probing && lp->lpsolstat == SCIP_LPSOLSTAT_ITERLIMIT) )
     ```
  6. In "soplex/src/" we need to add the following code to the corresponding files:
     ```markdown
     # Add this function to the SoPlexBase class in soplex.h (e.g. at line 664)
     double getPrimalRealIndex(int index);
     ```
     ```markdown
     # Add the declaration to soplex.hpp (e.g. at line 1157)
     template <class R>
     double SoPlexBase<R>::getPrimalRealIndex(int index){
        return _solReal._primal[index];
     }
     ```
     ```markdown
     # Add this function to soplex_interface.h (e.g. at line 71)
     double SoPlex_getPrimalRealIndex(void* soplex, int index);
     ```
     ```markdown
     # Add the declaration to soplex_interface.cpp (e.g. at line 227)
     double SoPlex_getPrimalRealIndex(void* soplex, int index){
        SoPlex* so = (SoPlex*)(soplex);
        return so->getPrimalRealIndex(index);
     }
     ```
  5. Now compile SCIP (see https://github.com/scipopt/scip/blob/master/INSTALL.md) and link the installation to this project by either setting the $ENV{SCIP_DIR} variable to the correct path or include -DSCIP_DIR="/path/to/scip/installation" in the compiling steps.

### Usage

The compilation of this project works as follows:
```markdown
$ mkdir build
$ cd build
$ cmake .. -DSCIP_DIR="/path/to/scip/installation"
$ make
```
Now, we can execute our algorithm. For example:
```markdown
$ ./rcfc ../../instances/pure_integer/irp.mps +useRCFC
```
