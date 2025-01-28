
#ifndef VRP_BRANCHINGRULE_ARCFLOW_H
#define VRP_BRANCHINGRULE_ARCFLOW_H

#include <iostream>

#include "objscip/objscip.h"
#include "model_data.h"
#include "probdata_vrp.h"
#include "scip/scip.h"
#include "pricer_vrp.h"
#include "vardata.h"
#include "var_tools.h"
#include "ConshdlrArcflow.h"
#include "random"

using namespace scip;

class ObjBranchruleArcflow : public ObjBranchrule{
public:
    int nProhibit_ = 0;
    int nEnforce_ = 0;
    /** default constructor */
    explicit ObjBranchruleArcflow(
            SCIP*   scip,
            int     priority,
            double  frac
    )
    : ObjBranchrule(scip, "ArcFlowBranching", "Branching rule for the arc flow variables", priority,
                    -1, 1.0),
                    branchingFrac_(frac)
    {
    }

    /** destructor */
    ~ObjBranchruleArcflow() override= default;

    /** branching execution method for fractional LP solutions */
    virtual SCIP_DECL_BRANCHEXECLP(scip_execlp);

    virtual SCIP_DECL_BRANCHINIT(scip_init);

private:
    std::vector<std::vector<SCIP_Real>> arcWeights_;
    double branchingFrac_;
};

#endif //VRP_BRANCHINGRULE_ARCFLOW_H
