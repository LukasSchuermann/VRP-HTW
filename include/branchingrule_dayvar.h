
#ifndef VRP_BRANCHINGRULE_DAYVAR_H
#define VRP_BRANCHINGRULE_DAYVAR_H

#include <iostream>

#include "objscip/objscip.h"
#include "model_data.h"
#include "probdata_vrp.h"
#include "scip/scip.h"
#include "pricer_vrp.h"
#include "ConshdlrDayVar.h"
#include "var_tools.h"
#include "random"

using namespace scip;

class ObjBranchruleDayVar : public ObjBranchrule{
public:
    int nProhibit_ = 0;
    int nEnforce_ = 0;
    /** default constructor */
    explicit ObjBranchruleDayVar(
            SCIP*   scip,
            int     priority,
            double  frac
    )
    : ObjBranchrule(scip, "DayVarBranching", "Branching rule for the day variables", priority,
                    -1, 1.0),
                    branchingFrac_(frac)
    {
    }

    /** destructor */
    ~ObjBranchruleDayVar() override= default;

    /** branching execution method for fractional LP solutions */
    virtual SCIP_DECL_BRANCHEXECLP(scip_execlp);

    virtual SCIP_DECL_BRANCHINIT(scip_init);

private:
    std::vector<std::vector<SCIP_Real>> valonDay_;
    double branchingFrac_;
};


#endif //VRP_BRANCHINGRULE_DAYVAR_H
