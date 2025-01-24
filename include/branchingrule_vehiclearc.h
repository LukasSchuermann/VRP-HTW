
#ifndef VRP_BRANCHINGRULE_VehicleArc_H
#define VRP_BRANCHINGRULE_VehicleArc_H

#include <iostream>

#include "objscip/objscip.h"
#include "model_data.h"
#include "probdata_vrp.h"
#include "scip/scip.h"
#include "pricer_vrp.h"
#include "vardata.h"
#include "var_tools.h"
#include "ConshdlrVehicleArc.h"
#include "random"

using namespace scip;

class ObjBranchruleVehicleArc : public ObjBranchrule{
public:
    int nProhibit_ = 0;
    int nEnforce_ = 0;
    /** default constructor */
    explicit ObjBranchruleVehicleArc(
            SCIP*   scip,
            int     priority,
            double  frac
    )
    : ObjBranchrule(scip, "VehicleArcBranching", "Branching rule for the arc flow variables", priority, -1, 1.0),
                    branchingFrac_(frac)
    {
    }

    /** destructor */
    ~ObjBranchruleVehicleArc() override= default;

    /** branching execution method for fractional LP solutions */
    virtual SCIP_DECL_BRANCHEXECLP(scip_execlp);

    virtual SCIP_DECL_BRANCHINIT(scip_init);

private:
    std::vector<std::vector<std::vector<SCIP_Real>>> vehicleArcVal_;
    double branchingFrac_;
};

#endif //VRP_BRANCHINGRULE_VehicleArc_H
