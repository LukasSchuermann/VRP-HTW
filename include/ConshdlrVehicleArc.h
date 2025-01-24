
#include "objscip/objscip.h"
#include "branchingrule_vehiclearc.h"
#include "model_data.h"
#include "probdata_vrp.h"
#include "scip/scip.h"
#include "vardata.h"
#include "iostream"


#define CONSHDLR_NAME3          "vehiclearc"
#define CONSHDLR_DESC3          "stores the local branching decisions"

#ifndef VRP_CONSHDLRVEHICLEARC_H
#define VRP_CONSHDLRVEHICLEARC_H

using namespace scip;

class ConshdlrVehicleArc : public ObjConshdlr
{
public:
    /** default constructor */
    explicit ConshdlrVehicleArc(
        SCIP*   scip
        )
        : ObjConshdlr(scip, CONSHDLR_NAME3, CONSHDLR_DESC3, 0, 0, 1,
                      -1, 1, 1, 0, FALSE, FALSE, TRUE,
                      SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_ALWAYS)
    {
    }

    /** destructor */
    ~ConshdlrVehicleArc() override= default;

    /** frees specific constraint data */
    virtual SCIP_DECL_CONSDELETE(scip_delete);

    /** transforms constraint data into data belonging to the transformed problem */
    virtual SCIP_DECL_CONSTRANS(scip_trans);

    /** constraint enforcing method of constraint handler for LP solutions */
    virtual SCIP_DECL_CONSENFOLP(scip_enfolp){
        return SCIP_OKAY;
    }

    /** constraint enforcing method of constraint handler for pseudo solutions */
    virtual SCIP_DECL_CONSENFOPS(scip_enfops){
        return SCIP_OKAY;
    }

    /** feasibility check method of constraint handler for integral solutions */
    virtual SCIP_DECL_CONSCHECK(scip_check){
        return SCIP_OKAY;
    }

    /** domain propagation method of constraint handler */
    virtual SCIP_DECL_CONSPROP(scip_prop);

    /** variable rounding lock method of constraint handler */
    virtual SCIP_DECL_CONSLOCK(scip_lock){
        return SCIP_OKAY;
    }

    /** constraint activation notification method of constraint handler */
    virtual SCIP_DECL_CONSACTIVE(scip_active);

    /** constraint deactivation notification method of constraint handler */
    virtual SCIP_DECL_CONSDEACTIVE(scip_deactive);

    /** constraint display method of constraint handler */
    virtual SCIP_DECL_CONSPRINT(scip_print);
};

SCIP_RETCODE SCIPcreateConsvehiclearc(
        SCIP*               scip,                /**< SCIP data structure */
        SCIP_CONS**         cons,                /**< pointer to hold the created constraint */
        const char*         name,                /**< name of the constraint */
        int                 tail,                /**< tail of the arc */
        int                 head,                /**< head of the arc */
        int                 day,
        CONSTYPE            type,                /**< stores whether arc gets enforced or prohibited */
        SCIP_NODE*          node,                /**< the node in the B&B-tree at which the cons is sticking */
        SCIP_Bool           local                /**< is constraint only valid locally? */
);

/** returns tail of the arc */
int SCIPgetTailvehiclearc(
        SCIP_CONS*            cons                /**< samediff constraint */
);

/** returns head of the arc */
int SCIPgetHeadvehiclearc(
        SCIP_CONS*            cons                /**< samediff constraint */
);

/** returns day of the arc */
int SCIPgetDayvehiclearc(
        SCIP_CONS*            cons                /**< samediff constraint */
);

/** return constraint type PROHIBIT or ENFORCE */
CONSTYPE SCIPgetTypevehiclearc(
        SCIP_CONS*            cons                /**< samediff constraint */
);


#endif //VRP_CONSHDLRVEHICLEARC_H
