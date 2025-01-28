
#ifndef __LABELING_ALGORITHMN_VRP__
#define __LABELING_ALGORITHMN_VRP__

#include "scip/scip.h"
#include "labellist.h"
#include "model_data.h"
#include "tourVRP.h"
#include "pricer_vrp.h"
#include "queue"

/* struct to pass arguments for labeling to worker threads */
typedef struct arg_struct {
    SCIP*               scip;
    model_data*         modelData;
    ObjPricerVRP*       pricerData;
    bool                isFarkas;        /**< TRUE for farkas-pricing, FALSE for redcost-pricing */
    bool                isHeuristic;
    bool                calcRedCosts;
    bool                noDomiance;
    int                 day;
    vector<tourVRP>*    bestTours;
} arg_struct;

SCIP_RETCODE labelingAlgorithmnParallel(
        SCIP*                       scip,
        vrp::ProbDataVRP*           probData,
        ObjPricerVRP*               pricerData,
        bool                        isFarkas,
        bool                        isHeuristic,
        bool                        calcRedCosts
);

SCIP_RETCODE generateLabelsBiDir(
        SCIP*               scip,
        model_data*         modelData,
        ObjPricerVRP*       pricerData,
        vector<tourVRP>&    bestTours,
        bool                isFarkas,
        bool                calcRedCosts,
        bool                isHeuristic,
        bool                noDominance,
        int                 day,
        bool                inCPC
);

SCIP_RETCODE propagateCustomer(
        SCIP*               scip,
        model_data          *modelData,
        ObjPricerVRP        *pricerData,
        int                 day,
        vector<LabelList*>  &propLabelLists,
        queue<int>          &indexQ,
        vector<LabelList*>  &labelLists,
        vector<bool>        &isInQ,
        bool                isFW,
        bool                calcRedCosts,
        bool                noDominance
);

SCIP_RETCODE concatenateLabels(
        SCIP*               scip,
        model_data*         modelData,
        ObjPricerVRP*       pricerData,
        int                 day,
        vector<LabelList*>  &fw_list,
        vector<LabelList*>  &bw_list,
        vector<pair<pair<LabelNode*, LabelNode*>, double>> &sol_pairs
);

SCIP_RETCODE concatenateLabelsFixing(
        SCIP*               scip,
        model_data*         modelData,
        ObjPricerVRP*       pricerData,
        int                 day,
        vector<LabelList*>  &fw_list,
        vector<LabelList*>  &bw_list
);

#endif