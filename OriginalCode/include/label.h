
#ifndef VRP_LABEL_H
#define VRP_LABEL_H

#include <utility>
#include "iostream"
#include "scip/scip.h"
#include "vector"
#include "bitset"
#include "model_data.h"
#include "tourVRP.h"
#include "pricer_vrp.h"

class Label {
public:
    SCIP_Real                   red_costs_;   /**< reduced costs */
    int                         current_;     /**< last customer */
    int                         cap_;         /**< current capacity */
    double                      time_;        /**< depature time */
    double                      obj_;         /**< current tour costs */
    bitset<neighborhood_size>   ng_memory_;   /**< active ng-neighborhood */
    bitset<neighborhood_size>   visitedEC_;   /**< keeps track of visited enforced customers */
    vector<double>              SRCstate_;    /**< states of subset row cuts */
    Label(
        SCIP_Real           red_costs,
        int                 current,
        int                 cap,
        int                 nSRC,
        double              time,
        double              obj,
        bitset<neighborhood_size>& ng_memory,
        bitset<neighborhood_size>& visitedEC
    ):
    red_costs_(red_costs),
    current_(current),
    cap_(cap),
    time_(time),
    obj_(obj),
    ng_memory_(ng_memory),
    visitedEC_(visitedEC)
    {
        SRCstate_.resize(nSRC, 0.0);
    }

    ~Label()
    = default;

    SCIP_Bool dominates(
            Label*              label,
            ObjPricerVRP*       pricerData,
            SCIP*               scip
    ) const
    {
        /** checks if THIS dominates label */
        assert(label != nullptr);

        /* capacity*/
        if(cap_ > label->cap_)
            return FALSE;
        /* departure time */
        if(SCIPisGT(scip, time_, label->time_))
            return FALSE;
        /* ng-neighborhood */
        if(!((ng_memory_ & label->ng_memory_) == ng_memory_))
            return FALSE;
        /* visited enforced customer */
        if(!(visitedEC_ == label->visitedEC_)) //TODO: only if day has EC
            return FALSE;
        /* reduced costs including potential influence of Subset row cuts */
        double rc_and_src = label->red_costs_;
        for(int c = 0; c < pricerData->nnonzSRC_; c++)
        {
            if(SCIPisGT(scip, SRCstate_[c], label->SRCstate_[c]))
                rc_and_src += pricerData->SRC_dualv[c];
        }
        if(SCIPisGT(scip, red_costs_, rc_and_src))
            return FALSE;

        return TRUE;
    }
    SCIP_Bool dominates_enumerate(
            Label*              label,
            SCIP*               scip
    ) const
    {
        /** checks if THIS dominates label in enumerate mode */
        assert(label != nullptr);

        /* route costs */
        if(SCIPisGT(scip, obj_, label->obj_))
            return FALSE;
        /* departure time */
        if(SCIPisGT(scip, time_, label->time_))
            return FALSE;
        /* same ng-neighborhood */
        if(!(ng_memory_ == label->ng_memory_))
            return FALSE;
        /* visited enforced customer */
        if(!(visitedEC_ == label->visitedEC_))
            return FALSE;

        return TRUE;
    }
};

Label* label_propagate2(
        model_data*     modelData,
        ObjPricerVRP*   pricerData,
        Label*         oldLabel,
        SCIP_Bool       isFW,
        int             end,
        int             day
);


#endif //VRP_LABEL_H
