
#include "prop_varfixing.h"
#include "pricer_vrp.h"
#include "vardata.h"
#include "tourVRP.h"

static
SCIP_RETCODE updateRootRedCosts(
    SCIP*               scip,
    ObjPricerVRP*       pricerData,
    model_data*         modelData,
    ObjPropVarFixing*   this_obj
){
    SCIP_Real best_redcosts;
    SCIP_Real best_objval;
    SCIP_Real current_redcosts;
    SCIP_Real current_objval = SCIPgetLocalLowerbound(scip);
    assert(SCIPisEQ(scip, current_objval, SCIPgetLPObjval(scip)));

    if(this_obj->withVeAss_){
        /* update vehicle assignment variables */
        for (int day = 0; day < modelData->nDays; day++) {
            for (auto cust: modelData->neighbors[0][day]) {
                if (!pricerData->global_timetable_[cust][day])
                    continue;
                /* best red costs of the vehicle assignment variable */
                best_redcosts = pricerData->root_dayVarRedCosts_[cust][day];
                best_objval = pricerData->root_dayVarLPObj_[cust][day];
                /* current red costs */
                current_redcosts = pricerData->dayVarRedCosts_[cust][day];
                /* update if current improves the best red costs */
                if(SCIPisGT(scip, current_redcosts + current_objval, best_redcosts + best_objval))
                {
                    pricerData->root_dayVarRedCosts_[cust][day] = current_redcosts;
                    pricerData->root_dayVarLPObj_[cust][day] = current_objval;
                }
            }
        }
    }

    if(this_obj->withArcFlow_){
        for(int day = 0; day < modelData->nDays; day++){
            modelData->neighbors[0][day].push_back(0);
            for(auto c1 : modelData->neighbors[0][day]){
                for(auto c2 : modelData->neighbors[0][day]){
                    if(c1 == c2 || pricerData->global_isForbidden_[day][c1][c2])
                        continue;
                    /* best red costs of the arc flow variable */
                    best_redcosts = pricerData->root_arcRedCosts_[day][c1][c2];
                    best_objval = pricerData->root_arcLPObj_[day][c1][c2];
                    /* current red costs */
                    current_redcosts = pricerData->arcRedCosts_[day][c1][c2];
                    /* update if current improves the best red costs */
                    if(SCIPisGT(scip, current_redcosts + current_objval, best_redcosts + best_objval))
                    {
                        pricerData->root_arcRedCosts_[day][c1][c2] = current_redcosts;
                        pricerData->root_arcLPObj_[day][c1][c2] = current_objval;
                    }
                }
            }
            modelData->neighbors[0][day].pop_back();
        }
    }

    return SCIP_OKAY;
}

static
SCIP_RETCODE varFixingRoot(
    SCIP*               scip,
    ObjPricerVRP*       pricerData,
    ObjPropVarFixing*   this_obj,
    model_data*         modelData,
    bool*               success
){
    int count_v = 0;
    int count_a = 0;
    SCIP_Real lhs;
    SCIP_Real cutoffbound = SCIPgetCutoffbound(scip);

    if(this_obj->withVeAss_){
        /* check vehicle assignment variables */
        for (int day = 0; day < modelData->nDays; day++) {
            for (auto cust: modelData->neighbors[0][day]) {
                if (!pricerData->global_timetable_[cust][day])
                    continue;
                /* check for day assignment variables with high reduced costs */
                lhs = pricerData->root_dayVarRedCosts_[cust][day] + pricerData->root_dayVarLPObj_[cust][day];
                if (SCIPisSumPositive(scip, lhs - cutoffbound)) {
                    count_v++;
                    *success = true;
                    pricerData->global_timetable_[cust][day] = false;
                    pricerData->timetable_[cust][day] = false;
                }
            }
        }
    }

    if(this_obj->withArcFlow_){
        /* check arc flow variables*/
        for(int day = 0; day < modelData->nDays; day++){
            modelData->neighbors[0][day].push_back(0);
            for(auto c1 : modelData->neighbors[0][day]){
                for(auto c2 : modelData->neighbors[0][day]){
                    if(c1 == c2 || pricerData->global_isForbidden_[day][c1][c2])
                        continue;
                    lhs = pricerData->root_arcRedCosts_[day][c1][c2] + pricerData->root_arcLPObj_[day][c1][c2];
                    if (SCIPisSumPositive(scip, lhs - cutoffbound)) {
                        count_a++;
                        *success = true;
                        pricerData->global_isForbidden_[day][c1][c2] = true;
                        pricerData->isForbidden_[day][c1][c2] = true;
                    }
                }
            }

            modelData->neighbors[0][day].pop_back();
        }
    }

    if(this_obj->verbose_){
        cout << "GlobalFixing: (" << SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) << ") cutoff: " << cutoffbound;
        cout << " --> DELETED VARS: (d) " << count_v << " (a) " << count_a << endl;
    }
    this_obj->arcflow_fixed_ += count_a;
    this_obj->vehicleass_fixed_ += count_v;

    return SCIP_OKAY;
}

static
SCIP_RETCODE varFixing(
        SCIP*               scip,
        ObjPricerVRP*       pricerData,
        ObjPropVarFixing*   this_obj,
        model_data*         modelData,
        bool*               success
){
    int count_v = 0;
    int count_a = 0;
    double gap = SCIPgetCutoffbound(scip) - SCIPgetLocalDualbound(scip);

    if(this_obj->withVeAss_){
        /* check vehicle assignment variables */
        for (int day = 0; day < modelData->nDays; day++) {
            for (auto cust: modelData->neighbors[0][day]) {
                if (!pricerData->timetable_[cust][day])
                    continue;
                /* check for day assignment variables that with high reduced costs */
                if (SCIPisSumPositive(scip, pricerData->dayVarRedCosts_[cust][day] - gap)) {
                    count_v++;
                    *success = true;
                    pricerData->timetable_[cust][day] = false;
                }
            }
        }
    }

    if(this_obj->withArcFlow_){
        /* check arc flow variables*/
        for(int day = 0; day < modelData->nDays; day++) {
            modelData->neighbors[0][day].push_back(0);
            for (auto c1: modelData->neighbors[0][day]) {
                if(!pricerData->timetable_[c1][day])
                    continue;
                for (auto c2: modelData->neighbors[0][day]) {
                    if(!pricerData->timetable_[c2][day])
                        continue;
                    if(c1 == c2 || pricerData->isForbidden_[day][c1][c2])
                        continue;
                    if (SCIPisSumPositive(scip, pricerData->arcRedCosts_[day][c1][c2] - gap)) {
                        count_a++;
                        *success = true;
                        pricerData->isForbidden_[day][c1][c2] = true;
                    }
                }
            }
            modelData->neighbors[0][day].pop_back();
        }
    }

    if(this_obj->verbose_){
        cout << "LocalFixing: (" << SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) << ") gap: " << gap;
        cout << " --> DELETED VARS: (d) " << count_v << " (a) " << count_a << endl;
    }
    this_obj->arcflow_fixed_ += count_a;
    this_obj->vehicleass_fixed_ += count_v;

    return SCIP_OKAY;
}

static
SCIP_RETCODE tourVarFixing(
    SCIP*               scip,
    vrp::ProbDataVRP*   probData,
    ObjPricerVRP*       pricerData,
    ObjPropVarFixing*   this_obj,
    bool                isGlobal
){
    vector<vector<bool>>& timetable = isGlobal ?  pricerData->global_timetable_ : pricerData->timetable_;
    vector<vector<vector<bool>>>& isForbidden = isGlobal ? pricerData->global_isForbidden_ : pricerData->isForbidden_;

    pricerData->fixed_nonzero_ = false;
    SCIP_Bool fixed;
    SCIP_Bool infeasible;
    int cnt = 0;
    for (auto var: probData->vars_) {
        fixed = false;
        if (!isGlobal && SCIPvarGetUbLocal(var) < 0.5)
            continue;
        if(isGlobal && SCIPvarGetUbGlobal(var) < 0.5)
            continue;

        auto *varData = dynamic_cast<ObjVarDataVRP*>(SCIPgetObjVardata(scip, var));
        tourVRP &tvrp = varData->tourVrp_;

        int day = tvrp.getDay();
        /* check first/last arc and first customer */
        if (isForbidden[day][tvrp.tour_[tvrp.length_ - 1]][0] ||
                isForbidden[day][0][tvrp.tour_[0]] || !timetable[tvrp.tour_[0]][day]) {
            fixed = TRUE;
        } else {
            /* check inner arcs and rest of customers */
            for (int v = 1; v < tvrp.length_; v++) {
                if (!timetable[tvrp.tour_[v]][day] ||
                    isForbidden[day][tvrp.tour_[v - 1]][tvrp.tour_[v]]) {
                    fixed = TRUE;
                    break;
                }
            }
        }
        if (fixed) {
            cnt++;
            if(!SCIPisZero(scip, SCIPvarGetLPSol(var)))
            {
//                cout << "var with val: " << SCIPvarGetLPSol(var) << endl;
//                SCIPprintVar(scip, var , nullptr);
                pricerData->fixed_nonzero_ = true;
            }
//                assert(SCIPisZero(scip, SCIPvarGetLPSol(var))); // TODO check if only due to EC
            if(isGlobal)
            {
                SCIP_CALL(SCIPtightenVarUbGlobal(scip, var, 0, TRUE, &infeasible, &fixed));
            }else
            {
                SCIP_CALL(SCIPfixVar(scip, var, 0.0, &infeasible, &fixed));
            }
            assert(!infeasible); // TODO - 60_0.25
            assert(fixed);
        }
    }

    this_obj->tourvar_fixed_ += cnt;

    return SCIP_OKAY;
}

/** execution method of propagator */
SCIP_DECL_PROPEXEC(ObjPropVarFixing::scip_exec) {
    bool success = false;
    bool success_root = false;
    auto *probData = dynamic_cast<vrp::ProbDataVRP *>(SCIPgetObjProbData(scip_));
    model_data *modelData = probData->getData();
    auto *pricerData = dynamic_cast<ObjPricerVRP *>(SCIPfindObjPricer(scip_, "VRP_Pricer"));

    *result = SCIP_DIDNOTRUN;
    /* check if new reduced costs have been calculated and if there is an optimal solution */
    if (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL) {
        return SCIP_OKAY;
    }
//    /* skip if no new reduced costs have been calculated, and we are not at the root */
//    if(!active_redcosts_ && SCIPgetDepth(scip) != 0){
//        return SCIP_OKAY;
//    }
    /* skip if there is no new cutoff bound and reduced costs at root */
    if(!active_redcosts_ && SCIPisEQ(scip, lastCutoff_, SCIPgetCutoffbound(scip)))
        return SCIP_OKAY;

    if(noRootFixing_ && !active_redcosts_)
        return SCIP_OKAY;

    if(verbose_){
        cout << "Start Propagator to FIX VARIABLES!!" << endl;
        cout << "active: " << active_redcosts_ << " lastco: " << lastCutoff_ << " currco: " << SCIPgetCutoffbound(scip) << endl;
    }
    *result = SCIP_DIDNOTFIND;

    if(SCIPgetDepth(scip) == 0)
    {
        if(active_redcosts_)
        {
            if(verbose_){
                cout << "UPDATE AT ROOT " << endl;
            }
            SCIP_CALL(updateRootRedCosts(scip, pricerData, modelData, this));
        }
        if(active_redcosts_ || SCIPisLT(scip, SCIPgetCutoffbound(scip), lastCutoff_))
        {
            if(verbose_){
                cout << "FIX AT ROOT" << endl;
            }

            SCIP_CALL(varFixingRoot(scip, pricerData, this, modelData, &success_root));
            lastCutoff_ = SCIPgetCutoffbound(scip);
        }
    }else{
        assert(SCIPgetDepth(scip) >= 1);
        if(SCIPisLT(scip, SCIPgetCutoffbound(scip), lastCutoff_))
        {
            assert(!active_redcosts_);
            if(verbose_){
                cout << "FIXforROOT!" << endl;
            }

            SCIP_CALL(varFixingRoot(scip, pricerData, this, modelData, &success_root));
            lastCutoff_ = SCIPgetCutoffbound(scip);
        }
        if(active_redcosts_)
        {
            if(verbose_){
                cout << "FIX AT TREE" << endl;
            }
            SCIP_CALL(varFixing(scip, pricerData, this, modelData, &success));
        }

    }

    if(success_root || success)
    {
        /* success XOR success_root */
        assert(success || success_root);
        assert(!success || !success_root);
        *result = SCIP_REDUCEDDOM;

        SCIP_CALL(tourVarFixing(scip, probData, pricerData, this, success_root));
        /* if no non-zero variable has been fixed, the next pricing iteration can be skipped  */
        pricerData->node_varsfixed_ = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
    }

    active_redcosts_ = false;
    return SCIP_OKAY;
}