
#include "ConshdlrRCFC.h"
#include "scip/cons_setppc.h"
#include "objscip/objscip.h"
#include "pricer_vrp.h"
#include "vardata.h"
#include "tourVRP.h"
#include "var_tools.h"
#include "ConshdlrRCFC.h"
#include <ctime>


/** Constraint data for "CPC" constraints */
struct SCIP_ConsData
{
    SCIP_Row*           cut;            /**< Corresponding row in the LP */
    SCIP_VAR**          vars;           /**< variables of the constraint */
    int                 nvars;          /**< number of variables in the constraint */
};

/** local methods */
static
SCIP_RETCODE consdataCreate(
        SCIP*           scip,
        SCIP_CONSDATA** consdata,
        vector<SCIP_VAR*>& vars
){
    SCIP_CALL(SCIPallocBlockMemory(scip, consdata));

    (*consdata)->cut = nullptr;
    (*consdata)->nvars = (int) vars.size();
    SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->nvars));
    for(int i = 0; i < (*consdata)->nvars; i++)
    {
        (*consdata)->vars[i] = vars[i];
    }

    return SCIP_OKAY;
}

static
SCIP_RETCODE addVarsToCons(
    SCIP*           scip,
    SCIP_CONSDATA*  consdata
){
    assert(consdata->cut != nullptr);
    for(int i = 0; i < consdata->nvars; i++)
    {
        SCIP_CALL(SCIPaddVarToRow(scip, consdata->cut, consdata->vars[i], 1));
    }

    return SCIP_OKAY;
}

static
SCIP_RETCODE createCPCCut(
        SCIP*               scip,
        vector<SCIP_VAR*>&  cut_vars
){
    SCIP_CONSDATA* consdata = nullptr;
    SCIP_CONSHDLR* conshdlr = nullptr;
    SCIP_Cons* cons = nullptr;
    SCIP_Bool infeasible;

    assert(scip != nullptr);

    conshdlr = SCIPfindConshdlr(scip, "RCFC");

    /* create constraint data */
    SCIP_CALL(consdataCreate(scip, &consdata, cut_vars));
    /* create constraint */
    SCIP_CALL(SCIPcreateCons(scip, &cons, "CPC-cut", conshdlr, consdata, FALSE, TRUE, TRUE, FALSE,
                             FALSE, TRUE, FALSE, FALSE, FALSE, TRUE));
    /* create lp row and add it */
    SCIP_CALL(SCIPcreateEmptyRowCons(scip, &(consdata->cut), cons, SCIPconsGetName(cons), 1,
                                     SCIPinfinity(scip), TRUE, FALSE, FALSE));
    SCIP_CALL(addVarsToCons(scip, consdata));
    SCIP_CALL(SCIPaddRow(scip, consdata->cut, true, &infeasible));
    assert(!infeasible);

    SCIPreleaseCons(scip, &cons);

    return SCIP_OKAY;
}

static
SCIP_RETCODE fixVariables(
        SCIP*               scip,
        const vector<SCIP_VAR*>&  fixableVars

){
    auto *probData = dynamic_cast<vrp::ProbDataVRP *>(SCIPgetObjProbData(scip));
    auto *pricerData = dynamic_cast<ObjPricerVRP *>(SCIPfindObjPricer(scip, "VRP_Pricer"));


    long long int current_node = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
    pricerData->tree_data_[current_node].gotFixed = true;
    for(auto var : fixableVars)
    {
        if(SCIPisLT(scip, SCIPvarGetLPSol(var), 1))
        {
            pricerData->tree_data_[current_node].gotFixed = false;
            pricerData->node_varsfixed_ = current_node;
        }
        SCIPchgVarLb(scip, var, 1.0);
        auto* varData = dynamic_cast<ObjVarDataVRP*>(SCIPgetObjVardata(scip, var));
        int day = varData->tourVrp_.getDay();

//        std::cout << varData->tourVrp_;
        if(probData->getData()->num_v[day] == 1)
            pricerData->fixedDay_[day] = true;
        for(auto u : varData->tourVrp_.tour_)
        {
            pricerData->eC_[u] = day;
            for(int d = 0; d < probData->getData()->nDays; d++)
            {
                if(d == day)
                    continue;
                if(current_node == 1)
                    pricerData->global_timetable_[u][d] = false;
                else
                    pricerData->timetable_[u][d] = false;
            }
        }
    }
    return SCIP_OKAY;
}


static
SCIP_RETCODE exactBranching(
        SCIP*                   scip,
        ConshdlrRCFC*            obj,
        SCIP_Var*               var,
        SCIP_Bool*              success
){
    SCIP_Bool cutoff, lperror;
    long long int old_iter = SCIPgetNLPIterations(scip);
    SCIPstartProbing(scip);
    SCIPchgVarUbProbing(scip, var, 0.0);

    SCIPsolveProbingLPWithPricing(scip, false, false, -1, &lperror, &cutoff);


    obj->LP_iters_exact_ += (SCIPgetNLPIterations(scip) - old_iter);
    if(SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)))
    {
        assert(cutoff);
        if(success != nullptr)
            *success = true;
        obj->nFixedExact_++;
    }

    SCIPendProbing(scip);

    return SCIP_OKAY;
}

static
SCIP_RETCODE heuristicNoMove(
        SCIP*                   scip,
        ConshdlrRCFC*            obj,
        SCIP_Col*               col,
        int                     index,
        double                  gap_abs,
        double                  varLPVal,
        SCIP_Bool*              success
){
    SCIP_Bool cutoff, lperror;
    SCIP_Var* var = SCIPcolGetVar(col);
    SCIP_LPI* lpi;
    int k;
    int nIters;

    double local_lpVal = varLPVal;

    SCIP_CALL( SCIPgetLPI(scip, &lpi) );

    obj->start_time_ = std::chrono::high_resolution_clock::now();
    SCIPstartProbing(scip);
    SCIPchgVarObjProbing(scip, var, SCIPvarGetObj(var) + gap_abs / local_lpVal);
    k = obj->iter_init_;
    nIters = 0;
    /* move from one optimum to the optimum for the next objective function */
    while (SCIPisEQ(scip, local_lpVal, varLPVal)) {
        if (SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL) {
            break;
        }
        SCIPsolveProbingLP(scip, k, &lperror, &cutoff);
        nIters += k;
        k = MAX((int) (obj->iter_perc_ * nIters), obj->iter_init_);
        local_lpVal = SCIPlpiColGetNewLPvalRCFC(lpi, index);
    }

    if(SCIPisEQ(scip, local_lpVal, varLPVal)){
        SCIPsolveProbingLPWithPricing(scip, false, false, -1, &lperror, &cutoff);
        if (SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip))){
            *success = true;
        }
    }

    SCIPendProbing(scip);


    return SCIP_OKAY;
}

static
SCIP_RETCODE heuristicExtended(
        SCIP*                   scip,
        ConshdlrRCFC*            obj,
        SCIP_Col*               col,
        int                     index,
        double                  gap_abs,
        double                  varLPVal,
        SCIP_Bool*              success
){
    SCIP_Bool cutoff, lperror;
    SCIP_Var* var = SCIPcolGetVar(col);
    SCIP_LPI* lpi;
    bool tooSteep = false;
    int k;
    int nIters;

    double local_gap = gap_abs;
    double local_lpVal = varLPVal;
    double orig_gap = local_gap;
    double addedObj = 0;
    double init_steepness;
    int iter = 0;

    double steepness = gap_abs / varLPVal;
    init_steepness = steepness;

    SCIP_CALL( SCIPgetLPI(scip, &lpi) );

    obj->start_time_ = std::chrono::high_resolution_clock::now();
    SCIPstartProbing(scip);
    /* move from one optimum to the optimum for the next objective function */
    while (SCIPisGT(scip, local_lpVal, 0)) {
        if(steepness > obj->lowest_fail_ * obj->factor_)
        {
            tooSteep = true;
            obj->steepness_skip_++;

            break;
        }

        addedObj += (local_gap / local_lpVal);

        k = obj->iter_init_;
        nIters = 0;
        SCIPchgVarObjProbing(scip, var, SCIPvarGetObj(var) + local_gap / local_lpVal);

        /* First iteration */
        SCIPsolveProbingLP(scip, k, &lperror, &cutoff);
        local_lpVal = SCIPlpiColGetNewLPvalRCFC(lpi, index);
        /* Abort Probing once x_i = 0 */
        while (SCIPisGT(scip, local_lpVal, 0)) {
            orig_gap = SCIPgetCutoffbound(scip) - (SCIPgetLPObjval(scip) - addedObj * local_lpVal);
            if(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL){
                break;
            }

            /* abort early if objective function is too steep relative to lp value of x_i */
            if((orig_gap/local_lpVal) > obj->lowest_fail_ * obj->factor_)
            {
                tooSteep = true;
                obj->steepness_skip_++;

                break;
            }

            nIters += k;
            k = MAX((int) (obj->iter_perc_ * nIters), obj->iter_init_);
            SCIPsolveProbingLP(scip, k, &lperror, &cutoff);

            local_lpVal = SCIPlpiColGetNewLPvalRCFC(lpi, index);
            if(SCIPisGT(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)))
            {
                tooSteep = true;
                break;
            }
            if(local_lpVal < varLPVal / 2)
            {
                tooSteep = true;
                break;
            }
        }

        iter++;
        if(SCIPisZero(scip, local_lpVal) || tooSteep)
        {
            auto end_t = std::chrono::high_resolution_clock::now();
            obj->lp_solving_time_ += std::chrono::duration_cast<std::chrono::milliseconds>(end_t - obj->start_time_).count();
            break;
        }
        assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);
        assert(SCIPvarGetLPSol(var) == SCIPlpiColGetNewLPvalRCFC(lpi, index));

        if (SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)))
        {
            SCIPsolveProbingLPWithPricing(scip, false, false, -1, &lperror, &cutoff);
            if (SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)))
            {
                obj->nFixedHeurExtv2_++;
                *success = true;
                break;
            }else if(SCIPisZero(scip, SCIPvarGetLPSol(var))){
                break;
            }
            /* potentially we want to stop after first pricing iteration -> too expensive */
            if(obj->used_Pricing_){
                break;
            }
        }
        obj->start_time_ = std::chrono::high_resolution_clock::now();
        local_gap = SCIPgetCutoffbound(scip) - SCIPgetLPObjval(scip);
        steepness += (local_gap / local_lpVal);
    }

    SCIPendProbing(scip);

    if(!*success)
    {
        if(obj->lowest_fail_ > init_steepness)
        {
            obj->lowest_fail_ = init_steepness;
        }
    }

    return SCIP_OKAY;
}


static
SCIP_RETCODE checkFixing(
        SCIP*                   scip,
        ConshdlrRCFC*           cons_obj,
        std::vector<int>        &checkLPIcols,
        SCIP_RESULT*            result,
        bool                    onlyFrac
){
    vector<SCIP_VAR*> fixableVars;
    vector<bool> violatedCut;
    std::vector<std::vector<SCIP_Var*>> cutVars;
    double cmp_low;

    bool fixed_frac_ = false;

    auto *probData = dynamic_cast<vrp::ProbDataVRP *>(SCIPgetObjProbData(scip));

    double gap_rel = 100 * (SCIPgetCutoffbound(scip) - SCIPgetLPObjval(scip)) / SCIPgetLPObjval(scip);

    int nSuccess = 0;

    for(auto index : checkLPIcols)
    {
        double gap_abs = SCIPgetCutoffbound(scip) - SCIPgetLPObjval(scip);
        if(SCIPisInfinity(scip, abs(gap_abs)))
            break;
        SCIP_COL** lpicols = SCIPGetlpiColsRCFC(scip);

        SCIP_Col* col = lpicols[index];
        if(col == nullptr || SCIPcolGetBasisStatus(col) != SCIP_BASESTAT_BASIC)
            continue;
        SCIP_Var* var = SCIPcolGetVar(col);
        auto varLPval = SCIPcolGetPrimsol(col);

        if(fixableVars.size() > 3)
           break;
        if(onlyFrac && fixableVars.size() >= 2)
            break;
        if(SCIPvarGetLbLocal(var) > 0.5)
            continue;
        if(SCIPvarGetLPSol(var) < cons_obj->minvalue_)
            break;
        bool isFrac = SCIPisLT(scip, varLPval, 1);

        cons_obj->nCheckVars_++;

        auto* vardata = dynamic_cast<ObjVarDataVRP*>(SCIPgetObjVardata(scip, var));
        auto* obj = dynamic_cast<ObjPropTourVarFixing*>(SCIPgetObjProp(scip, SCIPfindProp(scip, "tourVarFixing")));
        obj->tvrp_.copy(vardata->tourVrp_);
        cons_obj->obj_cmp_ = SCIPgetCutoffbound(scip);
        cons_obj->fixVar_ = var;

        SCIP_Bool success = false;
        SCIP_Bool lperror, cutoff;
        obj->up_ = true;
        obj->new_ = false;
        cons_obj->used_Pricing_ = false;
        if(probData->use_propagator_ == 1)
        {

            auto start = std::chrono::high_resolution_clock::now();
            SCIP_CALL( heuristicExtended(scip, cons_obj, col, index, gap_abs, varLPval, &success));

            auto end = std::chrono::high_resolution_clock::now();
            if(success){
                cons_obj->time_success_ += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            }else{
                cons_obj->time_fail_ += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            }

	    }else if(probData->use_propagator_ == 2)
        {
            auto start = std::chrono::high_resolution_clock::now();

            SCIP_CALL( heuristicNoMove(scip, cons_obj, col, index, gap_abs, varLPval, &success));

            auto end = std::chrono::high_resolution_clock::now();
            if(success){
                cons_obj->time_success_ += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            }else{
                cons_obj->time_fail_ += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            }
        }
        else{
            obj->new_ = true;

            SCIP_CALL( exactBranching(scip, cons_obj, var, &success) );
        }

        /* Can variable be fixed? */
        if(success)
        {
            nSuccess++;
            cons_obj->nFixed_++;
            if(isFrac)
            {
                fixed_frac_ = true;
                cons_obj->nFixed_frac_++;
            }

            /* Check if Node can be cut off - Currently deactivated */
            if(false && varLPval <= 1 - cons_obj->minvalue_)
            {
                cons_obj->n_try_cutoff_++;
                if(probData->use_propagator_ == 1)
                {
                    SCIPstartProbing(scip);
                    cons_obj->obj_cmp_ = cmp_low;
                    SCIPchgVarObjProbing(scip, var, SCIPvarGetObj(var) - gap_abs/(1 - varLPval));
                    SCIPsolveProbingLPWithPricing(scip, false, false, -1, &lperror, &cutoff);
                }else{
                    SCIPstartProbing(scip);
                    SCIPchgVarLb(scip, var, 1.0);
                    SCIPsolveProbingLPWithPricing(scip, false, false, -1, &lperror, &cutoff);
                }
                if((!obj->new_ && SCIPisLT(scip, cmp_low, SCIPgetLPObjval(scip))) || (obj->new_ && cutoff))
                {
                    *result = SCIP_CUTOFF;
                    SCIPendProbing(scip);
//                    std::cout << std::endl;
                    return SCIP_OKAY;
                }
                assert(!lperror);
                SCIPendProbing(scip);
            }
            fixableVars.push_back(var);
        }else
        {
            if(!cons_obj->withCuts_){
                if(SCIPisLT(scip, varLPval, std::min(2 * cons_obj->minvalue_, 1.0))){
                    break;
                }
            }else if(cons_obj->used_Pricing_)
            {
                cons_obj->cutSize_.push_back(cons_obj->nAdded_Vars_);
                cutVars.emplace_back();
                cutVars[cutVars.size()-1].push_back(var);
                violatedCut.push_back(isFrac);
                for(int i = 0; i < cons_obj->nAdded_Vars_; i++)
                {
                    cutVars[cutVars.size()-1].push_back(probData->vars_[probData->nVars_ - 1 - i]);
                }
            }else if(SCIPisLT(scip, varLPval, std::min(2 * cons_obj->minvalue_, 1.0)))
                break;
        }
    }

    if(!fixableVars.empty())
    {
        *result = SCIP_REDUCEDDOM;
        SCIP_CALL(fixVariables(scip, fixableVars));
    }
    if(!cutVars.empty())
    {
        cons_obj->nCuts_ += cutVars.size();
        assert(cutVars.size() == violatedCut.size());
        *result = SCIP_CONSADDED;
        cons_obj->nonViolatedCut_ = fixableVars.empty();
        for(int i = 0; i < cutVars.size(); i++)
        {
            createCPCCut(scip, cutVars[i]);
            if(violatedCut[i])
            {
                *result = SCIP_SEPARATED;
                cons_obj->nonViolatedCut_ = false;
            }
        }
    }
    if(*result != SCIP_SEPARATED && !fixed_frac_)
    {
        cons_obj->noChange_ = true;
    }

    if(nSuccess > 0)
    {
        cons_obj->depths_.push_back(SCIPgetDepth(scip));
        cons_obj->nVars_.push_back((int) checkLPIcols.size());
        cons_obj->success_.push_back(nSuccess);
        cons_obj->gaps_.push_back(gap_rel);
    }

    return SCIP_OKAY;
}


SCIP_DECL_CONSDELETE(ConshdlrRCFC::scip_delete)
{
    assert(conshdlr != nullptr);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), "CPC") == 0);
    assert(consdata != nullptr);
    assert(*consdata != nullptr);

    /* delete constraint data */
    SCIPreleaseRow(scip, &(*consdata)->cut);
    SCIPfreeBlockMemory(scip, consdata);
    SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->nvars);

    return SCIP_OKAY;
}

SCIP_DECL_CONSENFOLP(ConshdlrRCFC::scip_enfolp)
{
//    return SCIP_OKAY;
    if(stopfixing_)
        return SCIP_OKAY;
    if(SCIPinProbing(scip))
        return SCIP_OKAY;
    if(SCIPinDive(scip))
        return SCIP_OKAY;
    if(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_NOTSOLVED || SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_INFEASIBLE)
        return SCIP_OKAY;

    /* Our SCIP installation has issue if in probing mode when hitting the time limit */
    if(SCIPgetSolvingTime(scip) > 3590){
        return SCIP_OKAY;
    }
    auto *probData = dynamic_cast<vrp::ProbDataVRP *>(SCIPgetObjProbData(scip_));
    if(!probData->use_propagator_)
        return SCIP_OKAY;

    double gap_abs = SCIPgetCutoffbound(scip) - SCIPgetLPObjval(scip);
    double gap_rel = 100 * abs(gap_abs / SCIPgetLPObjval(scip));

    bool onlyFrac = false;

    if(gap_rel > 50)
        return SCIP_OKAY;
    if(SCIPgetDepth(scip) > std::max(depth_init_, SCIPgetMaxDepth(scip) / 2))
        return SCIP_OKAY;
    /** If expected sub tree is to small -> only try fixing fractional variables
     * since impact of the other variables comes too late */
    if(onlyFrac_ && gap_rel < root_gap_ / onlyFracGap_){
        onlyFrac = true;
    }


    auto* pricerData = dynamic_cast<ObjPricerVRP*>(SCIPfindObjPricer(scip_, "VRP_Pricer"));
    /* do not call this method too often at the same branching node */
    if(lastNode_ == SCIPnodeGetNumber(SCIPgetFocusNode(scip)))
    {
        /* if no fractional variable or violated cut got added in last iteration, skip */
        if(noChange_)
            return SCIP_OKAY;

        if((lastNode_ == 1 && nCalls_Node_ >= max_nCallsRoot_) ||
           (lastNode_ != 1 && nCalls_Node_ >= max_nCalls_))
        {
            nCalls_Node_ = 0;
            return SCIP_OKAY;
        }
        nCalls_Node_++;
    }else{
        nCalls_Node_ = 1;

        /* set lowest fail factor from parents data saved in tree data */
        lowest_fail_ = init_fail_factor_ * pricerData->tree_data_[SCIPnodeGetNumber(SCIPgetFocusNode(scip))].node_lowest_fail;
    }

    lastNode_ = SCIPnodeGetNumber(SCIPgetFocusNode(scip));

    assert(SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_OPTIMAL);

    SCIP_COL** lpicols = SCIPGetlpiColsRCFC(scip);
    std::vector<int> checkLPIcols;
    for(int i = 0; i < SCIPGetNlpiColsRCFC(scip); i++)
    {
        SCIP_Col* col = lpicols[i];
        if(SCIPcolGetLb(col) > 0.5)
            continue;
        if(onlyFrac && SCIPisGE(scip, SCIPcolGetPrimsol(col), 0.9))
            continue;
        if(SCIPisGE(scip, SCIPcolGetPrimsol(col), minvalue_))
        {
            checkLPIcols.push_back(i);
        }
    }

    std::sort(checkLPIcols.begin(), checkLPIcols.end(), [lpicols](auto &left, auto &right){
        return SCIPcolGetPrimsol(lpicols[left]) > SCIPcolGetPrimsol(lpicols[right]);
    });

    *result = SCIP_FEASIBLE;

    if(probData->nonzeroVars_.empty())
    {
        SCIP_CALL(SCIPsortNonzeroVars(scip, probData));
    }

    nonViolatedCut_ = false;
    noChange_ = false;

    SCIP_CALL(checkFixing(scip, this, checkLPIcols, result, onlyFrac));

    /* update tree data for lowest fail */
    pricerData->tree_data_[SCIPnodeGetNumber(SCIPgetFocusNode(scip))].node_lowest_fail = lowest_fail_;

    return SCIP_OKAY;
}

SCIP_DECL_CONSPRINT(ConshdlrRCFC::scip_print)
{
    return SCIP_OKAY;
}

/** creates and adds the row of a subset row constraint to the LP */
SCIP_RETCODE SCIPcreateAndAddRowCPC(
        SCIP*                   scip,
        SCIP_Cons*              cons
){
    SCIP_CONSDATA* consdata;
    assert(scip != nullptr);
    assert(cons != nullptr);
    SCIP_Bool infeasible;
    consdata = SCIPconsGetData(cons);
    assert(consdata->cut == nullptr);

    SCIP_CALL(SCIPcreateEmptyRowCons(scip, &(consdata->cut), cons, SCIPconsGetName(cons), 1,
                                     SCIPinfinity(scip), TRUE, FALSE, FALSE));

    assert(consdata->cut != nullptr);

    SCIP_CALL(addVarsToCons(scip, consdata));

    if(!SCIProwIsInLP(consdata->cut))
    {
        SCIP_CALL(SCIPaddRow(scip, consdata->cut, true, &infeasible));
        assert(!infeasible);
    }

    return SCIP_OKAY;
}

/** creates and captures a counterpart cut constraint */
SCIP_RETCODE SCIPcreateConsCPC(
        SCIP*                   scip,               /**< SCIP data structure */
        SCIP_CONS**             cons,               /**< pointer to hold the created constraint */
        const char*             name,               /**< name of constraint */
        vector<SCIP_VAR*>&      vars
){
    SCIP_CONSDATA* consdata = nullptr;
    SCIP_CONSHDLR* conshdlr = nullptr;

    assert(scip != nullptr);

    /* find the set partitioning constraint handler */
    conshdlr = SCIPfindConshdlr(scip, "CPC");
    assert(conshdlr != nullptr);
    /* create constraint data*/
    SCIP_CALL(consdataCreate(scip, &consdata, vars));
    assert(consdata != nullptr);
    /* create constraint */
    SCIP_CALL(SCIPcreateCons(scip, cons, name, conshdlr, consdata,FALSE, TRUE, TRUE, FALSE,
                             FALSE, TRUE, FALSE, FALSE, FALSE, FALSE));

    return SCIP_OKAY;
}
