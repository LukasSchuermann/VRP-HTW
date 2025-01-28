

#include "branchingrule_dayvar.h"
using namespace std;

static
SCIP_RETCODE setTreeData(
        SCIP*               scip,
        ObjPricerVRP*       pricerData
){
    long long int currNode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
    assert(!pricerData->tree_data_[currNode].visitedChild);
    if(currNode == 1)
    {
        pricerData->tree_data_[currNode].node_timetable_ = pricerData->global_timetable_;
        pricerData->tree_data_[currNode].node_isForbidden_ = pricerData->global_isForbidden_;
    }else
    {
        pricerData->tree_data_[currNode].node_timetable_ = pricerData->timetable_;
        pricerData->tree_data_[currNode].node_isForbidden_ = pricerData->isForbidden_;
    }
    pricerData->tree_data_[currNode].node_ng_DSSR_ = pricerData->ng_DSSR_;
    pricerData->tree_data_[currNode].num_vars = SCIPgetNVars(scip);
    pricerData->tree_data_[currNode].node_fixedDay_ = pricerData->fixedDay_;
    pricerData->tree_data_[currNode].node_varfixing_ = pricerData->varfixing_gap_;

    return SCIP_OKAY;
}

SCIP_DECL_BRANCHINIT(ObjBranchruleDayVar::scip_init){
    auto* probData = dynamic_cast<vrp::ProbDataVRP*>(SCIPgetObjProbData(scip_));
    model_data* modelData = probData->getData();
    valonDay_ = vector<vector<SCIP_Real>>(modelData->nC, vector<SCIP_Real>(modelData->nDays, 0.0));
    return SCIP_OKAY;
}

SCIP_DECL_BRANCHEXECPS(ObjBranchruleDayVar::scip_execlp)
{
    auto* probData = dynamic_cast<vrp::ProbDataVRP*>(SCIPgetObjProbData(scip_));
    auto* pricerData = dynamic_cast<ObjPricerVRP*>(SCIPfindObjPricer(scip_, "VRP_Pricer"));
    model_data* modelData = probData->getData();
    vector<int> numofDays(modelData->nC);
    int i, j;
    int best_cust = -1;
    int best_day = -1;
    SCIP_Real cur_val;
    SCIP_Real best_val = std::max(branchingFrac_, 1 - branchingFrac_);
    int bestnum = 0;
    bool useRandom = (branchingFrac_ == 0.0 || branchingFrac_ == 1.0);

    SCIP_NODE* childprohibit;
    SCIP_NODE* childenforce;
    SCIP_CONS* consprohibit;
    SCIP_CONS* consenforce;

    assert(probData != nullptr);
    assert(pricerData != nullptr);
    /* no Branching applied to node */
    assert(pricerData->tree_data_[SCIPnodeGetNumber(SCIPgetFocusNode(scip))].branchtype_ == BRANCHTYPE_NOTBRANCHED);

    /* update nodeData */
    SCIP_CALL(setTreeData(scip, pricerData));

    /* count branching decisions that did not lead to cutoff */
    if(SCIPgetDepth(scip) != 0){
        SCIP_Cons** addedCons = nullptr;
        SCIP_CALL(SCIPallocBlockMemoryArray(scip, &addedCons, 1));
        int nAddedCons = 0;
        SCIPnodeGetAddedConss(SCIPgetCurrentNode(scip), addedCons, &nAddedCons, 1);
        assert(nAddedCons == 1);
        assert(addedCons != nullptr);
        auto type = SCIPgetTypeOfDayCons(addedCons[0]);
        if(type == PROHIBIT){
            nProhibit_++;
        }else {
            assert(type == ENFORCE);
            nEnforce_++;
        }
        SCIPfreeBlockMemoryArray(scip, &addedCons, 1);
    }

    if(PRINT_BRANCHING_INFORMATION)
    {
        cout << "Node "<< SCIPnodeGetNumber(SCIPgetCurrentNode(scip))
        << " solved, but still fractional - start: DAY VAR BRANCHING!\n";
    }
    /* clean up valonDay */
    for(i = 0; i < modelData->nC; i++) {
        for (j = 0; j < modelData->nDays; j++) {
            valonDay_[i][j] = 0.0;
        }
    }

    /* get new vehicle assignment variable values */
    SCIP_CALL(getVehiAssValues(scip, valonDay_, numofDays));

    std::vector<std::pair<std::pair<int, int>, double>> candidates;

    for(i = 0; i < modelData->nC; i++)
    {
        for(j = 0; j < modelData->nDays; j++)
        {
            // TODO: Big numerical issues
            if(SCIPisSumLE(scip, valonDay_[i][j], 0 + 1e-5) || SCIPisSumGE(scip, valonDay_[i][j], 1 - 1e-5))
                continue;
            candidates.emplace_back(std::pair<int, int>(i, j), valonDay_[i][j]);
            cur_val = fabs(branchingFrac_ - valonDay_[i][j]);
            /* day variable must not be 'worse' than current best day variable */
            if(!SCIPisLE(scip, cur_val, best_val))
                continue;
            /* if equal, then compare number of days the customer gets served on */
            if(SCIPisEQ(scip, cur_val, best_val))
            {
                if(numofDays[i] <= bestnum)
                    continue;
            }
            /* new best variable found */
            best_cust = i;
            best_day = j;
            bestnum = numofDays[i];
            best_val = cur_val;
        }
    }

    if(candidates.empty()){
        *result = SCIP_DIDNOTFIND;
        return SCIP_OKAY;
    }

    if(useRandom){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dis(0, (int) candidates.size() - 1);

        int chosen_idx = dis(gen);
        best_cust = candidates[chosen_idx].first.first;
        best_day = candidates[chosen_idx].first.second;
        best_val = candidates[chosen_idx].second;
    }


    if(best_cust == -1)
    {
        return SCIP_OKAY;
    }


    assert(best_cust >= 1 && best_cust <= modelData->nC - 1);
    assert(best_day >= 0 && best_day <= modelData->nDays);
    if(PRINT_BRANCHING_INFORMATION)
    {
        cout << "\tChosen Var: (Cust: " << best_cust << ", day:" << best_day << ") with value: ";
        cout << valonDay_[best_cust][best_day] << " (num days: " << bestnum << ")\n";
    }

    /* create new nodes */
    SCIP_CALL( SCIPcreateChild(scip, &childprohibit, 0.0, SCIPgetLocalTransEstimate(scip)) );
    SCIP_CALL( SCIPcreateChild(scip, &childenforce, 0.0, SCIPgetLocalTransEstimate(scip)) );

    /* create corresponding constraints */
    SCIP_CALL( SCIPcreateConsDayVar(scip, &consprohibit, "prohibit", best_cust, best_day, PROHIBIT, childprohibit, TRUE) );
    SCIP_CALL( SCIPcreateConsDayVar(scip, &consenforce, "enforce", best_cust, best_day, ENFORCE, childenforce, TRUE) );

    /* add constraints to nodes */
    SCIP_CALL( SCIPaddConsNode(scip, childprohibit, consprohibit, nullptr) );
    SCIP_CALL( SCIPaddConsNode(scip, childenforce, consenforce, nullptr) );

    /* release constraints */
    SCIP_CALL( SCIPreleaseCons(scip, &consprohibit) );
    SCIP_CALL( SCIPreleaseCons(scip, &consenforce) );

    pricerData->tree_data_[SCIPnodeGetNumber(SCIPgetFocusNode(scip))].branchtype_ = BRANCHTYPE_DAYVAR;

    *result = SCIP_BRANCHED;

//    assert(FALSE);
    return SCIP_OKAY;
}