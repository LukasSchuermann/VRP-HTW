
#include "branchingrule_arcflow.h"

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

SCIP_DECL_BRANCHINIT(ObjBranchruleArcflow::scip_init){
    auto* probData = dynamic_cast<vrp::ProbDataVRP*>(SCIPgetObjProbData(scip_));
    model_data* modelData = probData->getData();
    arcWeights_ = vector<vector<SCIP_Real>>(modelData->nC, vector<SCIP_Real>(modelData->nC, 0.0));
    return SCIP_OKAY;
}

SCIP_DECL_BRANCHEXECPS(ObjBranchruleArcflow::scip_execlp)
{
    auto* probData = dynamic_cast<vrp::ProbDataVRP*>(SCIPgetObjProbData(scip_));
    auto* pricerData = dynamic_cast<ObjPricerVRP*>(SCIPfindObjPricer(scip_, "VRP_Pricer"));
    model_data* modelData = probData->getData();

    int narcs;
    int i, j;
    int start = 0;
    int end = 0;
    SCIP_Real bestval = std::max(branchingFrac_, 1 - branchingFrac_);
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
        auto type = SCIPgetTypeArcflow(addedCons[0]);
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
             << " solved, but still fractional - start: ARC FLOW BRANCHING!\n";
    }

    /* clean up arc weights */
    for(i = 0; i < modelData->nC; i++) {
        for (j = 0; j < modelData->nC; j++) {
            arcWeights_[i][j] = 0.0;
        }
    }

    SCIP_CALL(getArcFlowValues(scip, arcWeights_, &narcs));

    std::vector<std::pair<std::pair<int, int>, double>> candidates;

    for(i = 0; i < modelData->nC; i++)
    {
        for(j = 0; j < modelData->nC; j++)
        {
            if(SCIPisSumGT(scip, arcWeights_[i][j], 1e-5) && SCIPisSumLT(scip, arcWeights_[i][j], 1 - 1e-5))
            {
                candidates.emplace_back(std::pair<int, int>(i, j), arcWeights_[i][j]);
                if(fabs(branchingFrac_ - arcWeights_[i][j]) < bestval)
                {
                    bestval = fabs(branchingFrac_ - arcWeights_[i][j]);
                    start = i;
                    end = j;
                }
            }
        }
    }

    if(useRandom){
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dis(0, (int) candidates.size() - 1);

        int chosen_idx = dis(gen);
        start = candidates[chosen_idx].first.first;
        end = candidates[chosen_idx].first.second;
        bestval = candidates[chosen_idx].second;
    }


    assert(start >= 0 && start <= modelData->nC - 1);
    assert(end >= 0 && end <= modelData->nC - 1);
    if(PRINT_BRANCHING_INFORMATION)
        cout << "\tChosen arc: (" << start << "-" << end << ") with value: " << bestval << '\n';

    /* create new nodes */
    SCIP_CALL( SCIPcreateChild(scip, &childprohibit, 0.0, SCIPgetLocalTransEstimate(scip)) );
    SCIP_CALL( SCIPcreateChild(scip, &childenforce, 0.0, SCIPgetLocalTransEstimate(scip)) );

    /* create corresponding constraints */
    SCIP_CALL( SCIPcreateConsArcFlow(scip, &consprohibit, "prohibit", start, end, PROHIBIT, childprohibit, TRUE) );
    SCIP_CALL( SCIPcreateConsArcFlow(scip, &consenforce, "enforce", start, end, ENFORCE, childenforce, TRUE) );

    /* add constraints to nodes */
    SCIP_CALL( SCIPaddConsNode(scip, childprohibit, consprohibit, nullptr) );
    SCIP_CALL( SCIPaddConsNode(scip, childenforce, consenforce, nullptr) );

    /* release constraints */
    SCIP_CALL( SCIPreleaseCons(scip, &consprohibit) );
    SCIP_CALL( SCIPreleaseCons(scip, &consenforce) );

    pricerData->tree_data_[SCIPnodeGetNumber(SCIPgetFocusNode(scip))].branchtype_ = BRANCHTYPE_ARCFLOW;

    *result = SCIP_BRANCHED;

//    assert(FALSE);
    return SCIP_OKAY;

}