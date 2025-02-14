
#include "pricer_vrp.h"

#include <iostream>
#include <vector>

#include "model_data.h"
#include "scip/cons_setppc.h"
#include "scip/cons_linear.h"
#include "probdata_vrp.h"
#include "labeling_algorithm_vrp.h"
#include "local_search_pricing.h"
#include "vardata.h"
#include "tools_vrp.h"
#include "ConshdlrArcflow.h"
#include "ConshdlrDayVar.h"
#include "ConshdlrKPC.h"
#include "ConshdlrSRC.h"
#include "ConshdlrNVehicle.h"
#include "ConshdlrVehicle.h"
#include "heurDayVarRounding.h"
#include "printer.h"
#include "ctime"

using namespace std;
using namespace scip;

/**@name Local methods
 *
 * @{
 */

/** process the vehicle branching decisions to set enforced days */
static
SCIP_RETCODE setEnforcedDays(
        ObjPricerVRP*           pricerData
){
    SCIP_CONS** conss;
    SCIP_CONS* cons;
    int ncons;
    int day;
    int i;
    CONSTYPE type;

    /* reset enforced days */
    pricerData->eDays_.clear();


    conss = SCIPconshdlrGetConss(pricerData->cons_nVehicle_);
    ncons = SCIPconshdlrGetNConss(pricerData->cons_nVehicle_);

    for(i = 0; i < ncons; i++)
    {
        cons = conss[i];
        if (!SCIPconsIsActive(cons))
            continue;
        if (PRINT_BRANCHING_INFORMATION)
            cout << "\tNumber of Vehicles " << (SCIPgetTypeOfNVehicleCons(cons) == ENFORCE ? ">= " : "<= ") <<
                 SCIPgetnumVOfNVehicleCons(cons) << endl;
    }

    /* get all vehicle constraints */
    conss = SCIPconshdlrGetConss(pricerData->cons_vehicle_);
    ncons = SCIPconshdlrGetNConss(pricerData->cons_vehicle_);
    /* collect all branching decisions and update fixedDay */
    for(i = 0; i < ncons; i++)
    {
        cons = conss[i];

        /* Checks if current constraint is active. If not it can be ignored. */
        if(!SCIPconsIsActive(cons))
            continue;

        day = SCIPgetDayOfVehicleCons(cons);
        type = SCIPgetTypeOfVehicleCons(cons);

        assert(day >= 0);

        if(PRINT_BRANCHING_INFORMATION)
            cout << "\tday = " << day << ", type = " << (type == PROHIBIT ? "prohibit" : "enforce") << endl;

        if(type == PROHIBIT)
        {
            pricerData->fixedDay_[day] = true;
        }else
        {
            pricerData->eDays_[day] = cons;
        }
    }

    return SCIP_OKAY;
}

/** process the day var branching decisions */
static
SCIP_RETCODE setTimeTable(
        ObjPricerVRP*           pricerData,
        model_data*             modelData
){
    SCIP_CONS** conss;
    SCIP_CONS* cons;
    int ncons;
    int customer, day;
    CONSTYPE type;
    int i;

    /* get all day var constraints */
    conss = SCIPconshdlrGetConss(pricerData->cons_dayvar_);
    ncons = SCIPconshdlrGetNConss(pricerData->cons_dayvar_);
    /* collect all branching decisions and update timetable */
    for(i = 0; i < ncons; i++)
    {
        cons = conss[i];

        /* Checks if current constraint is active. If not it can be ignored. */
        if(!SCIPconsIsActive(cons))
            continue;

        customer = SCIPgetCustomerOfDayCons(cons);
        day = SCIPgetDayOfDayCons(cons);
        type = SCIPgetTypeOfDayCons(cons);

        assert(customer > 0);
        assert(day >= 0);

        if(PRINT_BRANCHING_INFORMATION)
            cout << "\tcustomer = " << customer << ", day = " << day << ", type = " << (type == PROHIBIT ? "prohibit" : "enforce") << endl;

        /* update customer's row in timetable */
        if( type == PROHIBIT )
        {
            /* A customer should not be prohibited twice */
//            assert(pricerData->timetable_[customer][day]);

            pricerData->timetable_[customer][day] = false;
        }
        else if( type == ENFORCE )
        {
            /* A customer should only be enforced, if he is available on the given day
             *
             * edit: this can happen due to root reduced cost fixing */
//            if(!pricerData->timetable_[customer][day] &&)
//            {
//                cout << "customer " << customer << " on day " << day << endl;
//                for(auto d1 : modelData->availableDays[customer])
//                {
//                    cout << "day " << d1 << " tw: (" << modelData->timeWindows[customer][d1].start << " " <<
//                    modelData->timeWindows[customer][d1].end << ") avail: " << pricerData->timetable_[customer][d1] << endl;
//                }
//            }
//            assert(pricerData->timetable_[customer][day]); //

            /* Set row except for the enforced day to FALSE */
            for (int j : modelData->availableDays[customer]) {
                assert(0 <= j && j < modelData->nDays);
                pricerData->timetable_[customer][j] = false;
            }
            pricerData->timetable_[customer][day] = true;
        }
        else
        {
            SCIPerrorMessage("unknow constraint type <%d>\n", type);
            return SCIP_INVALIDDATA;
        }
    }

    return SCIP_OKAY;
}

static
SCIP_RETCODE includeVehicleArcs(
        ObjPricerVRP*           pricerData,
        model_data*             modelData
){
    SCIP_CONS** conss;
    SCIP_CONS* cons;
    int ncons;
    int i;

    /* get all day var constraints */
    conss = SCIPconshdlrGetConss(pricerData->cons_veharc_);
    ncons = SCIPconshdlrGetNConss(pricerData->cons_veharc_);

    if(ncons == 0){
        return SCIP_OKAY;
    }

    /* collect all branching decisions and update timetable */
    for(i = 0; i < ncons; i++) {
        cons = conss[i];

        /* Checks if current constraint is active. If not it can be ignored. */
        if (!SCIPconsIsActive(cons))
            continue;

        int tail = SCIPgetTailvehiclearc(cons);
        int head = SCIPgetHeadvehiclearc(cons);
        int day = SCIPgetDayvehiclearc(cons);
        CONSTYPE type = SCIPgetTypevehiclearc(cons);

        assert(tail >= 0 && head >= 0);
        assert(day >= 0);

        if(PRINT_BRANCHING_INFORMATION)
            cout << "\tarc = (" << tail << ", " << head << ") on day = " << day << ", type = " << (type == PROHIBIT ? "prohibit" : "enforce") << endl;

        if( type == PROHIBIT ){
            pricerData->isForbidden_[day][tail][head] = TRUE;
            if(head == 0){
                pricerData->toDepot_[tail] = FALSE;
            }
        }
        else if( type == ENFORCE ){
            /* Set timetable row except for the enforced day to FALSE if tail/head are not the depot */
            if(head != 0){
                for (auto j : modelData->availableDays[head]) {
                    assert(0 <= j && j < modelData->nDays);
                    pricerData->timetable_[head][j] = false;
                }
                pricerData->timetable_[head][day] = true;
            }
            if(tail != 0){
                for (auto j : modelData->availableDays[tail]) {
                    assert(0 <= j && j < modelData->nDays);
                    pricerData->timetable_[tail][j] = false;
                }
                pricerData->timetable_[tail][day] = true;
            }

            /* Every other outgoing arc from tail and ingoing arc to head must be forbidden on that day */
            for(int j = 0; j < modelData->nC; j++){
                if(j != head) {
                    pricerData->isForbidden_[day][tail][j] = TRUE;
                }
                if(j != tail) {
                    pricerData->isForbidden_[day][j][head] = TRUE;
                }
            }
            if(head != 0){
                pricerData->toDepot_[tail] = FALSE;
            }
        }
        else
        {
            SCIPerrorMessage("unknow constraint type <%d>\n", type);
            return SCIP_INVALIDDATA;
        }

    }


    return SCIP_OKAY;
}

/**
 * Sets pricerdata->eC_ and pricerdata->nEC based on data in timetable
 * If customer is enforced on a certain day, the entry will be set to that day. Else -1
 */
static
SCIP_RETCODE setEnforcedCustomers(
    ObjPricerVRP*           pricerData,
    model_data*             modelData
){
    int dayAvailable;

    std::vector<int> countAvailable(modelData->nC, 0);
    std::vector<int> enforcedCustomers;

    /* reset count */
    for(int day = 0; day < modelData->nDays; day++)
    {
        pricerData->nEC_[day] = 0;
    }
    /* find enforced customers */
    for(int customer = 1; customer < modelData->nC; customer++)
    {
        dayAvailable = -1;
        for(int day : modelData->availableDays[customer])
        {
            /* customer is allowed on that day */
            if(pricerData->timetable_[customer][day])
            {
                dayAvailable = day;
                countAvailable[customer]++;
            }else
            {
//                cout << "customer " << customer << " not allowed on day " << day << endl;
            }
        }
        /* countAvailable == 1 indicates, that the customer is enforced on dayAvailable */
        if(countAvailable[customer] == 1 && modelData->num_v[dayAvailable] == 1)
        {
            pricerData->eC_[customer] = dayAvailable;
            pricerData->nEC_[dayAvailable]++;
            enforcedCustomers.push_back(customer);
        }else
        {
            pricerData->eC_[customer] = -1;
        }
    }
    if(PRINT_BRANCHING_INFORMATION){
        std::cout << "Enforced Customer: ";
        for(int i = 1; i < modelData->nC; i++){
            if(pricerData->eC_[i] != -1){
                std::cout << i << " (" << pricerData->eC_[i] << "), ";
            }
        }
        std::cout << std::endl;
    }

    /* Conflict analysis: Try to check for unreachable neighbors of enforced customers */
    unsigned int nEnforced = enforcedCustomers.size();
    for(unsigned int i = 0; i < nEnforced; i++){
        int u = enforcedCustomers[i];
        int day = pricerData->eC_[u];
        /* check every customer v available on that day if we can reach u from that or vice versa */
        for(auto v : modelData->neighbors[0][day]){
            if(!pricerData->timetable_[v][day]){
                continue;
            }
            if(u == v){
                continue;
            }
            /* if both directions not reachable in time */
            if(modelData->adjacency_k[day][u][v] || modelData->adjacency_k[day][v][u]) {
                continue;
            }
            pricerData->nConflicts_++;
            /* can not be an enforced customer */
            if(pricerData->eC_[v] != -1){
                std::cout << "EC on day: " << pricerData->eC_[v] << std::endl;
                throw std::runtime_error("Shall not happen!");
            }
            if(countAvailable[v] <= 1){
                throw std::runtime_error("Infeasible NODE! Cannot happen!");
            }
            pricerData->timetable_[v][day] = false;
            countAvailable[v]--;
            /* if there is only one available day for v left --> set as enforced customer */
            if(countAvailable[v] == 1){
                for(auto d : modelData->availableDays[v]){
                    if(pricerData->timetable_[v][d]){
                        pricerData->eC_[v] = d;
                        pricerData->nEC_[d]++;
                    }
                }
                enforcedCustomers.push_back(v);
                nEnforced++;
            }
        }
    }


    return SCIP_OKAY;
}


/** process the arc flow branching decisions */
static
SCIP_RETCODE setArcMatrix(
    ObjPricerVRP*           pricerData,
    model_data*             modelData
){
    SCIP_CONS** conss;
    SCIP_CONS* cons;
    int ncons;
    int tail, head;
    CONSTYPE type;
    int i, j;



    assert(pricerData->cons_arcflow_ != nullptr);

    /* safe prohibited and enforced arcs */
    for(i = 0; i < modelData->nC; i++){
        pricerData->toDepot_[i] = TRUE;
    }

    /* get all arc flow constraints */
    conss = SCIPconshdlrGetConss(pricerData->cons_arcflow_);
    ncons = SCIPconshdlrGetNConss(pricerData->cons_arcflow_);

    if(ncons == 0){
        return SCIP_OKAY;
    }
    /* safe prohibited and enforced arcs */
    vector<int> successor(modelData->nC);
    vector<int> predecessor(modelData->nC);
    for(i = 0; i < modelData->nC; i++){
        successor[i] = -1;
        predecessor[i] = -1;
    }

    for(i = 0; i < ncons; i++)
    {
        cons = conss[i];

        if(!SCIPconsIsActive(cons))
            continue;

        tail = SCIPgetTailArcflow(cons);
        head = SCIPgetHeadArcflow(cons);
        type = SCIPgetTypeArcflow(cons);

        assert(tail != head);
        if(PRINT_BRANCHING_INFORMATION)
            cout << "\tarc = (" << tail << ',' << head << "), type = " << (type == PROHIBIT ? "prohibit" : "enforce") << endl;

        if( type == PROHIBIT )
        {
            for(int day = 0; day < modelData->nDays; day++){
                pricerData->isForbidden_[day][tail][head] = TRUE;
            }
            if(head == 0)
            {
                pricerData->toDepot_[tail] = FALSE;
            }
        }
        else if( type == ENFORCE )
        {
            /* there can only be one enforced out/in-going arc for each customer */
            assert(tail == 0 || successor[tail] == -1);
            assert(head == 0 || predecessor[head] == -1);

            successor[tail] = head;
            predecessor[head] = tail;

            /* both are customers */
            if(head != 0 && tail != 0)
            {
                pricerData->toDepot_[tail] = FALSE;
                for(j = 0; j < modelData->nC; j++)
                {
                    /* prohibit every other arc out/in-going arc of tail/head */
                    if(j != head) {
                        for(int day = 0; day < modelData->nDays; day++){
                            pricerData->isForbidden_[day][tail][j] = TRUE;
                        }
                    }
                    if(j != tail) {
                        for(int day = 0; day < modelData->nDays; day++){
                            pricerData->isForbidden_[day][j][head] = TRUE;
                        }
                    }
                }
            }/* head is depot */
            else if(head == 0)
            {
                for(j = 1; j < modelData->nC; j++)
                {
                    for(int day = 0; day < modelData->nDays; day++){
                        pricerData->isForbidden_[day][tail][j] = TRUE;
                    }
                }
            }/* tail is depot */
            else
            {
                assert(tail == 0);
                for(j = 1; j < modelData->nC; j++)
                {
                    for(int day = 0; day < modelData->nDays; day++){
                        pricerData->isForbidden_[day][j][head] = TRUE;
                    }
                }
            }
        }
        else
        {
            SCIPerrorMessage("unknow constraint type <%d>\n", type);
            return SCIP_INVALIDDATA;
        }
    }

    return SCIP_OKAY;
}

/** Sets up neighborhood at the current branching node */
SCIP_RETCODE setCurrentNeighborhood(
    ObjPricerVRP*           pricerData,
    model_data*             modelData,
    bool                    atRoot
){
    assert(pricerData != nullptr);
    assert(modelData != nullptr);

    vector<vector<bool>>& timetable = atRoot ? pricerData->global_timetable_ : pricerData->timetable_;
    vector<vector<vector<bool>>>& isForbidden = atRoot ? pricerData->global_isForbidden_ : pricerData->isForbidden_;

    for(int u = 0; u < modelData->nC; u++)
    {
        for (int day = 0; day < modelData->nDays; day++)
        {
            /* clear old data */
            pricerData->neighbors_[u][day].clear();
            pricerData->predecessors_[u][day].clear();

            /* contine if customer is not available on day (possibly due to branching decisions) */
            if(!(timetable[u][day] && pricerData->global_timetable_[u][day]))
                continue;

            /* check for each neighbor if it is forbidden */
            for(auto nb : modelData->neighbors[u][day])
            {
                if(!(timetable[nb][day] && pricerData->global_timetable_[nb][day]))
                    continue;
                if(!(isForbidden[day][u][nb] || pricerData->global_isForbidden_[day][u][nb]))
                {
                    pricerData->neighbors_[u][day].push_back(nb);
                }
            }
            /* check for each predecessor if it is forbidden */
            for(auto pd : modelData->predecessors[u][day])
            {
                if(!(timetable[pd][day] && pricerData->global_timetable_[pd][day]))
                    continue;
                if(!(isForbidden[day][pd][u] || pricerData->global_isForbidden_[day][pd][u]))
                {
                    pricerData->predecessors_[u][day].push_back(pd);
                }
            }
        }
    }
    return SCIP_OKAY;
}

/** Sets the initial ng_dssr_data for the current branching node */
static
SCIP_RETCODE getParentsData(
    ObjPricerVRP*           pricerData,
    model_data*             modelData
){
    long long int parent_ID = SCIPnodeGetNumber(SCIPnodeGetParent(SCIPgetCurrentNode(pricerData->scip_)));
    node_data& nodeData = pricerData->tree_data_.at(parent_ID);

    /* initialize data */
    if(USE_DSSR)
        pricerData->ng_DSSR_ = nodeData.node_ng_DSSR_;
    pricerData->fixedDay_ = nodeData.node_fixedDay_;
    pricerData->varfixing_gap_ = nodeData.node_varfixing_;

    for(int u = 0; u < modelData->nC; u++)
    {
        for(int day = 0; day < modelData->nDays; day++)
        {
            pricerData->timetable_[u][day] = nodeData.node_timetable_[u][day] && pricerData->global_timetable_[u][day];
        }
        for(int v = u + 1; v < modelData->nC; v++)
        {
            for(int day = 0; day < modelData->nDays; day++){
                pricerData->isForbidden_[day][u][v] = nodeData.node_isForbidden_[day][u][v] || pricerData->global_isForbidden_[day][u][v];
                pricerData->isForbidden_[day][v][u] = nodeData.node_isForbidden_[day][v][u] || pricerData->global_isForbidden_[day][v][u];
            }
        }
    }
//        pricerData->timetable_ = nodeData.node_timetable_;
//        pricerData->isForbidden_[day] = nodeData.node_isForbidden_[day];

    /* check if the other child of parent node has already been processed
     * if true: delete parent node dssr data */
    if(nodeData.visitedChild)
    {
        pricerData->tree_data_.erase(parent_ID);
    }else
    {
        nodeData.visitedChild = true;
    }

    return SCIP_OKAY;
}

/** Handles information of cuts for the pricing problem */
static
SCIP_RETCODE includeCuts(
    ObjPricerVRP*       pricerData,
    model_data*         modelData,
    SCIP_Bool           isFarkas
){
    SCIP_CONS** conss;
    int ncons;
    int c;
    SCIP_Real dualvalue;

    /* robust constraints */
    conss = SCIPconshdlrGetConss(pricerData->cons_kpc_);
    ncons = SCIPconshdlrGetNConss(pricerData->cons_kpc_);
    for(c = 0; c < ncons; c++)
    {
        /* get dual values of corresponding cuts */
        if(isFarkas)
        {
            dualvalue = SCIPgetDualfarkasKPC(conss[c]);
        }else
        {
            dualvalue = SCIPgetDualsolKPC(conss[c]);
        }

        if(SCIPisEQ(pricerData->scip_, 0, dualvalue))
            continue;
        /* adjust prices for ingoing arcs of set S of conss[c] */
        vector<bool>& setS = SCIPgetSetOfKPC(conss[c]);
        for(int i = 0; i < modelData->nC; i++)
        {
            if(!setS[i]) /* i \in setS */
                continue;
            for(int j = 0; j < modelData->nC; j++)
            {
                if(setS[j]) /* j \notin setS */
                    continue;
                pricerData->arcPrices_[j][i] -= dualvalue;
            }
        }
    }
    return SCIP_OKAY;
}

/** include vehicle branching constrains dual values */
static
SCIP_RETCODE includeVehicleBranching(
    SCIP*               scip,
    vrp::ProbDataVRP*   probData,
    ObjPricerVRP*       pricerData,
    bool                isFarkas
){
    SCIP_CONS** conss;
    SCIP_CONS* cons;
    int c, ncons;
    int day;
    double val;

    conss = SCIPconshdlrGetConss(pricerData->cons_nVehicle_);
    ncons = SCIPconshdlrGetNConss(pricerData->cons_nVehicle_);
    pricerData->dual_nVehicle_ = 0.0;
    for(c = 0; c < ncons; c++)
    {
        cons = conss[c];
        if (!SCIPconsIsActive(cons))
            continue;
        if(isFarkas)
        {
            val = SCIPgetDualfarkasNVehicle(cons);
        }else
        {
            val = SCIPgetDualsolNVehicle(cons);
        }
        if(!SCIPisZero(scip, val))
        {
            assert(pricerData->dual_nVehicle_ == 0.0);
            pricerData->dual_nVehicle_ = val;
        }
    }

    conss = SCIPconshdlrGetConss(pricerData->cons_vehicle_);
    ncons = SCIPconshdlrGetNConss(pricerData->cons_vehicle_);
    for(c = 0; c < ncons; c++)
    {
        cons = conss[c];
        if(!SCIPconsIsActive(cons))
            continue;

        day = SCIPgetDayOfVehicleCons(cons);

        if(isFarkas)
        {
            val = SCIPgetDualfarkasVehicle(cons);
        }else
        {
            val = SCIPgetDualsolVehicle(cons);
        }
//        cout << "val: " << val << endl;
//        SCIPprintCons(scip, cons, nullptr);
        if(!SCIPisZero(scip, val))
        {
            assert(SCIPisZero(scip, pricerData->dualValues_[probData->getData()->nC + day]));
            pricerData->dualValues_[probData->getData()->nC + day] = val;
        }
    }

    return SCIP_OKAY;
}

/** check if primal heuristic found a new solution that shall be added */
static
SCIP_RETCODE checkPrimalHeuristic(
   SCIP*                scip,
   vrp::ProbDataVRP*    probData,
   bool*                success
){
    auto* heurData = dynamic_cast<HeurDayVarRounding *>(SCIPfindObjHeur(scip, "dayVarRounding"));
    if(SCIPisSumNegative(scip, heurData->best_ - SCIPgetUpperbound(scip)))
    {
        char algoName[] = "heurVarRounding";
        SCIP_Bool isfeasible;
        SCIP_SOL* sol;

        *success = true;

        SCIPcreateSol(scip, &sol, nullptr);
        SCIP_Real checksum = 0.0;
        for(auto& t : heurData->sol_tours_)
        {
            if(t.length_ > 0)
            {
                SCIP_Var* existingVar = nullptr;
                if(SCIPcontainsTourVar(scip, probData, &existingVar, t))
                {
                    assert(existingVar != nullptr);
                    SCIPsetSolVal(scip, sol, existingVar, 1.0);
                }else
                {
                    assert(existingVar == nullptr);
                    SCIP_CALL(add_tour_variable(scip, probData, FALSE, FALSE, algoName, t));
                    SCIPsetSolVal(scip, sol, probData->vars_[probData->nVars_ - 1], 1.0);
                }

                checksum += t.obj_;
            }
        }
        SCIP_CALL( SCIPtrySol(scip, sol, false, false, false, false, false, &isfeasible) );
        assert(isfeasible);
        assert(SCIPisEQ(scip, heurData->best_, checksum));
        assert(SCIPisEQ(scip, SCIPgetPrimalbound(scip), heurData->best_));
    }

    return SCIP_OKAY;
}

/** Process branching decisions of the new branching node to generate local instance */
SCIP_RETCODE ObjPricerVRP::set_current_graph(
    model_data*     modelData
)
{
    if(PRINT_BRANCHING_INFORMATION)
    {
        cout << "START PROCESSING OF NODE " <<  SCIPnodeGetNumber(SCIPgetCurrentNode(scip_)) <<
        " (parent-node = "<< SCIPnodeGetNumber(SCIPnodeGetParent(SCIPgetCurrentNode(scip_))) << ")" << endl;
        cout << "PB: " << SCIPgetPrimalbound(scip_) << ", DB: " << SCIPgetDualbound(scip_) << ", LP-Obj: "
        << SCIPgetLPObjval(scip_) << endl;
    }

    /* Set the initial parents data for the current branching node */
    SCIP_CALL(getParentsData(this, modelData));

    /* process vehicle branching decisions */
    SCIP_CALL(setEnforcedDays(this));

    /* process arc flow branching decisions */
    SCIP_CALL(setArcMatrix(this, modelData));

    /* process day var branching decisions */
    SCIP_CALL(setTimeTable(this, modelData));

    /* process vehicle arc branching decisions */
    SCIP_CALL(includeVehicleArcs(this, modelData));

    /* set up enforced customers based on timetable */
    SCIP_CALL(setEnforcedCustomers(this, modelData));

    /* get neighborhood based on branching decisions */
    SCIP_CALL(setCurrentNeighborhood(this, modelData, false));

    if(PRINT_BRANCHING_INFORMATION)
    {
        std::cout << "fixed days: ";
        for(int i = 0; i < modelData->nDays; i++)
        {
            std::cout << fixedDay_[i] << " ";
        }
        std::cout << std::endl;
    }
    return SCIP_OKAY;
}

/** Constructs the pricer object with the data needed */
ObjPricerVRP::ObjPricerVRP(
    SCIP*                                 scip,           /**< SCIP pointer */
    const char*                           p_name         /**< name of pricer */
    ):
    ObjPricer(scip, p_name, "Finds tour with negative reduced cost.", 0, TRUE),
    cons_arcflow_(nullptr),
    cons_dayvar_(nullptr),
    lastID_(1)
{}

/** Destructs the pricer object. */
ObjPricerVRP::~ObjPricerVRP()
= default;


/** initialization method of variable pricer (called after problem was transformed) */
SCIP_DECL_PRICERINIT(ObjPricerVRP::scip_init)
{
    auto* probData = dynamic_cast<vrp::ProbDataVRP*>(SCIPgetObjProbData(scip_));
    model_data* modelData = probData->getData();

    /* initialize graph for pricing */
    isForbidden_.resize(modelData->nDays, vector<vector<bool>>(modelData->nC, vector<bool>(modelData->nC, true)));
    global_isForbidden_.resize(modelData->nDays, vector<vector<bool>>(modelData->nC, vector<bool>(modelData->nC, true)));
    toDepot_.resize(modelData->nC);
    timetable_.resize(modelData->nC, vector<bool>(modelData->nDays, false));
    global_timetable_.resize(modelData->nC, vector<bool>(modelData->nDays, false));
    neighbors_.resize(modelData->nC, vector<vector<int>>(modelData->nDays));
    predecessors_.resize(modelData->nC, vector<vector<int>>(modelData->nDays));
    for(int u = 0; u < modelData->nC; u++)
    {
        toDepot_[u] = TRUE;
        for(int day = 0; day < modelData->nDays; day++)
        {
            /* init with global neighborhood */
            for(int nb : modelData->neighbors[u][day])
            {
                neighbors_[u][day].push_back(nb);
                isForbidden_[day][u][nb] = false;
                global_isForbidden_[day][u][nb] = false;
            }
            for(int pd: modelData->predecessors[u][day])
            {
                predecessors_[u][day].push_back(pd);
                isForbidden_[day][pd][u] = false;
                global_isForbidden_[day][pd][u] = false;
            }
            /* init based on time windows */
            if(modelData->timeWindows[u][day].end > 0)
            {
                timetable_[u][day] = TRUE;
                global_timetable_[u][day] = TRUE;
            }

            isForbidden_[day][u][0] = false;
            global_isForbidden_[day][u][0] = false;
        }
    }

    /* DSSR for ng-neighborhood */
    ng_DSSR_.resize(modelData->nC, bitset<neighborhood_size>());

    /* Arc prices */
    arcPrices_.resize(modelData->nC, vector<SCIP_Real>(modelData->nC));

    /* Day-Customer Assignment variables */
    dayVarRedCosts_.resize(modelData->nC, vector<SCIP_Real>(modelData->nDays, SCIP_DEFAULT_INFINITY));
    root_dayVarRedCosts_.resize(modelData->nC, vector<SCIP_Real>(modelData->nDays, 0.0));
    root_dayVarLPObj_.resize(modelData->nC, vector<SCIP_Real>(modelData->nDays, -SCIP_DEFAULT_INFINITY));
    arcRedCosts_.resize(modelData->nDays, vector<vector<SCIP_Real>>(modelData->nC, vector<SCIP_Real>(modelData->nC, SCIP_DEFAULT_INFINITY)));
    root_arcRedCosts_.resize(modelData->nDays, vector<vector<SCIP_Real>>(modelData->nC, vector<SCIP_Real>(modelData->nC, 0.0)));
    root_arcLPObj_.resize(modelData->nDays, vector<vector<SCIP_Real>>(modelData->nC, vector<SCIP_Real>(modelData->nC, -SCIP_DEFAULT_INFINITY)));

    /* Fixed days */
    fixedDay_.resize(modelData->nDays, false);
    eDays_.resize(modelData->nDays, nullptr);

    /* enforced customers */
    eC_.resize(modelData->nC);
    nEC_.resize(modelData->nDays);
    for(int u = 0; u < modelData->nC; u++)
    {
        if(modelData->availableDays[u].size() == 1 && modelData->num_v[modelData->availableDays[u][0]] == 1)
        {
            nEC_[modelData->availableDays[u][0]]++;
            eC_[u] = modelData->availableDays[u][0];
        }else
        {
            eC_[u] = -1;
        }
    }
    /* tree data */
    tree_data_[1] = node_data();
    tree_data_[1].node_lowest_fail = 10000;

    varfixing_gap_ = 0.1;
    dual_nVehicle_ = 0.0;

    dualValues_.resize(modelData->nC + modelData->nDays);

    /* constraint handler */
    cons_arcflow_ = SCIPfindConshdlr(scip, "arcflow");
    cons_dayvar_ = SCIPfindConshdlr(scip, "dayVar");
    cons_veharc_ = SCIPfindConshdlr(scip, "vehiclearc");
    cons_kpc_ = SCIPfindConshdlr(scip, "KPC");
    cons_src_ = SCIPfindConshdlr(scip, "SRC");
    cons_vehicle_ = SCIPfindConshdlr(scip, "Vehicle");
    cons_nVehicle_ = SCIPfindConshdlr(scip, "nVehicle");

    /* propagator */
    prop_varfixing_ = dynamic_cast<ObjPropVarFixing*>(SCIPgetObjProp(scip, SCIPfindProp(scip, "varFixing")));

    node_varsfixed_ = -1;
    fixed_nonzero_ = false;

    nSRC_ = 0;
    nnonzSRC_ = 0;

    atRoot_ = true;

    enumerate_ = false;

    return SCIP_OKAY;
} /*lint !e715*/

/** perform pricing */
SCIP_RETCODE ObjPricerVRP::pricing(
    SCIP*                 scip,               /**< SCIP data structure */
    bool                  isFarkas            /**< whether we perform Farkas pricing */
    )
{
    int nvars;
    bool calcRedCosts = false;
    auto* probData = dynamic_cast<vrp::ProbDataVRP*>(SCIPgetObjProbData(scip_));
    long long int currNode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

    /* include current branching node information */
    if(lastID_ != currNode)
    {
        atRoot_ = false;
        long long int parent_ID = SCIPnodeGetNumber(SCIPnodeGetParent(SCIPgetCurrentNode(scip)));
        bool sepfirst = false;
        if(tree_data_[parent_ID].branchtype_ >= 4)
        {
            sepfirst = true;
        }

        lastID_ = currNode;

        tree_data_[currNode] = node_data();

        tree_data_[currNode].node_lowest_fail = tree_data_[parent_ID].node_lowest_fail;

        set_current_graph(probData->getData());

        if(sepfirst)
        {
            if(SCIPisLE(scip, SCIPgetLPObjval(scip), SCIPgetPrimalbound(scip)))
            {
                return SCIP_OKAY;
            }
        }
//        assert(false);
    }

    /* check for primal heuristic solution to add */
    bool success = false;
    SCIP_CALL( checkPrimalHeuristic(scip, probData, &success));

    if(tree_data_[currNode].gotFixed)
    {
//        std::cout << "got fixed at " << currNode << " focus: " << SCIPnodeGetNumber(SCIPgetFocusNode(scip)) << std::endl;
        SCIP_CALL(setCurrentNeighborhood(this, probData->getData(), currNode == 1));
        tree_data_[currNode].gotFixed = false;
        return SCIP_OKAY;
    }
    if(node_varsfixed_ == currNode)
    {
//        std::cout << "got vixed at " << currNode << " focus: " << SCIPnodeGetNumber(SCIPgetFocusNode(scip))  << std::endl;
        SCIP_CALL(setCurrentNeighborhood(this, probData->getData(), currNode == 1));
        node_varsfixed_ = -1;
    }

    /* set dual values */
    SCIP_CALL( setDualValues(scip, probData, isFarkas));

    /* set arc prices */
    SCIP_CALL(setArcPrices(probData->getData(), isFarkas));

    /* check for subset row cuts */
    if(SCIPconshdlrGetNConss(cons_src_) > 0)
    {
        if(nSRC_ < SCIPconshdlrGetNConss(cons_src_))
        {
            nSRC_ = SCIPconshdlrGetNConss(cons_src_);
            SRC_dualv.resize(nSRC_);
            SRC_para_.resize(nSRC_);
            SRC_Set_.resize(nSRC_);
        }
        SCIP_CALL(includeSRCpricingData(scip, isFarkas));
//        return SCIP_OKAY;
    }

    /** START PRICING */
    /* Heuristic local search pricing */
    nvars = probData->nVars_;
    if(SCIPconshdlrGetNConss(cons_src_) == 0)
    {
        heuristic_pricing(scip, probData, isFarkas);
        if(nvars != probData->nVars_)
            return SCIP_OKAY;
    }

    /* Exact pricing */
    auto start = std::chrono::high_resolution_clock::now();
    nvars = probData->nVars_;
    bool stopEarly = false;
    SCIP_CALL( generate_tours(scip, probData, isFarkas, calcRedCosts, false, &stopEarly));
    auto end = std::chrono::high_resolution_clock::now();
    auto duration_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    labelTime_ += duration_ms.count();
    nLabelCalls_++;

    /* In case the lagrangian bound is higher than the primal bound we stop */
    if(stopEarly){
        return SCIP_OKAY;
    }

    /* perform reduced cost fixing runs */
    if(nvars == probData->nVars_)
    {
        probData->nonzeroVars_.clear();


        if(!SCIPisInfinity(scip, SCIPgetLPObjval(scip)) && SCIPisPositive(scip, (SCIPgetPrimalbound(scip) - SCIPgetLPObjval(scip))))
        {
            assert(SCIPisEQ(scip, SCIPgetCutoffbound(scip), SCIPgetPrimalbound(scip)));
            /* do not calculate if there is a new cutoffbound, and we are not at the root */
            if(currNode == 1 || SCIPisEQ(scip, SCIPgetCutoffbound(scip), prop_varfixing_->lastCutoff_))
            {
                /* if the gap decreased by at least 20% since the last variable fixing, repeat it */
                double local_gap = (SCIPgetCutoffbound(scip) - SCIPgetLPObjval(scip)) / SCIPgetLPObjval(scip);

//                std::cout << local_gap << " vs " << varfixing_gap_ * gap_decay_rc_ << " max src " << maxSRC_rc_ << std::endl;
                if(local_gap <= varfixing_gap_ * gap_decay_rc_ && SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) <= max_depth_rc_ &&
                   SCIPconshdlrGetNConss(cons_src_) <= maxSRC_rc_)
                {
                    if(prop_varfixing_->verbose_){
                        std::cout << "(Node: " << currNode << "): ";
                        std::cout << "Check for reduced costs of arc flow and vehicle assignment variables at gap: ";
                        std::cout << local_gap << " (abs: " << (SCIPgetCutoffbound(scip) - SCIPgetLPObjval(scip)) << ")" << std::endl;
                    }

                    start = std::chrono::high_resolution_clock::now();
                    varfixing_gap_ = local_gap;

                    calcRedCosts = true;
                    bool earlyStop = false;
                    SCIP_CALL( generate_tours(scip, probData, isFarkas, calcRedCosts, false, &earlyStop));

                    prop_varfixing_->active_redcosts_ = true;
                    for(int i = 0; i < probData->getData()->nDays; i++){
                        if(SCIPisNegative(scip, dayVarRedCosts_[0][i]))
                            cout << "NEGATIVE: " << dayVarRedCosts_[0][i] << " i: " << i << " nvar " << SCIPgetNVars(scip) <<  endl;
                        assert(!SCIPisNegative(scip, dayVarRedCosts_[0][i]));
                    }

                    end = std::chrono::high_resolution_clock::now();

                    duration_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                    labelRCTime_ += duration_ms.count();
                    nRCLabelCalls_++;
                }
            }
        }
    }

    return SCIP_OKAY;
}

/** Pricing of additional variables if LP is feasible. */
SCIP_DECL_PRICERREDCOST(ObjPricerVRP::scip_redcost)
{
   SCIPdebugMsg(scip, "call scip_redcost ...\n");

   /* set result pointer, see above */
   *result = SCIP_SUCCESS;

   /* call pricing routine */
   SCIP_CALL( pricing(scip, false) );

//   /* abort if too long at current node */
//   if(SCIPgetDepth(scip) > 0 && SCIPgetNPriceRounds(scip) > 20 && SCIPgetLPObjval(scip) < SCIPgetPrimalbound(scip)){
//       *stopearly = true;
//   }

   return SCIP_OKAY;
} /*lint !e715*/


/** Pricing of additional variables if LP is infeasible. */
SCIP_DECL_PRICERFARKAS(ObjPricerVRP::scip_farkas)
{
   SCIPdebugMsg(scip, "call scip_farkas ...\n");

    /* set result pointer, see above */
    *result = SCIP_SUCCESS;

   /* call pricing routine */
   SCIP_CALL( pricing(scip, true) );

   return SCIP_OKAY;
} /*lint !e715*/

/** receives dual values of current LP */
SCIP_RETCODE ObjPricerVRP::setDualValues(
        SCIP*               scip,
        vrp::ProbDataVRP*   probData,
        bool                isFarkas
){
    model_data* modelData = probData->getData();
    dualValues_[0] = 0.0;

    if(isFarkas)
    {
        for(int i = 1; i < modelData->nC; i++)
        {
            dualValues_[i] = SCIPgetDualfarkasLinear(scip, probData->cons_[i-1]);
        }
        for(int i = 0; i < modelData->nDays; i++)
        {
            if(modelData->num_v[i] == 1)
                dualValues_[i + modelData->nC] = SCIPgetDualfarkasSetppc(scip, probData->cons_[modelData->nC + i - 1]);
            else
                dualValues_[i + modelData->nC] = SCIPgetDualfarkasLinear(scip, probData->cons_[modelData->nC + i - 1]);
        }
        if(SCIPconshdlrGetNConss(cons_vehicle_) > 0 || SCIPconshdlrGetNConss(cons_nVehicle_) > 0)
        {
            SCIP_CALL(includeVehicleBranching(scip, probData, this, true));
        }
    }else
    {
        for(int i = 1; i < modelData->nC; i++)
        {
            dualValues_[i] = SCIPgetDualsolLinear(scip, probData->cons_[i-1]);
        }
        for(int i = 0; i < modelData->nDays; i++)
        {
            if(modelData->num_v[i] == 1)
                dualValues_[i + modelData->nC] = SCIPgetDualsolSetppc(scip, probData->cons_[modelData->nC + i - 1]);
            else
                dualValues_[i + modelData->nC] = SCIPgetDualsolLinear(scip, probData->cons_[modelData->nC + i - 1]);
        }
        if(SCIPconshdlrGetNConss(cons_vehicle_) > 0 || SCIPconshdlrGetNConss(cons_nVehicle_) > 0)
        {
            SCIP_CALL(includeVehicleBranching(scip, probData, this, false));
        }
    }
    return SCIP_OKAY;
}

/** sets the arc prices for current pricing iteration */
SCIP_RETCODE ObjPricerVRP::setArcPrices(
        model_data*         modelData,
        bool                isFarkas
){
    /* initialize for new dualvalues and isFarkas-state */
    for(int i = 0; i < modelData->nC; i++)
    {
        for(int j = i+1; j < modelData->nC; j++)
        {
            arcPrices_[i][j] = (!isFarkas)*modelData->travel[i][j] - dualValues_[j];
            arcPrices_[j][i] = (!isFarkas)*modelData->travel[j][i] - dualValues_[i];
        }
    }
    /* if there are robust cuts, include their dual information */
    if(SCIPconshdlrGetNConss(cons_kpc_) > 0)
    {
        SCIP_CALL(includeCuts(this, modelData, isFarkas));
    }
    return SCIP_OKAY;
}

/** includes subset row cut data into pricing */
SCIP_RETCODE ObjPricerVRP::includeSRCpricingData(
        SCIP*               scip,
        bool                isFarkas
)
{
    SCIP_CONS** conss;
    int ncons;
    int c;
    SCIP_Real dualvalue = 0;
    nnonzSRC_ = 0;

    /* robust constraints */
    conss = SCIPconshdlrGetConss(cons_src_);
    ncons = SCIPconshdlrGetNConss(cons_src_);
    for(c = 0; c < ncons; c++) // only consider new cuts
    {
        if(!SCIPconsIsActive(conss[c]))
            continue;
        /* get dual values of corresponding cuts */
        if (isFarkas) {
            dualvalue = SCIPgetDualfarkasSRC(conss[c]);
        } else {
            dualvalue = SCIPgetDualsolSRC(conss[c]);
        }
        /* we only need to consider cuts with negative dualvalue */
        if(SCIPisNegative(scip, dualvalue))
        {
            /* receive cut data */
            SRC_dualv[nnonzSRC_] = dualvalue;
            SRC_para_[nnonzSRC_] = SCIPgetpOfSRC(conss[c]);
            SRC_Set_[nnonzSRC_] = SCIPgetSetOfSRC(conss[c]);
            nnonzSRC_++;
        }
    }
    return SCIP_OKAY;
}

SCIP_Real ObjPricerVRP::getTourVRPredcosts(
        model_data*     modelData,
        tourVRP&        tvrp
){
    int day = tvrp.getDay();
    SCIP_Real srcVal;
    SCIP_Real redcosts = -dualValues_[modelData->nC + day] - dual_nVehicle_;

    for(int c = 0; c < nnonzSRC_; c++)
    {
        srcVal = 0;
        assert(SCIPisNegative(scip_, SRC_dualv[c]));
        assert(SRC_para_[c] == 0.5);
        for(auto u : tvrp.tour_)
        {
            if((*SRC_Set_[c])[u])
                srcVal += SRC_para_[c];
        }
        redcosts -= (floor(srcVal) * SRC_dualv[c]);
    }

    vector<int>& tour = tvrp.tour_;
    redcosts += arcPrices_[0][tour[0]];
    for(int i = 1; i < tvrp.length_; i++)
    {
        redcosts += arcPrices_[tour[i-1]][tour[i]];
    }
    redcosts += arcPrices_[tour[tvrp.length_-1]][0];

    return redcosts;
}

/** generates negative reduced cost tours (uses restricted shortest path dynamic programming algorithm) */
SCIP_RETCODE ObjPricerVRP::generate_tours(
    SCIP*                   scip,
    vrp::ProbDataVRP*       probData,
    bool                    isFarkas,
    bool                    getDayVarRed,
    bool                    isHeuristic,
    bool                    *stopEarly
)
{
    char algoName[] = "pricingLabel";

    if(!PARALLEL_LABELING)
    {
        for(int i = 0; i < probData->getData()->nDays; i++)
        {
            if(fixedDay_[i])
                continue;
            vector<tourVRP> bestTours;
            double rc = 0;
            std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;
            start_time_ = std::chrono::high_resolution_clock::now();

            SCIP_CALL(generateLabelsBiDir(scip, probData->getData(), this, bestTours, isFarkas, getDayVarRed, false, false, i, false, &rc));

            int count = 0;

            for(auto& tvrp : bestTours)
            {
                /* check for ng-path violations */
                if(USE_DSSR)
                {
                    if(violatesNGProperty(probData->getData()->ng_set, ng_DSSR_, tvrp.tour_))
                        continue;
                }
                assert(tvrp.length_ > 0);
                SCIP_CALL(add_tour_variable(scip, probData, isFarkas, FALSE, algoName, tvrp));


                count++;
                if(count > probData->getData()->nC / 4)
                {
                    break;
                }
            }
        }
    }else
    {
        SCIP_CALL( labelingAlgorithmnParallel(scip, probData, this, isFarkas, isHeuristic, getDayVarRed, stopEarly) );
    }
    return SCIP_OKAY;
}

/** generates negative reduced cost tour (uses local search heuristics) */
bool ObjPricerVRP::heuristic_pricing(
        SCIP*                       scip,
        vrp::ProbDataVRP*           probData,
        bool                   isFarkas
){
    localSearchPricing(scip, this, probData, dualValues_, isFarkas);

    return FALSE;
}
