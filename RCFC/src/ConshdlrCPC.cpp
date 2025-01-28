//
// Created by lukasschuermann on 29.09.23.
//

#include <iostream>

#include "scip/cons_linear.h"
#include "ConshdlrCPC.h"
#include "objscip/objscip.h"
#include "pricer_vrp.h"
#include "vardata.h"
#include "labeling_algorithm_vrp.h"

#define DEFAULT_MAXROUNDS             3 /**< maximal number of GMI separation rounds per node (-1: unlimited) */
#define DEFAULT_MAXROUNDSROOT        15 /**< maximal number of GMI separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXSEPACUTS          10 /**< maximal number of GMI cuts separated per separation round */
#define DEFAULT_MAXSEPACUTSROOT      25 /**< maximal number of GMI cuts separated per separation round in root node */

#define DEFAULT_AWAY              0.01 /**< minimal fractionality of a basic variable in order to try GMI cut - default */
#define DEFAULT_MAX_DYN          1.0e+6 /**< maximal valid range max(|weights|)/min(|weights|) of cut coefficients - default */
#define DEFAULT_MAX_SUPP_ABS        100 /**< maximum cut support - absolute value in the formula - default */
#define DEFAULT_MAX_SUPP_REL       0.05 /**< maximum cut support - relative value in the formula - default */
#define CONSHDLR_NAME   "CPC"


/**
 * UPWARDS = What if x = 1
 * DOWNWARDS = What if x = 0
 **/
enum CutDir
{
    UPWARDS = -1,
    DOWNWARDS = 1
};

static
bool getCPCFromRowAggr(
        SCIP*                 scip,               /**< pointer to the SCIP environment */
        ConshdlrCPC*          objcpc,             /**< pointer to the conshdlr object */
        int                   ncols,              /**< number of columns in the LP */
        int                   nrows,              /**< number of rows in the LP */
        SCIP_COL**            cols,               /**< columns of the LP */
        SCIP_ROW**            rows,               /**< rows of the LP */
        const double*         binvarow,
        int                   basisind,           /**< index of basic variable */
        const SCIP_Real*      binvrow,            /**< row of the basis inverse */
        SCIP_Real*            cutcoefs,           /**< array for cut elements in sparse format - must be of size ncols */
        int*                  cutind,             /**< array for indices of nonzero cut coefficients - must be of size ncols */
        int*                  cutnz,              /**< pointer to store number of nonzero elements in the cut */
        SCIP_Real*            cutrhs,             /**< pointer to store cut rhs */
        SCIP_Real*            workcoefs,          /**< working array of size ncols, allocated by caller for efficiency */
        const SCIP_Real       factor              /**< factor for the reduced cost comparison */
){
    SCIP_COL* col;
    SCIP_ROW* row;
    SCIP_Var* var;
    double rowElem;

    assert(scip != nullptr);
    assert(cols != nullptr);
    assert(rows != nullptr);
    assert(binvrow != nullptr);
    assert(cutcoefs != nullptr);
    assert(cutind != nullptr);
    assert(cutnz != nullptr);
    assert(cutrhs != nullptr);
    assert(workcoefs != nullptr);

    /* clear cutcoefs and initialize cutcoefs */
    BMSclearMemoryArray(workcoefs, ncols);
    *cutrhs = 1;
    workcoefs[SCIPcolGetLPPos(cols[basisind])] = 1;

    std::vector<double> max_factors_vars = std::vector<double>();
    std::vector<double> max_factors_slack = std::vector<double>();

    double red_costs = 0;
    int cnt = 0;
    /* check which original variables join the cut */
    for(int i = 0; i < ncols; i++)
    {
        col = cols[i];
        int c = SCIPcolGetLPPos(col);
        /* get variable and check if not fixed */
        var = SCIPcolGetVar(col);
        if(SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)))
            continue;

        col = SCIPvarGetCol(var);

        switch (SCIPcolGetBasisStatus(col))
        {
            case SCIP_BASESTAT_LOWER:
                /* Check if element in simplex row has correct sign */
                rowElem = binvarow[c];
                if(SCIPisPositive(scip, rowElem))
                {
                    /* check if column can be excluded due to optimality condition */
                    red_costs = SCIPgetVarRedcost(scip, var);
                    if(SCIPisLE(scip, factor * rowElem, red_costs))
                        continue;

                    /* set coefficient */
                    cnt++;
                    workcoefs[SCIPcolGetLPPos(col)] = 1;
                    if(!SCIPvarIsBinary(var))
                    {
                        *cutrhs += SCIPvarGetLbLocal(var);
                    }
                }
                continue;
            case SCIP_BASESTAT_UPPER:
                assert(false); // Should not happen for this VRP
                /* Check if element in simplex row has correct sign */
                rowElem = binvarow[c];
                if(SCIPisPositive(scip, -rowElem))
                {
                    /* check if column can be excluded due to optimality condition */
                    red_costs = SCIPgetVarRedcost(scip, var);
                    max_factors_vars.push_back(red_costs / rowElem);
                    if(SCIPisGE(scip, factor * rowElem, red_costs))
                        continue;
                    /* set coefficient and rhs */
                    cnt++;
                    workcoefs[SCIPcolGetLPPos(col)] = -1;
                    if(!SCIPvarIsBinary(var))
                    {
                        *cutrhs -= SCIPvarGetUbLocal(var);
                    }else
                        *cutrhs -= 1;
                }
                continue;
            default:
                /* In any other case: skip */
                continue;
        }
    }


    int cnt_s = 0;
    std::vector<SCIP_Row*> vec_rows;
    /* check which slack variables join the cut */
    for(int c = 0; c < nrows; ++c)
    {
        row = rows[c];
        assert(row != nullptr);
//        assert(SCIPisEQ(scip, 0, SCIProwGetConstant(row)));
        /* If the slack variable is fixed, we can ignore this cut coefficient. */
        if( SCIPisFeasZero(scip, SCIProwGetRhs(row) - SCIProwGetLhs(row)) )
            continue;
        /* Get simplex tableau element. */ //TODO: RANGED CONSTRAINTS?
        switch ( SCIProwGetBasisStatus(row) )
        {
            case SCIP_BASESTAT_LOWER:
                /* check if coefficient has correct sign */
                if(SCIPisPositive(scip, -binvrow[SCIProwGetLPPos(row)]))
                {
                    /* Take element if nonbasic at lower bound.
                    * Then Slack Variable is at upper bound (\leq 0) -> flip it */

                    /* we might exclude the variable due to the gap */
                    /* would not have negative reduced cost with the changed objective function */
                    if(SCIPisGE(scip, factor * binvrow[SCIProwGetLPPos(row)], -SCIProwGetDualsol(row)))
                    {
                        continue;
                    }
                    if (cnt_s >= 0)
                        return false;
                    /* we do not actually add slack variable to cut but its row */
                    *cutrhs += (SCIProwGetLhs(row) - SCIProwGetConstant(row));
                    SCIP_COL** rowcols = SCIProwGetCols(row);
                    SCIP_Real* rowvals = SCIProwGetVals(row);
                    cnt_s++;
		            for(int i = 0; i < SCIProwGetNLPNonz(row); ++i )
                    {
                        if(SCIPisLT(scip, SCIPcolGetLb(rowcols[i]), SCIPcolGetUb(rowcols[i])))
                        {
                            workcoefs[SCIPcolGetLPPos(rowcols[i])] += rowvals[i];
                        }else// if(SCIPisEQ(scip, SCIPcolGetLb(rowcols[i]), 1))
                        {/* if fixed on 1, we have to adjust the rhs */
                            *cutrhs -= rowvals[i] * SCIPcolGetPrimsol(rowcols[i]);
                        }
                    }
                }
                continue;
            case SCIP_BASESTAT_UPPER:
                if(SCIPisPositive(scip, binvrow[SCIProwGetLPPos(row)]))
                {
                    /* we might exclude the variable due to the gap */
                    /* would not have negative reduced cost with the changed objective function */
                    if(SCIPisLE(scip, factor * binvrow[SCIProwGetLPPos(row)], -SCIProwGetDualsol(row)))
                    {
                        continue;
                    }

                    if (cnt_s >= 0)
                        return false;

                    /* we do not actually add slack variable to cut but its row */
                    *cutrhs -= (SCIProwGetRhs(row) - SCIProwGetConstant(row));
                    SCIP_COL** rowcols = SCIProwGetCols(row);
                    SCIP_Real* rowvals = SCIProwGetVals(row);
                    cnt_s++;
                    for(int i = 0; i < SCIProwGetNLPNonz(row); ++i )
                    {
                        if(SCIPisLT(scip, SCIPcolGetLb(rowcols[i]), SCIPcolGetUb(rowcols[i])))
                        {
                            workcoefs[SCIPcolGetLPPos(rowcols[i])] -= rowvals[i];
                        }else// if(SCIPisEQ(scip, SCIPcolGetLb(rowcols[i]), 1))
                        {/* if fixed on 1, we have to adjust the rhs */
                            *cutrhs += rowvals[i] * SCIPcolGetPrimsol(rowcols[i]);
                        }
                    }
                }
                continue;
            default:
                /* Basic variable: skip */
                continue;
        }
    }
    std::cout << "nSlack: " << cnt_s << " orig_var: " << cnt << std::endl;

    auto* probData = dynamic_cast<vrp::ProbDataVRP*>(SCIPgetObjProbData(scip));
    auto* pricerData = dynamic_cast<ObjPricerVRP*>(SCIPfindObjPricer(scip, "VRP_Pricer"));


    std::cout << "FACTOR: " << factor << std::endl;

    for(int i = 0; i < probData->nCons_; i++)
    {
        int pos = SCIProwGetLPPos(SCIPconsGetRow(scip, probData->cons_[i]));
        pricerData->dualValues_[i] += (factor * binvrow[pos]);
    }
    pricerData->setArcPrices(probData->getData(), false);
    setCurrentNeighborhood(pricerData, probData->getData(), false);
    for(int i = 0; i < probData->getData()->nDays; i++)
    {
        vector<tourVRP> bestTours;
        SCIP_CALL(generateLabelsBiDir(scip, probData->getData(), pricerData, bestTours, false, false, false, false, i, true));
        std::cout << "FOUND " << bestTours.size() << " tours on day " << i << std::endl;
//        std::cout << bestTours[0];
    }

    assert(false);

    /* set cut coefficients */
    *cutnz = 0;
    for(int c = 0; c < ncols; ++c )
    {
        col = cols[c];
        assert(col != nullptr);
        auto i = SCIPcolGetLPPos(col);
        assert(0 <= i);

        if(!SCIPisZero(scip, workcoefs[i]))
        {
            cutcoefs[*cutnz] = workcoefs[i];
            cutind[*cutnz] = c;
            (*cutnz)++;
        }
    }
    std::cout << "TOTAL nVars: " << *cutnz << std::endl;


    if(objcpc->max_supp_abs_ + objcpc->max_supp_rel_ * ncols < *cutnz)
        return false;

    return true;
}


static
SCIP_RETCODE addCPCToModel(
        SCIP*           scip,
        SCIP_Conshdlr*  conshdlr,
        SCIP_Result*    result,
        SCIP_Col**      cols,
        SCIP_Real       cutlhs,
        CutDir          cutdir,
        bool            cutislocal,
        SCIP_Real*      cutcoefs,
        const int*      cutind,
        int             cutnz,
        int*            ncuts,
        int             c,
        double*         effica
){
    char cutname[SCIP_MAXSTRLEN];
    SCIP_ROW* cut;
    int j;

    /* construct cut name */
    if( c >= 0)
        (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "cpc%d_x%d_d%d", SCIPgetNLPs(scip), c, cutdir);
    else
        (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "cpc%d_s%d_d%d", SCIPgetNLPs(scip), -c-1, cutdir);

    /* create empty cut */
    SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &cut, conshdlr, cutname, cutlhs, SCIPinfinity(scip), cutislocal, FALSE, FALSE) );

    /* cache the row extension and only flush them if the cut gets added */
    SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

    /* collect all non-zero coefficients */
    for( j = 0; j < cutnz; ++j )
    {
        SCIP_CALL( SCIPaddVarToRow(scip, cut, SCIPcolGetVar(cols[cutind[j]]), cutcoefs[j]) );
    }

    if( SCIProwGetNNonz(cut) == 0 )
    {
        assert(SCIPisFeasNegative(scip, cutlhs));
        SCIPdebugMsg(scip, " -> GMI cut detected infeasibility with cut 0 <= %f.\n", cutlhs);
        *result = SCIP_CUTOFF;
//                    break;
    }
    *effica = SCIPgetCutEfficacy(scip, nullptr, cut);

    /* Only take efficacious cuts, except for cuts with one non-zero coefficient (= bound
     * changes); the latter cuts will be handeled internally in sepastore. */
    if( SCIProwGetNNonz(cut) == 1 || SCIPisCutEfficacious(scip, nullptr, cut) )
    {
        SCIP_Bool infeasible;

        SCIPdebugMsg(scip, " -> found GMI cut <%s>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f).\n",
                     cutname, SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                     SCIPgetCutEfficacy(scip, nullptr, cut),
                     SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                     SCIPgetRowMaxCoef(scip, cut)/SCIPgetRowMinCoef(scip, cut));

        /* flush all changes before adding the cut */
        SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

        SCIP_CALL( SCIPaddRow(scip, cut, FALSE, &infeasible) );

        /* add global cuts that are not implicit bound changes to the cut pool */
        if( ! cutislocal && SCIProwGetNNonz(cut) > 1 )
        {
            SCIP_CALL( SCIPaddPoolCut(scip, cut) );
        }

        if ( infeasible )
            *result = SCIP_CUTOFF;
        else
            *result = SCIP_SEPARATED;
        (*ncuts)++;
    }

    /* release the row */
    SCIP_CALL( SCIPreleaseRow(scip, &cut) );

    return SCIP_OKAY;
}

/**@name Callback methods
 *
 * @{
 */
/** frees specific constraint data */

/** transforms constraint data into data belonging to the transformed problem */
SCIP_DECL_CONSTRANS(ConshdlrCPC::scip_trans)
{   /*lint --e{715}*/
    SCIP_CONSDATA* sourcedata;
    SCIP_CONSDATA* targetdata;

    assert(conshdlr != nullptr);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), "CPC") == 0);
    assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
    assert(sourcecons != nullptr);
    assert(targetcons != nullptr);

    sourcedata = SCIPconsGetData(sourcecons);
    assert(sourcedata != nullptr);

//    /* create constraint data for target constraint */
//    SCIP_CALL( consdataCreate(scip, &targetdata,
//                              sourcedata->customer, sourcedata->day, sourcedata->type, sourcedata->node) );

    /* create target constraint */
    SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
                              SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
                              SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
                              SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
                              SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

    return SCIP_OKAY;
}


SCIP_DECL_CONSSEPALP(ConshdlrCPC::scip_sepalp){
    SCIP_VAR** vars;
    SCIP_COL** cols;
    SCIP_ROW** rows;
    SCIP_Real* binvrow;
    SCIP_Real* cutcoefs;
    SCIP_Real* workcoefs;
    SCIP_Real cutlhs;
    SCIP_Real absgap;
    SCIP_Real relgap;
    int* cutind;
    int* basisind;
    int nvars;
    int ncols;
    int nrows;
    int maxsepacuts;
    int ncuts;
    int cutnz;
    int i, c;
    return SCIP_OKAY;

    double rel_gap = fabs(SCIPgetLPObjval(scip) - SCIPgetCutoffbound(scip))/ SCIPgetLPObjval(scip);
    if(rel_gap > 0.1)
        return SCIP_OKAY;
    *result = SCIP_DIDNOTRUN;
    if(SCIPgetGap(scip) > rel_gap_high_)
        return SCIP_OKAY;

    /* Only call separator, if an optimal LP solution is at hand. */
    if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
        return SCIP_OKAY;

    if(SCIPisInfinity(scip, SCIPgetLocalDualbound(scip)))
        return SCIP_OKAY;

    /* Only call separator, if the LP solution is basic. */
    if( ! SCIPisLPSolBasic(scip) )
        return SCIP_OKAY;

    /* Only call separator, if there are fractional variables. */
    if( SCIPgetNLPBranchCands(scip) == 0 )
        return SCIP_OKAY;

//    if(SCIPgetDepth(scip) > 0)
//        return SCIP_OKAY;
    std::cout << "++++++++++++++++++++++++ ITERATION " << ncallsnode_ << " ++++++++++++++++++++++++" << std::endl;

//    if(ncallsnode_ > 2)
//        assert(false);
    if(SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == lastnode_)
    {
        ncallsnode_++;
        if((SCIPgetDepth(scip) == 0 && maxroundsroot_ >= 0 && ncallsnode_ > maxroundsroot_) ||
           (SCIPgetDepth(scip) > 0 && maxrounds_ >= 0 && ncallsnode_ > maxrounds_))
            return SCIP_OKAY;
    }else{
        lastnode_ = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
        ncallsnode_ = 1;
    }

    /* get variables data */
    SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, nullptr, nullptr, nullptr, nullptr) );

    /* get LP data */
    SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
    SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

    /* exit if LP is trivial */
    if( ncols == 0 || nrows == 0 )
        return SCIP_OKAY;

//    assert(!SCIPisGE(scip, SCIPgetLocalDualbound(scip), primalb_))
//        return SCIP_OKAY;

//    assert(SCIPgetPrimalbound(scip) < primalb_);
    assert(SCIPisLE(scip, SCIPgetLocalDualbound(scip), SCIPgetPrimalbound(scip)));
    /* check local absolute gap */
    if(SCIPgetPrimalbound(scip) - SCIPgetLocalDualbound(scip) < SCIPinfinity(scip))
    {
//        absgap = SCIPgetPrimalbound(scip) - SCIPgetLocalDualbound(scip);
        absgap = SCIPgetPrimalbound(scip) - SCIPgetLocalDualbound(scip);
        if(SCIPisZero(scip, SCIPgetLocalDualbound(scip)))
            relgap = -1;
        else
            relgap = absgap / abs(SCIPgetLocalDualbound(scip));

        assert(SCIPgetPrimalbound(scip) < SCIPinfinity(scip));
        assert(0 < absgap < SCIPinfinity(scip));
    }else
    {
        return SCIP_OKAY;
    }

//    dualvals_.push_back(SCIPgetDualbound(scip));

    /* we do not expect good cuts for high gaps */
    int max_nofind = max_noFind_high_;
    int max_ntries = max_nTries_high_;

    if(relgap > rel_gap_low_)
    {
        max_nofind = max_noFind_low_;
        max_ntries = max_nTries_low_;
    }

    *result = SCIP_DIDNOTFIND;

    double* binvarow;
    /* allocate temporary memory */
    SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, ncols) );
    SCIP_CALL( SCIPallocBufferArray(scip, &workcoefs, ncols) );
    SCIP_CALL( SCIPallocBufferArray(scip, &cutind, ncols) );
    SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );
    SCIP_CALL( SCIPallocBufferArray(scip, &binvarow, ncols) );
    SCIP_CALL( SCIPallocBufferArray(scip, &binvrow, nrows) );

    /* get basis indices */
    SCIP_CALL( SCIPgetLPBasisInd(scip, basisind) );

    /* get the maximal number of cuts allowed in a separation round */
    if(SCIPgetDepth(scip) == 0 )
        maxsepacuts = maxsepacutsroot_;
    else
        maxsepacuts = maxsepacuts_;

    if( maxsepacuts == -1 )
        maxsepacuts = INT_MAX;

    ncuts = 0;
    double lowestsupp = SCIPinfinity(scip);
    int ntries = 0;

    double effica = 0;
    double maxeffica = 0;
    nRuns_++;

//    std::cout << "------------- Start Propa at gap " << relgap*100 << " abs: " << absgap << " ---------------" << std::endl;

    /* For all basic columns belonging to integer variables, try to generate a CPC cut. */
    for( i = 0; i < nrows && ncuts < maxsepacuts && ! SCIPisStopped(scip) && *result != SCIP_CUTOFF; ++i )
    {
        if(ncuts == 0 && ntries > max_nofind)
        {
            break;
        }
        if(ntries > max_ntries)
            break;
        SCIP_Bool tryrow;
        SCIP_Real primsol;

        tryrow = FALSE;
        c = basisind[i];
        primsol = SCIP_INVALID;
        if( c >= 0 )
        {
            if(!SCIPvarIsBinary(SCIPcolGetVar(cols[c])))
            {
                continue;
            }
            ntries++;
            assert(c < ncols);
            assert(cols[c] != nullptr);

            primsol = SCIPcolGetPrimsol(cols[c]);

            if( (SCIPfeasFrac(scip, primsol) >= away_) && (SCIPfeasFrac(scip, primsol) <= 1.0 - away_) )
            {
                tryrow = TRUE;
            }
        }
        if (tryrow)
        {
            SCIP_Real factor;
            bool cutislocal;
            bool success = false;

            /* get the row of B^-1 for this basic integer variable with fractional solution value */
            SCIP_CALL( SCIPgetLPBInvRow(scip, i, binvrow, nullptr, nullptr) );

            cutislocal = (SCIPgetDepth(scip) != 0);

            /* aggregate the simplex tableau row */
            SCIP_CALL (SCIPgetLPBInvARow(scip, i, binvrow, binvarow, nullptr, nullptr));

            /* create a CPC out of the simplex tableau row */
            /* downwards */
            /* for constraint type DOWNWARDS the LP-value must not be too low */
            if(primsol >= 0.5)
            {
                auto* varData = dynamic_cast<ObjVarDataVRP*>(SCIPgetObjVardata(scip, SCIPcolGetVar(cols[c])));
//                std::cout << varData->tourVrp_;
//                for(int o = 0; o < nrows; o++)
//                {
//                    std::cout << binvrow[o] << " (" << o << ") ";
//                }
//                std::cout << std::endl;
                factor = absgap / primsol;
                success = getCPCFromRowAggr(scip, this, ncols, nrows, cols, rows, binvarow, c, binvrow, cutcoefs, cutind,
                                        &cutnz, &cutlhs, workcoefs, factor);
                if(success)
                {
                    std::cout << "-------------------------------------------------------" << std::endl;
                    std::cout << "Success for variable with LP-val: " << primsol << std::endl;
                    if(cutnz < lowestsupp)
                        lowestsupp = cutnz;
                    SCIP_CALL(addCPCToModel(scip, conshdlr, result, cols, cutlhs, DOWNWARDS, cutislocal, cutcoefs,
                                            cutind, cutnz, &ncuts, c, &effica));
                    if(effica > maxeffica)
                        maxeffica = effica;
                    if(*result == SCIP_SEPARATED)
                        ncuts_++;
//                    std::cout << " -> down - Effica: " << effica << " soldist: " << 1 - primsol;
//                    std::cout << " nnz: " << cutnz << " factor: " << factor << std::endl;
                }

//                success = getGOMIRFromRow(scip, this, ncols, nrows, cols, rows, c, binvrow, primsol, cutcoefs, cutind,
//                                                &cutnz, &cutlhs, &cutact, workcoefs, bigR, smallG, aggrrow, &mir_effica);
            }

        }
    }


    /* free temporary memory */
    SCIPfreeBufferArray(scip, &binvrow);
    SCIPfreeBufferArray(scip, &basisind);
    SCIPfreeBufferArray(scip, &workcoefs);
    SCIPfreeBufferArray(scip, &cutcoefs);
    SCIPfreeBufferArray(scip, &cutind);
//    assert(false);
    return SCIP_OKAY;

}

SCIP_DECL_CONSPRINT(ConshdlrCPC::scip_print)
{
    return SCIP_OKAY;
}


ConshdlrCPC::ConshdlrCPC(SCIP *scip)
: ObjConshdlr(scip, "CPC", "counterpart cuts conshdlr", 50000000, 400, 0,
                                                      1, 0, 0, 0, false, FALSE, FALSE,
                                                      SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_ALWAYS)
{
    checkedIfInt_ = false;
    modelIsInt_ = false;
    lastnode_ = -1;
    ncallsnode_ = 0;
    skipped_slack_ = 0;
    checked_slack_ = 0;
    nRuns_ = 0;
    ncuts_ = 0;
    start_gap_ = 0.0;

    maxrounds_ = DEFAULT_MAXROUNDS;
    maxroundsroot_ = DEFAULT_MAXROUNDSROOT;
    maxsepacuts_ = DEFAULT_MAXSEPACUTS;
    maxsepacutsroot_ = DEFAULT_MAXSEPACUTSROOT;
    away_ = DEFAULT_AWAY;
    maxdynamism_ = DEFAULT_MAX_DYN;
    max_supp_abs_ = DEFAULT_MAX_SUPP_ABS;
    max_supp_rel_ = DEFAULT_MAX_SUPP_REL;
}

