
#ifndef VRP_CONSHDLRRCFC_H
#define VRP_CONSHDLRRCFC_H



#include "objscip/objscip.h"
#include "model_data.h"
#include "probdata_vrp.h"
#include "scip/scip.h"
#include "vardata.h"

using namespace scip;

class ConshdlrRCFC : public ObjConshdlr
{
public:
    bool stopfixing_;
    /** default constructor */
    explicit ConshdlrRCFC(
            SCIP*   scip,
            double  factor,
            double  init_fail,
            bool    withCuts
    )
            : ObjConshdlr(scip, "RCFC", "counterpart cuts conshdlr", 0, 5000, 0,
                          0, 0, 1, 0, TRUE, FALSE, FALSE,
                          SCIP_PROPTIMING_BEFORELP, SCIP_PRESOLTIMING_ALWAYS),
                          stopfixing_(false)
    {
        LP_iters_old_ = 0;
        LP_iters_exact_ = 0;
        nFixedExact_ = 0;
        nFixedHeurExtv2_ = 0;
        time_exact_ = 0;
        iter_init_ = 10;
        iter_perc_ = 0.5;
        validGap_ = 0.1;
        obj_cmp_ = -1;
        lowest_fail_ = 10000;
        fixVar_ = nullptr;
        factor_ = factor;
        steepness_skip_ = 0;

        init_fail_factor_ = init_fail;
        lastNode_ = -1;
        nCalls_Node_ = 0;
        max_nCalls_ = 1;
        max_nCallsRoot_ = 5;
        nonViolatedCut_ = false;
        noChange_ = false;

        nCuts_ = 0;

        nCheckVars_ = 0;
        pricing_time_ = 0;
        lp_solving_time_ = 0;
        lpIter_ = 0;
        n_try_cutoff_ = 0;
        withCuts_ = withCuts;
    }

    /** destructor */
    ~ConshdlrRCFC() override= default;

    /** frees specific constraint data */
    virtual SCIP_DECL_CONSDELETE(scip_delete);

    /** transforms constraint data into data belonging to the transformed problem */
    virtual SCIP_DECL_CONSTRANS(scip_trans){
        return SCIP_OKAY;
    }

    /** separation method of constraint handler for LP solution */
    virtual SCIP_DECL_CONSSEPALP(scip_sepalp){
        return SCIP_OKAY;
    }

    /** constraint enforcing method of constraint handler for LP solutions */
    virtual SCIP_DECL_CONSENFOLP(scip_enfolp);

    /** constraint enforcing method of constraint handler for pseudo solutions */
    virtual SCIP_DECL_CONSENFOPS(scip_enfops){
        return SCIP_OKAY;
    }

    /** feasibility check method of constraint handler for integral solutions */
    virtual SCIP_DECL_CONSCHECK(scip_check){
        return SCIP_OKAY;
    }

    /** variable rounding lock method of constraint handler */
    virtual SCIP_DECL_CONSLOCK(scip_lock){
        return SCIP_OKAY;
    }

    /** constraint display method of constraint handler */
    virtual SCIP_DECL_CONSPRINT(scip_print);


    long long int LP_iters_old_;
    long long int LP_iters_exact_;
    int nFixedExact_;
    long time_exact_;
    int iter_init_;
    double iter_perc_;
    double validGap_;
    int nFixedHeurExtv2_;
    double root_gap_ = 100;
    long long int time_success_ = 0;
    long long int time_fail_ = 0;
    double minvalue_ = 0.1;

    int depth_init_ = 4;
    int onlyFracGap_ = 5;
    bool onlyFrac_ = false;

    double init_fail_factor_;
    int steepness_skip_;
    double obj_cmp_;
    double lowest_fail_;
    double factor_;
    SCIP_Var* fixVar_;

    int max_nCalls_;
    int max_nCallsRoot_;
    int nCalls_Node_;
    long long int lastNode_;

    bool used_Pricing_;
    int nAdded_Vars_;
    bool nonViolatedCut_;
    bool noChange_;
    bool withCuts_;

    int nCuts_;
    int nFixed_ = 0;
    int nFixed_frac_ = 0;
    int nCheckVars_;

    std::vector<int> cutSize_;

    std::chrono::time_point<std::chrono::high_resolution_clock> start_time_;
    long long int pricing_time_;
    long long int lp_solving_time_;
    long long int lpIter_;
    vector<double> fail_factors_;
    int n_try_cutoff_;

    std::vector<int> depths_;
    std::vector<int> nVars_;
    std::vector<int> success_;
    std::vector<double> gaps_;
};

/** creates and adds row for a subset row constraint */
SCIP_RETCODE SCIPcreateAndAddRowCPC(
        SCIP*                   scip,
        SCIP_Cons*              cons
);

/** creates and captures a subset row constraint based on Set S */
SCIP_RETCODE SCIPcreateConsCPC(
        SCIP*                   scip,               /**< SCIP data structure */
        SCIP_CONS**             cons,               /**< pointer to hold the created constraint */
        const char*             name,               /**< name of constraint */
        vector<SCIP_VAR*>&      vars
);



#endif //VRP_CONSHDLRRCFC_H
