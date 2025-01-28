//
// Created by lukasschuermann on 27.09.23.
//

#ifndef VRP_ConshdlrCPC_H
#define VRP_ConshdlrCPC_H

#include "objscip/objscip.h"
#include "scip/scip.h"
#include "vector"

class ConshdlrCPC : public scip::ObjConshdlr
{
public:
    std::vector<double> dualvals_;
    long skipped_slack_;
    long checked_slack_;
    int nRuns_;
    double start_gap_;
    int ncuts_;
    int             max_supp_abs_;         /**< maximum cut support - absolute value in the formula */
    SCIP_Real       max_supp_rel_;         /**< maximum cut support - relative value in the formula */
    int max_noFind_high_;
    int max_noFind_low_;
    int max_nTries_high_;
    int max_nTries_low_;
    double rel_gap_high_;
    double rel_gap_low_;
    /** default constructor */
    explicit ConshdlrCPC(
            SCIP*       scip
    );

    /** destructor */
    ~ConshdlrCPC() override= default;

    /** frees specific constraint data */
    virtual SCIP_DECL_CONSDELETE(scip_delete){
        return SCIP_OKAY;
    };

    /** transforms constraint data into data belonging to the transformed problem */
    virtual SCIP_DECL_CONSTRANS(scip_trans);

    /** separation method of constraint handler for LP solution */
    virtual SCIP_DECL_CONSSEPALP(scip_sepalp);

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

    /** variable rounding lock method of constraint handler */
    virtual SCIP_DECL_CONSLOCK(scip_lock){
            return SCIP_OKAY;
    }

    virtual SCIP_DECL_CONSCOPY(scip_copy){
        return SCIP_OKAY;
    }

    /** constraint display method of constraint handler */
    virtual SCIP_DECL_CONSPRINT(scip_print);


    SCIP_Real getMaxDynamism() const{return maxdynamism_;}
    SCIP_Real getAway() const{return away_;}
    bool isint() const{return modelIsInt_;}

private:
    bool            checkedIfInt_;
    bool            modelIsInt_;
    long long int   lastnode_;
    int             ncallsnode_;
    int             maxrounds_;          /**< maximal number of GMI separation rounds per node (-1: unlimited) */
    int             maxroundsroot_;      /**< maximal number of GMI separation rounds in the root node (-1: unlimited) */
    int             maxsepacuts_;        /**< maximal number of GMI cuts separated per separation round */
    int             maxsepacutsroot_;    /**< maximal number of GMI cuts separated per separation round in root node */
    SCIP_Real       away_;               /**< minimal fractionality of a basis variable in order to try GMI cut */
    SCIP_Real       maxdynamism_;        /**< maximal valid range max(|weights|)/min(|weights|) of cut coefficients */
};


#endif //VRP_ConshdlrCPC_H
