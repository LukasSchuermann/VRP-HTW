
#ifndef __SCIP_VRPPROBDATA_H__
#define __SCIP_VRPPROBDATA_H__


#include <cassert>
#include <utility>

#include "stdio.h"
#include "model_data.h"
#include "vector"
#include "scip/scip.h"
#include "objscip/objscip.h"
#include "tourVRP.h"

using namespace std;
namespace vrp
{

/** SCIP user problem data for VRP */
class ProbDataVRP : public scip::ObjProbData
{
    private:
        model_data*     modelData_;
    public:
        int                     nVars_;
        int                     nCons_;
        vector<SCIP_VAR*>       vars_;
        vector<SCIP_CONS*>      cons_;
        vector<SCIP_VAR*>       emptyVars_;
        vector<SCIP_VAR*>       nonzeroVars_;
        int                     use_propagator_;
        double                  fw_time_;
        double                  bd_time_;
        vector<vector<SCIP_Real>> prices;
        std::vector<std::vector<double>> depth_gaps_;

        /** default constructor */
        ProbDataVRP(
            model_data*             modelData,
            int                     nVars,
            int                     nCons,
            int                     activate_propagator
        ):
        modelData_(modelData),
        nVars_(nVars),
        nCons_(nCons),
        use_propagator_(activate_propagator)
        {
            cons_ = std::vector<SCIP_CONS*>(nCons);
            vars_ = std::vector<SCIP_VAR*>(nVars);
            fw_time_ = 0;
            bd_time_ = 0;
        }

        /** destructor */
        ~ProbDataVRP() override
        = default;

        /** destructor of user problem data to free original user data (called when original problem is freed) */
        SCIP_RETCODE scip_delorig(
                SCIP*              scip                /**< SCIP data structure */
        ) override;

        /** destructor of user problem data to free transformed user data (called when transformed problem is freed) */
        SCIP_RETCODE scip_deltrans(
                SCIP*              scip                /**< SCIP data structure */
        ) override;

        /** creates user data of transformed problem by transforming the original user problem data
         *  (called after problem was transformed)
         */
        SCIP_RETCODE scip_trans(
                SCIP*              scip,               /**< SCIP data structure */
                ObjProbData**      objprobdata,        /**< pointer to store the transformed problem data object */
                SCIP_Bool*         deleteobject        /**< pointer to store whether SCIP should delete the object after solving */
        ) override;

        /** solving process initialization method of transformed data (called before the branch and bound process begins) */
        SCIP_RETCODE scip_initsol(
                SCIP*              scip                /**< SCIP data structure */
        ) override;

        /** solving process deinitialization method of transformed data (called before the branch and bound data is freed) */
        SCIP_RETCODE scip_exitsol(
                SCIP*              scip,                /**< SCIP data structure */
                SCIP_Bool          restart              /**< was this exit solve call triggered by a restart? */
        ) override;

        model_data* getData()
        {
            return modelData_;
        }
    };

} /* namespace vrp */

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(
        SCIP*                 scip,               /**< SCIP data structure */
        model_data*           modeldata,           /**< model data */
        std::vector<tourVRP>& sol_tvrps,
        int                   activate_propagator
);

/** adds given variable to the problem data */
SCIP_RETCODE SCIPprobdataAddVar(
        SCIP*                   scip,                   /**< SCIP data structure */
        vrp::ProbDataVRP*       objprobdata,            /**< problem data */
        SCIP_VAR*               var                     /**< variables to add */
);

#endif
