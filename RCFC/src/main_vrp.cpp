
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"
#include "objscip/objscip.h"

#include "model_data.h"
#include "probdata_vrp.h"
#include "branchingrule_arcflow.h"
#include "branchingrule_dayvar.h"
#include "branchingrule_vehicle.h"
#include "ConshdlrArcflow.h"
#include "ConshdlrDayVar.h"
#include "ConshdlrKPC.h"
#include "ConshdlrSRC.h"
#include "ConshdlrRCFC.h"
#include "ConshdlrCPC.h"
#include "ConshdlrVehicle.h"
#include "ConshdlrNVehicle.h"
#include "tourVRP.h"
#include "tools_vrp.h"
#include "eventhdlr_nodeInit.hpp"
#include "printer.h"
#include "heurDayVarRounding.h"
#include "prop_varfixing.h"
#include "prop_tourvarfixing.h"


#include "iostream"
#include "fstream"
#include "json.hpp"

using json = nlohmann::json;

static
SCIP_RETCODE readArguments(
    int                 argc,               /**< number of shell parameters */
    char**              argv,               /**< array with shell parameters */
    std::string&        input_file,
    int*                dayVarBranching,
    int*                activate_propagator,
    int*                seed,
    double*             fail_factor,
    double*             init_factor,
    bool*               withCuts,
    std::string&             sol_file,
    std::string&             out_file,
    bool*               withSol,
    bool*               useSRC,
    int*                nMaxSRC,
    bool*               useKPC,
    double*             gap_decay,
    int*                max_depth,
    int*                maxSRCrcf,
    bool*               useArcRC,
    bool*               useVeAssRC,
    bool*               noRootFixing,
    double*             branchingFrac,
    double*             minvalue,
    int*                depth_init,
    int*                onlyFracGap,
    bool*               onlyFrac
)
{
    char usage[SCIP_MAXSTRLEN];
    int status;
    char* locstr;

    assert( argc >= 1 );
    assert( argv != nullptr );
    assert( input_file.empty() );

    /* init usage text */
    status = snprintf(usage, SCIP_MAXSTRLEN - 1, "usage: %s <path of inputfile> ", argv[0]);
    assert( 0 <= status && status < SCIP_MAXSTRLEN );


    /* mandatory argument: inputfile */
    input_file = argv[1];
    if ( input_file.empty() )
    {
        fprintf(stderr, "No path of data supplied.\n");
        fprintf(stderr, "%s\n", usage);
        return SCIP_ERROR;
    }
    /* check for branching strategy */
    for (int i = 2; i < argc; i++)
    {
        if ( ! strcmp(argv[i], "-b"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing branching choice. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *dayVarBranching = atoi(locstr);
            if(*dayVarBranching != 0 && *dayVarBranching != 1)
            {
                fprintf(stderr, "Invalid branching strategy -> Choose 1 for day var and 0 for arc flow. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            free(locstr);
        }else if(! strcmp(argv[i], "-p"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing propagation info. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *activate_propagator = atoi(locstr);
            if(*activate_propagator < 0 || *activate_propagator > 3)
            {
                fprintf(stderr, "Invalid prop! -> Choose propagation setting in [0 (noprop), 1 (it-RCFC), 2 (RCFC), "
                                "3 (strong branching)]!");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            free(locstr);
            if(*activate_propagator == 1)
                std::cout << "USE it-RCFC PROPAGATOR!" << std::endl;
            else if(*activate_propagator == 2)
                std::cout << "USE RCFC PROPAGATOR!" << std::endl;
            else if(*activate_propagator == 3)
                std::cout << "USE STRONG BRANCHING PROPAGATOR!" << std::endl;
        }else if ( ! strcmp(argv[i], "-s"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing seed number. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *seed = atoi(locstr);
            if(*seed < 0)
            {
                fprintf(stderr, "Invalid seed! -> Choose seed >= 0!");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            free(locstr);
        }else if ( ! strcmp(argv[i], "-f"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing fail factor. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *fail_factor = atof(locstr);
            if(*fail_factor <= 0)
            {
                fprintf(stderr, "Invalid factor! -> Choose fail factor > 0 !");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }else
            {
                std::cout << "Chosen fail factor of " << *fail_factor << std::endl;
            }
            free(locstr);
        }else if ( ! strcmp(argv[i], "-i"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing init factor. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *init_factor = atof(locstr);
            if(*init_factor <= 0)
            {
                fprintf(stderr, "Invalid init factor! -> Choose init factor > 0 !");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }else
            {
                std::cout << "Chosen fail factor of " << *init_factor << std::endl;
            }
            free(locstr);
        }else if ( ! strcmp(argv[i], "-c"))
        {
            *withCuts = true;
        }else if ( ! strcmp(argv[i], "-sol"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing solution file. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            *withSol = true;
            sol_file = argv[i];
        }else if(! strcmp(argv[i], "-out"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing output file. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            out_file = argv[i];
        }else if ( ! strcmp(argv[i], "-src"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing max SRC number. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            *nMaxSRC = atoi(argv[i]);
            if(*nMaxSRC == 0){
                *useSRC = false;
            }
        }else if ( ! strcmp(argv[i], "-nkpc"))
        {
            *useKPC = false;
        }else if ( ! strcmp(argv[i], "-bf"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing branching fractionality. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            *branchingFrac = atof(argv[i]);
        }else if ( ! strcmp(argv[i], "-rcf"))
        {
            if( i + 3 > argc - 1 || (! strncmp(argv[i+1], "-",1)) ||
                (! strncmp(argv[i+2], "-",1)) || (! strncmp(argv[i+3], "-",1)))
            {
                printf("i %d, argc %d\n",i,argc);
                fprintf(stderr, "reduced cost fixing parameters. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            *gap_decay = std::atof(argv[i++]);
            *max_depth = std::atoi(argv[i++]);
            *maxSRCrcf = std::atoi(argv[i]);
        }else if ( ! strcmp(argv[i], "-nArcRC"))
        {
            *useArcRC = false;
        }else if ( ! strcmp(argv[i], "-nVeAssRC")){
            *useVeAssRC = false;
        }else if ( ! strcmp(argv[i], "-nR")){
            *noRootFixing = true;
        }else if ( ! strcmp(argv[i], "-mv"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Min value for rcfc missing. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            *minvalue = atof(argv[i]);
        }else if ( ! strcmp(argv[i], "-di"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Depth init for rcfc missing. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            *depth_init = atoi(argv[i]);
        }else if ( ! strcmp(argv[i], "-of"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "only frac activation value missing. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            *onlyFrac = true;
            *onlyFracGap = atoi(argv[i]);
        }
    }

    return SCIP_OKAY;
}

static
SCIP_RETCODE setUpScip(
        SCIP**              scip,
        int                 dayVarBranching,
        bool                useRCFC,
        bool                useCPC,
        double              fail_factor,
        double              init_fail,
        bool                withCuts,
        bool                withSol,
        bool                useKPC,
        bool                useSRC,
        int                 nMaxSRC,
        double              branchingFrac
) {
    /* initialize SCIP */
    SCIP_CALL(SCIPcreate(scip));

    /* include default SCIP plugins */
    SCIP_CALL( SCIPincludeDefaultPlugins(*scip) );


    /* include branching rules */
    SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleVehicle(*scip, 5000), TRUE));
    if(dayVarBranching)
    {
        SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleArcflow(*scip, 50000, branchingFrac), TRUE));
        SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleDayVar(*scip, 50001, branchingFrac), TRUE));
    }else{
        SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleArcflow(*scip, 50001, branchingFrac), TRUE));
        SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleDayVar(*scip, 50000, branchingFrac), TRUE));
    }

    /* turn off all separation algorithms */
    SCIP_CALL( SCIPsetSeparating(*scip, SCIP_PARAMSETTING_OFF, TRUE) );

    /* turn off primal heuristics if optimal solution is read */
    if(withSol){
        SCIP_CALL( SCIPsetHeuristics(*scip, SCIP_PARAMSETTING_OFF, TRUE) );
    }

    /* include robust cuts separator */
//    SCIP_CALL(SCIPincludeObjSepa(*scip, new ObjCvrpSep(*scip), TRUE));

    /* include propagator */
    SCIP_CALL(SCIPincludeObjProp(*scip, new ObjPropVarFixing(*scip, false), TRUE));
    SCIP_CALL(SCIPincludeObjProp(*scip, new ObjPropTourVarFixing(*scip), TRUE));

    /* include primal heuristics */
    SCIP_CALL(SCIPincludeObjHeur(*scip, new HeurDayVarRounding(*scip), TRUE));

    /* include event handler */
    SCIP_CALL(SCIPincludeObjEventhdlr(*scip, new EventhdlrNodeInit(*scip), TRUE));

    /* include constraint handler */
    SCIP_CALL( SCIPincludeObjConshdlr(*scip, new ConshdlrArcflow(*scip), TRUE) );
    SCIP_CALL( SCIPincludeObjConshdlr(*scip, new ConshdlrDayVar(*scip), TRUE) );
    SCIP_CALL( SCIPincludeObjConshdlr(*scip, new ConshdlrKPC(*scip, useKPC), TRUE) );
    SCIP_CALL( SCIPincludeObjConshdlr(*scip, new ConshdlrSRC(*scip, useSRC, nMaxSRC), TRUE) );
    SCIP_CALL( SCIPincludeObjConshdlr(*scip, new ConshdlrVehicle(*scip), TRUE) );
    SCIP_CALL( SCIPincludeObjConshdlr(*scip, new ConshdlrNVehicle(*scip), TRUE) );

    if(useCPC){
        SCIP_CALL( SCIPincludeObjConshdlr(*scip, new ConshdlrCPC(*scip), TRUE) );
    }
    if(useRCFC){
        SCIP_CALL( SCIPincludeObjConshdlr(*scip, new ConshdlrRCFC(*scip, fail_factor, init_fail, withCuts), TRUE) );
    }

    /** change parameters */

    /* Constraint handler */
    SCIP_CALL( SCIPsetIntParam(*scip,"constraints/Vehicle/maxrounds",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"constraints/NVehicle/maxrounds",2) );

    /* change display columns */
    SCIP_CALL( SCIPsetIntParam(*scip,"display/nfrac/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/curcols/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/cuts/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/lpobj/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/primalgap/active",0) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/gap/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/lpavgiterations/active",0) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/lpiterations/active",0) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/vars/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/conflicts/active",0) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/strongbranchs/active",0) );
    /* Branching */
    SCIP_CALL( SCIPsetIntParam(*scip,"display/nnodes/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/nodesleft/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip,"display/maxdepth/active",2) );
    SCIP_CALL( SCIPsetIntParam(*scip, "display/freq", 50));

    /* for column generation instances, disable restarts */
    SCIP_CALL( SCIPsetIntParam(*scip,"presolving/maxrestarts",0) );

    /* activate reduced costs propagator */
//    SCIP_CALL(SCIPsetBoolParam(*scip, "propagating/redcost/force", true));

    /* set the limit for the relative gap after which the solving process is stopped */
//    SCIP_CALL (SCIPsetRealParam(*scip, "limits/gap", 0.0));

    SCIP_CALL( SCIPsetIntParam(*scip,"display/verblevel", 4) );

    /* set a time limit */
    SCIP_CALL( SCIPsetRealParam(*scip, "limits/time", 3600));


    return SCIP_OKAY;
}

static
SCIP_RETCODE addInitTours(
    SCIP*           scip,
    vector<int>&    hash_day,
    SCIP*           scip_new
)
{
    char algoName[] = "initTour";
    auto* probData = dynamic_cast<vrp::ProbDataVRP*>(SCIPgetObjProbData(scip));
    for(auto var : probData->vars_)
    {
        if(SCIPgetSolVal(scip, SCIPgetBestSol(scip), var) < 0.5)
            continue;
        tourVRP tvrp;
        auto* vardata = dynamic_cast<ObjVarDataVRP*>(SCIPgetObjVardata(scip, var));
        tvrp.copy(vardata->tourVrp_);
        tvrp.setDay(hash_day[vardata->getDay()]);
        SCIP_CALL(add_tour_variable(scip_new, dynamic_cast<vrp::ProbDataVRP*>(SCIPgetObjProbData(scip_new)),
                                    FALSE, TRUE, algoName, tvrp));
    }
    return SCIP_OKAY;
}

static
SCIP_RETCODE readSolution(
    std::string&             sol_file,
    vector<tourVRP>&    tvrps
){
    int day, length, cap;
    double obj;
    int count = 0;
    std::string line;
    ifstream solstream;
    solstream.open(sol_file);
    if(solstream.is_open())
    {
        while (getline(solstream, line))
        {
            day = stoi(line.substr(0, line.find(' ')));
            line = line.substr(line.find(' ') + 1);
            obj = stod(line.substr(0, line.find(' ')));
            line = line.substr(line.find(' ') + 1);
            cap = stoi(line.substr(0, line.find(' ')));
            line = line.substr(line.find(' ') + 1);
            length = stoi(line.substr(0, line.find(' ')));
            line = line.substr(line.find(' ') + 1);
            tvrps.emplace_back(length, day);
            tvrps[count].obj_ = obj;
            tvrps[count].capacity_ = cap;
            for(int i = 0; i < length; i++)
            {
                tvrps[count].tour_[i] = stoi(line.substr(0, line.find(' ')));
                line = line.substr(line.find(' ') + 1);
            }

            count++;
        }
        solstream.close();
    }
    return SCIP_OKAY;
}

/** output solution */
static
SCIP_RETCODE outputSolution(
    SCIP*           scip,
    std::string&    path,
    std::string&    out_file,
    bool            vehAssBr,
    double          frac
)
{
    size_t lastSlashPos = path.find_last_of('/');
    std::string filename = path.substr(lastSlashPos + 1);

    size_t underscorePos = filename.find('_');
    std::string beforeUnderscore = filename.substr(0, underscorePos);

    size_t nPos = filename.find('n');
    size_t pPos = filename.find('p');

    int nValue = std::stoi(filename.substr(nPos + 1, pPos - nPos - 1));
    float pValue = std::stof(filename.substr(pPos + 1));


    std::ofstream sol_file;
    sol_file.open(out_file, std::ios_base::app);
    sol_file << "Base;N;P;SolTime;nNodes;Gap;propTime;pricingTime;successTime;failTime;nFixed;nFrac;nCuts;avgSize" << std::endl;
    sol_file << beforeUnderscore << ";" << nValue << ";" << pValue << ";" << SCIPgetSolvingTime(scip) << ";";
    sol_file << SCIPgetNNodes(scip) << ";" << SCIPgetGap(scip) * 100 << ";";

    auto* prop = dynamic_cast<ConshdlrRCFC*>(SCIPfindObjConshdlr(scip, "RCFC"));


    int sum = std::accumulate(prop->cutSize_.begin(), prop->cutSize_.end(), 0);
    // Return the mean
    auto avg = sum > 0 ? static_cast<double>(sum) / prop->cutSize_.size() : 0.0;


    if(prop != nullptr){
        sol_file << SCIPconshdlrGetEnfoLPTime(SCIPfindConshdlr(scip, "RCFC")) << ";" << prop->pricing_time_ << ";";
        sol_file << prop->time_success_ << ";" << prop->time_fail_ << ";" << prop->nFixed_ << ";" << prop->nFixed_frac_;
        sol_file << ";" << prop->nCuts_ << ";" << avg << std::endl;
    }else
        sol_file << std::endl;

    return SCIP_OKAY;
}


static
SCIP_RETCODE findInfeasiblityReason(
    SCIP*   scip
){
    auto *pricerData = dynamic_cast<ObjPricerVRP *>(SCIPfindObjPricer(scip, "VRP_Pricer"));
    auto *probData = dynamic_cast<vrp::ProbDataVRP *>(SCIPgetObjProbData(scip));
    model_data* modelData = probData->getData();
    std::cout << "Infeasibility might come from the following days / customers" << std::endl;
    for(int i = 0; i < modelData->nDays; i++)
    {
        if(!SCIPisZero(scip, pricerData->dualValues_[modelData->nC + i]))
            std::cout << "day " << i << std::endl;
    }
    for(int i = 0; i < modelData->nC; i++)
    {
        if(!SCIPisZero(scip, pricerData->dualValues_[i]))
        {
            assert(SCIPisEQ(scip, pricerData->dualValues_[i], SCIPgetDualfarkasLinear(scip, probData->cons_[i-1])));
            std::cout << "cust " << i << " available on days: ";
            for(auto day : modelData->availableDays[i])
                std::cout << day << " (" << modelData->timeWindows[i][day].start << "-"
                << modelData->timeWindows[i][day].end << ")  ";
            std::cout << std::endl;
        }
    }
    return SCIP_OKAY;
}

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runColumnGenerationModel(
    int                   argc,               /**< number of shell parameters */
    char**                argv                /**< array with shell parameters */
)
{
    SCIP* scip = nullptr;
    std::string input_file;
    auto* modelData = new model_data;
    int ng_parameter = 8;
    int seed = 0;
    double fail_factor = 10;
    double init_factor = 10;
    bool withCuts = false;
    bool useRCFC = true;
    bool useCPC = false;
    std::string sol_file;
    std::string out_file;
    bool withSol = false;
    /* rcfc */
    int activate_propagator = 0;
    double minvalue = 0.1;
    int depth_init = 4;
    int onlyFracGap = 5;
    bool onlyFrac = false;
    /* cutting settings */
    bool useSRC = true;
    int nMaxSRC = 200;
    bool useKPC = true;
    /* reduced costs fixing parameters */
    double gap_decay_rc = 0.8;
    int max_depth_rc = 0;
    int maxSRC_rc = 75;
    bool useArcRC = true;
    bool useVeAssRC = true;
    bool noRootFixing = false;
    /* branching settings */
    double branchingFrac = 0.5;
    int dayVarBranching = 1;

    SCIP_CALL(readArguments(argc, argv, input_file, &dayVarBranching, &activate_propagator, &seed, &fail_factor,
                            &init_factor, &withCuts, sol_file, out_file, &withSol, &useSRC, &nMaxSRC, &useKPC,
                            &gap_decay_rc, &max_depth_rc, &maxSRC_rc, &useArcRC, &useVeAssRC, &noRootFixing,
                            &branchingFrac, &minvalue, &depth_init, &onlyFracGap, &onlyFrac));
    useRCFC = (activate_propagator > 0);
    if(withSol){
        std::cout << "Include solution from file " << sol_file << std::endl;
    }
    std::cout << "Activate: " << (dayVarBranching ? "Vehicle Assignment Branching" : "Arc Flow Branching") << std::endl;
    if(useRCFC)
        std::cout << "Activate Reduced Cost Fix and Cut with fail_factor " << fail_factor << " and init factor " <<
        init_factor << " with cuts?: " << withCuts << std::endl;
    if(useCPC)
        std::cout << "Activate GO-cuts" << std::endl;
    if(!useKPC)
        std::cout << "Deactivate KPC cut separation" << std::endl;
    if(!useSRC)
        std::cout << "Deactivate SRC cut separation" << std::endl;
    else{
        std::cout << "Max number of SRC cuts was set to " << nMaxSRC << (nMaxSRC == 200 ? " (default)" : "") << std::endl;
    }
    if(branchingFrac == 0.0 || branchingFrac == 1.0){
        std::cout << "Activate Random Branching" << std::endl;
    }else{
        std::cout << "Set Branching Fractional-choice to " << branchingFrac << (branchingFrac == 0.5 ? " (default)" : "") << std::endl;
    }

    std::cout << "Set SCIP seed to " << seed << std::endl;


    SCIP_CALL(setUpScip(&scip, dayVarBranching, useRCFC, useCPC, fail_factor, init_factor, withCuts, withSol,
                        useKPC, useSRC, nMaxSRC, branchingFrac));

    SCIP_CALL(getModelDataFromJson(modelData, input_file, ng_parameter));
    modelData->minTravel = true;


    vector<tourVRP> sol_tvrps;
    /** deactivate for no warm start */
    if(!sol_file.empty()){
        SCIP_CALL(readSolution(sol_file, sol_tvrps));
    }

    SCIP_CALL(SCIPprobdataCreate(scip, modelData, sol_tvrps, activate_propagator));

    SCIPsetIntParam(scip, "randomization/randomseedshift", seed);

    /* set RCF parameters */
    auto* pricer = dynamic_cast<ObjPricerVRP*>(SCIPfindObjPricer(scip, "VRP_Pricer"));
    pricer->gap_decay_rc_ = gap_decay_rc;
    pricer->max_depth_rc_ = max_depth_rc;
    pricer->maxSRC_rc_ = maxSRC_rc;
    auto *prop = dynamic_cast<ObjPropVarFixing*>(SCIPfindObjProp(scip, "varFixing"));
    prop->withVeAss_ = useVeAssRC;
    prop->withArcFlow_ = useArcRC;
    prop->noRootFixing_ = noRootFixing;
    auto* cons = dynamic_cast<ConshdlrRCFC*>(SCIPfindObjConshdlr(scip, "RCFC"));
    if(cons != nullptr){
        cons->minvalue_ = minvalue;
        cons->depth_init_ = depth_init;
        cons->onlyFrac_ = onlyFrac;
        cons->onlyFracGap_ = onlyFracGap;
    }


    SCIP_CALL(SCIPsolve(scip));
    /********************
    * Print Solution *
    ********************/

    if(SCIPgetBestSol(scip) == nullptr)
    {
        SCIP_CALL(findInfeasiblityReason(scip));
    }else{
        SCIPprintBestSol(scip, nullptr, false);
    }

    /** Output */
    SCIP_CALL(outputSolution(scip, input_file, out_file, dayVarBranching, branchingFrac));

    /********************
    * Deinitialization *
    ********************/
    delete modelData;
    SCIP_CALL( SCIPfree(&scip) );

    return SCIP_OKAY;
}


int
main(
    int                        argc,
    char**                     argv
)
{
    SCIP_RETCODE retcode;

    retcode = runColumnGenerationModel(argc, argv);
    if( retcode != SCIP_OKAY )
    {
      SCIPprintError(retcode);
      return -1;
    }

    return 0;
}

