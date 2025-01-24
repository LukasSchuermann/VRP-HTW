
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
#include "branchingrule_vehiclearc.h"
#include "ConshdlrArcflow.h"
#include "ConshdlrDayVar.h"
#include "ConshdlrVehicleArc.h"
#include "ConshdlrKPC.h"
#include "ConshdlrSRC.h"
#include "ConshdlrVehicle.h"
#include "ConshdlrNVehicle.h"
#include "tourVRP.h"
#include "tools_vrp.h"
#include "eventhdlr_nodeInit.hpp"
#include "printer.h"
#include "heurDayVarRounding.h"
#include "prop_varfixing.h"


#include "iostream"
#include "fstream"
#include "json.hpp"

using json = nlohmann::json;

static
SCIP_RETCODE readArguments(
    int                 argc,               /**< number of shell parameters */
    char**              argv,               /**< array with shell parameters */
    std::string&        input_file,
    int*                branchingRule,
    int*                seed,
    std::string&        sol_file,
    std::string&        out_file,
    bool*               withSol,
    bool*               useSRC,
    int*                nMaxSRC,
    bool*               useKPC,
    double*             gap_decay,
    int*                max_depth,
    int*                maxSRCrcf,
    bool*               noRootFixing,
    double*             branchingFrac,
    int*                timeout
)
{
    char usage[SCIP_MAXSTRLEN];
    char* locstr;

    assert( argc >= 1 );
    assert( argv != nullptr );
    assert( input_file.empty() );

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
            *branchingRule = atoi(locstr);
            if(*branchingRule != 1 && *branchingRule != 2 && *branchingRule != 3)
            {
                fprintf(stderr, "Invalid branching strategy -> Choose 1 for vehicle assignment and 2 for arc flow"
                                " and 3 vehicle arc branching. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            free(locstr);
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
        }else if ( ! strcmp(argv[i], "-nR")){
            *noRootFixing = true;
        }else if ( ! strcmp(argv[i], "-timeout")){
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing timeout. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            *timeout = atoi(argv[i]);
        }
    }

    return SCIP_OKAY;
}

static
SCIP_RETCODE setUpScip(
        SCIP**              scip,
        int                 branchingRule,
        bool                withSol,
        bool                useKPC,
        bool                useSRC,
        int                 nMaxSRC,
        double              branchingFrac,
        int                 timeout
) {
    /* initialize SCIP */
    SCIP_CALL(SCIPcreate(scip));

    /* include default SCIP plugins */
    SCIP_CALL( SCIPincludeDefaultPlugins(*scip) );


    /* include branching rules */
    SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleVehicle(*scip, 5000), TRUE));
    switch (branchingRule) {
        case BRANCHTYPE_ARCFLOW:
            SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleArcflow(*scip, 50002, branchingFrac), TRUE));
            SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleDayVar(*scip, 50001, branchingFrac), TRUE));
            SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleVehicleArc(*scip, 50000, branchingFrac), TRUE));
            break;
        case BRANCHTYPE_DAYVAR:
            SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleArcflow(*scip, 50001, branchingFrac), TRUE));
            SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleDayVar(*scip, 50002, branchingFrac), TRUE));
            SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleVehicleArc(*scip, 50000, branchingFrac), TRUE));
            break;
        case BRANCHTYPE_VEHICLEARC:
            SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleArcflow(*scip, 50000, branchingFrac), TRUE));
            SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleDayVar(*scip, 50001, branchingFrac), TRUE));
            SCIP_CALL(SCIPincludeObjBranchrule(*scip, new ObjBranchruleVehicleArc(*scip, 50002, branchingFrac), TRUE));
            break;
        default:
            throw std::runtime_error("No feasible branching rule.");
    }

    /* turn off all separation algorithms */
    SCIP_CALL( SCIPsetSeparating(*scip, SCIP_PARAMSETTING_OFF, TRUE) );

    /* turn off primal heuristics if optimal solution is read */
    if(withSol){
        SCIP_CALL( SCIPsetHeuristics(*scip, SCIP_PARAMSETTING_OFF, TRUE) );
    }

    /* include propagator */
    SCIP_CALL(SCIPincludeObjProp(*scip, new ObjPropVarFixing(*scip, false), TRUE));

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
    SCIP_CALL( SCIPincludeObjConshdlr(*scip, new ConshdlrVehicleArc(*scip), TRUE) );

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
    SCIP_CALL( SCIPsetRealParam(*scip, "limits/time", timeout));


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
    int             branchingRule,
    double          frac
)
{
    size_t lastSlashPos = path.find_last_of('/');
    std::string filename = path.substr(lastSlashPos + 1);

    size_t underscorePos = filename.find('_');
    std::string beforeUnderscore = filename.substr(0, underscorePos);

    size_t nPos = filename.find('n');
    size_t pPos = filename.find('p');
    int nValue;
    float pValue;
    if(nPos < pPos){
        nValue = std::stoi(filename.substr(nPos + 1, pPos - nPos - 1));
        pValue = std::stof(filename.substr(pPos + 1));
    }else{
        nValue = std::stoi(filename.substr(nPos + 1, pPos - nPos - 1));
        pValue = -1;
    }
    auto* pricer_obj = dynamic_cast<ObjPricerVRP*>(SCIPfindObjPricer(scip, "VRP_Pricer"));
    auto* pricer = SCIPfindPricer(scip, "VRP_Pricer");
    std::ofstream sol_file;
    sol_file.open(out_file, std::ios_base::app);
    sol_file << "Base;N;P;SolTime;nNodes;Gap;maxDepth;nVars;nPricerCalls;pricerTime;conflicts;nProhibit;nEnforce;brTime;brFrac;nOther" << std::endl;
    sol_file << beforeUnderscore << ";" << nValue << ";" << pValue << ";" << SCIPgetSolvingTime(scip) << ";";
    sol_file << SCIPgetNNodes(scip) << ";" << SCIPgetGap(scip) * 100 << ";" << SCIPgetMaxDepth(scip);
    sol_file << ";" << SCIPgetNVars(scip) << ";" << SCIPpricerGetNCalls(pricer) << ";" << SCIPpricerGetTime(pricer);
    sol_file << ";" << pricer_obj->nConflicts_;

    if(branchingRule == BRANCHTYPE_DAYVAR){
        auto* bRule = dynamic_cast<ObjBranchruleDayVar*>(SCIPfindObjBranchrule(scip, "DayVarBranching"));
        sol_file << ";" << bRule->nProhibit_ << ";" << bRule->nEnforce_ << ";";
        sol_file << SCIPbranchruleGetTime(SCIPfindBranchrule(scip, "DayVarBranching")) << ";" << frac;
        sol_file << ";" << SCIPbranchruleGetNChildren(SCIPfindBranchrule(scip, "ArcFlowBranching")) << std::endl;
    }else if(branchingRule == BRANCHTYPE_ARCFLOW){
        auto* bRule = dynamic_cast<ObjBranchruleArcflow*>(SCIPfindObjBranchrule(scip, "ArcFlowBranching"));
        sol_file << ";" << bRule->nProhibit_ << ";" << bRule->nEnforce_ << ";";
        sol_file << SCIPbranchruleGetTime(SCIPfindBranchrule(scip, "ArcFlowBranching")) << ";" << frac;
        sol_file << ";" << SCIPbranchruleGetNChildren(SCIPfindBranchrule(scip, "DayVarBranching")) << std::endl;
    }else{
        assert(branchingRule == BRANCHTYPE_VEHICLEARC);
        auto* bRule = dynamic_cast<ObjBranchruleVehicleArc*>(SCIPfindObjBranchrule(scip, "VehicleArcBranching"));
        sol_file << ";" << bRule->nProhibit_ << ";" << bRule->nEnforce_ << ";";
        sol_file << SCIPbranchruleGetTime(SCIPfindBranchrule(scip, "VehicleArcBranching")) << ";" << frac;
        sol_file << ";" << SCIPbranchruleGetNChildren(SCIPfindBranchrule(scip, "DayVarBranching")) << std::endl;
    }


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
    std::string sol_file;
    std::string out_file;
    bool withSol = false;
    bool useSRC = true;
    int nMaxSRC = 200;
    bool useKPC = true;
    int timeout = 1800;
    /* reduced costs fixing parameters */
    double gap_decay_rc = 0.8;
    int max_depth_rc = 0;
    int maxSRC_rc = 75;
    bool noRootFixing = false;
    /* branching parameters */
    int branchingRule = BRANCHTYPE_DAYVAR;
    double branchingFrac = 0.5;

    SCIP_CALL(readArguments(argc, argv, input_file, &branchingRule, &seed, sol_file, out_file, &withSol, &useSRC,
                            &nMaxSRC, &useKPC, &gap_decay_rc, &max_depth_rc, &maxSRC_rc, &noRootFixing, &branchingFrac,
                            &timeout));

    if(withSol){
        std::cout << "Include solution from file " << sol_file << std::endl;
    }
    if(branchingRule == BRANCHTYPE_DAYVAR){
        std::cout << "Activate: Vehicle Assignment Branching" << std::endl;
    }else if(branchingRule == BRANCHTYPE_ARCFLOW){
        std::cout << "Activate: Arc Flow Branching" << std::endl;
    }else{
        assert(branchingRule == BRANCHTYPE_VEHICLEARC);
        std::cout << "Activate: Vehicle Arc Branching" << std::endl;
    }

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


    SCIP_CALL(setUpScip(&scip, branchingRule, withSol, useKPC, useSRC, nMaxSRC, branchingFrac, timeout));

    SCIP_CALL(getModelDataFromJson(modelData, input_file, ng_parameter));
    modelData->minTravel = true;


    vector<tourVRP> sol_tvrps;
    /** deactivate for no warm start */
    if(!sol_file.empty()){
        SCIP_CALL(readSolution(sol_file, sol_tvrps));
    }

    SCIP_CALL(SCIPprobdataCreate(scip, modelData, sol_tvrps));

    SCIPsetIntParam(scip, "randomization/randomseedshift", seed);

    /* set RCF parameters */
    auto* pricer = dynamic_cast<ObjPricerVRP*>(SCIPfindObjPricer(scip, "VRP_Pricer"));
    pricer->gap_decay_rc_ = gap_decay_rc;
    pricer->max_depth_rc_ = max_depth_rc;
    pricer->maxSRC_rc_ = maxSRC_rc;
    auto *prop = dynamic_cast<ObjPropVarFixing*>(SCIPfindObjProp(scip, "varFixing"));
    prop->noRootFixing_ = noRootFixing;

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

//    SCIPprintStatistics(scip, nullptr);
    /** Output */
    if(!out_file.empty()){
        SCIP_CALL(outputSolution(scip, input_file, out_file, branchingRule, branchingFrac));
    }


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

