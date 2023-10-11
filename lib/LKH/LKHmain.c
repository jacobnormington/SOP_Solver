#include "LKHmain.h"

/*
 * This file contains the main function of the program.
 */

void Clean_Mem() {
    if(CostMatrix) free(CostMatrix);
    if(Depot) free(Depot);
    if(FirstActive) free(FirstActive);
    if(LastActive) free(LastActive);
    if(FirstConstraint) free(FirstConstraint);
    if(FirstNode) free(FirstNode);
    if(FirstSegment) free(FirstSegment);
    if(FirstSSegment) free(FirstSSegment);
}

void Reset_Parameter() {
    AscentCandidates = 0;
    Asymmetric = 0;
    BackboneTrials = 0;     /* Number of backbone trials in each run */
    Backtracking = 0;       /* Specifies whether backtracking is used for*/
    BestCost = 0;      /* Cost of the tour in BestTour */
    BestPenalty = 0;   /* Penalty of the tour in BestTour */
    BetterCost = 0;    /* Cost of the tour stored in BetterTour */
    BetterPenalty = 0; /* Penalty of the tour stored in BetterTour */
    BWTSP_B = 0;            /* Number of black nodes in a BWTSP instance */ 
    BWTSP_Q = 0;            /* Maximum number of subsequent white nodes in a
                           BWTSP instance */ 
    BWTSP_L = 0;            /* Maximum length of any path between two black
                           nodes in a BTWSP instance */
    CacheMask = 0;  /* Mask for indexing the cache */
    CandidateFiles = 0;     /* Number of CANDIDATE_FILEs */
    EdgeFiles = 0;          /* Number of EDGE_FILEs */
    CostMatrix = 0;        /* Cost matrix */
    CurrentGain = 0;
    CurrentPenalty = 0;
    DemandDimension = 0;    /* Number of commodities in a M-PDTSP instance */
    Depot = 0;
    Dimension = 0;          /* Number of nodes in the problem */
    DimensionSaved = 0;     /* Saved value of Dimension */
    Dim = 0;                /* DimensionSaved - Salesmen + 1 */
    DistanceLimit = 0;   /* Maxixim route distance for a CVRP instance */
    Excess = 0;          /* Maximum alpha-value allowed for any 
                        candidate edge is set to Excess times the 
                        absolute value of the lower bound of a 
                        solution tour */
    ExtraCandidates = 0;    /* Number of extra neighbors to be added to 
                        the candidate set of each node */
    FirstActive = 0, LastActive = 0; /* First and last node in the list 
                                of "active" nodes */
    FirstConstraint = 0; /* First constraint in the list of SOP 
                             precedence constraints */
    FirstNode = 0;        /* First node in the list of nodes */
    FirstSegment = 0;  /* A pointer to the first segment in the cyclic 
                        list of segments */
    FirstSSegment = 0;   /* A pointer to the first super segment in
                          the cyclic list of segments */
    Gain23Used = 0; /* Specifies whether Gain23 is used */
    GainCriterionUsed = 0;  /* Specifies whether L&K's gain criterion is 
                        used */
    GridSize = 0;        /* Grid size used in toroidal instances */
    GroupSize = 0;  /* Desired initial size of each segment */
    SGroupSize = 0; /* Desired initial size of each super segment */
    Groups = 0;     /* Current number of segments */
    SGroups = 0;    /* Current number of super segments */
    Hash = 0;  /* Hash value corresponding to the current tour */
    InitialPeriod = 0;      /* Length of the first period in the ascent */
    InitialStepSize = 0;    /* Initial step size used in the ascent */
    InitialTourFraction = 0; /* Fraction of the initial tour to be 
                            constructed by INITIAL_TOUR_FILE edges */
    Kicks = 0;      /* Specifies the number of K-swap-kicks */
    KickType = 0;   /* Specifies K for a K-swap-kick */
    LastLine = 0; /* Last input line */
    LowerBound = 0;      /* Lower bound found by the ascent */
    M = 0;          /* The M-value is used when solving an ATSP-
                instance by transforming it to a STSP-instance */
    MaxBreadth = 0; /* The maximum number of candidate edges 
                considered at each level of the search for a move */
    MaxCandidates = 0;      /* Maximum number of candidate edges to be 
                        associated with each node */
    DistanceLimit = 0; /* Maxixim route distance for a CVRP instance */
    MaxMatrixDimension = 0; /* Maximum dimension for an explicit
                        cost matrix */
    MaxSwaps = 0;   /* Maximum number of swaps made during the 
                search for a move */
    MaxTrials = 0;  /* Maximum number of trials in each run */
    MergeTourFiles = 0;     /* Number of MERGE_TOUR_FILEs */
    MoveType = 0;   /* Specifies the sequantial move type to be used 
                in local search. A value K >= 2 signifies 
                that a k-opt moves are tried for k <= K */
    MoveTypeSpecial = 0; /* A special (3- or 5-opt) move is used */
    NodeSet = 0;  /* Array of all nodes */
    Norm = 0;       /* Measure of a 1-tree's discrepancy from a tour */
    NonsequentialMoveType = 0; /* Specifies the nonsequential move type to
                            be used in local search. A value 
                            L >= 4 signifies that nonsequential
                            l-opt moves are tried for l <= L */
    Optimum = 0;       /* Known optimal tour length. 
                        If StopAtOptimum is 1, a run will be 
                        terminated as soon as a tour length 
                        becomes equal this value */
    PatchingA = 0;  /* Specifies the maximum number of alternating
                cycles to be used for patching disjunct cycles */
    PatchingC = 0;  /* Specifies the maximum number of disjoint cycles to be 
                patched (by one or more alternating cycles) */
    PenaltyGain = 0;
    Precision = 0;  /* Internal precision in the representation of 
                transformed distances */
    PredSucCostAvailable = 0; /* PredCost and SucCost are available */
    POPMUSIC_InitialTour = 0;  /* Specifies whether the first POPMUSIC tour
                           is used as initial tour for LK */
    POPMUSIC_MaxNeighbors = 0; /* Maximum number of nearest neighbors used 
                           as candidates in iterated 3-opt */
    POPMUSIC_SampleSize = 0;   /* The sample size */
    POPMUSIC_Solutions = 0;    /* Number of solutions to generate */
    POPMUSIC_Trials = 0;       /* Maximum trials used for iterated 3-opt */
    Rand = 0;           /* Table of random values */
    Recombination = 0;        /* IPT or GPX2 */
    RestrictedSearch = 0;     /* Specifies whether the choice of the first
                          edge to be broken is restricted */
    Reversed = 0; /* Boolean used to indicate whether a tour has 
                been reversed */
    Run = 0;        /* Current run number */
    Runs = 0;       /* Total number of runs */
    Scale = 0;      /* Scale factor for Euclidean and ATT instances */ 
    ServiceTime = 0; /* Service time for a CVRP instance */
    Serial = 0;
    Seed = 0;  /* Initial seed for random number generation */
    StartTime = 0;       /* Time when execution starts */
    StopAtOptimum = 0;      /* Specifies whether a run will be terminated if 
                        the tour length becomes equal to Optimum */
    Subgradient = 0;        /* Specifies whether the Pi-values should be 
                        determined by subgradient optimization */
    SubproblemSize = 0;     /* Number of nodes in a subproblem */
    SubsequentMoveType = 0; /* Specifies the move type to be used for all 
                        moves following the first move in a sequence 
                        of moves. The value K >= 2 signifies that a 
                        K-opt move is to be used */
    SubsequentMoveTypeSpecial = 0; /* A special (3- or 5-opt) subsequent move
                              is used */
    SubsequentPatching = 0; /* Species whether patching is used for 
                        subsequent moves */
    Swaps = 0;      /* Number of swaps made during a tentative move */
    OldSwaps = 0;   /* Saved number of swaps */
    TimeLimit = 0;       /* The time limit in seconds */
    TotalDemand = 0;        /* Sum of demands for a CVRP instance */
    TraceLevel = 0; /* Specifies the level of detail of the output 
                given during the solution process. 
                The value 0 signifies a minimum amount of 
                output. The higher the value is the more 
                information is given */
    Trial = 0;      /* Ordinal number of the current trial */
    TSPTW_CurrentMakespanCost = 0;
    TSPTW_Makespan = 0;

    return;
}

int LKH(char *problem_file, bool initial_LKHRun)
{
    Reset_Parameter();
    GainType Cost, OldOptimum;
    double Time, LastTime;
    Node *N;
    int i;

    ParameterFileName = problem_file;
    ReadParameters();
    ProblemFileName = problem_file;
    
    Gain23Used = 0;
    KickType = 4;
    MaxSwaps = 0;
    MoveType = 5;
    MoveTypeSpecial = 1;
    MaxPopulationSize = 10;
    TraceLevel = 0;

    StartTime = LastTime = GetTime();
    MaxMatrixDimension = 20000;
    MergeWithTour = Recombination == IPT ? MergeWithTourIPT :
        MergeWithTourGPX2;
    ReadProblem();

    if (!initial_LKHRun) {
        pthread_mutex_lock(&Sol_lock);
        Read_BBCost();
        pthread_mutex_unlock(&Sol_lock);
    }

    if (SubproblemSize > 0) {
        if (DelaunayPartitioning)
            SolveDelaunaySubproblems();
        else if (KarpPartitioning)
            SolveKarpSubproblems();
        else if (KCenterPartitioning)
            SolveKCenterSubproblems();
        else if (KMeansPartitioning)
            SolveKMeansSubproblems();
        else if (RohePartitioning)
            SolveRoheSubproblems();
        else if (MoorePartitioning || SierpinskiPartitioning)
            SolveSFCSubproblems();
        else
            SolveTourSegmentSubproblems();
        return EXIT_SUCCESS;
    }
    AllocateStructures();
    if (ProblemType == TSPTW)
        TSPTW_Reduce();
    if (ProblemType == VRPB || ProblemType == VRPBTW)
        VRPB_Reduce();
    if (ProblemType == PDPTW)
        PDPTW_Reduce();

    CreateCandidateSet();
    InitializeStatistics();

    if (Norm != 0 || Penalty) {
        Norm = 9999;
        BestCost = PLUS_INFINITY;
        BestPenalty = CurrentPenalty = PLUS_INFINITY;
    } else {
        /* The ascent has solved the problem! */
        Optimum = BestCost = (GainType) LowerBound;
        UpdateStatistics(Optimum, GetTime() - LastTime);
        RecordBetterTour();
        RecordBestTour();
        CurrentPenalty = PLUS_INFINITY;
        BestPenalty = CurrentPenalty = Penalty ? Penalty() : 0;
        WriteTour(OutputTourFileName, BestTour, BestCost);
        WriteTour(TourFileName, BestTour, BestCost);
        Runs = 0;
    }

    /* Find a specified number (Runs) of local optima */

    MaxTrials = 10000;

    while (true) {
        LastTime = GetTime();
        if (BB_Complete || BB_SolFound || LastTime - StartTime >= TimeLimit) {
            /*
            if (TraceLevel >= 1)
                printff("*** Time limit exceeded ***\n");
            */
            break;
        }
        Cost = FindTour();      /* using the Lin-Kernighan heuristic */
        if (MaxPopulationSize > 1 && !TSPTW_Makespan) {
            /* Genetic algorithm */
            int i;
            for (i = 0; i < PopulationSize; i++) {
                GainType OldPenalty = CurrentPenalty;
                GainType OldCost = Cost;
                Cost = MergeTourWithIndividual(i);
                if (TraceLevel >= 1 &&
                    (CurrentPenalty < OldPenalty ||
                     (CurrentPenalty == OldPenalty && Cost < OldCost))) {
                    
                    /*
                    if (CurrentPenalty)
                        printff("  Merged with %d: Cost = " GainFormat,
                                i + 1, Cost);
                    else
                        printff("  Merged with %d: Cost = " GainFormat "_"
                                GainFormat, i + 1, CurrentPenalty, Cost);
                    */

                    if (Optimum != MINUS_INFINITY && Optimum != 0) {
                        if (ProblemType != CCVRP && ProblemType != TRP &&
                            ProblemType != MLP &&
                            MTSPObjective != MINMAX &&
                            MTSPObjective != MINMAX_SIZE)
                            printff(", Gap = %0.4f%%",
                                    100.0 * (Cost - Optimum) / Optimum);
                        else
                            printff(", Gap = %0.4f%%",
                                    100.0 * (CurrentPenalty - Optimum) /
                                    Optimum);
                    }
                    //printff("\n");
                }
            }
            if (!HasFitness(CurrentPenalty, Cost)) {
                if (PopulationSize < MaxPopulationSize) {
                    AddToPopulation(CurrentPenalty, Cost);
                    /*
                    if (TraceLevel >= 1)
                        PrintPopulation();
                    */
                } else if (SmallerFitness(CurrentPenalty, Cost,
                                          PopulationSize - 1)) {
                    i = ReplacementIndividual(CurrentPenalty, Cost);
                    ReplaceIndividualWithTour(i, CurrentPenalty, Cost);
                    /*
                    if (TraceLevel >= 1)
                        PrintPopulation();
                    */
                }
            }
        } else if (Run > 1 && !TSPTW_Makespan)
            Cost = MergeTourWithBestTour();
        if (CurrentPenalty < BestPenalty ||
            (CurrentPenalty == BestPenalty && Cost < BestCost)) {
            BestPenalty = CurrentPenalty;
            BestCost = Cost;
            RecordBetterTour();
            RecordBestTour();
            WriteTour(TourFileName, BestTour, BestCost);
        }
        OldOptimum = Optimum;
        if (!Penalty ||
            (MTSPObjective != MINMAX && MTSPObjective != MINMAX_SIZE)) {
            if (CurrentPenalty == 0 && Cost < Optimum)
                Optimum = Cost;
        } else if (CurrentPenalty < Optimum)
            Optimum = CurrentPenalty;
        if (Optimum < OldOptimum) {
            printff("*** New OPTIMUM = " GainFormat " ***\n", Optimum);
            if (FirstNode->InputSuc) {
                Node *N = FirstNode;
                while ((N = N->InputSuc = N->Suc) != FirstNode);
            }
        }
        Time = fabs(GetTime() - LastTime);
        UpdateStatistics(Cost, Time);
        if (TraceLevel >= 1 && Cost != PLUS_INFINITY) {
            //printff("Run %d: ", Run);
            //StatusReport(Cost, LastTime, "");
            //printff("\n");
        }
        if (StopAtOptimum && MaxPopulationSize >= 1) {
            if (ProblemType != CCVRP && ProblemType != TRP &&
                ProblemType != MLP &&
                MTSPObjective != MINMAX &&
                MTSPObjective != MINMAX_SIZE ?
                CurrentPenalty == 0 && Cost == Optimum :
                CurrentPenalty == Optimum) {
                Runs = Run;
                break;
            }
        }
        if (PopulationSize >= 2 &&
            (PopulationSize == MaxPopulationSize ||
             Run >= 2 * MaxPopulationSize) && Run < Runs) {
            Node *N;
            int Parent1, Parent2;
            Parent1 = LinearSelection(PopulationSize, 1.25);
            do
                Parent2 = LinearSelection(PopulationSize, 1.25);
            while (Parent2 == Parent1);
            ApplyCrossover(Parent1, Parent2);
            N = FirstNode;
            do {
                if (ProblemType != HCP && ProblemType != HPP) {
                    int d = C(N, N->Suc);
                    AddCandidate(N, N->Suc, d, INT_MAX);
                    AddCandidate(N->Suc, N, d, INT_MAX);
                }
                N = N->InitialSuc = N->Suc;
            }
            while (N != FirstNode);
        }
        SRandom(++Seed);
    }
    PrintStatistics();
    if (Salesmen > 1) {
        if (Dimension == DimensionSaved) {
            for (i = 1; i <= Dimension; i++) {
                N = &NodeSet[BestTour[i - 1]];
                (N->Suc = &NodeSet[BestTour[i]])->Pred = N;
            }
        } else {
            for (i = 1; i <= DimensionSaved; i++) {
                Node *N1 = &NodeSet[BestTour[i - 1]];
                Node *N2 = &NodeSet[BestTour[i]];
                Node *M1 = &NodeSet[N1->Id + DimensionSaved];
                Node *M2 = &NodeSet[N2->Id + DimensionSaved];
                (M1->Suc = N1)->Pred = M1;
                (N1->Suc = M2)->Pred = N1;
                (M2->Suc = N2)->Pred = M2;
            }
        }
        CurrentPenalty = BestPenalty;
        MTSP_Report(BestPenalty, BestCost);
        MTSP_WriteSolution(MTSPSolutionFileName, BestPenalty, BestCost);
        SINTEF_WriteSolution(SINTEFSolutionFileName, BestCost);
    }
    if (ProblemType == ACVRP ||
        ProblemType == BWTSP ||
        ProblemType == CCVRP ||
        ProblemType == CTSP ||
        ProblemType == CVRP ||
        ProblemType == CVRPTW ||
        ProblemType == MLP ||
        ProblemType == M_PDTSP ||
        ProblemType == M1_PDTSP ||
        MTSPObjective != -1 ||
        ProblemType == ONE_PDTSP ||
        ProblemType == OVRP ||
        ProblemType == PDTSP ||
        ProblemType == PDTSPL ||
        ProblemType == PDPTW ||
        ProblemType == RCTVRP ||
        ProblemType == RCTVRPTW ||
        ProblemType == SOP ||
        ProblemType == TRP ||
        ProblemType == TSPTW ||
        ProblemType == VRPB ||
        ProblemType == VRPBTW || ProblemType == VRPPD) {
        //printff("Best %s solution:\n", Type);
        CurrentPenalty = BestPenalty;
        //SOP_Report(BestCost);
    }
    //printff("\n");
    //Clean_Mem();
    return BestCost;
}
