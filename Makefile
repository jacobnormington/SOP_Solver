CC = g++
CG = gcc
LINK = -o
CXX_VERSION = -std=c++14 -pthread
OPTIMIZATION = -O3 -c
OPTIMIZATION_LINK = -O3
CXXFLAG = -Wall $(CXX_VERSION)
SRC = ./src
PBB_LIB = ./lib
LKH_LIB = ./lib/LKH

OBJDIR = obj

LKH_OBJS = $(patsubst %,$(OBJDIR)/%,$(_LKH_OBJS))
PBB_OBJS = $(patsubst %,$(OBJDIR)/%,$(_PBB_OBJS))

_PBB_OBJS = main.o solver.o hash.o hungarian.o active_tree.o precedence.o active_allocate.o
_LKH_OBJS = Activate.o AddCandidate.o AddExtraCandidates.o                  \
		AddTourCandidates.o AdjustCandidateSet.o AdjustClusters.o       \
		AllocateStructures.o Ascent.o                                   \
		Best2OptMove.o Best3OptMove.o Best4OptMove.o Best5OptMove.o     \
		BestKOptMove.o BestSpecialOptMove.o                             \
		Between.o Between_SL.o Between_SSL.o                            \
		BridgeGain.o BuildKDTree.o C.o CandidateReport.o                \
		ChooseInitialTour.o Connect.o CreateCandidateSet.o              \
		CreateDelaunayCandidateSet.o CreateNNCandidateSet.o             \
		Create_POPMUSIC_CandidateSet.o CreateQuadrantCandidateSet.o     \
		CTSP_InitialTour.o CVRP_InitialTour.o Delaunay.o                \
		Distance.o Distance_MTSP.o Distance_SOP.o Distance_SPECIAL.o    \
		eprintf.o ERXT.o Excludable.o Exclude.o FindTour.o              \
		FixedOrCommonCandidates.o Flip.o Flip_SL.o Flip_SSL.o           \
		Forbidden.o FreeStructures.o                                    \
		fscanint.o Gain23.o GenerateCandidates.o Genetic.o              \
		GeoConversion.o GetTime.o GreedyTour.o Hashing.o Heap.o         \
		Improvement.o IsBackboneCandidate.o IsCandidate.o               \
		IsCommonEdge.o IsPossibleCandidate.o KSwapKick.o LinKernighan.o \
		LKHmain.o                                                       \
		Make2OptMove.o Make3OptMove.o Make4OptMove.o Make5OptMove.o     \
		MakeKOptMove.o MergeTourWithBestTour.o MergeWithTourIPT.o       \
		Minimum1TreeCost.o MinimumSpanningTree.o                        \
		MTSP2TSP.o MTSP_InitialTour.o MTSP_Report.o                     \
		MTSP_WriteSolution.o                                            \
		NormalizeNodeList.o NormalizeSegmentList.o                      \
		OrderCandidateSet.o PatchCycles.o                               \
		Penalty_ACVRP.o Penalty_BWTSP.o Penalty_CCVRP.o                 \
		Penalty_CVRP.o Penalty_CVRPTW.o Penalty_CTSP.o                  \
		Penalty_1_PDTSP.o Penalty_MLP.o Penalty_M_PDTSP.o               \
		Penalty_M1_PDTSP.o Penalty_MTSP.o Penalty_OVRP.o                \
		Penalty_PDPTW.o Penalty_PDTSP.o Penalty_PDTSPF.o                \
		Penalty_PDTSPL.o Penalty_RCTVRP.o                               \
		Penalty_SOP.o Penalty_TRP.o Penalty_TSPDL.o Penalty_TSPPD.o     \
		Penalty_TSPTW.o Penalty_VRPB.o Penalty_VRPBTW.o Penalty_VRPPD.o \
		PDPTW_Reduce.o printff.o PrintParameters.o                      \
		Random.o ReadCandidates.o ReadEdges.o ReadLine.o                \
		ReadParameters.o ReadPenalties.o ReadProblem.o RecordBestTour.o \
		RecordBetterTour.o RemoveFirstActive.o                          \
		ResetCandidateSet.o RestoreTour.o                               \
		SegmentSize.o Sequence.o SFCTour.o SolveCompressedSubproblem.o  \
		SINTEF_WriteSolution.o SOP_RepairTour.o STTSP2TSP.o             \
		SolveDelaunaySubproblems.o SolveKarpSubproblems.o               \
		SolveKCenterSubproblems.o SolveKMeansSubproblems.o              \
		SolveRoheSubproblems.o SolveSFCSubproblems.o SolveSubproblem.o  \
		SolveSubproblemBorderProblems.o SolveTourSegmentSubproblems.o   \
		SOP_InitialTour.o SOP_Report.o StatusReport.o                   \
		Statistics.o StoreTour.o SymmetrizeCandidateSet.o               \
		TrimCandidateSet.o TSPDL_InitialTour.o TSPTW_MakespanCost.o     \
		TSPTW_Reduce.o VRPB_Reduce.o BIT.o                              \
		WriteCandidates.o WritePenalties.o WriteTour.o                  \
		MergeWithTourGPX2.o gpx.o LKH.o

PROG = sop_solver

$(PROG): $(PBB_OBJS) $(LKH_OBJS)
	$(CC) $(CXXFLAG) $(OPTIMIZATION_LINK) $(LINK) $(PROG) $^

$(OBJDIR)/%.o: $(SRC)/%.cpp
	mkdir -p $(@D)
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $< -o $@ 

$(OBJDIR)/%.o: $(PBB_LIB)/%.cpp
	$(CC) $(CXXFLAG) $(OPTIMIZATION) $< -o $@

$(OBJDIR)/%.o: $(LKH_LIB)/%.c
	$(CG) $(OPTIMIZATION) $< -o $@

clean:
	rm -rf $(OBJDIR)/*.o sop_solver obj
