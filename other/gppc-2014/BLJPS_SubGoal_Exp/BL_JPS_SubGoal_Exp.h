#pragma once
#define DIAG_UNBLOCKED
//#include "Node.h"
#include "PathFindingAlgorithm.h"
//#include "binHeap2.h"
//#include "PairingHeap.h"
//#include "JumpPointNode.h"
#include "Coordinate.h"
#include <limits>
#include "SubgoalDefinitions.h"
#include <deque>
#include <queue>
#include <vector>
//#include <ext/hash_map>
#include "Timer.h"
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
//#include <stdio.h>
//#include <stdlib.h>
using namespace std;


//#define PREPROCESS_STARTING_NODES
//#define PREPROCESS_END_NODES
//#define USE_PAIR_WISE_DISTANCES
struct DiagonalJumpEntry2
{
	DiagonalJumpEntry2(Coordinate _first, Coordinate _second, subgoalId _id) :from(_first), to(_second), id(_id)
	{}
	Coordinate from, to;
	subgoalId id;
};
#ifdef DIAG_UNBLOCKED
#define BLJPS_SUBGOAL_ALG_EXP_BLOCKED_NAME "UNBLOCKED"
#else
#define BLJPS_SUBGOAL_ALG_EXP_BLOCKED_NAME "BLOCKED"
#endif
#ifdef PREPROCESS_END_NODES
#define BLJPS_SUBGOAL_ALG_EXP_ENDS_NAME "-ENDS"
#else
#define BLJPS_SUBGOAL_ALG_EXP_ENDS_NAME ""
#endif

#ifdef PREPROCESS_STARTING_NODES
#define BLJPS_SUBGOAL_ALG_EXP_START_NAME "-STARTS"
#else
#define BLJPS_SUBGOAL_ALG_EXP_START_NAME ""
#endif
#ifdef PAIRWISE_DISTANCES
#define BLJPS_SUBGOAL_ALG_EXP_DIST_NAME "-DIST"
#else
#define BLJPS_SUBGOAL_ALG_EXP_DIST_NAME ""
#endif

#define BLJPS_SUBGOAL_ALG_EXP_NAME (string("BLJPS_SUBGOAL-EXP-")+string(BLJPS_SUBGOAL_ALG_EXP_BLOCKED_NAME)+string(BLJPS_SUBGOAL_ALG_EXP_ENDS_NAME)+string(BLJPS_SUBGOAL_ALG_EXP_START_NAME)+string(BLJPS_SUBGOAL_ALG_EXP_DIST_NAME))

//Directions
// N, NE, E, SE, S, SW, W, NW , StartPosition
#define NO_DIRECTION 8






static std::vector<subgoalId> defaultAbstractPath;	// For default parameters
static std::vector<xyLoc> defaultXYPath;
static std::vector<mapLoc> defaultPath;

class BL_JPS_SUBGOAL_EXP : public PathFindingAlgorithm
{
public:
	//SubgoalGraph();

	void backupPreProcess(){}
	void useBackupData(){}
	bool isCoordinateBlocked(const Coordinate &c){ return false; };
	int getGridWidth(){ return width; };
	int getGridHeight(){ return height; };


#ifdef SG_RUNNING_IN_HOG
	SubgoalGraph(Map* map);
	void LoadMap(Map* map);
#else	// For competition
	BL_JPS_SUBGOAL_EXP(std::vector<bool> &bits, int width, int height, int memoryLimit = MEMORY_LIMIT, int timeLimit = TIME_LIMIT);
	void LoadMap(std::vector<bool> &bits, int width, int height);
#endif
	BL_JPS_SUBGOAL_EXP(const char *filename);	// Read a saved graph from the file
	~BL_JPS_SUBGOAL_EXP();

	void InitializeValues();

	/// Main preprocessing functions
	void SetDirections();		// Initialize the deltaMapLoc array (after we know the map dimensions)
	void IdentifySubgoals();	// Place a subgoal at every corner
	void ComputeClearances();	// Compute clearances (from an obstacle or subgoal) for the 4 cardinal directions
	void LinkSubgoals();		// Add edges so that every subgoal is optimally reachable from each other in the graph
	void PruneSubgoals();		// Prune subgoals that are unnecessary for global travel (make them local subgoals)
	void FinalizeGraph();		// Now that the number of edges are fixed, move them from vector to array to save memory
	void MemoryAnalysis(int memoryLimit = MEMORY_LIMIT);	// Decide if we have extra space for pairwise distances
	// (or if we even have space for the actual subgoals)
	// OUTDATED
	void CalculatePairwiseDistances();		// Only between global subgoals
	void SaveGraph(const char *filename);	// Save all the relevant data to the provided file
	void LoadGraph(const char *filename);	// Load all the relevant data from the provided file

	/// New functions (6/10/2013)
	double UnprunedPairwiseMemory();	// Returns the memory requirement (in bytes) of the pairwise distance matrix with the unpruned graph
	double PrunedPairwiseMemory();	// Returns the memory requirement (in bytes) of the pairwise distance matrix with the pruned graph
	int GetDirectedEdgeCount();	// Returns the total number of directed edges (Only works if the graph is not finalized)

	/// Functions for finding low-level paths / areas
	void GetMoveDirectionDetails	// Decide how many diagonal and cardinal steps to make and in which direction
		(xyLoc & from, xyLoc & to, direction & c, direction & d, int & nTotalMoves, int & nDiagMoves);
	bool IsHReachable		// Check if two locations are h-reachable, save the path to the path
		(xyLoc & from, xyLoc & to, std::vector<mapLoc> & path = defaultPath, bool append = false);
	bool IsQuickHReachable	// Quickly check for a path between the two locations (no guarantee of finding an existing path)
		(xyLoc & from, xyLoc & to, std::vector<mapLoc> & path = defaultPath, bool append = false);
	bool IsLookaheadHReachable	// Similar to QuickHReachable, better chance of finding a path (some extra cost)
		(xyLoc & from, xyLoc & to, std::vector<mapLoc> & path = defaultPath, bool append = false);
	void GetDirectHReachableSubgoals	// Given a location, find all the subgoals it should be connected to
		(xyLoc & from, std::vector<subgoalId> & subgoals);
	//Copy of the above function but will also store the direction in which the subGoal was found
	void GetDirectHReachableSubgoals(xyLoc & from, std::vector<subGoalIdDir> & subgoals,bool jumpPoint=false);

	/// Functions related to pruning
	cost CostOtherPath	// Look for the best path between sg1 and sg2, that does not go through sg
		(subgoalId & sg, subgoalId & sg1, subgoalId & sg2, cost limit = INFINITE_COST);
	bool IsNecessaryToConnect	// Check if we can find the optimal path between sg1 and sg2 without using sg
		(subgoalId sg, subgoalId sg1, subgoalId sg2);
	void AddEdge(subGoalIdDir sg1, subGoalIdDir sg2);	// Add the edge sg1->sg2 to the neighborhood vector (also to the edge vector if sg2 is not pruned)
	void RemoveEdge(subgoalId sg1, subgoalId sg2,char dir);	// Remove the edge sg1->sg2 from the edge vector
	void PruneSubgoal(subgoalId sg, char dir, char checkDirs,char checkDirsReversed);	// Mark the subgoal as local and make the relevant edge changes

	/// Functions for managing clearances
	void SetClearance(mapLoc loc, direction d, int clearance)
	{
		cellInfo[loc] = (cellInfo[loc] & ~(CLEARANCE_MASK << clearanceShift[d])) | (((clearance <= CLEARANCE_LIMIT) ? clearance : 0) << clearanceShift[d]);
	}
	int GetClearance(mapLoc loc, direction d)
	{
		return (cellInfo[loc] >> clearanceShift[d]) & CLEARANCE_MASK;
	}
	int GetTotalClearance(mapLoc loc, direction d)
	{
		int totalClearance = 0;
		while (GetClearance(loc, d) == 0)
		{
			loc += deltaMapLoc[d] * (CLEARANCE_LIMIT);
			totalClearance += (CLEARANCE_LIMIT);
		}
		return totalClearance + GetClearance(loc, d);
	}

#ifdef USE_FLOAT_DIST
	/// Heuristic functions
	inline cost HCost(const xyLoc & l1, const xyLoc &l2)	// The octile distance between l1 and l2
	{
		//int dx = (l1.x>l2.x) ? (l1.x - l2.x) : (l2.x - l1.x);	int dy = (l1.y>l2.y) ? (l1.y - l2.y) : (l2.y - l1.y);
		return Node::estimateDistance(l1,l2);
	}
	inline cost HCost(subgoalId sg1, subgoalId sg2) { return Node::estimateDistance(location[sg1], location[sg2]); }
#else
	/// Heuristic functions
	cost HCost(xyLoc l1, xyLoc l2)	// The octile distance between l1 and l2
	{
		int dx = (l1.x>l2.x) ? (l1.x - l2.x) : (l2.x - l1.x);	int dy = (l1.y>l2.y) ? (l1.y - l2.y) : (l2.y - l1.y);
		return (dx>dy) ? (dx*CARD_COST + dy*DIAG_DIFF) : (dy*CARD_COST + dx*DIAG_DIFF);
	}
	cost HCost(subgoalId sg1, subgoalId sg2) { return HCost(location[sg1], location[sg2]); }

#endif

	void jump(const Coordinate &c, const char dir,const char origDir, std::vector<subGoalIdDir> & subgoals);

	/// Functions for managing the heap
	inline void AddToOpen(subgoalId sg, cost fVal, char dir, unsigned char followUpDir);

	inline void PopMin();
#ifdef REPLACE_POPPED
	void PopReplacableTop(); // Force the heap to remove the top element, without waiting for a replacement
#endif
	inline heapElement GetMin();
	inline void PercolateUp(int index);
	inline void PercolateDown(int index);
	/* The stack and heap work together as follows:
	* If an element is expanded from the heap, it is not immediately removed.
	* All its successors with the same fVal are inserted into the stack and
	* all its successors with a higher fVal are inserted into the heap
	* If at the end of the expansion the stack is empty, only then we remove the element from the heap.
	* Essentially, as long as all its same-fVal successors are around, the original state stays in the heap,
	* mainly as a representative of the stack (for instance, theHeap[0].fVal is used for termination)
	*/

	/// Conversion functions
	/* When generating our version of the map, we add a padding. Any xyLoc refers
	* to a location in the original map, any mapLoc refers to the index of a cell
	* in our padded and linearized map. The following functions do the conversion
	* (height and width take the padding into account, that is, height = 2 + original height)
	*/
	xyLoc ToXYLoc(mapLoc loc) { return xyLoc((loc%width) - 1, (loc / width) - 1); }
	mapLoc ToMapLoc(xyLoc loc) { return (loc.y + 1)*width + loc.x + 1; }
	subgoalId ToSubgoalId(mapLoc loc) { return cellInfo[loc]; }
	void ToXYLocPath	// Convert a mapLocPath to xyLocPath
		(std::vector<mapLoc> & mapLocPath, std::vector<xyLoc> & xyLocPath);

	/// Functions for doing an A* search, using the subgoal graph
	void ResetSearch();	// Reset the search counter and the 'generated' values
	void ConnectStartAndGoalToGraph	// As preparation for a search over the subgoal graph
		(xyLoc & start, xyLoc & goal, subgoalId & sgStart, subgoalId & sgGoal, std::vector<subGoalIdDir> & startSubgoals, std::vector<subGoalIdDir> & goalSubgoals, std::vector<subgoalId> & addedEdges);
	cost SubgoalAStarSearch	// Do a search over the subgoal graph, between two subgoals
		(subgoalId & start, subgoalId & goal, cost limit = INFINITE_COST, std::vector<subgoalId> & abstractPath = defaultAbstractPath);
	cost TryLocalPair	// Try to find a local subgoal pair connecting start and goal (not guaranteed)
		(subgoalId & sgStart, subgoalId & sgGoal, std::vector<subgoalId> & startDirectHReachableSubgoals, std::vector<subgoalId> & goalDirectHReachableSubgoals, std::vector<subgoalId> & localPath);

	/// Functions about using pairwise distances (and covering the cases where pairwise distances are insufficient)
	cost CheckCommonLocal	// Check if there is a single local subgoal connecting start and goal
		(subgoalId & sgStart, subgoalId & sgGoal, std::vector<subGoalIdDir> & startDirectHReachableSubgoals, std::vector<subGoalIdDir> & goalDirectHReachableSubgoals, std::vector<subgoalId> & localPath);
	cost CheckAllLocalPairs	// Try to find a local subgoal pair connecting start and goal (guaranteed)
		(subgoalId & sgStart, subgoalId & sgGoal, std::vector<subgoalId> & startDirectHReachableSubgoals, std::vector<subgoalId> & goalDirectHReachableSubgoals, std::vector<subgoalId> & localPath);
	void AppendOptimalPath		// Given two global subgoals, construct the abstract path between them
		(subgoalId sg1, subgoalId sg2, std::vector<subgoalId> & path);
	void GetGlobalConnections	// Find which global subgoals should be used for lookups (and their distances)
		(bool start, xyLoc & loc, subgoalId & locSg, std::vector<subGoalIdDir> &safeConnections, std::vector<cost> & goalGlobalCosts);
	cost LookupOptimalPath	// Main function for looking up pairwise distances
		(xyLoc & startLoc, xyLoc & goalLoc, subgoalId & start, subgoalId & goal, std::vector<subGoalIdDir> & startDirectConnections, std::vector<subGoalIdDir> & goalDirectConnections, std::vector<subgoalId> & path);

	/// Main functions for getting a low level path between two locations
	cost SubgoalAStarSearch	// Either use search or lookup, convert abstract path to low level path
		(xyLoc  start, xyLoc  goal, cost limit = INFINITE_COST, std::vector<xyLoc> & thePath = defaultXYPath);
	//cost findSolution(xyLoc & start, xyLoc & goal, std::vector<xyLoc> & thePath = defaultXYPath);
	void findSolution(int sX, int sY, int _eX, int _eY, std::vector<xyLoc> & sol = defaultXYPath);

	void preprocessBoundaryLookupTables();
	int getSpaceIdY(short spaceX, short spaceY);

	bool isSpaceIdY(int spaceId, short spaceX, short spaceY);
	int getSpaceIdX(short spaceX, short spaceY);
	bool isSpaceIdX(int spaceId, short spaceX, short spaceY);
	bool directSolution(short sX, short sY, short eX, short eY, vector<Coordinate> & sol);
	void buildDiagEndJumps();

    int binarySearchL(const vector<DiagonalJumpEntry2> & v, short val);
	Coordinate diagonalDirection(Coordinate c, int dir);
	void getAllEndNodes(short eX, short eY, vector<subGoalIdDir> &endNodes);
	int jump(const Coordinate &c, const char dir);

	/// Statistics/output functions
	double initializeTime, searchTime, finalizeTime;
	void SearchStatistics()
	{
#ifdef SG_STATISTICS
		std::cout << "Subgoal statistics:" << std::endl;
		std::cout << "Total initialize time: " << initializeTime*1000.0 << "ms" << std::endl;
		std::cout << "Total search time: " << searchTime*1000.0 << "ms" << std::endl;
		std::cout << "Total finalize time: " << finalizeTime*1000.0 << "ms" << std::endl;
#endif
	}
	void PrintGraphStatistics();
	void GetHeapCapacity()
	{
#ifdef SG_STATISTICS
		std::cout << "The heap capacity is: " << theHeap.capacity() << ", which is about " << (theHeap.capacity()*sizeof(heapElement)) / 1024 << "kb." << std::endl;
#endif
	}
	bool UsingSubgoals(){ return useSubgoals; }
#ifdef SG_RUNNING_IN_HOG
	void OpenGLDraw(const MapEnvironment *env);
#endif

	/// Auxiliary (mostly inline) functions:
#ifdef USE_BOOL
	bool IsTraversable(mapLoc loc){ return traversable[loc]; }
	void SetTraversable(mapLoc loc){ traversable[loc] = true; }

	bool IsSubgoal(mapLoc loc){ return subgoal[loc]; }
	void SetSubgoal(mapLoc loc){ subgoal[loc] = true; }

	bool IsPruned(subgoalId sg){ return pruned[sg]; }
	void SetPruned(subgoalId sg){ pruned[sg] = true; }
	void SetUnpruned(subgoalId sg){ pruned[sg] = false; }

	bool HasExtraEdge(subgoalId sg){ return hasExtraEdge[sg]; }
	void SetExtraEdgeFlag(subgoalId sg){ hasExtraEdge[sg] = true; }
	void RemoveExtraEdgeFlag(subgoalId sg){ hasExtraEdge[sg] = false; }

	bool IsOpen(subgoalId sg){ return open[sg]; }
	void SetOpen(subgoalId sg){ open[sg] = true; }
	void SetClosed(subgoalId sg){ open[sg] = false; }
#else
	bool IsTraversable(mapLoc loc){ return traversable[(loc >> 3)] & (1 << (loc & 7)); }
	bool IsTraversable(xyLoc loc){ return IsTraversable(ToMapLoc(loc)); }
	void SetTraversable(mapLoc loc){ traversable[loc >> 3] |= (1 << (loc & 7)); }



	bool IsSubgoal(mapLoc loc){ return subgoal[(loc >> 3)] & (1 << (loc & 7)); }
	void SetSubgoal(mapLoc loc){ subgoal[loc >> 3] |= (1 << (loc & 7)); }

	bool isCompletlyPruned(subgoalId sg){ return pruned[sg] == 0xFF; };
	bool IsPruned(subgoalId sg){ return pruned[sg]!=0 ; }
	bool IsPruned(subgoalId sg, char dir){ return pruned[sg] & (1 << (dir & 7)); }
	void SetPruned(subgoalId sg, char dir){ pruned[sg ] |= (1 << (dir & 7)); }
	//void SetUnpruned(subgoalId sg){ pruned[sg] =0; }

	bool HasExtraEdge(subgoalId sg){ return hasExtraEdge[(sg >> 3)] & (1 << (sg & 7)); }
	void SetExtraEdgeFlag(subgoalId sg){ hasExtraEdge[sg >> 3] |= (1 << (sg & 7)); }
	void RemoveExtraEdgeFlag(subgoalId sg){ hasExtraEdge[sg >> 3] &= ~(1 << (sg & 7)); }

	bool IsOpen(subgoalId sg){ return open[(sg >> 3)] & (1 << (sg & 7)); }
	void SetOpen(subgoalId sg){ open[sg >> 3] |= (1 << (sg & 7)); }
	void SetClosed(subgoalId sg){ open[sg >> 3] &= ~(1 << (sg & 7)); }
	char getDir(int from, int to);
	char getDir2(int from, int to);
	void jumpSpecial(const Coordinate &c, const char dir, const char origDir, std::vector<subGoalIdDir> & subgoals);
	int jumpNew(const Coordinate &c, const char dir, std::vector<subGoalIdDir> & subgoals);
	bool getJumpPointNew(Coordinate s, const char direction, Coordinate & jp);
	short binarySearchL(const vector<pair<short, short> > & v, short val);
	short binarySearchR(const vector<pair<short, short> > & v, short val);
	void preProcessSubGoalGraph();
	bool directSolution(Coordinate s, Coordinate e, bool &firstConnection);
	void getDiagCoordinate2(const Coordinate& nxt, vector<Coordinate> & solution);
	char getDir2(const Coordinate& c, const Coordinate& nxt);
	char getDir(const Coordinate& c, const Coordinate& nxt);
	inline int mag(const int i)
	{
		return i < 1 ? -1 : 1;
	}
#endif
	void preProcessGrid();

	unsigned char forcedNeighbours(const xyLoc &coord, const int dir);
	void dumpPreprocessedDataToFile(const char * fileName) { SaveGraph(fileName); }
	void readPreprocessedDataToFile(const char * fileName) { LoadGraph(fileName); }

private:
	/// Variables for controling the behaviour of the subgoal graph
	bool useSubgoals;	// False if the subgoal graph exceeds the memory limit (signal to use buckets instead)
	char usePairwise;	// 0: No pairwise distances 1: Pairwise distances 2: Also pairwise connections (not used)
	bool keepLocalEdges;

	unsigned int height, width, mapSize, nSubgoals, nGlobalSubgoals, nLocalSubgoals;
	/* A note about directions:
	* There are 8 directions. We label them 0-7 as follows: N = 0, NE, E, SE, S, SW, W, NW = 7
	* Even numbers are cardinal directions and odd numbers are diagonal directions.
	* For a diagonal direction d, its associated cardinal directions are d+1 and d-1
	*/
	int deltaMapLoc[24];	/* Given a mapLoc loc and a direction d, this array can be used to determine the
							* loc' we end up at if we follow d from loc (loc' = loc + deltaMapLoc[d]). Note that,
							* even there is only 8 directions, we keep 3 copies for each direction for a total of
							* 24 values. This is to minimize the number of mod operations. This array should be
							* initialized after we know the width and height of the (padded) map.
							*/

	/// Information about each location on the map
#ifdef USE_BOOL
	bool* traversable;		// 1: Traversable, 0: Not traversable (size = mapSize) // std::vector<bool>
	bool* subgoal;			// 1: Subgoal, 0: Not a subgoal (size = mapSize) // std::vector<bool>
#else
	char* traversable;
	char* subgoal;
#endif
	subgoalId* cellInfo;	/* For subgoals, keep the subgoalId to directly link a cell to the subgoal graph;
							* For traversable cells, keep the obstacle-subgoal clearances for the four
							* cardinal directions (to speed up the linking to nearby subgoals)
							* Currently not used for obstacles (maybe use it to label different chunks of obstacles?)
							* (size = mapSize)
							*/

	/// Information about each subgoal
	xyLoc* location;		// Keep the location of the subgoal (size = nSubgoals)
#ifdef USE_BOOL
	bool* pruned;			// Keep track of whether a subgoal is pruned or not (size = nSubgoals) // std::vector<bool>
#else
	char* pruned;
#endif

	/// Variables for storing edges
	bool finalized;	// Whether to use the vectors or edges to lookup neighbors
	std::vector<std::vector<subGoalIdDir> > edgeVector;	// For storing the edges that will be used for search (size = nSubgoals)
	std::vector<std::vector<subGoalIdDir> > neighborhoodVector;	// For storing the subgoals that should be considered to see if the subgoal can be pruned
	std::vector<std::vector<subGoalIdDir> > incomingConnectionsVector;
	vector<vector<DiagonalJumpEntry2>> diagonalJumps[4];
	vector<vector<pair<short, short> > > jumpLookup[4];

	subGoalIdDir** neighbors;	// Once the subgoal graph is pruned and the number of edges are fixed, use this instead of vectors (for memory)
	uint16_t* nNeighbors;	// Keep the number of neighbors for each subgoal

	subGoalIdDir** localNeighbors;	// Once the subgoal graph is pruned and the number of edges are fixed, use this instead of vectors (for memory)
	uint16_t* nLocalNeighbors;	// Keep the number of local neighbors for each subgoal

#ifdef USE_BOOL
	bool* hasExtraEdge;		/* Even after the graph is finalized, the number of neighbors of a subgoal
							* can change when a new start and goal state are linked to the graph.
							* However, handled carefully, the change in this number can be limited
							* to be at most 1. Therefore, we allocate space for an extra neighbor
							* for each subgoal and keep track of whether there is an extra neighbor
							* with this boolean vector
							*/
#else
	char* hasExtraEdge;
#endif

	/// Search related variables
	/* Instead of keeping the list of generated states for each search and use hash tables
	* to access that information, preallocate the memory for search info for each subgoal
	* and directly access the info with the subgoalId. This is manageable since we only keep
	* the information for subgoals, not any arbitrary cell in the map.
	*/
	uint16_t* generated;	// Keep track of the last time a subgoal has been generated
	uint16_t search;	/* The current search number. When a search terminates, only increment this
						* to reset the previous search. search == generated[subgoalX] iff
						* subgoalX is generated for the current search
						*/
	cost* gCost;		// g-cost of each subgoal
	subgoalId* parent;	// Parent of each subgoal

#ifdef USE_BOOL
	bool* open;			// 0: In closed list; 1: In open list. (generated[x] < search: x not yet generated) // std::vector<bool>
#else
	char* open;
#endif

#ifdef REPLACE_POPPED
	bool canReplaceTop;
#endif
	std::vector<heapElement> theHeap;
	std::vector<heapElement> theStack;

	/// Pairwise-distance variables
	cost** dist;

	//Boundary Lookup Table
	vector<vector<short>> xBoundaryPoints, yBoundaryPoints;

};

