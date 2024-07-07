#ifndef TOOLFUNCTIONS_H
#define TOLLFUNCTIONS_H

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
//#include <bitset>
#include <ctime>
#include <unordered_set>
#include <boost/dynamic_bitset.hpp>


#include "Object.h"
#include "SetNumberPointFloat.h"
#include "Patterns.h"

typedef boost::dynamic_bitset<> BitSet;



// Statement some global variables
extern std::vector<ObjWithoutCoordinate> verticesVec;
extern std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc> neiOneBlock;

// Definition global variables
extern clock_t timeConstCoLHashmap;
extern double totalTimeConstCoLHashmap;

extern clock_t timeComputeUPIs;
extern double totalTimeComputeUPIs;


extern std::unordered_map<std::string, std::unordered_map<char, std::unordered_set<ObjWithoutCoordinate, myHashFunc>>> CoLHM;
extern std::unordered_map<std::string, std::unordered_map<char, std::unordered_set <ObjWithoutCoordinate, myHashFunc>>>::iterator it;

extern std::unordered_map<char, float> uf; // Save the weight of each feature


//extern std::unordered_map<std::string, int> allMaxCl;
//extern int numOfMaxCl;

/*
* @brief: This function count the number of instances of each feature.
* @param: dataList: an input dataset;
* @retval: wf: (a global variable) a hash map that stores each feature name and its instance number.
*/
void getFeatureWeight(std::vector<std::vector<std::string>> dataList);


/**
* @brief This function make a grid on an input dataset
* @param dataList: an input dataset; dist_thres: a distance threshold
* @retval a hash map that stores cell id with instances fall in it.
*/
std::map<std::pair<int, int>, std::vector<Objects>> makeGrid(
	std::vector<std::vector<std::string>> &dataList, float dist_thres);


/**
* @bref This function calculates the distances of instances in the current block.
* @param alll instance in the current block
* @retval A vector that save the distance of instances.
*/
float calculateDistanceTwoInstances(Objects currentInst, Objects checkInst);


/**
* @brief: This function generates star neighborhoods of instances.
* @param: grid: a grid posing of the input dataset
*         dist_thres: a distance threshold
* @retval: SN: a hash map that stores as <instance, <<neighbors>,<neighbors>>. This struture is different with star neighbors in Join-less
*/
void genStarNeighborhoods(
	std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc>& neighbors,
	std::unordered_map<ObjWithoutCoordinate, int, myHashFunc> & vertices,
	std::map<std::pair<int, int>, std::vector<Objects>>& grid,
	float dist_thres);


/**
* @brief: This function gets the neighbors of instances in the current block
* @param: neighbors: neighboring instances
* @retval: size_t: the maximal degeree
*/
void getNeiOneBlock(std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc>& neighbors);


/**
* @brief: This function sorts vertices by their degeneracy.
* @param: vertices: a set of indexes of all vertices
*         neighbors: star neighbors of all vertices
* @retval: a vector of oredered vertices sorted by degeneracy
*/
void OrderedDegeneracy(std::vector<ObjWithoutCoordinate>& orderDeg,
	std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc>& neiOneBlock);


/**
* @brief: This function get maximal cliques
* @param: vertices: a set of indexes of all vertices
*         neighbors: star neighbors of all vertices
* @retval: a vector of oredered vertices sorted by degeneracy
*/
void BronKerboschDeg(std::vector<ObjWithoutCoordinate> P,
	std::vector<ObjWithoutCoordinate> X,
	std::vector<ObjWithoutCoordinate>& orderDeg);


/**
* @brief: This function enumerates all maximal cliques based on bron-kerbosch pivot.
* @param: vertices: a set of all vertices (instances)
*         neighbors: a hashmap stores neighbors of each instance
* @retval: allMaxCliques: all maximal cliques.
*/
void BronKerboschPivot(std::vector<ObjWithoutCoordinate> P, std::vector<ObjWithoutCoordinate> R, std::vector<ObjWithoutCoordinate> X);



/**
* @brief: This function buids the co-location hashmap structure.
* @param: allMaxCliques: all maximal cliques.
* @retval: CoLoHM: a co-location hashmap.
*/
void constructCoLoHashmap(std::vector<ObjWithoutCoordinate> R);


/**
* @brief: This function call PIs and filter prevalent patterns.
* @param: CoLHM: the co-location hash map.
*         pre_thress: a minimum prevalent threshold.
* @retval: allPrevPats: a co-location hashmap.
*/
void calPIby2Phases(std::map<string, float>& allPats);


/*
*@brief: This function calculate the pi of the current pattern
*@param: pattern: the current pattern which need to calculate its pi.
*        superPatterns:
*        cpHashMap: a co-location pattern hash map
*@retval: PI: the participation index of the pattern
*         superPatterns: all super patterns of the current pattern
*/
float calculateWeightWithSuperPatterns2Phases(std::string pattern, int sizeofPatttern);


/*
*@brief: This function generate all combination of a char vector
*@param: c: a vector of charl; combo: the size of combo; C(n, m) = n!/(m!(n-m)!)
*@retval: a vector of sub vectors of the vector
*/
template<typename T>
std::string getCombination(const T& c, int combo);

template<typename T>
std::vector<std::string> combo(const T& c, int k);


/**
* @brief: This function filters prevalent patterns
* @param: CoLHM: the co-location hash map.
*         pre_thress: a minimum prevalent threshold.
* @retval: allPrevPats: a co-location hashmap.
*/
void filterPrevPatts(
	std::map<string, float>& allPrevPats,
	std::map<std::string, float>& allPats,
	float prev_thres);


// Comparator function to sort pairs
// according to second value
bool cmp(std::pair<std::string, float>& a,
	std::pair<std::string, float>& b);


// to value in a (key-value) pairs
void sorthashmap(std::unordered_map<std::string, float>& M,
	std::vector<std::pair<std::string, float>>& A);



/**
* @brief: This function classify patterns by sizes
* @param: allPrevPats: all patterns
*         classPattBySize: classified pattens by sizes
* @retval:
*/
void classifyPattsBySize(std::unordered_map<string, float>& allPrevPats, std::map<int, std::vector<Patterns>>& classPattBySize);



/**
* @brief: This function calculate the avg error of delta closed
* @param: allPrevPats: all patterns
*
* @retval:
*/
void calculateAvgError(std::unordered_map<std::string, float>& allPCP,
	float delta_thres,
	std::pair<int, float>& results);


/*
* * @brief: This function calculates utility index and filter delta-closed PCPs
* @param: CoLHM: a co-location hash map
*	utility_thres: the utility threshold
*	delta_thres: the delta threshold
* Method: delCPCP: a set of all delta closed patterns
*	first: get all key and sort by size put into Ck
*	second: get the largest size of key
*	third: calculate upi of this key c, if upi is larger > threshold, check its super set in delCPCP,
* if having any supers that |upi(super) - upi(c)| <= delta, c is deleted, else put c into delCPCP.
	fourth: generate subsets of c and put into new candidate set.
* @retval: a set of delta-closed PCPs
*/
void filterDeltaClosedPCPs(std::unordered_map<std::string, float>& deltaCPCP,
	std::unordered_map<std::string, float>& allCPCP,
	std::unordered_map<std::string, float>& closedPCP,
	float utility_thres, float delta_thres, float alpha, float beta);



/*
*@brief: This function calculates UPI of each pattern
*@param: pattern: the current pattern which need to calculate its pi.
*        superPatterns:
*        cpHashMap: a co-location pattern hash map
*			uf: utility of feature
*@retval: UPI
*         superPatterns: all super patterns of the current pattern
*/
float calculatePIWithSuperPatterns(std::string pattern, int sizeofPatttern, float alpha, float beta);


/**
* @brief: This function filters prevalent patterns
* @param: CoLHM: the co-location hash map.
*         pre_thress: a minimum prevalent threshold.
* @retval: allPrevPats: a co-location hashmap.
*/
int findMaxSizeCand(std::map<std::string, int>& candPatts);



/*
*@brief: This function filters closed PCPS
*@param:
*@retval: none
*/
bool findClosedPCPs(std::unordered_map<std::string, float>& closedPCP,
	std::string oneLCand, float UPI);



/*
*@brief: This function filters delta PCPS
*@param: deltaCPCP: a set of all delta CPCP
*        oneLCand: the current pattern
*        delta_thres
*@retval: none
*/
bool findDeltaPats(std::unordered_map<std::string, float>& deltaCPCP,
	std::string oneLCand, float UPI, float delta_thres);









#endif