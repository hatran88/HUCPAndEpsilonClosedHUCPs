
#ifndef TOOLFUNCTIONS_H
#define TOOLFUNCTIONS_H


#include <iostream>
#include <map>
#include <vector>
#include <unordered_map>
#include <set>


#include "Object.h"
#include "fptree.h"


#include"tree.hh"



/**
* @brief: This strut sorts prevalent patterns by size of .
* @param: dataList: an input dataset;
* @retval: a hash map that stores each feature name and its instance number.
*/
struct cmpMapBySizePatt {
	bool operator()(const std::string& a, const std::string& b) const
	{
		return a.size() <= b.size();
	}
};


/**
* @brief: This function count the number of instances of each feature.
* @param: dataList: an input dataset;
* @retval: a hash map that stores each feature name and its instance number.
*/
std::map<char, int> countNumberInstance(std::vector<std::vector<std::string>> dataList);


/**
* @brief This function make a grid on an input dataset
* @param dataList: an input dataset; dist_thres: a distance threshold
* @retval a hash map that stores cell id with instances fall in it.
*/
std::map<std::pair<int, int>, std::vector<ObjWithCoord>> makeGrid(
	std::vector<std::vector<std::string> > dataList, 
	float dist_thres);


/**
* @bref This function calculates the distances of instances in the current block.
* @param alll instance in the current block
* @retval A vector that save the distance of instances.
*/
float calculateDistanceTwoInstances(ObjWithCoord currentInst, ObjWithCoord checkInst);


/**
* @brief: This function generates star neighborhoods of instances.
* @param: grid: a grid posing of the input dataset
*         dist_thres: a distance threshold
* @retval: SN: a hash map that stores as <instance, <<neighbors>,<neighbors>>. This struture is different with star neighbors in Join-less
*/
void genStarNeighborhoods(
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid,
	float dist_thres,
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>& SN);


/**
* @brief: This function generates star neighborhoods of features.
* @param: SN: a map of star neighbor instances
* @retval: a hash map that stores as <feat, <feat>>
*/
std::map<char, std::vector<std::set<char>>> genFeatNeibTrans(
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> & SN);


/**
* @brief: This function builds candidate pattern FP_tree of each feature type.
* @param: featNeibTrans: a map of neighbor transaction
*			prev_thres: prevalence threshold
* @retval: a hash map that stores candidates
*/
std::map<char, std::map<std::string, float>> buildCandPattTrees(
	std::map<char, std::vector<std::set<char>>> featNeibTrans);



/**
* @brief: This function gets candidate patterns if these patterns are in all feature trees.
* @param: Treei: a set of FP-Growth tree
* @retval: a set of candidates
*/
std::map<std::string, float> genCands(
	std::map<char, std::map<std::string, float>> Treei, 
	float prev_thres);




/**
* @brief: This function gets all candidates with thier size equal to k
* @param: C: a set of candidates
*			k: size k candidate
* @retval: the set of size k candidates
*/
std::map<std::string, float> getSizeKCands(std::map<std::string, float> C, int k);


/**
* @brief: This function gets instances of size k candidates
* @param: Ck: a set of size k candidates
*			SN: the star neighborhoods
* @retval: a map that stores instances of size k candidates
*/
void findStarInsts(
	std::map<std::string, float> Ck,
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>& SN,
	std::map<std::string, std::vector<std::vector<std::vector<ObjWithoutCoord>>>>& SIk);


/**
* @brief: This function gets table instances of size-2 patterns
* @param: SIk: a set of size k candidates
*			Ck: a set of candidates
* @retval: table instances of size 2 patterns
*/
std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>> collectSize2TableInst(
	std::map<std::string, std::vector<std::vector<std::vector<ObjWithoutCoord>>>>& SIk,
	int k,
	std::map<std::string, float>& Ck);


/**
* @brief: This function generate cartesian product of vector<vector<int>>.
* @param: v: a vector<int>.
* @retval: resultCartesian a vector<vector<int>>.
*/
std::vector<std::vector<ObjWithoutCoord>> cartesianProduct(
	std::vector<std::vector<ObjWithoutCoord>> v);


/**
* @brief: This function find real row instances of candidate patterns.
* @param: CIKplus: a hash map of table instances of size k+1 candidate patterns.
*         CIk: a hash map of table instances of size k candiate patterns.
* @retval: tabInstkplus: a hash map of real table instances of candidate patterns.
*/
std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>> filterCliqueInstances2(
	std::map<std::string, std::vector<std::vector<std::vector<ObjWithoutCoord>>>>& SIkplus,
	std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>> CIk);



/**
* @brief: This function calculate PIs of candidate patterns and filter prevalent patterns.
* @param: CIKplus: a hash map of table instances of size k+1 candidate patterns.
*         totalInstEachFeat: the number of instances of each features.
*         prev_thres: a minimum prevalence threshold
*         k: the size of patterns
*			Pk: prevalent patterns
*			nonPk: non-prevalent patterns
* @retval: Pk: a hash map of prevalent patterns.
*/
void selectPrevalentPatterns(
	std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>>& CIkplus,
	std::map<char, int>& totalInstEachFeat,
	float prev_thres,
	int k,
	std::map<std::string, float, cmpMapBySizePatt>& Pk,
	std::map < std::string, float>& nonPk);


/**
* @brief: This function adds new prevalent patterns into the result
* @param: Pk : the size k prevalent patterns
*			Rk: all prevalent patterns
* @retval: Nope
*/
void addPrePatts(
	std::map<std::string, float, cmpMapBySizePatt>& Pk,
	std::map<std::string, float, cmpMapBySizePatt>& Rk);


/**
* @brief: This function removes checked candidates from the set of candidates
* @param: C : a set of orignial candidates
*		Ck: a set of checked candidates
* @retval: Nope
*/
void delCkFromC(std::map<std::string, float>& C, std::map<std::string, float>& Ck);


/**
*@brief: This function removes super candidates from the set of original candidates by a non-prevalent pattern
* @param : nonPk : a set of prevalent maximal patterns
* C : a set of orignial candidates
* @retval : C : a set of prunned candidates
*/
void superPrun(
	std::map<std::string, float>& nonPk, 
	std::map<std::string, float>& C);

/**
* @brief: This function select coarse prevalent patterns.
* @param: SIK: a hash map of table instances
* @retval: CIk: a hash map of table instances which satify coarse min_prev
*/
void selectCoarsePrevalentPatterns(
	std::map<std::string, std::vector<std::vector<std::vector<ObjWithoutCoord>>>>& SIkplus,
	std::map<char, int>& totalInstEachFeat,
	float prev_thres,
	int k);


/**
* @brief: This function selects closed patterns
* @param: Pk: a hash map of size k prevalent co-location patterns
* @retval: CloseP: a hash map of closed patterns
*/
void filterClosedPatts(std::map<std::string, float, cmpMapBySizePatt>& Pk,
	std::map<float, std::vector<std::string>>& CloseP);


#endif