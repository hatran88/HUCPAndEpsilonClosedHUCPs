#include <iostream>
#include <map>>
#include <vector>
#include <set>

#include "Object.h"
#include "toolFunctions.h"


#include "tree.hh"


#ifndef PRINTFUNCTIONS_H
#define PRINTFUNCTIONS_H


using namespace std;


/**
* @brief: This function print the instance number of each feature.
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
void printInstanceNumber(std::map<char, int> totalInstEachFeat);

/**
* @brief: This function print the star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printGrid(std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid);


/**
* @brief: This function print candidate patterns and their table instances
* @param: SIk: a hash map of candidate pateterns and their table instasnces
* @retval: Nope
*/
void printStarInstance(std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN);


/**
* @brief: This function prints feature neighborhood transactions.
* @param: featNeibTrans: a hash map of feature transactions
* @retval: Nope
*/
void printFeatNeibTrans(std::map<char, std::vector<std::set<char>>> featNeibTrans);


/**
* @brief: This function prints the star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printStarNeighborhood(std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> SN);


/**
* @brief: This function prints neighbor pairs
* @param: NP: the neighbor pairs
* @retval: Nope
*/
void printNeiborPairs(std::map<std::pair<ObjWithoutCoord, ObjWithoutCoord>, int>& NP);


/**
* @brief: This function prints size k candidates
* @param: Ck: a set of size k candidates
* @retval: Nope
*/
void printCands(std::map<std::string, float> Ck);


/**
* @brief: This function prints size k candidates
* @param: Ck: a set of size k candidates
* @retval: Nope
*/
void printSizekCands(std::map<std::string, float>  Ck);


/**
* @brief: This function prints all row instances of candidates.
* @param: SIk: a hash map of row instances
* @retval: Nope
*/
void printOneRowInst(std::set<ObjWithoutCoord> rowInst);


/**
* @brief: This function prints all row instances of candidates.
* @param: SIk: a hash map of row instances
* @retval: Nope
*/
void printCoarseTableInsts(std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>> SIk);


void printCoarseTableInsts_2(std::map<ObjWithoutCoord, std::vector<std::vector<ObjWithoutCoord>>> coarseTable);

/**
* @brief: This function prints the candidate got by FP-growth tree
* @param: Treei: a set of candidates grouped by features*			
* @retval: Nope
*/
void printTreei(std::map<char, std::map<std::string, float>> Treei);


/**
* @brief: This function prints the features in the prunning tree
* @param: Treei: a set of candidates grouped by features*
* @retval: Nope
*/
void printPrunTree(tree<char> treeC);


/**
* @brief: This function prints all true row instances of candidates.
* @param: CIc: the true row instances of the current pattern
* @retval: Nope
*/
void printTrueTableInsts(std::vector<std::vector<ObjWithoutCoord>> CIc);


/**
* @brief: This function print all prevalent patterns.
* @param: PK: a hash map of prevalent patterns
* @retval: Nope
*/
void printPrevalentPattern(std::map<std::string, float, cmpMapBySizePatt> Pk);


/**
* @brief: This function prints execution time of construct tree of candidate by size
* @param: NonPk: a vector of non prevalent patterns
* @retval: Nope
*/
void printExectionTimeConstTreeBySizeCandidate(std::map<int, double>& constrTreeTime);


/**
* @brief: This function prints closed patterns
* @param: CloseP: a map of closed patterns
* @retval: Nope
*/
void printClosedPatts(std::map<float, std::vector<std::string>>& CloseP);




#endif
