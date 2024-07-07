
#ifndef PRINTFUNCTIONS_H
#define PRINTFUNCTIONS_H

#include <iostream>
#include "Object.h"
#include <string>
#include <map>
#include <vector>
#include <set>
#include <unordered_map>

#include "ToolFunctions.h"


using namespace std;


/**
* @brief: This function print the instance number of each feature.
* @param: wf: weight of featrues
* @retval: Nope
*/
void printFeatureWeight();


/**
* @brief: This function print the star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printStarNeighborhood(
	std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc> & SN);


/**
* @brief: This function prints all instances of the input dataset.
* @param: vertices: a set of all instances.
* @retval: Nope
*/
void printVetices(std::unordered_map<ObjWithoutCoordinate, int, myHashFunc> & vertices);


/**
* @brief: This function prints neighbors of an instances of in one block.
* @param: vertices: a set of all instances.
* @retval: Nope
*/
void printNeiOneBlock(std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc> oneNeiBlock);


/**
* @brief: This function prints all instances of the input dataset.
* @param: vertices: a set of all instances.
* @retval: Nope
*/
void printVeticeVec(std::vector<ObjWithoutCoordinate> verticesVec);




/**
* @brief: This function prints all indexes of a vector.
* @param: P: a set of indexes.
* @retval: Nope
*/
void printIndexVetices(std::set<unsigned int> P);


/**
* @brief: This function prints all maximal cliques.
* @param: allMaximal: a set of all maximal cliques
* @retval: Nope
*/
void printAllMaximalClique(std::vector<std::set<unsigned int>> allMaxCl);



/**
* @brief: This function prints all maximal cliques.
* @param: allMaximal: a set of all maximal cliques
* @retval: Nope
*/
void printOneMaximalClique(std::set<unsigned int> oneMaxCl);



/**
* @brief: This function prints colohashmap
* @param: CoLHM: a hash map of all maximal cliques
* @retval: Nope
*/
void printCoLHM();



/**
* @brief: This function prints all prevelant co-location patterns.
* @param: allPrevPats: all prevelant patterns.
* @retval: Nope
*/
void printPrevPatts(std::unordered_map<std::string, float>& allPrevPats);


/**
* @brief: This function prints sorted prevelant co-location patterns by sizes of patterns
* @param: allPrevPats: all prevelant patterns.
* @retval: Nope
*/
void printSortPrevPats(std::vector<std::pair<std::string, float>> sortPrevPats);


/**
* @brief: This function prints patterns by sizes
* @param: classPattBySize: all prevelant patterns are classified by size
* @retval: Nope
*/
void printNumberPattsBySize(std::map<int, std::vector<Patterns>>& classPattBySize);


/**
* @brief: This function prints top-k pattern of size 3, 4, 5
* @param: classPattBySize: all prevelant patterns are classified by size
*		topk: top k patterns
* @retval: Nope
*/
void printTopKPatts(std::map<int, std::vector<Patterns>>& classPattBySize, int tokp);




#endif
