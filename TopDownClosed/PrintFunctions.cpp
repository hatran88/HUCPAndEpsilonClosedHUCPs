#include <iostream>
#include "Object.h"
#include <string>
#include <map>
#include <vector>
#include <set>
#include "PrintFunctions.h"
#include <unordered_map>

using namespace std;



/**
* @brief: This function print the instance number of each feature.
* @param: wf: weight of featrues
* @retval: Nope
*/
void printFeatureWeight()
{
	std::cout << "The weight of each feature is: " << endl;
	for (auto const& inst : uf)
	{
		std::cout << inst.first << ":" << inst.second << ",";
	}
	std::cout << endl;
}


/**
* @brief: This function print the star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printStarNeighborhood(
	std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc> & SN)
{
	std::cout << "The star neighborhoods are: " << std::endl;
	for (auto const& feat : SN)
	{
		std::cout << feat.first.feature << "." << feat.first.instance << " : {";
		for (auto const& inst: feat.second)
		{
			std::cout << inst.feature << "." << inst.instance << ", ";
		}
		std::cout << " } \n";
	}
	std::cout << endl;
}


/**
* @brief: This function prints neighbors of an instances of in one block.
* @param: vertices: a set of all instances.
* @retval: Nope
*/
void printNeiOneBlock(std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc> oneNeiBlock)
{
	std::cout << "Instances and their neighbors: \n";
	for (auto const& inst : oneNeiBlock)
	{
		std::cout << inst.first.feature << "." << inst.first.instance << " : "; 
		for (auto const& nei : inst.second)
		{
			std::cout << nei.feature << "." << nei.instance << ", ";
		}
		std::cout << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function prints all indexes of a vector.
* @param: P: a set of indexes.
* @retval: Nope
*/
void printIndexVetices(std::set<unsigned int> P)
{
	for (auto const& index : P)
	{
		std::cout << index << ',';
	}
	std::cout << endl;
}


/**
* @brief: This function prints all instances of the input dataset.
* @param: vertices: a set of all instances.
* @retval: Nope
*/
void printVetices(std::unordered_map<ObjWithoutCoordinate, int, myHashFunc> vertices)
{
	for (auto const& vetex : vertices)
	{
		std::cout << vetex.first.feature << "." << vetex.first.instance<<", ";
	}
	std::cout << endl;
}



/**
* @brief: This function prints all instances of the input dataset.
* @param: vertices: a set of all instances.
* @retval: Nope
*/
void printVeticeVec(std::vector<ObjWithoutCoordinate> verticesVec)
{
	for (auto const& inst : verticesVec)
	{
		std::cout << inst.feature << "." << inst.instance << ", ";
	}
	std::cout << endl;
}


/**
* @brief: This function prints all maximal cliques.
* @param: allMaximal: a set of all maximal cliques
* @retval: Nope
*/
void printAllMaximalClique(std::vector<std::set<unsigned int>> allMaxCl)
{
	for (auto const& cl : allMaxCl)
	{
		std::cout << "{ ";
		for (auto const& inst : cl)
		{
			std::cout << inst << ", ";
		}
		std::cout << " } ";
		std::cout << endl;
	}
	std::cout << endl;
}



/**
* @brief: This function prints all maximal cliques.
* @param: allMaximal: a set of all maximal cliques
* @retval: Nope
*/
void printOneMaximalClique(std::set<unsigned int> oneMaxCl)
{
	for (auto const& inst : oneMaxCl)
	{
		std::cout << inst << ", ";
	}
	std::cout << endl;
}



/**
* @brief: This function prints colohashmap
* @param: CoLHM: a hash map of all maximal cliques
* @retval: Nope
*/
void printCoLHM()
{
	std::cout << "The CoLoHM is: \n";
	for (auto const& item : CoLHM)
	{
		// Print key
		std::cout << item.first << " : \n";
		// Print value
		for (auto const& featinst : item.second)
		{
			// Print char or feature
			std::cout << featinst.first << " : ";
			// Print instances
			for (auto const& inst : featinst.second)
			{
				std::cout << inst.feature<<"."<<inst.instance << ", ";
			}
			std::cout << endl;
		}
		std::cout << endl;
	}
	std::cout << endl;
}

/**
* @brief: This function prints all prevelant co-location patterns.
* @param: allPrevPats: all prevelant patterns.
* @retval: Nope
*/
void printPrevPatts(std::unordered_map<std::string, float> & allPrevPats)
{
	std::cout << "All prevalent patterns: " << endl;
	for (auto const& item : allPrevPats)
	{
		std::cout << item.first << " : " << item.second << endl;
	}
	std::cout << endl;
}

/**
* @brief: This function prints sorted prevelant co-location patterns by sizes of patterns
* @param: allPrevPats: all prevelant patterns.
* @retval: Nope
*/
void printSortPrevPats(std::vector<std::pair<std::string, float>> sortPrevPats)
{
	for (auto const& item : sortPrevPats)
	{
		std::cout << item.first << " : " << item.second << "\n ";
	}
	std::cout << endl;
}



/**
* @brief: This function prints patterns by sizes
* @param: classPattBySize: all prevelant patterns are classified by size
* @retval: Nope
*/
void printNumberPattsBySize(std::map<int, std::vector<Patterns>>& classPattBySize)
{
	std::cout << "Patterns are classified by sizes: \n";
	for (auto const& item : classPattBySize)
	{
		std::cout << item.first << " : " << item.second.size();
		std::cout << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function prints top-k pattern of size 3, 4, 5
* @param: classPattBySize: all prevelant patterns are classified by size
*		topk: top k patterns
* @retval: Nope
*/
void printTopKPatts(std::map<int, std::vector<Patterns>>& classPattBySize, int topk)
{
	for (auto const& item : classPattBySize)
	{
		if (item.first == 3)
		{
			int k = 1;
			for (auto const& pat : item.second)
			{
				if (k <= topk)
				{
					
				}
				else
				{

				}
			}
		}
	}


}



