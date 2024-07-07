#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <ctime>


#include "printFunctions.h"
#include "Object.h"
#include "toolFunctions.h"


#include "tree.hh"


using namespace std;

/**
* @brief: This function print the instance number of each feature.
* @param: totalInstEachFeat: instance number of each feature
* @retval: Nope
*/
void printInstanceNumber(std::map<char, int> totalInstEachFeat)
{
	std::cout << "The instance number of each feature is: " << std::endl;
	for (auto const& inst : totalInstEachFeat)
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
void printGrid(std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid)
{
	for (auto const& cell : grid)
	{
		std::cout << "The cell is: " << "(" << cell.first.first << cell.first.second << "): ";
		std::cout << "Instances: ";
		for (auto const& obj : cell.second)
		{
			std::cout << obj.feature << "." << obj.instance << ", ";
		}
		std::cout << endl;
	}

}


/**
* @brief: This function print candidate patterns and their table instances
* @param: SIk: a hash map of candidate pateterns and their table instasnces
* @retval: Nope
*/
void printStarInstance(std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN)
{
	std::cout << "The star instances: " << endl;
	for (auto const& si : SN)
	{
		std::cout << si.first.feature << "." << si.first.instance << " : ";
		for (auto const& tbist : si.second)
		{
			std::cout << tbist.feature << "." << tbist.instance << ", ";
		}
		std::cout << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function print the star neighborhoods.
* @param: SN: the star neighborhoods of each instances.
* @retval: Nope
*/
void printStarNeighborhood(std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> SN)
{
	std::cout << "The star neighborhoods are: " << endl;
	for (auto const& feat : SN)
	{
		std::cout << "Feature: " << feat.first << endl;
		for (auto const& inst : feat.second)
		{
			std::cout << "Instance: " << inst.first.feature << "." << inst.first.instance << ": ";
			for (auto const& neigh : inst.second)
			{
				std::cout << neigh.feature << "." << neigh.instance << ", ";
			}
			std::cout << endl;
		}
		std::cout << endl;
	}
}


/**
* @brief: This function prints neighbor pairs
* @param: NP: the neighbor pairs
* @retval: Nope
*/
void printNeiborPairs(std::map<std::pair<ObjWithoutCoord, ObjWithoutCoord>, int>& NP)
{
	for (auto const& item : NP)
	{
		std::cout << "<" << item.first.first.feature << "." << item.first.first.instance
			<< ", " << item.first.second.feature << "." << item.first.second.instance << ">" << endl;
	}
}


/**
* @brief: This function prints feature neighborhood transactions.
* @param: featNeibTrans: a hash map of feature transactions
* @retval: Nope
*/
void printFeatNeibTrans(std::map<char, std::vector<std::set<char>>> featNeibTrans)
{
	std::cout << "Feature neighborhood transactions. \n";
	for (auto const& item : featNeibTrans)
	{
		std::cout << item.first << " : "<< endl;
		for (auto const& trans : item.second)
		{
			for (auto const& inst : trans)
			{
				std::cout << inst << ", ";
			}
			std::cout << endl;
		}		
		std::cout << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function prints all candidates
* @param: C: a set of all candidates
* @retval: Nope
*/
void printCands(std::map<std::string, float> C)
{
	std::cout << "All candidates are: \n";
	for (auto const& cand : C)
	{
		std::cout << cand.first<<":"<<cand.second << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function prints size k candidates
* @param: Ck: a set of size k candidates
* @retval: Nope
*/
void printSizekCands(std::map<std::string, float>  Ck)
{
	std::cout << "Size k candidates are: \n";
	for (auto const& cand : Ck)
	{
		std::cout << cand.first << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function prints all row instances of candidates.
* @param: SIk: a hash map of row instances
* @retval: Nope
*/
void printOneRowInst(std::set<ObjWithoutCoord> rowInst)
{
	for (auto const& inst : rowInst)
	{
		std::cout << inst.feature << "." << inst.instance << ", ";
	}
	std::cout << endl;
}


/**
* @brief: This function prints all row instances of candidates.
* @param: SIk: a hash map of row instances
* @retval: Nope
*/
void printCoarseTableInsts(std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>> SIk)
{
	std::cout << "All row instances of patterns are: \n";
	for (auto const& item : SIk)
	{
		std::cout << item.first << " : " << endl;
		for (auto const& row : item.second)
		{
			std::cout << "{ ";
			for (auto const& inst : row)
			{
				std::cout << inst.feature << "." << inst.instance << ", ";
			}
			std::cout << "}" << endl;
		}
		std::cout << endl;
	}
	std::cout << endl;
}


void printCoarseTableInsts_2(std::map<ObjWithoutCoord, std::vector<std::vector<ObjWithoutCoord>>> coarseTable)
{
	for (auto const& item : coarseTable)
	{
		std::cout << item.first.feature << "." << item.first.instance << ": " << endl;
		
		for (auto const& feat : item.second)
		{
			std::cout << "{ ";
			for (auto const& inst : feat)
			{
				std::cout << inst.feature << "." << inst.instance << ", ";
			}
			std::cout << " } "<<endl;
		}
		std::cout << endl;
	}
	cout << endl;
}

/**
* @brief: This function prints the candidate got by FP-growth tree
* @param: Treei: a set of candidates grouped by features*
* @retval: Nope
*/
void printTreei(std::map<char, std::map<std::string, float>> Treei)
{
	std::cout << "Candidates got by FP-growth trees: \n";
	for (auto const& item : Treei)
	{
		std::cout << item.first << " : ";
		for (auto const& cand : item.second)
		{
			std::cout << cand.first << ": " << cand.second << endl;
		}
		std::cout << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function prints the features in the prunning tree
* @param: Treei: a set of candidates grouped by features*
* @retval: Nope
*/
void printPrunTree(tree<char> treeC)
{
	tree<char>::iterator itpt = treeC.begin();
	while (itpt != treeC.end())
	{
		for (int i = 0; i < treeC.depth(itpt); ++i)
			cout << "   ";
		cout << (*itpt) << endl;
		++itpt;
	}
	std::cout << endl << endl;
}



/**
* @brief: This function prints all true row instances of candidates.
* @param: CIc: the true row instances of the current pattern
* @retval: Nope
*/
void printTrueTableInsts(std::vector<std::vector<ObjWithoutCoord>> CIc)
{
	std::cout << "The true table instance of the pattern is: \n";
	for (auto const& item : CIc)
	{
		std::cout << "{ ";
		for (auto const& row : item)
		{			
			std::cout << row.feature << "." << row.instance << ", ";			
		}
		std::cout << "}" << endl;
		std::cout << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function print all prevalent patterns.
* @param: PK: a hash map of prevalent patterns
* @retval: Nope
*/
void printPrevalentPattern(std::map<std::string, float, cmpMapBySizePatt> Pk)
{
	std::cout << "The prevalent patterns are: " << endl;
	for (auto const& pattern : Pk)
	{
		std::cout << pattern.first << " : " << pattern.second << endl;	
	}
	std::cout << endl;
}




/**
* @brief: This function prints execution time of construct tree of candidate by size
* @param: NonPk: a vector of non prevalent patterns
* @retval: Nope
*/
void printExectionTimeConstTreeBySizeCandidate(std::map<int, double>& constrTreeTime)
{
	std::cout << "The execution of collecting co-location instances by sizes: \n";
	for (auto const& item : constrTreeTime)
	{
		std::cout << item.first << " : " << item.second / CLOCKS_PER_SEC << endl;
	}
	std::cout << endl;
}


/**
* @brief: This function prints closed patterns
* @param: CloseP: a map of closed patterns
* @retval: Nope
*/
void printClosedPatts(std::map<float, std::vector<std::string>>& CloseP)
{
	std::cout << "All the closed patterns are: \n";
	for (auto const& item : CloseP)
	{
		std::cout << item.first << " : ";
		for (auto const& c : item.second)
		{
			std::cout << c << ", ";
		}
		std::cout << endl;
	}
	std::cout << endl;
}