#include <iostream>
#include <map>
#include <vector>
#include <unordered_set>
#include <string>
#include <algorithm>
#include <math.h>
#include <cstdlib>
#include <numeric>

#include "toolFunctions.h"
#include "Object.h"
#include "printFunctions.h"
#include "fptree.h"
#include "SetNumberPointFloat.h"


#include "tree.hh"


#include <boost/algorithm/string.hpp>


using namespace std;




/**
* @brief: This function count the number of instances of each feature.
* @param: dataList: an input dataset;
* @retval: a hash map that stores each feature name and its instance number.
*/
std::map<char, int> countNumberInstance(std::vector<std::vector<std::string>> dataList)
{
	std::map<char, int> totalInstNumEachFeat;

	char feature;
	int instance;

	for (std::vector<std::string> vec : dataList)
	{
		// calc the cell id
		feature = vec[0][0]; // get the feature	

		if (totalInstNumEachFeat.empty())
		{
			totalInstNumEachFeat.insert({ feature, 1 });
		}
		else if (totalInstNumEachFeat.find(feature) == totalInstNumEachFeat.end())
		{
			totalInstNumEachFeat.insert({ feature, 1 });
		}
		else
		{
			instance = totalInstNumEachFeat.find(feature)->second;
			instance += 1;
			totalInstNumEachFeat[feature] = instance;
		}
	}

	return totalInstNumEachFeat;
}


/**
* @brief This function make a grid on an input dataset
* @param dataList: an input dataset; dist_thres: a distance threshold
* @retval a hash map that stores cell id with instances fall in it.
*/
std::map<std::pair<int, int>, std::vector<ObjWithCoord>> makeGrid(
	std::vector<std::vector<std::string> > dataList, 
	float dist_thres)
{
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid;

	int cell_x, cell_y;
	std::pair<int, int> cell_key;

	for (std::vector<std::string> vec : dataList)
	{
		// calc the cell id
		cell_x = ceil(std::stof(vec[2]) / dist_thres); // x coordinate
		cell_y = ceil(std::stof(vec[3]) / dist_thres); // y coordinate		
													   // make keys
		cell_key = make_pair(cell_x, cell_y);
		//cout << "One cell key: " << cell_key.first << ":" << cell_key.second << endl;

		// package the value
		std::vector<ObjWithCoord> value;
		ObjWithCoord instance = { vec[0][0], std::stoi(vec[1]), std::stof(vec[2]), std::stof(vec[3]) };
		value.push_back(instance);

		// check if this key is exsiting?
		if (grid.empty()) //The grid is empty
		{
			grid[cell_key] = value;
		}
		else if (grid.find(cell_key) == grid.end()) // if the key is not exist
		{
			grid[cell_key] = value;
		}
		else // the key has already existed
		{
			std::vector<ObjWithCoord> old_value = grid.find(cell_key)->second;
			old_value.push_back(instance);
			grid[cell_key] = old_value;
		}
	}

	return grid;
}


/**
* @bref This function calculates the distances of instances in the current block.
* @param alll instance in the current block
* @retval A vector that save the distance of instances.
*/
float calculateDistanceTwoInstances(ObjWithCoord currentInst, ObjWithCoord checkInst)
{
	return  sqrt((currentInst.x - checkInst.x)*(currentInst.x - checkInst.x)
		+ (currentInst.y - checkInst.y)*(currentInst.y - checkInst.y));
}


/**
* @brief: This function generates star neighborhoods of instances.
* @param: grid: a grid posing of the input dataset
*         dist_thres: a distance threshold
* @retval: SN: a hash map that stores as <instance, <<neighbors>,<neighbors>>. This struture is different with star neighbors in Join-less
*/
void genStarNeighborhoods(
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid, 
	float dist_thres,
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> & SN)
{	
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>::iterator itSN;

	int i, j; // the index of cells in x and y
	std::vector<std::pair<int, int>> fiveCells;
	float dist;
	std::vector<ObjWithCoord> fiveCellInst; // save all instance in the five cells and two neighboring instances.
		
	// Loop each cell in grid and computation neighbors in five cells
	for (std::map<std::pair<int, int>, std::vector<ObjWithCoord>>::iterator it = grid.begin(); it != grid.end(); ++it)
	{
		// Create the five cells
		i = it->first.first;
		j = it->first.second;

		fiveCells.push_back(it->first);
		fiveCells.push_back(make_pair(i, j + 1));
		fiveCells.push_back(make_pair(i + 1, j + 1));
		fiveCells.push_back(make_pair(i + 1, j));
		fiveCells.push_back(make_pair(i + 1, j - 1));

		// Get all instances in the current cell and five cells
		for (auto const& cell : fiveCells)
		{
			if (grid.count(cell))
			{
				fiveCellInst.insert(fiveCellInst.end(), grid.find(cell)->second.begin(), grid.find(cell)->second.end());
			}
		}

		// Sort all instances in the five cells
		std::sort(fiveCellInst.begin(), fiveCellInst.end());

		// Iterator each instance in the currentCellInst and check its neighbors		
		for (auto const& currentInst : it->second)
		{
			for (auto const& checkInst : fiveCellInst) // check with each instance in the five cells
			{
				if (currentInst.feature != checkInst.feature) // only check two instances belong to different features.
				{
					dist = calculateDistanceTwoInstances(currentInst, checkInst);
					if (dist <= dist_thres) // the two instances have neighbor relationship
					{
						// Put them into SN
						ObjWithoutCoord curInt = ObjWithoutCoord{ currentInst.feature, currentInst.instance };
						ObjWithoutCoord cheInt = ObjWithoutCoord{ checkInst.feature, checkInst.instance };
						// Put current instance
						itSN = SN.find(curInt);
						if (itSN != SN.end()) // this instance has already existed
						{
							// update value							
							itSN->second.insert(cheInt);
						}
						else // This feature has not existed in SN, directly put into SN	
						{
							// Put into SN							
							SN.insert({ curInt, std::set<ObjWithoutCoord> {cheInt} });
						}
						// Put cheInt instance
						itSN = SN.find(cheInt);
						if (itSN != SN.end()) // this instance has already existed
						{
							// update value							
							itSN->second.insert(curInt);
						}
						else // This feature has not existed in SN, directly put into SN	
						{
							// Put into SN							
							SN.insert({ cheInt, std::set<ObjWithoutCoord> {curInt} });
						}
					}
				}
			}
		}
		// clear all element for the next iterator
		fiveCells.clear();
		fiveCellInst.clear();
	}

}


/**
* @brief: This function generates star neighborhoods of features.
* @param: SN: a map of star neighbor instances
* @retval: a hash map that stores as <feat, <feat>>
*/
std::map<char, std::vector<std::set<char>>> genFeatNeibTrans(
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> & SN)
{
	std::map<char, std::vector<std::set<char>>> featNeibTrans;
	std::map<char, std::vector<std::set<char>>>::iterator itFeat;

	std::set<char> newValue;
	std::vector<std::set<char>> inner;

	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>::iterator itSN = SN.begin();
	while (itSN != SN.end())
	{
		itFeat = featNeibTrans.find(itSN->first.feature);
		if (itFeat != featNeibTrans.end()) // find out this key
		{
			// Build new value
			for (auto const& instnei : itSN->second)
			{
				newValue.insert(instnei.feature);
			}
			// Put into the key
			itFeat->second.push_back(newValue);
		}
		else
		{	// this key is not in the map, put new item			
			for (auto const& instnei : itSN->second)
			{
				newValue.insert(instnei.feature);
			}			
			inner.push_back(newValue);
			featNeibTrans.insert({ itSN->first.feature,  inner });
		}
		// Clear
		inner.clear();
		newValue.clear();
		// Terminal
		++itSN;
	}

	return featNeibTrans;
}



/**
* @brief: This function builds candidate pattern FP_tree of each feature type.
* @param: featNeibTrans: a map of neighbor transaction
*			prev_thres: prevalence threshold 
* @retval: a hash map that stores candidates
*/
std::map<char, std::map<std::string, float>> buildCandPattTrees(
	std::map<char, std::vector<std::set<char>>> featNeibTrans)
{
	std::map<char, std::map<std::string, float>> Treei;

	std::map<std::string, float> candvec;

	// Save each neighborhood set in transaction
	std::vector<Transaction> trans;		

	// Loop each item of featNeibTrans to buil FP-growth for each feature
	for (auto const& item : featNeibTrans)
	{				
		// Build FP-growth tree		
		Transaction td;
		// 2. Convert feature to itemset		
		for (auto const& rowI : item.second)
		{						
			// Loop each feature to create a transaction
			for (auto const& neiFeat : rowI)
			{
				Item nf{neiFeat};
				td.push_back(nf);				
			}
			// Put into transaction
			trans.push_back(td);		
			// Next iterator			
			td.clear();
		}			
		// 2. Create a FP-growth tree with the transaction
		FPTree fptree{trans, 0.0f };
		std::set<Pattern> patterns = fptree_growth(fptree);	
		// Add the first feature to the patterns to create a candidate		
		for (auto const& pat : patterns)
		{			
			// Put the first feature
			std::string onePat{item.first};
			// Put the other 
			for (auto const& feat : pat.first)
			{
				onePat += feat;
			}

			// Sort
			std::sort(onePat.begin(), onePat.end());
			
			// Call upper PIs
			float PI = (float)pat.second / (float)item.second.size();			

			// Put into the result
			candvec.insert({onePat, PI});
		}

		// Put the result
		Treei.insert({ item.first, candvec });
		// Go to the next feature
		trans.clear();
		candvec.clear();
	}

	return Treei;
}




/**
* @brief: This function gets candidate patterns if these patterns are in all feature trees.
* @param: Treei: a set of FP-Growth tree
* @retval: a set of candidates
*/
std::map<std::string, float> genCands(
	std::map<char, std::map<std::string, float>>  Treei, 
	float prev_thres)
{
	std::map<std::string, float> C; // Final candidate

	std::set<std::string> allC;

	// 1. Get all candidate
	std::map<char, std::map<std::string, float>>::iterator itTreei = Treei.begin();
	while (itTreei != Treei.end())
	{
		std::map<std::string, float>::iterator itinner = itTreei->second.begin();
		while (itinner != itTreei->second.end())
		{
			allC.insert(itinner->first);
			++itinner;
		}
		++itTreei;
	}
			
	// 2. Check a feature in a candidate in C must be in all FP-Growth tree	
	std::set<std::string>::iterator itallC = allC.begin();
	while (itallC != allC.end())
	{		
		// Check each candidate in C is in each FP-growth tree		
		std::vector<float> PRs; // Save upper PIs of the candidate in each feature of the tree

		for (auto const & feat : *itallC)
		{
			itTreei = Treei.find(feat);
			
			if (itTreei != Treei.end()) // Found out this feature of the candidate in the tree
			{
				// Find the current candidate is in the key of this feature
				std::map<std::string, float>::iterator itFindTree = itTreei->second.find(*itallC);
				if (itFindTree != itTreei->second.end())
				{					
					// Save the PI of this pattern
					PRs.push_back(itFindTree->second);					
				}
			}			
		}

		// Check this candidate is not in here
		if (PRs.size() == (*itallC).size()) // This is a true candidate
		{
			// Get the upper PI
			float PI = *std::min_element(std::begin(PRs), std::end(PRs));
			
			// Check if upper PI smaller than prev, delete it
			if (PI >= prev_thres)
			{
				C.insert({ *itallC, PI });
			}			
		}
		// Terminate
		++itallC;
	}
	return C;
}



/**
* @brief: This function gets all candidates with thier size equal to k
* @param: C: a set of candidates
*			k: size k candidate
* @retval: the set of size k candidates
*/
std::map<std::string, float> getSizeKCands(std::map<std::string, float> C, int k)
{
	std::map<std::string, float> Ck;

	for (auto const& cand : C)
	{
		if (cand.first.size() == k)
		{
			Ck.insert({ cand.first, cand.second });
		}
	}

	return Ck;
}


/**
* @brief: This function gets instances of size k candidates
* @param: Ck: a set of size k candidates
*			SN: the star neighborhoods
* @retval: a map that stores instances of size k candidates
*/
void findStarInsts(
	std::map<std::string, float> Ck,
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>& SN,
	std::map<std::string, std::vector<std::vector<std::vector<ObjWithoutCoord>>>>& SIk)
{
	for (auto const& cand : Ck)
	{
		int n = cand.first.size();
		// Save the table instance for one candidate
		std::vector<std::vector<std::vector<ObjWithoutCoord>>> Tk;
		
		// Loop element in SN to find star neighbors
		for (auto const& item : SN)
		{
			std::vector<std::vector<ObjWithoutCoord>> oneTk;

			// Only get item that feature of key is same as the first feature in cand
			if (item.first.feature == cand.first[0])
			{
				// Put the center star neighbor
				oneTk.push_back(std::vector<ObjWithoutCoord> {item.first});

				// Get instances in a star neighbor item
				std::map<char, std::vector<ObjWithoutCoord>> groupInst;
				std::map<char, std::vector<ObjWithoutCoord>>::iterator itgroupInst;
				for (auto const& inst : item.second)
				{
					// Only take the instances that their feature type is in candidate
					if (cand.first.find(inst.feature) != std::string::npos)
					{
						itgroupInst = groupInst.find(inst.feature);
						if (itgroupInst != groupInst.end()) // the fearture has already existed
						{
							itgroupInst->second.push_back(inst);
						}
						else
						{
							// Directly put into gp				
							groupInst.insert({ inst.feature, std::vector<ObjWithoutCoord> {inst} });
						}
					}
				}
				// Check to put into table instance				
				if (groupInst.size() == n - 1) // If the current star neighbor includes all features of the cand
				{
					for (auto const& its : groupInst)
					{
						oneTk.push_back(its.second);
					}
				}
			}
			// Gather to table instance
			if (oneTk.size() == n)
			{
				Tk.push_back(oneTk);
			}
		}
		// Put the table instance belongs to one candidate
		if (!Tk.empty())
		{
			SIk.insert({ cand.first, Tk });
		}
	}
}




/**
* @brief: This function generate cartesian product of vector<vector<int>>.
* @param: v: a vector<int>.
* @retval: resultCartesian a vector<vector<int>>.
*/
std::vector<std::vector<ObjWithoutCoord>> cartesianProduct(
	std::vector<std::vector<ObjWithoutCoord>> v)
{
	std::vector<std::vector<ObjWithoutCoord>> resultCartesian;

	auto product = [](long long a, std::vector<ObjWithoutCoord>& b)
	{
		return a * b.size();
	};

	const long long N = accumulate(v.begin(), v.end(), 1LL, product);

	std::vector<ObjWithoutCoord> u(v.size());

	for (long long n = 0; n < N; ++n)
	{
		lldiv_t q{ n, 0 };
		for (long long i = v.size() - 1; 0 <= i; --i) {
			q = div(q.quot, v[i].size());
			u[i] = v[i][q.rem];
		}

		resultCartesian.push_back(u);
	}

	return resultCartesian;
}




/**
* @brief: This function gets table instances of size-2 patterns
* @param: SIk: a set of size k candidates
*			Ck: a set of candidates
* @retval: table instances of size 2 patterns
*/
std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>> collectSize2TableInst(
	std::map<std::string, std::vector<std::vector<std::vector<ObjWithoutCoord>>>>& SIk,
	int k,
	std::map<std::string, float>& Ck)
{
	// Save size 2 table instances
	std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>> T2;
	// Get table instances
	for (auto const& item : SIk)
	{
		std::vector<std::vector<ObjWithoutCoord>> oneT2;
		for (auto const& rows : item.second)
		{
			std::vector<std::vector<ObjWithoutCoord>> innerValue = cartesianProduct(rows);
			oneT2.insert(oneT2.end(), innerValue.begin(), innerValue.end());
		}
		// Save 
		T2.insert({item.first, oneT2});
	}
	// Terminate
	return T2;
}




/**
* @brief: This function find real row instances of candidate patterns.
* @param: SIKplus: a hash map of table instances of size k+1 candidate patterns.
*         CIk: a hash map of table instances of size k candiate patterns.
* @retval: tabInstkplus: a hash map of real table instances of candidate patterns.
*/
std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>> filterCliqueInstances2(
	std::map<std::string, std::vector<std::vector<std::vector<ObjWithoutCoord>>>>& SIkplus,
	std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>> CIk)
{
	// Save the size (k+1) table instances
	std::map <std::string, std::vector<std::vector<ObjWithoutCoord>>> tabInstkplus;

	// Save the sub pattern of pattern of size (k+1). Example (A, B, C) -> patternkSub = (B, C)
	std::string patternkSub;
	std::vector<std::vector<ObjWithoutCoord>> innerValue, realTableInstkplus;

	std::vector<ObjWithoutCoord> checkRowkPlus;

	// Loop each element of CIkplus
	for (auto const& item : SIkplus)
	{
		// Get sub pattern (execpt the first element)
		patternkSub.insert(patternkSub.end(), item.first.begin() + 1, item.first.end());
		// Check each row instance in coarseTabInstkPlus
		if (CIk.find(patternkSub) != CIk.end())
		{
			for (auto const& rows : item.second)
			{
				// Product coarse row instances			
				innerValue = cartesianProduct(rows);
				// Check real row instances				
				for (auto const& rowInstKPlus : innerValue)
				{
					checkRowkPlus.insert(checkRowkPlus.end(), rowInstKPlus.begin() + 1, rowInstKPlus.end());
					if (std::find(CIk.find(patternkSub)->second.begin(),
						CIk.find(patternkSub)->second.end(), checkRowkPlus)
						!= CIk.find(patternkSub)->second.end())
					{
						// The current row instance is a clique
						realTableInstkplus.push_back(rowInstKPlus);
					}
					checkRowkPlus.clear();
				}
			}

			// Save the table instance of the current candidate
			if (realTableInstkplus.size() > 0)
			{
				// Put into the final result
				tabInstkplus.insert({ item.first, realTableInstkplus });
				realTableInstkplus.clear();
			}
		}

		// Clear all
		patternkSub.clear();

	}
	return tabInstkplus;
}



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
	std::map < std::string, float> & nonPk)
{
	std::vector<float> PRs;
	float PI;
	
	std::vector<std::set<ObjWithoutCoord>> uniqueTableInstance;
	int numberInst, totalInst;

	for (std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>>::iterator iter = CIkplus.begin();
		iter != CIkplus.end(); ++iter)
	{
		// Retrieve patterns and table instances
		std::string pattern = iter->first;
		// Loop each row instance in tableInstance
		for (auto const& rowInst : iter->second)
		{
			if (uniqueTableInstance.empty())
			{
				for (int i = 0; i < k; ++i)
				{
					std::set<ObjWithoutCoord> temp;
					temp.insert(rowInst[i]);
					uniqueTableInstance.push_back(temp);
				}
			}
			else
			{
				for (int i = 0; i < k; ++i)
				{
					uniqueTableInstance[i].insert(rowInst[i]);
				}
			}
		}
		// Loop each element
		for (int i = 0; i < k; ++i)
		{
			totalInst = totalInstEachFeat.find(pattern[i])->second;
			//std::cout << "Total number instances: " << totalInst << endl;
			numberInst = uniqueTableInstance[i].size();
			//std::cout << "Number of instances: " << numberInst<<endl;

			PRs.push_back((float)numberInst / (float)totalInst);
		}

		// find the minimum elemnet
		PI = *std::min_element(std::begin(PRs), std::end(PRs));
		// Set the number of float point
		PI = Precision(PI, 3);
		
		// Check prevalent patterns
		if (PI >= prev_thres)
		{
			Pk.insert({ pattern, PI });
		}
		else
		{
			nonPk.insert({ pattern, PI });
		}

		// Clear all temporary varibles
		uniqueTableInstance.clear();
		PRs.clear();
		pattern.clear();
	}
}


/**
* @brief: This function adds new prevalent patterns into the result
* @param: Pk : the size k prevalent patterns
*			Rk: all prevalent patterns
* @retval: Nope
*/
void addPrePatts(
	std::map<std::string, float, cmpMapBySizePatt>& Pk, 
	std::map<std::string, float, cmpMapBySizePatt>& Rk)
{
	for (auto const& patt : Pk)
	{
		Rk.insert({ patt.first, patt.second });
	}
}



/**
* @brief: This function removes checked candidates from the set of candidates
* @param: C : a set of all orignial candidates
*		Ck: a set of checked size k candidates
* @retval: Nope
*/
void delCkFromC(
	std::map<std::string, float> & C, 
	std::map<std::string, float> & Ck)
{
	std::map<std::string, float>::iterator itC;
	for (auto const & checkedCand : Ck)
	{
		itC = C.begin();
		while (itC != C.end())
		{
			if (checkedCand.first == itC->first) // Find out, delete
			{
				itC = C.erase(itC); // Delete from C
				break;
			}
			else
			{
				++itC;
			}			
		}
	}
}



/**
*@brief: This function removes super candidates from the set of original candidates by a non-prevalent pattern
* @param : nonPk : a set of prevalent maximal patterns
* C : a set of orignial candidates
* @retval : C : a set of prunned candidates
*/
void superPrun(
	std::map<std::string, float>& nonPk,
	std::map<std::string, float>& C)
{
	for (auto const& patt : nonPk)
	{
		std::map<std::string, float>::iterator itC = C.begin();
		while (itC != C.end())
		{
			if (std::includes(itC->first.begin(), itC->first.end(), 
				patt.first.begin(), patt.first.end()))
			{
				itC = C.erase(itC);
			}
			else
			{
				++itC;
			}
		}
	}
}




/**
* @brief: This function select coarse prevalent patterns.
* @param: SIK: a hash map of table instances
* @retval: CIk: a hash map of table instances which satify coarse min_prev
*/
void selectCoarsePrevalentPatterns(
	std::map<std::string, std::vector<std::vector<std::vector<ObjWithoutCoord>>>>& SIkplus,
	std::map<char, int>& totalInstEachFeat,
	float prev_thres,
	int k)
{
	std::map<std::string, std::vector<std::vector<std::vector<ObjWithoutCoord>>>>::iterator iter = SIkplus.begin();
	while (iter != SIkplus.end())
	{
		// Loop each row instance to get tableInstance
		std::vector<std::set<ObjWithoutCoord>> uniqueTableInstance(k);
		for (auto const& rows : iter->second)
		{
			for (size_t t = 0; t < k; ++t)
			{
				uniqueTableInstance[t].insert(rows[t].begin(), rows[t].end());
			}
		}

		// Calculate PRs
		std::vector<float> PRs;
		for (size_t t = 0; t < k; ++t)
		{
			int totalInst = totalInstEachFeat.find(iter->first[t])->second;
			int numberInst = uniqueTableInstance[t].size();
			PRs.push_back((float)numberInst / (float)totalInst);
		}

		// Find the minimum elemnet
		float PI = *std::min_element(std::begin(PRs), std::end(PRs));
		//std::cout << " : " << PI << endl;

		// Check prevalent patterns
		if (PI < prev_thres)
		{
			iter = SIkplus.erase(iter);
		}
		else
		{
			++iter;
		}
	}
}



/**
* @brief: This function selects closed patterns
* @param: Pk: a hash map of size k prevalent co-location patterns
* @retval: CloseP: a hash map of closed patterns
*/
void filterClosedPatts(std::map<std::string, float, cmpMapBySizePatt>& Pk,
	std::map<float, std::vector<std::string>>& CloseP)
{
	// 7. Filter closed patterns
	// Since Rk is sort by sizes, the first element is allway the smallest	
	std::map<float, std::vector<std::string>>::iterator itCloseP;
	bool issuper;
	int s;
	std::map<std::string, float, cmpMapBySizePatt>::iterator itPk = Pk.begin();
	std::vector<std::string>::iterator itchecksuperset;

	while (itPk != Pk.end())
	{
		// Check this pattrn in CloseP
		itCloseP = CloseP.find(itPk->second);
		if (itCloseP == CloseP.end()) // not find out this key
		{			
			// Put into the result directly
			//std::vector<std::string> key{ itPk->first };			
			CloseP.insert({ itPk->second, std::vector<std::string> { itPk->first } });
		}
		else // closeP has this pattern, check it suppersets
		{
			itchecksuperset = itCloseP->second.begin();
			while (itchecksuperset != itCloseP->second.end())
			{
				// check superset relationship
				issuper = std::includes(itPk->first.begin(), itPk->first.end(),
					(*itchecksuperset).begin(), (*itchecksuperset).end());
				if (issuper)
				{
					itchecksuperset = itCloseP->second.erase(itchecksuperset);
				}
				else
				{
					++itchecksuperset;
				}
			}
			// after check superset, put the current pattern into the result
			itCloseP->second.push_back(itPk->first);			
		}	
		// Terminate
		++itPk;		
	}
}



