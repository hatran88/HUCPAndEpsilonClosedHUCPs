
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <ctime>
#include "boost/dynamic_bitset.hpp"
#include <unordered_set>
#include <numeric>      // std::accumulate


#include "Object.h"
#include "ReaderCSV.h"
#include "PrintFunctions.h"
#include "SetNumberPointFloat.h"
#include "Patterns.h"





using namespace std;

// Definition global variables
std::vector<ObjWithoutCoordinate> verticesVec; // the instances in a block
std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc> neiOneBlock; // neighbor instances in a block

// Definition global variables
clock_t timeConstCoLHashmap;
double totalTimeConstCoLHashmap = 0.0f;
clock_t timeFilterHUCPs;
double totalTimeFilHUCPs = 0.0f;

std::unordered_map<std::string, std::unordered_map<char, std::unordered_set<ObjWithoutCoordinate, myHashFunc>>> CoLHM;
std::unordered_map<std::string, std::unordered_map<char, std::unordered_set <ObjWithoutCoordinate, myHashFunc>>>::iterator it;

// Save the weight of each feature
std::unordered_map<char, float> uf;


/*
* @brief: This function count the number of instances of each feature.
* @param: dataList: an input dataset;
* @retval: wf: (a global variable) a hash map that stores each feature name and its instance number.
*/
void getFeatureWeight(std::vector<std::vector<std::string>> dataList)
{
	std::unordered_map<char, float>::iterator iw;

	for (auto const& vec : dataList)
	{	
		//std::cout << vec[0][0] << ", " << std::stof(vec[4]) << endl;
		if (uf.empty())
		{
			uf.insert({ vec[0][0], std::stof(vec[4]) });
		}
		else
		{
			iw = uf.find(vec[0][0]);
			if (iw != uf.end()) // if the current feature exists
			{
				iw->second = iw->second + std::stof(vec[4]);
			}
			else // the current feature does not exist
			{
				uf.insert({ vec[0][0], std::stof(vec[4]) });
			}
		}		
	}	
}


/**
* @brief This function make a grid on an input dataset
* @param dataList: an input dataset; dist_thres: a distance threshold
* @retval a hash map that stores cell id with instances fall in it.
*/
std::map<std::pair<int, int>, std::vector<Objects>> makeGrid(
	std::vector<std::vector<std::string>> & dataList, float dist_thres)
{
	std::map<std::pair<int, int>, std::vector<Objects>> grid;

	int cell_x, cell_y;
	pair<int, int> cell_key;

	std::vector<Objects> value;

	for (auto const& vec : dataList)
	{
		// calc the cell id
		cell_x = ceil(std::stof(vec[2]) / dist_thres); // x coordinate
		cell_y = ceil(std::stof(vec[3]) / dist_thres); // y coordinate		
		cell_key = make_pair(cell_x, cell_y);

		Objects instance = { vec[0][0], std::stoi(vec[1]), std::stof(vec[4]), std::stof(vec[2]), std::stof(vec[3])};
		//Objects instance = { vec[0][0], std::stoi(vec[1]), std::stof(vec[4]), std::stof(vec[2]), std::stof(vec[3]) };
		value.insert(value.end(), instance);

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
			grid.find(cell_key)->second.push_back(instance);
		}
		// Clear
		value.clear();
	}

	return grid;
}


/**
* @bref This function calculates the distances of instances in the current block.
* @param alll instance in the current block
* @retval A vector that save the distance of instances.
*/
float calculateDistanceTwoInstances(Objects currentInst, Objects checkInst)
{
	return sqrt((currentInst.x - checkInst.x)*(currentInst.x - checkInst.x)
		+ (currentInst.y - checkInst.y)*(currentInst.y - checkInst.y));
}


/**
* @brief: This function generates star neighborhoods of instances.
* @param: grid: a grid posing of the input dataset
*         dist_thres: a distance threshold
* @retval: SN: a hash map that stores as <instance, <<neighbors>,<neighbors>>. This struture is different with star neighbors in Join-less
*/
void genStarNeighborhoods(
	std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc> & neighbors,
	std::unordered_map<ObjWithoutCoordinate, int, myHashFunc>& vertices,
	std::map<std::pair<int, int>, std::vector<Objects>>& grid,
	float dist_thres)
{
	std::unordered_map<ObjWithoutCoordinate, std::unordered_set<ObjWithoutCoordinate, myHashFunc>, myHashFunc> tempNei;

	int i, j; // the index of cells in x and y
	std::vector<std::pair<int, int>> fiveCells;
	float dist;
	std::vector<Objects> fiveCellInst; // save all instance in the five cells and two neighboring instances.	
	// For get index of elements in vectors
	std::vector<Objects>::iterator itvert;
	unsigned int indexCurrentInst, indexCheckInst;

	// Loop each cell in grid and computation neighbors in five cells
	for (std::map<std::pair<int, int>, std::vector<Objects>>::iterator it = grid.begin(); it != grid.end(); ++it)
	{
		// Create the five cells
		i = it->first.first;
		j = it->first.second;

		fiveCells.push_back(it->first);
		fiveCells.push_back(std::make_pair(i, j + 1));
		fiveCells.push_back(std::make_pair(i + 1, j + 1));
		fiveCells.push_back(std::make_pair(i + 1, j));
		fiveCells.push_back(std::make_pair(i + 1, j - 1));

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
		for (Objects const& currentInst : it->second)
		{
			ObjWithoutCoordinate currentInstNoCoor{ currentInst.feature, currentInst.instance, currentInst.utility };
			// Save all instances
			vertices.insert({ currentInstNoCoor, 1 });

			// Put thi point in the neighbor first
			//indexCurrentInst = lower_bound(vertices.begin(), vertices.end(), currentInst) - vertices.begin();
			if (tempNei.find(currentInstNoCoor) == tempNei.end())
			{
				tempNei.insert({ currentInstNoCoor, std::unordered_set<ObjWithoutCoordinate, myHashFunc> {} });
			}

			for (Objects const& checkInst : fiveCellInst) // check with each instance in the five cells
			{
				if (currentInst.feature != checkInst.feature) // only check two instances belong to different features.
				{
					dist = calculateDistanceTwoInstances(currentInst, checkInst);

					if (dist <= dist_thres) // the two instances have neighbor relationship
					{
						// Find index of instances						
						//indexCheckInst = lower_bound(vertices.begin(), vertices.end(), checkInst) - vertices.begin();
						ObjWithoutCoordinate checkInstNoCoor{ checkInst.feature, checkInst.instance, checkInst.utility };

						// Put the current instance into neighbors
						if (tempNei.find(currentInstNoCoor) != tempNei.end()) // this instance has already existed
						{
							// update value							
							tempNei.find(currentInstNoCoor)->second.insert(checkInstNoCoor);
						}
						else // This feature has not existed in SN, directly put into SN	
						{
							// Put into SN							
							tempNei.insert({ currentInstNoCoor, std::unordered_set<ObjWithoutCoordinate, myHashFunc>{checkInstNoCoor} });
						}

						// Put the current instance into neighbors
						if (tempNei.find(checkInstNoCoor) != tempNei.end()) // this instance has already existed
						{
							// Update value							
							tempNei.find(checkInstNoCoor)->second.insert(currentInstNoCoor);
						}
						else // This feature has not existed in SN, directly put into SN	
						{
							tempNei.insert({ checkInstNoCoor, std::unordered_set<ObjWithoutCoordinate, myHashFunc>{currentInstNoCoor} });
						}
					}
				}
			}
		}
		// clear all element for the next iterator
		fiveCells.clear();
		fiveCellInst.clear();
	}

	// Convert to neighbors (using vector)
	std::unordered_map<ObjWithoutCoordinate, std::unordered_set<ObjWithoutCoordinate, myHashFunc>, myHashFunc>::iterator itN = tempNei.begin();
	while (itN != tempNei.end())
	{
		std::vector<ObjWithoutCoordinate> value(itN->second.begin(), itN->second.end());
		std::sort(value.begin(), value.end());
		
		neighbors.insert({ itN->first, value });

		itN = tempNei.erase(itN);
	}

}



/**
* @brief: This function gets the neighbors of instances in the current block
* @param: neighbors: neighboring instances
* @retval: size_t: the maximal degeree
*/
void getNeiOneBlock(std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc>& neighbors)
{
	std::vector<ObjWithoutCoordinate> commonInst;

	std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc>::iterator itneighbors;
	for (auto const& inst : verticesVec)
	{
		itneighbors = neighbors.find(inst);
		if (itneighbors != neighbors.end())
		{
			// Get the intersection elements
			std::set_intersection(verticesVec.begin(), verticesVec.end(), 
				itneighbors->second.begin(), itneighbors->second.end(),
				std::back_inserter(commonInst));
			// Put into the result
			neiOneBlock.insert({inst, commonInst});
			// Clear
			commonInst.clear();
		}
	}
}



/**
* @brief: This function sorts vertices by their degeneracy.
* @param: vertices: a set of indexes of all vertices
*         neighbors: star neighbors of all vertices
* @retval: a vector of oredered vertices sorted by degeneracy
*/
void OrderedDegeneracy(std::vector<ObjWithoutCoordinate>& orderDeg,
	std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc>& neiOneBlock)
{
	// Store vertices and their neighbor number in to a map
	std::map<ObjWithoutCoordinate, int> vertDeg;

	std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc>::iterator itNei = neiOneBlock.begin();
	while (itNei != neiOneBlock.end())
	{
		vertDeg.insert({ itNei->first, itNei->second.size() });
		++itNei;
	}

	std::map<ObjWithoutCoordinate, int>::iterator itVertDeg;
	// Loop each vertices and calc deg
	while (!vertDeg.empty())
	{
		// Assum the first item is the minimum
		ObjWithoutCoordinate minV{ vertDeg.begin()->first };
		int minDeg = vertDeg.begin()->second;

		// Loop to find the real minimum
		itVertDeg = vertDeg.begin();
		// start index number 1
		++itVertDeg;

		while (itVertDeg != vertDeg.end())
		{
			if (itVertDeg->second < minDeg)
			{
				minV = itVertDeg->first;
				minDeg = itVertDeg->second;
			}
			++itVertDeg;
		}
		//Put the vertex to the result
		orderDeg.push_back(minV);
		// Delete this minV from vertDeg
		vertDeg.erase(minV);
		// Minus one deg form vertices which are neighbor of minV
		for (auto const& item : neiOneBlock.find(minV)->second)
		{
			// get value of item and minus 1
			itVertDeg = vertDeg.find(item);
			if (itVertDeg != vertDeg.end())
			{
				itVertDeg->second--;
			}
		}
	}
}


/**
* @brief: This function get maximal cliques
* @param: vertices: a set of indexes of all vertices
*         neighbors: star neighbors of all vertices
* @retval: a vector of oredered vertices sorted by degeneracy
*/
void BronKerboschDeg(std::vector<ObjWithoutCoordinate> P,
	std::vector<ObjWithoutCoordinate> X,
	std::vector<ObjWithoutCoordinate>& orderDeg)
{
	for (auto const& v : orderDeg)
	{
		std::vector<ObjWithoutCoordinate> newP, newX, newR;

		newR.push_back(v);

		std::set_intersection(P.begin(), P.end(),
			neiOneBlock.find(v)->second.begin(), neiOneBlock.find(v)->second.end(),
			std::back_inserter(newP));

		std::set_intersection(X.begin(), X.end(),
			neiOneBlock.find(v)->second.begin(), neiOneBlock.find(v)->second.end(),
			std::back_inserter(newX));		

		BronKerboschPivot(newP, newR, newX);

		// Delete
		P.erase(std::remove(P.begin(), P.end(), v), P.end());
		//std::sort(P.begin(), P.end());
		// and add v to X
		X.push_back(v);
		std::sort(X.begin(), X.end());
	}
}



/**
* @brief: This function enumerates all maximal cliques based on bron-kerbosch pivot.
* @param: vertices: a set of all vertices (instances)
*         neighbors: a hashmap stores neighbors of each instance
* @retval: allMaxCliques: all maximal cliques.
*/
void BronKerboschPivot(std::vector<ObjWithoutCoordinate> P, std::vector<ObjWithoutCoordinate> R, std::vector<ObjWithoutCoordinate> X)
{
	if (P.empty() && X.empty())
	{
		if (R.size() > 1)
		{	
			/*std::cout << "get one : ";
			for (auto const& inst : R)
			{
				std::cout << inst.feature << "." << inst.instance << ", ";
			}
			std::cout << endl;*/
			// Get one maximal cliques and construct hashmap
			timeConstCoLHashmap = clock();
			constructCoLoHashmap(R);
			totalTimeConstCoLHashmap = totalTimeConstCoLHashmap + double(clock() - timeConstCoLHashmap);
		}
	}
	else
	{
		// Find the pivot		
		std::vector<ObjWithoutCoordinate> unionPX;
		std::set_union(P.begin(), P.end(), X.begin(), X.end(), std::back_inserter(unionPX));
		
		// Find maximum
		std::vector<ObjWithoutCoordinate> intersecPNu;
		std::vector<ObjWithoutCoordinate> maxu{ unionPX[0] };
		int maxnei = neiOneBlock.find(unionPX[0])->second.size();

		for (size_t j = 1; j < unionPX.size(); ++j)
		{
			std::set_intersection(P.begin(), P.end(),
				neiOneBlock.find(unionPX[j])->second.begin(), neiOneBlock.find(unionPX[j])->second.end(),
				std::back_inserter(intersecPNu));

			if (intersecPNu.size() > maxnei)
			{
				maxnei = intersecPNu.size();
				maxu.clear();
				maxu.push_back(unionPX[j]);
			}
		}

		std::vector<ObjWithoutCoordinate> deffPNu;
		std::set_difference(P.begin(), P.end(), neiOneBlock.find(maxu[0])->second.begin(), neiOneBlock.find(maxu[0])->second.end(),
			std::back_inserter(deffPNu));

		// Loop and recursion
		int psize = deffPNu.size();

		for (size_t i = 0; i < psize; ++i)
		{
			ObjWithoutCoordinate v = deffPNu[i];
			std::vector<ObjWithoutCoordinate> newP, newX;

			std::set_intersection(P.begin(), P.end(),
					neiOneBlock.find(v)->second.begin(), neiOneBlock.find(v)->second.end(),
					back_inserter(newP));
				std::sort(newP.begin(), newP.end());
			
			std::set_intersection(X.begin(), X.end(),
				neiOneBlock.find(v)->second.begin(), neiOneBlock.find(v)->second.end(),
				back_inserter(newX));
			std::sort(newX.begin(), newX.end());			

			std::vector<ObjWithoutCoordinate> newR = R;
			newR.push_back(v);
			std::sort(newR.begin(), newR.end());

			BronKerboschPivot(newP, newR, newX);

			// Delete v from P
			P.erase(std::remove(P.begin(), P.end(), v), P.end());
			// and add v to X
			X.push_back(v);
			std::sort(X.begin(), X.end());
		}
	}
}



/**
* @brief: This function buids the co-location hashmap structure.
* @param: allMaxCliques: all maximal cliques.
* @retval: CoLoHM: a co-location hashmap.
*/
void constructCoLoHashmap(std::vector<ObjWithoutCoordinate> R)
{	
	// Temporary varibales
	std::string pattern;

	// Loop each instance in the maximal clique to build key
	for (auto const& inst : R)
	{
		pattern += inst.feature;			
	}
	// Check this key has not already in CoLHM
	it = CoLHM.find(pattern);
	if (it != CoLHM.end()) // This key has existed, update old values
	{
		// update other instances in the clique
		for (auto const& inst : R)
		{
			it->second.find(inst.feature)->second.insert(inst);
		}
	}
	else // not exist
	{
		// build value and add new item
		std::unordered_map<char, std::unordered_set<ObjWithoutCoordinate, myHashFunc>> outValue;
		for (auto const& inst : R) {

			outValue.insert({ inst.feature, std::unordered_set<ObjWithoutCoordinate, myHashFunc> {inst} });
		}

		CoLHM.insert({ pattern, outValue });
	}
}




/*
*@brief: This function generate all combination of a char vector
*@param: c: a vector of charl; combo: the size of combo; C(n, m) = n!/(m!(n-m)!)
*@retval: a vector of sub vectors of the vector
*/
template<typename T>
std::string getCombination(const T& c, int combo)
{
	std::string result;
	int n = c.size();
	for (int i = 0; i < n; ++i) {
		if ((combo >> i) & 1)
			result.push_back(c[i]);
	}
	return result;
}


template<typename T>
std::vector<std::string> combo(const T& c, int k)
{
	std::vector<std::string> combination;

	int n = c.size();
	int combo = (1 << k) - 1;       // k bit sets

	while (combo < 1 << n)
	{
		combination.insert(combination.end(), getCombination(c, combo));
		int x = combo & -combo;
		int y = combo + x;
		int z = (combo & ~y);
		combo = z / x;
		combo >>= 1;
		combo |= y;
	}

	return combination;
}



/*
*@brief: This function calculate the pi of the current pattern
*@param: pattern: the current pattern which need to calculate its pi.
*        superPatterns:
*        cpHashMap: a co-location pattern hash map
*@retval: PI: the participation index of the pattern
*         superPatterns: all super patterns of the current pattern
*/
float calculateWeightWithSuperPatterns2Phases(std::string pattern, int sizeofPatttern)
{
	// 1. Get instances of the super patterns
	//std::cout << pattern << ": ";
	std::vector<std::unordered_set<ObjWithoutCoordinate, myHashFunc>> numInstPatt(sizeofPatttern);
	bool isInHashMap;		
	for (it = CoLHM.begin(); it != CoLHM.end(); ++it)
	{
		// Check if the key in cpHashMap is super set of the current pattern
		isInHashMap = std::includes(it->first.begin(), it->first.end(), pattern.begin(), pattern.end()); // method 1: use std		
		// isInHashMap = a & (~b) ;// method 3: use bitwise difference
		if (isInHashMap)
		{
			// Get Instances
			for (int t = 0; t < sizeofPatttern; t++)
			{
				numInstPatt[t].insert(it->second.find(pattern[t])->second.begin(), it->second.find(pattern[t])->second.end());
			}
		}
	}

	// 2. Call the weight of the pattern
	std::vector<float> wc(sizeofPatttern, 0.0f);
	
	for (int i = 0; i < sizeofPatttern; i++)
	{
		for (auto const& insta : numInstPatt[i])
		{
			wc[i] = wc[i] + insta.utility;
			//show instance
			//std::cout << insta.feature << "." << insta.instance << ", ";
		}		
		wc[i] = wc[i] / (float)uf.find(pattern[i])->second;
		//std::cout << pattern[i] << " : w[" << i << "] = " << wc[i] << endl;
		//std::cout << endl;
	}
	//std::cout << "sum: " << std::accumulate(wc.begin(), wc.end(), 0.0f) / sizeofPatttern << endl;
	float wpat = Precision(std::accumulate(wc.begin(), wc.end(), 0.0f)/ sizeofPatttern, 3);
	
	return wpat;
}


/**
* @brief: This function call PIs and filter prevalent patterns.
* @param: CoLHM: the co-location hash map.
*         pre_thress: a minimum prevalent threshold.
* @retval: allPrevPats: a co-location hashmap.
*/
void calPIby2Phases(std::map<std::string, float>& allPats)
{
	// Phase 1. Get all keyset from CHash and calculate the PI
	// 1.1 Get all keys
	std::vector<std::string> candPatts;
	for (it = CoLHM.begin(); it != CoLHM.end(); ++it)
	{		
		candPatts.push_back(it->first);
	}
	//1.2 Calculate the weight of patterns
	float wc;
	for (auto oneCand : candPatts)
	{
		wc = calculateWeightWithSuperPatterns2Phases(oneCand,oneCand.size());
		//std::cout << "PI : " << PI << endl;
		allPats.insert({ oneCand, wc });
	}

	// Phase 2. Calculate the sub-patterns
	// 2.1 Sort by size	
	std::sort(candPatts.begin(), candPatts.end(), []
	(const std::string& first, const std::string& second) {
			return first.size() > second.size();
		});
	// 2.2 Process each pattern from large to small	
	std::map<string, float>::iterator itAllPatts;
		
	// Checking 
	std::vector<std::string>::iterator itCands = candPatts.begin();
	while (itCands != candPatts.end())
	{
		if ((*itCands).size() > 2)
		{
			// 2.3 Gen direct sub-pattern
			std::vector<std::string> subPatts = combo(*itCands, (*itCands).size() - 1);
			// 2.4 Check this subPatts has allready calcualated
			for (auto const& subPatt : subPatts)
			{
				itAllPatts = allPats.find(subPatt);
				if (itAllPatts == allPats.end()) // this candidate has not calculated
				{
					wc = calculateWeightWithSuperPatterns2Phases(subPatt, subPatt.size());					
					allPats.insert({ subPatt, wc });
				}
			}
		}
		// Terminate
		++itCands;
	}
}




/**
* @brief: This function filters prevalent patterns
* @param: CoLHM: the co-location hash map.
*         pre_thress: a minimum prevalent threshold.
* @retval: allPrevPats: a co-location hashmap.
*/
void filterPrevPatts(
	std::map<string, float>& allPrevPats,
	std::map<std::string, float>& allPats,
	float weight_thres)
{
	std::map<std::string, float>::iterator itAllPatts = allPats.begin();
	while (itAllPatts != allPats.end())
	{
		//std::cout << itAllPatts->first << " : " << itAllPatts->second << endl;
		if (itAllPatts->second >= weight_thres)
		{
			allPrevPats.insert({ itAllPatts->first, itAllPatts->second });
		}

		++itAllPatts;
	}
}


// Comparator function to sort pairs
// according to second value
bool cmp(std::pair<std::string, float>& a,
	std::pair<std::string, float>& b)
{
	return a.first.size() <= b.first.size();
}


// to value in a (key-value) pairs
void sorthashmap(std::unordered_map<std::string, float>& M,
	std::vector<std::pair<std::string, float>>& A)
{
	// Copy key-value pair from Map
	// to vector of pairs
	for (auto& it : M) {
		A.push_back(it);
	}
	// Sort using comparator function
	sort(A.begin(), A.end(), cmp);
}


/**
* @brief: This function classify patterns by sizes
* @param: allPrevPats: all patterns
*         classPattBySize: classified pattens by sizes
* @retval:
*/
void classifyPattsBySize(std::unordered_map<string, float>& allPrevPats, std::map<int, std::vector<Patterns>>& classPattBySize)
{
	int sizeP;
	std::map<int, std::vector<Patterns>>::iterator itc;
	for (auto iter = allPrevPats.begin(); iter != allPrevPats.end(); ++iter)
	{
		Patterns onep{ iter->first, iter->second };

		sizeP = iter->first.size();
		itc = classPattBySize.find(sizeP);
		if (itc != classPattBySize.end())
		{
			itc->second.push_back(onep);
		}
		else
		{
			classPattBySize.insert({ sizeP, std::vector<Patterns>{onep} });
		}
	}
}



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
void filterHPCPs(std::unordered_map<std::string, float>& allCPCP,
	float utility_thres, float alpha, float beta)
{
	// Phase 1. Get all keyset from CHash and calculate the PI
	// 1.1 Get all keys
	std::map<std::string, int> candPatts;

	std::map<std::string, int>::iterator itcand;

	for (it = CoLHM.begin(); it != CoLHM.end(); ++it)
	{
		//candPatts.insert({ itCoLHM->first, 1 });
		candPatts.insert({ it->first, 1 });
	}

	//1.3 Get the largest size of candidate
	int l = findMaxSizeCand(candPatts);
	int numC = candPatts.size();
	float UPI;
	std::map<std::string, int>::iterator itcan;

	// 1.4 Loop each cand
	while (l >= 2 || numC > 0)
	{
		std::vector<std::string> Cl; // Save size l candidates
		std::vector<std::string> remainCand; // Save the delete cand

		// Find size l candidates
		itcan = candPatts.begin();
		while (itcan != candPatts.end())
		{
			if (itcan->first.size() == l)
			{
				// Save it
				Cl.push_back(itcan->first);
				// Delete it
				itcan = candPatts.erase(itcan);
			}
			else
			{
				++itcan;
			}
		}
		// 3. Calculate the UPI of each size l candidate		
		for (auto oneLCand : Cl)
		{
			UPI = calculatePIWithSuperPatterns(oneLCand, oneLCand.size(), alpha, beta);
			
			// filter for HUCPs
			timeFilterHUCPs = clock();
			if (UPI >= utility_thres)
			{
				// Save all PCPS
				allCPCP.insert({ oneLCand, UPI });				
			}
			// Get the direct subset of the current candidate
			if (oneLCand.size() > 2)			{
				
				std::vector<std::string> subPatts = combo(oneLCand, oneLCand.size() - 1);
				// check these candidates are not in candPatts
				for (auto const& subP : subPatts)
				{		
					itcand = candPatts.find(subP);
					if (itcand == candPatts.end()) // this candidate has not calculated
					{
						candPatts.insert({ subP, 1 });
					}
				}
			}
			totalTimeFilHUCPs = totalTimeFilHUCPs + double(clock() - timeFilterHUCPs);
		}

		// Next iterator		
		numC = candPatts.size();
		l = findMaxSizeCand(candPatts);
	}
}


/*
*@brief: This function filters delta PCPS
*@param: deltaCPCP: a set of all delta CPCP
*        oneLCand: the current pattern
*        delta_thres
*@retval: none
*/
bool findDeltaPats(std::unordered_map<std::string, float>& deltaCPCP,
	std::string oneLCand, float UPI, float delta_thres)
{
	// Check the direct superset
	bool flagDel = true; // Assume itSortPats is a meta patterns
	std::unordered_map<std::string, float>::iterator itdelta = deltaCPCP.begin();
	while (itdelta != deltaCPCP.end())
	{
		if (itdelta->first.size() > oneLCand.size())
		{
			bool checkinclude = std::includes(itdelta->first.begin(), itdelta->first.end(),
				oneLCand.begin(), oneLCand.end());

			if (checkinclude)
			{
				// check upi difference				
				if (Precision(abs(itdelta->second - UPI), 3) <= delta_thres)
				{
					// itSortPats is needed to deleted
					flagDel = false;

					break;
				}
			}
		}

		++itdelta;
	}

	return flagDel;
}



/*
*@brief: This function filters closed PCPS
*@param:
*@retval: none
*/
bool findClosedPCPs(std::unordered_map<std::string, float>& closedPCP,
	std::string oneLCand, float UPI)
{
	// Check the direct superset
	bool flagDel = true; // Assume itSortPats is a meta patterns
	std::unordered_map<std::string, float>::iterator itdelta = closedPCP.begin();
	while (itdelta != closedPCP.end())
	{
		if (itdelta->first.size() > oneLCand.size())
		{
			bool checkinclude = std::includes(itdelta->first.begin(), itdelta->first.end(),
				oneLCand.begin(), oneLCand.end());

			if (checkinclude)
			{
				// check upi difference				
				if (itdelta->second == UPI)
				{
					// itSortPats is needed to deleted
					flagDel = false;

					break;
				}
			}
		}

		++itdelta;
	}

	return flagDel;
}




/**
* @brief: This function filters prevalent patterns
* @param: CoLHM: the co-location hash map.
*         pre_thress: a minimum prevalent threshold.
* @retval: allPrevPats: a co-location hashmap.
*/
int findMaxSizeCand(std::map<std::string, int>& candPatts)
{
	std::map<std::string, int>::iterator it = candPatts.begin();
	int t = 1;
	while (it != candPatts.end())
	{
		if (it->first.size() > t)
		{
			t = it->first.size();
		}

		++it;
	}
	return t;
}



/*
*@brief: This function calculates UPI of each pattern
*@param: pattern: the current pattern which need to calculate its pi.
*        superPatterns:
*        cpHashMap: a co-location pattern hash map
*			uf: utility of feature
*@retval: UPI
*         superPatterns: all super patterns of the current pattern
*/
float calculatePIWithSuperPatterns(std::string pattern, int sizeofPatttern, float alpha, float beta)
{
	// 1. Get all participating instances
	std::vector<std::unordered_set<ObjWithoutCoordinate, myHashFunc>> numInstPatt(sizeofPatttern);
	bool isInHashMap;

	// 1. Get instances of the super patterns	
	for (it = CoLHM.begin(); it != CoLHM.end(); ++it)
	{
		// Check if the key in cpHashMap is super set of the current pattern
		isInHashMap = std::includes(it->first.begin(), it->first.end(), pattern.begin(), pattern.end()); // method 1: use std		
		// isInHashMap = a & (~b) ;// method 3: use bitwise difference
		if (isInHashMap)
		{
			// Get Instances
			for (int t = 0; t < sizeofPatttern; t++)
			{
				numInstPatt[t].insert(it->second.find(pattern[t])->second.begin(), it->second.find(pattern[t])->second.end());
			}
		}
	}
	// 2. Calculate inter utility
	std::vector<float> upr(sizeofPatttern, 0.0f);
	for (size_t i = 0; i < sizeofPatttern; ++i)
	{
		for (auto const& pinst : numInstPatt[i])
		{
			upr[i] += pinst.utility;
		}

		upr[i] = upr[i] / (float)uf.find(pattern[i])->second;
	}

	// 3. Calculate extern utility
	std::vector<float> extupr(sizeofPatttern, 0.0f);
	float sumReFeats = 0.0f;
	for (int i = 0; i < sizeofPatttern; i++)
	{
		//std::cout << "i: " << i << endl;
		for (int j = 0; j < sizeofPatttern; ++j)
		{
			if (j != i)
			{
				//std::cout << "j: " << j << endl;
				for (auto const& parinst : numInstPatt[j])
				{
					//std::cout << "ca \n";
					extupr[i] += parinst.utility;
				}
				// Calculate utility of the remain features
				//std::cout << "feature: " << pattern[i] << endl;
				sumReFeats += (float)uf.find(pattern[j])->second;
			}
		}
		extupr[i] = extupr[i] / sumReFeats;
		//std::cout << "sum uf: " << sumReFeats;
		sumReFeats = 0.0f;
	}

	// 4. Calculate UPI
	for (int t = 0; t < sizeofPatttern; ++t)
	{
		upr[t] = alpha * upr[t] + beta * extupr[t];
	}

	float UPI = *std::min_element(std::begin(upr), std::end(upr));

	UPI = Precision(UPI, 3);
	//std::cout << pattern << " : " << UPI << endl;

	return UPI;
}



/**
* @brief: This function calculate the avg error of delta closed
* @param: allPrevPats: all patterns
*
* @retval:
*/
void calculateAvgError(std::unordered_map<std::string, float>& allPCP,
	float delta_thres,
	std::pair<int, float>& results)
{
	// 1. Group into meta set patterns

	std::map<std::string, std::vector<std::string>> Y; // Save Y and it's covered patterns
	std::map<std::string, std::vector<std::string>>::iterator itY;

	std::unordered_map < std::string, float>::iterator itall = allPCP.begin();
	std::unordered_map < std::string, float>::iterator itinner;

	while (itall != allPCP.end())
	{
		bool flagDel = true; // Assume itSortPats is a meta patterns

		for (itinner = allPCP.begin(); itinner != allPCP.end(); ++itinner)
		{
			// Check if other patterns are super sets of this pattern and if their PI difference is smaller than delta_thres
			if (itinner->first.size() > itall->first.size())
			{
				bool checkinclude = std::includes(itinner->first.begin(), itinner->first.end(),
					itall->first.begin(), itall->first.end());

				if (checkinclude)
				{
					float temp = Precision(abs(itinner->second - itall->second), 3);

					if (temp < delta_thres)
					{
						flagDel = false; // not a delta-closed PCP
						//std::cout << "insert to Y /n";
						// Check this key exsit
						if (Y.size() == 0)
						{
							Y.insert({ itinner->first, std::vector<std::string> {itall->first} });
						}
						itY = Y.find(itinner->first);
						if (itY != Y.end())
						{
							itY->second.push_back(itall->first);
						}
						else
						{
							Y.insert({ itinner->first, std::vector<std::string> {itall->first} });
						}
						// Check in Y to delete RK
						itY = Y.find(itall->first);
						if (itY != Y.end())
						{
							itY = Y.erase(itY);
						}
					}
				}
			}
		}
		// check for the current check pattern is or not deleted
		if (flagDel)
		{
			Y.insert({ itall->first, std::vector<std::string> {itall->first} });
		}
		// Terminate	
		++itall;
	}

	// 2. Calcualte avg error
	float avger = 0.0;
	float subv = 0.0;
	float PIY = 0.0;

	//printPrevalentPattern(Rk);
	for (itY = Y.begin(); itY != Y.end(); ++itY)
	{
		// Get PI of itY
		for (itall = allPCP.begin(); itall != allPCP.end(); ++itall)
		{
			if (itall->first == itY->first)
			{
				PIY = itall->second;
				break;
			}
		}
		// Calculate error
		for (auto const& c : itY->second)
		{
			for (itall = allPCP.begin(); itall != allPCP.end(); ++itall)
			{
				if (itall->first == c)
				{
					subv = subv + std::abs(itall->second - PIY);
					break;
				}
			}
		}
	}

	results.first = Y.size();
	results.second = subv / allPCP.size();

}



