/* 
// General + MaxColoc : This code implies the maximal co-location patterns
// A framework for generating condensed co-location sets from spatial databases
// Proposed by: Jin Soung Yoo and Mark Bow
// 
//
// Writtern by: Hatran
// Date: 25/0/2021
*/



#include <iostream>
#include <ctime>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <string>
#include <unordered_map>

#include "tree.hh"

#include "Object.h"
#include "printFunctions.h"
#include "readerCSV.h"
#include "toolFunctions.h"

// For collect memory usage
#include <Windows.h>
#include <stdio.h>
#include <Psapi.h>
#pragma comment(lib, "psapi.lib")


using namespace std;


int main()
{
	// Set time
	clock_t startTime; // start and end time
	clock_t timeMaterializeNeighbor; // time for materializing star neighborhoods
	clock_t timeGenCandidate; // time for generating candidates
	clock_t timeFindCliqueInsts;	
	clock_t timeCalPIAndFlter;
	clock_t timeStarWhile;

	double time_taken; // count total execution time
	double totalTimeMatStarNei = 0.0f; // total time for materializing star neighborhoods
	double totalTimeGenCand = 0.0f; // total time for generating candidates	
	double totalTimeFilCliqueInst = 0.0f; // total time for collecting candidates	
	double totalTimeFilPrevCoLoPat = 0.0f; // total time for filgering candidates

	startTime = clock(); // Begin count time

	// 1. Set parameters
	float dist_thres = 2.50f;  // set the distance threshold
	float prev_thres = 0.2f; // set the prevalence

	// 2. Load data sets
	string file_name = "./Data/testDataset.csv";
	freopen("./Data/testDataset_distance_2.5_PI_03.txt", "w", stdout);

	// 1.1 Load the dataset
	CSVReader reader(file_name);
	std::vector<std::vector<std::string>> dataList = reader.getData();
	// Delete the first line
	dataList.erase(dataList.begin());

	// 3. Count the instance number of each feature
	std::map<char, int> totalInstEachFeat = countNumberInstance(dataList);
	//printInstanceNumber(totalInstEachFeat);

	// 4. Make grid and find neighbors
	// 4.1 Make grid
	std::map<std::pair<int, int>, std::vector<ObjWithCoord>> grid = makeGrid(dataList, dist_thres);
	//printGrid(grid);
	// 4.2 Find neighbors, generate A1: <B.1, C.1>
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN;
	genStarNeighborhoods(grid, dist_thres, SN);
	//printStarInstance(SN);	

	// 4.4 Find event neighborhood transactions
	std::map<char, std::vector<std::set<char>>> featNeibTrans = genFeatNeibTrans(SN);
	//printFeatNeibTrans(featNeibTrans);	

	timeMaterializeNeighbor = clock();
	totalTimeMatStarNei = double(timeMaterializeNeighbor - startTime);

	// 5. Generate candidates
	// 5.1 Build FP-growth for each feature to get a set of candidates (Step 3-5)
	std::map<char, std::map<std::string, float>> Treei = buildCandPattTrees(featNeibTrans);	

	//printTreei(Treei);	
	// Collect candidate from tree (Step 6)
	std::map<std::string, float> C = genCands(Treei, prev_thres);
	//std::cout << "Num of cand: " << C.size() << endl;
	//printCands(C);
	// 5.2 Starting mining from size 2 candidates (Step 7)
	int k = 2;

	timeGenCandidate = clock();
	totalTimeGenCand = totalTimeGenCand + double(timeGenCandidate - timeMaterializeNeighbor);
	
	// save size k prevalent patterns 
	std::map<std::string, float, cmpMapBySizePatt> Pk;
	// and non - prevalent of size k
	std::map < std::string, float> nonPk;
	// Save all prevalent patterns
	std::map<std::string, float, cmpMapBySizePatt> Rk;
	// Save size k candidates
	std::map<std::string, float> Ck;
	// Save table instanes of size k and k+1
	std::map<std::string, std::vector<std::vector<ObjWithoutCoord>>> Tk;

	int pass = 0;
	std::map<int, int> numCand;
	std::map<int, double> timeCollectInstSizek;	
	std::map<float, std::vector<std::string>> CloseP;

	// 6. Process each candidate to determine prevalence (Steps 8-19)
	while (!C.empty())
	{
		// total time for collecting co-location instances of candidates
		double totalTimeFilCliqueInstSizek = 0.0f; 
		// Count the number of iterators
		pass += 1;		

		timeStarWhile = clock();
		// Step 9: Get all size k candidates
		Ck = getSizeKCands(C, k);
	//	printCands(Ck);
		// Count the number of size k patterns
		numCand.insert({ k, Ck.size() });

		//printSizekCands(Ck);		
		timeGenCandidate = clock();
		totalTimeGenCand = totalTimeGenCand + double(timeGenCandidate - timeStarWhile);

		// Step 10 Collect star instances
		// SIK ->  ABC : <<instances A>, <instances B>, <instances C>>
		std::map<std::string, std::vector<std::vector<std::vector<ObjWithoutCoord>>>> SIk;		

		findStarInsts(Ck, SN, SIk);
		//printCoarseTableInsts(SIk);	
				
		// If k == 2, directly construct table instances and calcualate PI
		if (k == 2)
		{
			//Step 11:  Collect table instances
			Tk = collectSize2TableInst(SIk, k, Ck);
			
			timeFindCliqueInsts = clock();
			totalTimeFilCliqueInst = totalTimeFilCliqueInst + double(timeFindCliqueInsts - timeGenCandidate);

			totalTimeFilCliqueInstSizek = totalTimeFilCliqueInstSizek 
				+ double(timeFindCliqueInsts - timeGenCandidate);

			// Step 12: Calculate PIs and filter prevalent patterns
			selectPrevalentPatterns(Tk, totalInstEachFeat, prev_thres, k, Pk, nonPk);

			timeCalPIAndFlter = clock();
			totalTimeFilPrevCoLoPat = totalTimeFilPrevCoLoPat + double(timeCalPIAndFlter - timeFindCliqueInsts);
		}
		else
		{
			// 4.3 Use coarse prun		
			selectCoarsePrevalentPatterns(SIk, totalInstEachFeat, prev_thres, k);

			// Step 11: Use size k table instances to validate size k+1 table instances
			Tk = filterCliqueInstances2(SIk, Tk);

			timeFindCliqueInsts = clock();
			totalTimeFilCliqueInst = totalTimeFilCliqueInst + double(timeFindCliqueInsts - timeGenCandidate);

			totalTimeFilCliqueInstSizek = totalTimeFilCliqueInstSizek
				+ double(timeFindCliqueInsts - timeGenCandidate);

			// Step 12: Calculate PI and filter prevalent patterns
			selectPrevalentPatterns(Tk, totalInstEachFeat, prev_thres, k, Pk, nonPk);

			timeCalPIAndFlter = clock();
			totalTimeFilPrevCoLoPat = totalTimeFilPrevCoLoPat + double(timeCalPIAndFlter - timeFindCliqueInsts);
		}

		//std::cout << "size k prevalent patterns: \n";
		//printPrevalentPattern(Pk);

		// Step 13: Filter colosed patterns
		filterClosedPatts(Pk, CloseP);
		Pk.clear();
		//printClosedPatts(CloseP);
		// Step 14: Delete already checked candidates
		delCkFromC(C, Ck);
		Ck.clear();

		// Step 15: Candidate prune
		superPrun(nonPk, C);
		nonPk.clear();

		// Add execution time of collecting table instance
		timeCollectInstSizek.insert({k, totalTimeFilCliqueInstSizek});

		// Update 		
		k = k + 1;	
	}	

	// 8. Print results
	int numClosedP = 0;
	for (auto const& c : CloseP)
	{
		numClosedP += c.second.size();
	}
	std::cout << "The number of closed prevalent pattern is: " << numClosedP << endl;

	std::cout << "Execution time for materizaling starneighbors is: " << totalTimeMatStarNei / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Execution time for generating candidates: " << totalTimeGenCand / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for finding clique instances: " << totalTimeFilCliqueInst / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for filtering prevalent co-locations: " << totalTimeFilPrevCoLoPat / CLOCKS_PER_SEC << "s." << endl;

	std::cout << "Time taken by the program is (No neighboring) : " << 
		(totalTimeGenCand + totalTimeFilCliqueInst + totalTimeFilPrevCoLoPat) / CLOCKS_PER_SEC << " s." << endl;
	time_taken = double(clock() - startTime);
	std::cout << "Time taken by the program is (all) : " << time_taken / CLOCKS_PER_SEC << " s." << endl;

	// 9. Show memory usage
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS memCounter;
	GetProcessMemoryInfo(handle, &memCounter, sizeof(memCounter));
	SIZE_T physMemUsedByMe = memCounter.WorkingSetSize;
	SIZE_T peakMemUsedByMe = memCounter.PeakWorkingSetSize;
	std::cout << "Memory usage: " << physMemUsedByMe / 1024/1024 << "(MB)" << std::endl;
	std::cout << "Peak memory usage: " << peakMemUsedByMe / 1024 /2024 << "(MB)" << std::endl;

	
	// print the result
	//printClosedPatts(CloseP);
	

	// Other codes here




	return 0;

}



// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
