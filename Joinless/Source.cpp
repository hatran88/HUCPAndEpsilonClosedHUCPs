
/* This version take class Object and ObjectWithoutCoordinate as individual .h and .cpp
* 
* This version: 1) size k table instances are validate by size (k-1) table instances
*				2) to reduce memory, do not product coarse row instances, 
*				each generate one coarse row instances, use size (k-1) to check the row instance is or not a real clique
* 
* Writer: Hatran
* Date: 1.8.2020
*/



#include <iostream>
#include <ctime>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

// For collect memory usage
#include <Windows.h>
#include <stdio.h>
#include <Psapi.h>
#pragma comment(lib, "psapi.lib")

// Import my function
#include "Object.h"
#include "printFunctions.h"
#include "readerCSV.h"
#include "toolFunctions.h"



using namespace std;


int main()
{
	// Set time
	clock_t startTime, endTime; // start and end time
	clock_t timeMaterializeNeighbor; // time for materializing star neighborhoods	
	clock_t timeGenerateCandidate; // time for generating candidates
	clock_t timeFilStarInstance;
	clock_t timeFilCoarseColoPat;
	clock_t timeFilCliqueInstance; // time for collecting co-location instances
	clock_t timeFilPrevCoLoPat; // time for filtering co-location instances
	clock_t timeWhile;

	double time_taken; // count total execution time
	double totalTimeMatStarNei = 0.0f; // total time for materializing star neighborhoods	
	double totalTimeGenCand = 0.0f; // total time for generating candidates
	double totalTimeFilStarInst = 0.0f;
	double totalTimeFilCoarseCoLoPat = 0.0f;
	double totalTimeFilCliqueInst = 0.0f; // total time for collecting candidates
	double totalTimeFilPrevCoLoPat = 0.0f; // total time for filgering candidates


	startTime = clock(); // Begin count time



	// 1. Set parameters
	float dist_thres = 280.0f;  // set the distance threshold
	float prev_thres = 0.15f; // set the prevalence

	// 2. Load data sets
	string file_name = "./Data/Exp_1_2/Shenzhen_type_1_sample_04_20_feature.csv";
	freopen("./Data/Exp_1_2/Shenzhen_type_1_sample_04_20_feature_distance_280_PI_015.txt", "w", stdout);

	// 1.1 Load the dataset																	  
	CSVReader reader(file_name); // Creating an object of CSVReader/CSVWriter	
	std::vector<std::vector<std::string>> dataList = reader.getData(); // Get the data from CSV File
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
	std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>> SN = genStarNeighborhoods(grid, dist_thres);	
	//printSN(SN);	

	std::map<char, std::map<ObjWithoutCoord, std::set<ObjWithoutCoord>>> groupSNByFeat = groupStarNeighByFeatures(SN);
	//printStarNeighborhood(groupSNByFeat);

	timeMaterializeNeighbor = clock();
	totalTimeMatStarNei = double(timeMaterializeNeighbor - startTime);
	// Free memory
	SN.clear();

	// 5. Initialation
	// 5.1 Set Pk, k, Ck
	std::map<vector<char>, float> Pk, PkAll; // Store all prevalent patterns
	std::vector<char> pattern;
	std::vector<vector<char>> Ckplus; // Store all candidate patterns
	// Store table isntances of size k patterns
	std::map <std::vector<char>, std::vector<std::vector<ObjWithoutCoord>>> CIk; 
	
	int k = 2;
	// 5.2 Get size 1 pattern
	for (auto const& pat : totalInstEachFeat)
	{
		pattern.push_back(pat.first);
		Pk.insert({ pattern, 1.0f });
		pattern.clear();
	}	

	timeFilPrevCoLoPat = clock();
	totalTimeFilPrevCoLoPat = totalTimeFilPrevCoLoPat + double(timeFilPrevCoLoPat - timeMaterializeNeighbor);

	// Step 4. Find size k patterns
	while (Pk.size())
	{		
		timeWhile = clock();

		// 4.1 Generate candidate patterns
		Ckplus = genCandidatePatterns(Pk, k);
		//printCandidates(Ckplus);

		timeGenerateCandidate = clock();
		totalTimeGenCand = totalTimeGenCand + double(timeGenerateCandidate - timeWhile);
		
		// 4.3 Select coarse prevalent patterns
		if (k <= 2) // directly get table instances of candidate patterns and calculate their PIs
		{
			// 4.3 Filter star instances			
			CIk.clear();
			CIk = filterStarInstancesSize2(Ckplus, groupSNByFeat, k);
						
			timeFilCliqueInstance = clock();
			totalTimeFilCliqueInst = totalTimeFilCliqueInst + double(timeFilCliqueInstance - timeGenerateCandidate);
			//printTableInstances(CIkplus);
			// 4.5 Select prevalent patterns			
			Pk = selectPrevalentPatterns(CIk, totalInstEachFeat, prev_thres, k);	
			//printPrevalentPattern(Pk);
			
			// 4.6 Add Pk to the result PkAll
			PkAll.insert(Pk.begin(), Pk.end());
			timeFilPrevCoLoPat = clock();
			totalTimeFilPrevCoLoPat = totalTimeFilPrevCoLoPat + double(timeFilPrevCoLoPat - timeFilCliqueInstance);
		}
		else
		{
			// 4.2 Filter star instances : to store table isntances of size k+1 patterns
			std::map<std::vector<char>, std::vector<std::vector<std::vector<ObjWithoutCoord>>>> SIkplus 
				= filterStarInstances(Ckplus, groupSNByFeat, k);

			timeFilStarInstance = clock();
			totalTimeFilStarInst = totalTimeFilStarInst + double(timeFilStarInstance - timeGenerateCandidate);			
			//printStarInstance(SIkplus);
			
			// 4.3 Use coarse prun			
			selectCoarsePrevalentPatterns(SIkplus, totalInstEachFeat, prev_thres, k);

			timeFilCoarseColoPat = clock();
			totalTimeFilCoarseCoLoPat = totalTimeFilCoarseCoLoPat + double(timeFilCoarseColoPat - timeFilStarInstance);
			
			// 4.4 Filter clique instances			
			CIk = filterCliqueInstances2(SIkplus, CIk);
			//printTableInstances(CIkplus);					
			
			timeFilCliqueInstance = clock();
			totalTimeFilCliqueInst = totalTimeFilCliqueInst + double(timeFilCliqueInstance - timeFilCoarseColoPat);

			// 4.5 Select prevalent patterns			
			Pk = selectPrevalentPatterns(CIk, totalInstEachFeat, prev_thres, k);
			//printPrevalentPattern(Pk);
			
			// 4.6 Add Pk to the result PkAll
			PkAll.insert(Pk.begin(), Pk.end());
			
			timeFilPrevCoLoPat = clock();
			totalTimeFilPrevCoLoPat = totalTimeFilPrevCoLoPat + double(timeFilPrevCoLoPat - timeFilCliqueInstance);
		}

		// 4.6 Go to next iterator
		k += 1;
		Ckplus.clear();
	}

	// Print prevalent patterns	
	std::cout << "The number of prevalent pattern is: " << PkAll.size() << endl;
	//printPrevalentPattern(PkAll);

	// Print running time
	std::cout << "Execution time for materizaling starneighbors is: " << totalTimeMatStarNei / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Execution time for generating candidates: " << totalTimeGenCand / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for filtering star instances: " << totalTimeFilStarInst / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for filtering coarse co-locations: " << totalTimeFilCoarseCoLoPat / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for filtering clique instances: " << totalTimeFilCliqueInst / CLOCKS_PER_SEC << "s." << endl;
	std::cout << "Execution time for filtering prevalent co-locations: " << totalTimeFilPrevCoLoPat / CLOCKS_PER_SEC << "s." << endl;

	// Calculate time
	endTime = clock();
	time_taken = double(endTime - startTime);
	std::cout << "Time taken by the program is (No neighboring) : " << 
		(totalTimeGenCand + totalTimeFilStarInst + totalTimeFilCoarseCoLoPat + totalTimeFilCliqueInst + totalTimeFilPrevCoLoPat) / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time taken by the program is (all) : " << time_taken / CLOCKS_PER_SEC << " s." << endl;

	// Show memory usage
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS memCounter;
	GetProcessMemoryInfo(handle, &memCounter, sizeof(memCounter));
	SIZE_T physMemUsedByMe = memCounter.WorkingSetSize;
	SIZE_T physPeackMemUsedByMe = memCounter.PeakWorkingSetSize;
	std::cout << "Memory usage: " << physMemUsedByMe / 1024 << "(kB)" << std::endl;
	std::cout << "Peak memory usage: " << physPeackMemUsedByMe / 1024 << "(kB)" << std::endl;
	
	// 10. Count the number of each size
	std::vector<int> maxsize;
	int onek;
	std::map<int, int> sizePats;
	std::map<int, int>::iterator itsizePats;

	std::map<vector<char>, float>::iterator itAllPats = PkAll.begin();
	while (itAllPats != PkAll.end())
	{
		onek = itAllPats->first.size();
		maxsize.push_back(onek);

		itsizePats = sizePats.find(onek);
		if (itsizePats != sizePats.end())
		{
			itsizePats->second += 1;
		}
		else
		{
			sizePats.insert({ onek, 1 });
		}
		// Terminate
		++itAllPats;
	}

	int maxSize = *std::max_element(maxsize.begin(), maxsize.end());
	std::cout << "The maximal size of patterns is: " << maxSize << endl;

	std::cout << "Patters by sizes: \n";
	printPattbySize(sizePats);




	// The next code in here



	return 0;

}








