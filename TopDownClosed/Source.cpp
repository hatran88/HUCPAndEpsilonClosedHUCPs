/* Source.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Date: 19/3/2022
// Writtern by: Hatran
*/


#include <iostream> 
#include <fstream> 
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include <map>
#include <set>
#include <ctime>
#include <unordered_map>
#include <unordered_set>


#include "Object.h"
#include "ReaderCSV.h"
#include "PrintFunctions.h"
#include "ToolFunctions.h"
#include "findMetaPatts.h"


// For collect memory usage
#include <Windows.h>
#include <stdio.h>
#include <Psapi.h>
#pragma comment(lib, "psapi.lib")


//#include "single_include/nlohmann/json.hpp" // using json
//using json = nlohmann::json;


using namespace std;


int main()
{

	// Set time
	clock_t startTime, endTime; // start and end time
	clock_t timeMaterializeNeighbor; 
	clock_t timeFindMaxCliques;
	clock_t timeQueryAndFilter; 
	
	double time_taken; // count total execution time
	double totalTimeMatStarNei = 0.0f; 
	double totalTimeFinMaxCliques = 0.0f;
	double totalTimeQueryPartInsts = 0.0f;

	startTime = clock(); // Begin count time

	// 1. Set parameters
	float dist_thres = 350.0f;  // set the distance threshold
	float utility_thres = 0.15f; // set the prevalence

	float alpha = 0.6f;
	float beta = 0.4f;
	float delta_thres = 0.05f;
	
	// 2. Load data sets
	string file_name = "./Data/Ex_1/Shanghai_POI_sample_035_checkin_exp_dis.csv";
	freopen("./Data/Ex_1/Shanghai_POI_sample_035_checkin_exp_dis_distance_350_PI_015_delta_005.txt", "w", stdout);

	// Creating an object of CSVReader/CSVWriter
	CSVReader reader(file_name); 
	std::vector<std::vector<std::string>> dataList = reader.getData(); // Get the data from CSV File																   
	// Delete the first line
	dataList.erase(dataList.begin());

	// 1.2 Count the instance number of each feature
	// the result is saved into a global variable wf
	getFeatureWeight(dataList);
	//printFeatureWeight();	
	int n = dataList.size(); // the number of instances

	// Step 2. Make grid and find neighbors
	// 2.1 Make grid and get all instances (vertices)	
	std::map<std::pair<int, int>, std::vector<Objects>> grid = makeGrid(dataList, dist_thres);
	dataList.clear();
	//printGrid(grid);

	std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc> neighbors;

	// Save all instances
	std::unordered_map<ObjWithoutCoordinate, int, myHashFunc> vertices(dataList.size());
	// Get neighbors
	genStarNeighborhoods(neighbors, vertices, grid, dist_thres);
	//printStarNeighborhood(neighbors);	

	// Get the size neighboring instance of each instance
	//getSizeNeiOfInst(neighbors);

	timeMaterializeNeighbor = clock();
	totalTimeMatStarNei = totalTimeMatStarNei + double(timeMaterializeNeighbor - startTime);

	timeFindMaxCliques = clock();
	// Step 3. Finding maximal cliques by iterator each instance and its neighbors
	//std::vector<ObjWithoutCoordinate> verticesOneBlock;
	bool checkdel;
	std::unordered_map<ObjWithoutCoordinate, std::vector<ObjWithoutCoordinate>, myHashFunc>::iterator itSN, itIn;
	// Process each instance and its neighboring instances as S
	while (vertices.size())
	{
		// Get one instance and put into the verticesOneBlock
		verticesVec.push_back(vertices.begin()->first);		

		// Get the neighbor instances of this instance
		itSN = neighbors.find(vertices.begin()->first);
		if (itSN != neighbors.end())
		{
			// 1. Put all the neighbor into verticesOneBlock to build G+
			verticesVec.insert(verticesVec.end(), itSN->second.begin(), itSN->second.end());
			// Del the first vertice
			vertices.erase(vertices.begin());
			//std::cout << "After sort: \n";
			std::sort(verticesVec.begin(), verticesVec.end());
			//printVeticeVec(verticesVec);

			// 2. Build G-
			// If the neighboring instances of an instance oi in G+ are included in G+, oi can be deleted from the data set
			for (auto const& inst : itSN->second)
			{
				// Get all neighboring instances of the instance
				itIn = neighbors.find(inst);
				//if (itIn != neighbors.end() && itIn->second.size())
				if (itIn != neighbors.end())
				{
					checkdel = std::includes(verticesVec.begin(), verticesVec.end(), itIn->second.begin(), itIn->second.begin());				
					if (checkdel)
					{
						// Delete this instance from vertices
						vertices.erase(inst);
						//std::cout << "Del vertices: " << vertices.size() << endl;
					}			
				}				
			}		
		}
		// Find maximal cliques in G+
		if (verticesVec.size() > 1)
		{
			//std::sort(verticesVec.begin(), verticesVec.end());
			// Get neighboring instances			
			getNeiOneBlock(neighbors);
			// Print neighbor of one block
			//printNeiOneBlock(neiOneBlock);
			// Compute the degency			
			std::vector<ObjWithoutCoordinate> orderDeg;
			OrderedDegeneracy(orderDeg, neiOneBlock);			
			// Find maximal cliques
			std::vector<ObjWithoutCoordinate> X;
			BronKerboschDeg(verticesVec, X, orderDeg);
			
		}
		// Clear
		neiOneBlock.clear();
		verticesVec.clear();
	}

	totalTimeFinMaxCliques = totalTimeFinMaxCliques + double(clock() - timeFindMaxCliques) - totalTimeConstCoLHashmap;

	//std::cout << "Number of maximal cliques: " << numOfMaxCl << endl;
	
	// Ste 4. Construct co-location hashmap	
	// This step is included in the listing maximal cliques step.
	//printCoLHM();

	// This step is included in the listing maximal cliques step.
	timeQueryAndFilter = clock();
	std::unordered_map<std::string, float> deltaCPCP, allPCP, closedPCP;
	filterDeltaClosedPCPs(deltaCPCP, allPCP, closedPCP, utility_thres, delta_thres, alpha, beta);
	totalTimeQueryPartInsts = double(clock() - timeQueryAndFilter) - totalTimeComputeUPIs;
		
	std::cout << "Delta closed PCPs: " << deltaCPCP.size() << endl;
	//printPrevPatts(deltaCPCP);
	
	// Get the maximal size of patterns
	std::vector<int> sizePats;
	std::unordered_map<string, float>::iterator itAllPats = deltaCPCP.begin();
	while (itAllPats != deltaCPCP.end())
	{
		sizePats.push_back(itAllPats->first.size());
		++itAllPats;
	}
	int maxSize = *std::max_element(sizePats.begin(), sizePats.end());
	std::cout << "The maximal size of patterns is: " << maxSize << endl;

	// Calculate the tolerance
	// allPrevPats: all prevalent patterns
	// 1. sort all pattern by size
	std::vector<std::pair<std::string, float>> sortPrevPats;
	sorthashmap(deltaCPCP, sortPrevPats);
	//printSortPrevPats(sortPrevPats);

	//2. Classify pattern by sizes
	std::map<int, std::vector<Patterns>> classPattBySize;
	classifyPattsBySize(deltaCPCP, classPattBySize);
	printNumberPattsBySize(classPattBySize);

	//printTopKPatts(classPattBySize, topk);

	// 3. Calculate the avg error
	//std::pair<int, float> results; // int is the number of meta patterns and float is the error 
	//calculateAvgError(allPCP, delta_thres, results);

	////std::cout << "The number of meta patterns: " << results.first << endl;
	//std::cout << "The average error: " << results.second << endl;

	//// 3. Compute compression
	//std::cout << "The compression rate of closed: " << float(allPCP.size() - closedPCP.size()) / (float)allPCP.size() * 100.0f << endl;
	//std::cout << "The compression rate of delta-closed: " << float(allPCP.size() - deltaCPCP.size()) / (float)allPCP.size() * 100.0f << endl;

	// Step 4. Count time and print time
	endTime = clock();
	time_taken = double(endTime - startTime);

	std::cout << "Time for meterializing neighbors: " << totalTimeMatStarNei / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time for finding maximal cliques: " << totalTimeFinMaxCliques / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time for constructingco-location hashmap: " << totalTimeConstCoLHashmap / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time for querying participating instances: " << totalTimeQueryPartInsts / CLOCKS_PER_SEC << " s." << endl;
	std::cout << "Time for filtering HUCPs: " << totalTimeComputeUPIs / CLOCKS_PER_SEC << " s." << endl;	
	std::cout << "Time taken by the program is : " << time_taken / CLOCKS_PER_SEC << " s." << endl;

	// Show memory usage
	HANDLE handle = GetCurrentProcess();
	PROCESS_MEMORY_COUNTERS memCounter;
	GetProcessMemoryInfo(handle, &memCounter, sizeof(memCounter));
	SIZE_T physMemUsedByMe = memCounter.WorkingSetSize;
	SIZE_T peakMemUsedByMe = memCounter.PeakWorkingSetSize;
	std::cout << "Memory usage: " << (float)physMemUsedByMe / 1024 /1024 << "(MB)" << std::endl;
	std::cout << "Peak memory usage: " << (float)peakMemUsedByMe / 1024 / 1024 << "(MB)" << std::endl;

	
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
