#include <iostream>
#include "Object.h"
#include <string>
#include <map>
#include <vector>
#include <set>
#include <unordered_map>

#include "ToolFunctions.h"


void findMetaPats(std::vector<std::pair<std::string, float>> &sortPrevPats, 
	std::vector<std::pair<std::string, float>>& metaPrevPats,
	float prev_thres, float theta) 
{
	std::vector<std::pair<std::string, float>>::iterator itSortPats = sortPrevPats.begin();
	int numPat = sortPrevPats.size(); // the number of all prevalent patterns
	int t = 0;

	while (itSortPats != sortPrevPats.end())
	{
		bool flagDel = true; // Assume itSortPats is a meta patterns		
		for (int i = t; i < numPat; i++) {

			// Check if other patterns are super sets of this pattern and if the PI (1-theta)
			if (sortPrevPats[i].first.size() > itSortPats->first.size())
			{
				bool checkinclude = std::includes(sortPrevPats[i].first.begin(), sortPrevPats[i].first.end(),
					itSortPats->first.begin(), itSortPats->first.end());

				if (checkinclude)
				{
					float temp = Precision((1 - theta) * prev_thres + theta * itSortPats->second, 3);
					if (sortPrevPats[i].second >= temp)
					{
						// itSortPats is needed to deleted
						flagDel = false;

						break;
					}
				}
			}
		}
		// Check itSortPats is a meta patterns
		if (flagDel)
		{
			metaPrevPats.push_back(*itSortPats);
		}
		// Terminate
		++t;
		++itSortPats;
	}
}

