

#ifndef FINDMETAPATTS_H
#define FINDMETAPATTS_H

#include <iostream>
#include <string>
#include <vector>



using namespace std;


/**
* @brief: This function filters meta patterns.
* @param: sortPrevPats: sorted by sizes from low to high
* metaPrevPats: the result
* prev_thres: prevalence threshold
* theta: tolerance
* @retval: Nope
*/
void findMetaPats(std::vector<std::pair<std::string, float>>& sortPrevPats,
	std::vector<std::pair<std::string, float>>& metaPrevPats,
	float prev_thres, float theta);




#endif
