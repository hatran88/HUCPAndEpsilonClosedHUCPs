#pragma once
#ifndef PATTERNS_H
#define PATTERNS_H

#include <iostream>


using namespace std;


// Define Objects without coordinate information
struct Patterns
{
	string patts;
	float PI;
	
	// Match both features and instances for hash_map 
	bool operator==(const Patterns& obj) const;
	// compare two instances
	bool operator<(const Patterns& obj) const;
};



#endif
