#include <iostream>
#include "Patterns.h"

using namespace std;


// Match both features and instances for hash_map 
bool Patterns::operator==(const Patterns& obj) const
{
	return PI == obj.PI;
}


bool Patterns::operator<(const Patterns& obj) const
{
	if (PI < obj.PI) return true;	
	else return false;
}


//// Hash function for unordered_map
//std::size_t myHashFunc::operator()(const ObjWithoutCoordinate& inst) const
//{
//	std::size_t h1 = std::hash<char>()(inst.feature);
//	std::size_t h2 = std::hash<int>()(inst.instance);
//
//	//return (std::hash<char>()(inst.feature)) ^ (std::hash<int>()(inst.instance));
//	return h1 ^ h2;
//}

