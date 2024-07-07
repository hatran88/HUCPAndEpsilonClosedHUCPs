#include <iostream>
#include "Object.h"

using namespace std;


// Match both features and instances for hash_map 
bool Objects::operator==(const Objects& obj) const
{
	return feature == obj.feature && instance == obj.instance;
}
// compare two instances
bool Objects::operator<(const Objects& obj) const
{
	if (feature < obj.feature) return true;
	else if (feature == obj.feature && instance < obj.instance) return true;
	else return false;
}




// Match both features and instances for hash_map 
bool ObjWithoutCoordinate::operator==(const ObjWithoutCoordinate& obj) const
{
	return feature == obj.feature && instance == obj.instance;
}


bool ObjWithoutCoordinate::operator<(const ObjWithoutCoordinate& obj) const
{
	if (feature < obj.feature) return true;
	else if (feature == obj.feature && instance < obj.instance) return true;
	else return false;
}


// Hash function for unordered_map
std::size_t myHashFunc::operator()(const ObjWithoutCoordinate& inst) const
{
	std::size_t h1 = std::hash<char>()(inst.feature);
	std::size_t h2 = std::hash<int>()(inst.instance);

	//return (std::hash<char>()(inst.feature)) ^ (std::hash<int>()(inst.instance));
	return h1 ^ h2;
}

