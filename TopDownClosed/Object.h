#ifndef OBJECT_H
#define OBJECT_H

#include <iostream>


using namespace std;

// Define an instance structure
struct Objects
{
	char feature;
	int instance;
	float utility; // save weight of instance
	float x;
	float y;

	Objects(char feat, int inst, float u, float xcor, float ycor)
	{
		feature = feat;
		instance = inst;
		utility = u;
		x = xcor;
		y = ycor;
	}


	// Match both features and instances for hash_map 
	bool operator==(const Objects& obj) const;
	// compare two instances
	bool operator<(const Objects& obj) const;
};

// Define Objects without coordinate information
struct ObjWithoutCoordinate
{
	char feature;
	int instance;
	float utility; // save weight of instance

	ObjWithoutCoordinate(char feat, int inst, float u)
	{
		feature = feat;
		instance = inst;
		utility = u;
	}

	// Match both features and instances for hash_map 
	bool operator==(const ObjWithoutCoordinate& obj) const;

	// compare two instances
	bool operator<(const ObjWithoutCoordinate& obj) const;
};


// Hash function for unordered_map
struct myHashFunc
{
	std::size_t operator()(const ObjWithoutCoordinate& inst) const;
};



#endif
