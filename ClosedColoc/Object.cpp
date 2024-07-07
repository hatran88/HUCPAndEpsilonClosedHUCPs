#include <iostream>
#include "Object.h"


// For ObjWithCoord
bool ObjWithCoord::operator==(const ObjWithCoord& obj) const
{
	return feature == obj.feature && instance == obj.instance;
}


bool ObjWithCoord::operator<(const ObjWithCoord& obj) const
{
	if (feature < obj.feature) return true;
	else if (feature == obj.feature && instance < obj.instance) return true;
	else return false;
}


// For ObjWithoutCoord
bool ObjWithoutCoord::operator==(const ObjWithoutCoord& obj) const
{
	return feature == obj.feature && instance == obj.instance;
}

bool ObjWithoutCoord::operator<(const ObjWithoutCoord& obj) const
{
	if (feature < obj.feature) return true;
	else if (feature == obj.feature && instance < obj.instance) return true;
	else return false;
}

