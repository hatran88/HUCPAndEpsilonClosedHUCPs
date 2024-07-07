
#ifndef READERCSV_H
#define READERCSV_H

#include <iostream>
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>



/*
* A class to read data from a csv file.
*/
struct CSVReader
{
	std::string fileName;
	std::string delimeter;

	CSVReader(std::string filename, std::string delm = ",") :
		fileName(filename), delimeter(delm) { }  // construct function

	// Function to fetch data from a CSV File
	std::vector<std::vector<std::string> > getData();
};

#endif
