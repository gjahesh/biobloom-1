/*
 * BloomFilterCategorizerTests.cpp
 *
 *  Created on: Feb 15, 2013
 *      Author: cjustin
 */

#include "Common/Uncompress.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char **argv) {
	string filename = "test.txt.gz";

	//TEST IS NOT WORKING, couldn't get gzip to print to file
	//Test functionality of Uncompress.h
	//write out file with extention .gz as compressed file
	ofstream testfile(filename.c_str(), ios::out);
	testfile << "test" << endl;
	testfile.close();
	//read in file with extension .gz as compressed file
//	ifstream testfile2(filename.c_str(), ios::in);
//	string temp;
//	testfile2 >> temp;
//	cout << temp << endl;
	//remove file
//	remove(filename.c_str());
}
