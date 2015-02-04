/*
 * WindowedFileParser.h
 *
 *  Created on: Jul 18, 2012
 *      Author: cjustin
 */

#ifndef WINDOWEDFILEPARSER_H_
#define WINDOWEDFILEPARSER_H_
#include <vector>
#include <boost/unordered/unordered_map.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include "DataLayer/FastaReader.h"
#include "Common/ReadsProcessor.h"
#include <deque>

using namespace std;
using namespace boost;

class WindowedFileParser {
public:
	//constructor/destructor
	explicit WindowedFileParser(const string &fileName, unsigned windowSize, unsigned blocksize = 100000);
	const vector<string> getHeaders() const;
	void setLocationByHeader(const string &header);
	size_t getSequenceSize(const string &header) const;
	string getNextSeq();

	virtual ~WindowedFileParser();

private:
	struct FastaIndexValue {
		size_t index;
		size_t size;
		size_t start;
		size_t bpPerLine;
		size_t charsPerLine;
		size_t numLines;
	};

	unordered_map<string, FastaIndexValue> m_fastaIndex;
	ifstream m_fastaFileHandle;
	unsigned m_windowSize;
	vector<string> m_headers;
	string m_currentHeader;
	size_t m_currentLineNum; //current sequence location from header
	string m_currentString;
	bool m_sequenceNotEnd;

	string m_bufferString; //so reallocation does not need to occur

	unsigned m_blocksize; //minimum size to distribute sequences as

	//helper methods
	void initializeIndex(string const &fileName);

};

#endif /* WINDOWEDFILEPARSER_H_ */
