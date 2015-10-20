/*
 * ReadsProcessor.h
 *
 *  Created on: Aug 8, 2012
 *      Author: cjustin
 */

#ifndef READSPROCESSOR_H_
#define READSPROCESSOR_H_
#include <string>
#include <stdint.h>

using namespace std;

class ReadsProcessor {
public:

	ReadsProcessor(unsigned windowSize);
	const char* prepSeq(string const &sequence, size_t position);

	virtual ~ReadsProcessor();
private:
	char * m_output; //containers preventing reallocation of mem
	unsigned m_kmerSize;
};

#endif /* READSPROCESSOR_H_ */
