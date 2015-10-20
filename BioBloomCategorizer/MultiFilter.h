/*
 * MultiFilter.h
 *
 *  Created on: Sep 7, 2012
 *      Author: cjustin
 */

#ifndef MULTIFILTER_H_
#define MULTIFILTER_H_
#include <string>
#include <vector>
#include "bloomfilter/BloomFilter.hpp"
#include "boost/unordered/unordered_map.hpp"
#include "boost/shared_ptr.hpp"
using namespace std;

class MultiFilter {
public:
	MultiFilter(uint16_t hashNum, uint16_t kmerSize);
	void addFilter(string const &filterID, boost::shared_ptr<BloomFilter> filter);
	const boost::unordered_map<string, bool> multiContains(const char* kmer);
	const boost::unordered_map<string, bool> multiContains(const char* kmer,
			vector<string> const &tempFilters);
	const BloomFilter &getFilter(const string &filterID);
	const vector<string> &getFilterIds() const;
	virtual ~MultiFilter();
private:
	boost::unordered_map<string, boost::shared_ptr<BloomFilter> > m_filters;
	uint16_t m_hashNum;
	uint16_t m_kmerSize;
	vector<string> m_filterIDs;
};

#endif /* MULTIFILTER_H_ */
