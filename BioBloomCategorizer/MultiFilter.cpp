/*
 * MultiFilter.cpp
 *
 *  Created on: Sep 7, 2012
 *      Author: cjustin
 */

#include "MultiFilter.h"

MultiFilter::MultiFilter(uint16_t hashNum, uint16_t kmerSize) :
		m_hashNum(hashNum), m_kmerSize(kmerSize) {

}

void MultiFilter::addFilter(string const &filterID,
		boost::shared_ptr<BloomFilter> filter) {
	m_filters[filterID] = filter;
	m_filterIDs.push_back(filterID);
}

//todo: implement partial hash function hashing (ie. Only half the number of hashing values for one filter)
/*
 * checks filters for kmer, hashing only single time
 */
const boost::unordered_map<string, bool> MultiFilter::multiContains(
		const char* kmer) {
	const vector<size_t> &hashResults = m_filters.begin()->second->multiHash(kmer);
	boost::unordered_map<string, bool> tempResults;
	for (boost::unordered_map<string, boost::shared_ptr<BloomFilter> >::iterator it =
			m_filters.begin(); it != m_filters.end(); ++it) {
		tempResults[(*it).first] = ((*it).second)->contains(hashResults);
	}
	return tempResults;
}

/*
 * checks filters for k-mer, given a list of filterIDs, hashing only single time
 */
const boost::unordered_map<string, bool> MultiFilter::multiContains(
		const char* kmer, vector<string> const &tempFilters) {
	const vector<size_t> &hashResults = m_filters.begin()->second->multiHash(kmer);
	boost::unordered_map<string, bool> tempResults;
	for (vector<string>::const_iterator it = tempFilters.begin();
			it != tempFilters.end(); ++it) {
		tempResults[*it] = m_filters.at(*it)->contains(hashResults);
	}

	return tempResults;
}

const vector<string> &MultiFilter::getFilterIds() const {
	return m_filterIDs;
}

MultiFilter::~MultiFilter() {
}

