/*
 * ResultsManager.cpp
 *
 *  Created on: Jun 25, 2013
 *      Author: cjustin
 */

#include "ResultsManager.h"
#include <sstream>
#include <iostream>
#include <BioBloomClassifier.h>
#if _OPENMP
# include <omp.h>
#endif

ResultsManager::ResultsManager(const vector<string> &filterOrderRef,
		bool inclusive) :
		m_filterOrder(filterOrderRef), m_multiMatch(
				0), m_noMatch(0), m_inclusive(inclusive)
{
	//initialize variables and print filter ids
	for (vector<string>::const_iterator i = m_filterOrder.begin();
			i != m_filterOrder.end(); ++i)
	{
		m_aboveThreshold[*i] = 0;
		m_unique[*i] = 0;
	}
}

/*
 * Records data for read summary based on thresholds
 * Returns filter ID that this read equals
 */
const string ResultsManager::updateSummaryData(
		const unordered_map<string, bool> &hits)
{
	string filterID;
	bool noMatchFlag = true;
	bool multiMatchFlag = false;

	for (vector<string>::const_iterator i = m_filterOrder.begin();
			i != m_filterOrder.end(); ++i)
	{
		if (hits.at(*i)) {
#pragma omp atomic
			++m_aboveThreshold[*i];
			if (noMatchFlag) {
				noMatchFlag = false;
				filterID = *i;
			} else {
				multiMatchFlag = true;
			}
		}
	}
	if (noMatchFlag) {
		filterID = NO_MATCH;
#pragma omp atomic
		++m_noMatch;
	} else {
		if (multiMatchFlag) {
			filterID = MULTI_MATCH;
#pragma omp atomic
			++m_multiMatch;
		} else {
		//TODO : USE LOCKS, # of locks == # of filterIDs
#pragma omp atomic
			++m_unique[filterID];
		}
	}
	return filterID;
}

/*
 * Records data for read summary based on thresholds
 * Returns filter ID that this read pair equals
 */
const string ResultsManager::updateSummaryData(
		const unordered_map<string, bool> &hits1,
		const unordered_map<string, bool> &hits2)
{
	string filterID;
	bool noMatchFlag = true;
	bool multiMatchFlag = false;

	for (vector<string>::const_iterator i = m_filterOrder.begin();
			i != m_filterOrder.end(); ++i)
	{
		if (m_inclusive) {
			if (hits1.at(*i) || hits2.at(*i)) {
#pragma omp atomic
				++m_aboveThreshold[*i];
				if (noMatchFlag) {
					noMatchFlag = false;
					filterID = *i;
				} else {
					multiMatchFlag = true;
				}
			}
		} else {
			if (hits1.at(*i) && hits2.at(*i)) {
#pragma omp atomic
				++m_aboveThreshold[*i];
				if (noMatchFlag) {
					noMatchFlag = false;
					filterID = *i;
				} else {
					multiMatchFlag = true;
				}
			}
		}
	}
	if (noMatchFlag) {
		filterID = NO_MATCH;
#pragma omp atomic
		++m_noMatch;
	} else {
		if (multiMatchFlag) {
			filterID = MULTI_MATCH;
#pragma omp atomic
			++m_multiMatch;
		} else {
#pragma omp atomic
			++m_unique[filterID];
		}
	}
	return filterID;
}

const string ResultsManager::getResultsSummary(size_t readCount) const
{

	stringstream summaryOutput;

	//print header
	summaryOutput
			<< "filter_id\thits\tmisses\tshared\trate_hit\trate_miss\trate_shared\n";

	for (vector<string>::const_iterator i = m_filterOrder.begin();
			i != m_filterOrder.end(); ++i)
	{
		summaryOutput << *i;
		summaryOutput << "\t" << m_aboveThreshold.at(*i);
		summaryOutput << "\t" << readCount - m_aboveThreshold.at(*i);
		summaryOutput << "\t" << (m_aboveThreshold.at(*i) - m_unique.at(*i));
		summaryOutput << "\t"
				<< double(m_aboveThreshold.at(*i)) / double(readCount);
		summaryOutput << "\t"
				<< double(readCount - m_aboveThreshold.at(*i))
						/ double(readCount);
		summaryOutput << "\t"
				<< double(m_aboveThreshold.at(*i) - m_unique.at(*i))
						/ double(readCount);
		summaryOutput << "\n";
	}

	summaryOutput << MULTI_MATCH;
	summaryOutput << "\t" << m_multiMatch;
	summaryOutput << "\t" << readCount - m_multiMatch;
	summaryOutput << "\t" << 0;
	summaryOutput << "\t" << double(m_multiMatch) / double(readCount);
	summaryOutput << "\t"
			<< double(readCount - m_multiMatch) / double(readCount);
	summaryOutput << "\t" << 0.0;
	summaryOutput << "\n";

	summaryOutput << NO_MATCH;
	summaryOutput << "\t" << m_noMatch;
	summaryOutput << "\t" << readCount - m_noMatch;
	summaryOutput << "\t" << 0;
	summaryOutput << "\t" << double(m_noMatch) / double(readCount);
	summaryOutput << "\t" << double(readCount - m_noMatch) / double(readCount);
	summaryOutput << "\t" << 0.0;
	summaryOutput << "\n";

	cerr << summaryOutput.str() << endl;
	return summaryOutput.str();
}

ResultsManager::~ResultsManager()
{
}

