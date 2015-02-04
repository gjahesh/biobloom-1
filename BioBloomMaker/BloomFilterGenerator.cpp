/*
 * BloomFilterGenerator.cpp
 *
 * Using an input file (currently fasta), generates a bloom filter
 * The filter is returned. Intended to create a filter information
 * object as well (not implemented).
 *
 *  Created on: Jun 20, 2012
 *      Author: cjustin
 */
#include "BloomFilterGenerator.h"
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include "WindowedFileParser.h"
#include "Common/BloomFilter.h"
#include "Common/BloomFilterInfo.h"
#include <cassert>
#include <cmath>
#if _OPENMP
# include <omp.h>
#endif

/*
 * Constructor:
 * User must specify kmer size used in sliding window and give it a list of
 * filenames with corresponding headers to make filter from.
 */
BloomFilterGenerator::BloomFilterGenerator(vector<string> const &filenames,
		unsigned kmerSize, unsigned hashNum) :
		m_kmerSize(kmerSize), m_hashNum(hashNum), m_expectedEntries(0), m_filterSize(
				0), m_totalEntries(0)
{

	//for each file loop over all headers and obtain max number of elements
	for (vector<string>::const_iterator i = filenames.begin();
			i != filenames.end(); ++i)
	{
		WindowedFileParser parser(*i, kmerSize);
		m_fileNamesAndHeaders[*i] = parser.getHeaders();
		for (vector<string>::iterator j = m_fileNamesAndHeaders[*i].begin();
				j != m_fileNamesAndHeaders[*i].end(); ++j)
		{
			//subtract kmer size for max number of possible kmers
			m_expectedEntries += parser.getSequenceSize(*j) - kmerSize + 1;
		}
	}
}

/*
 * Constructor:
 * User must specify kmer size used in sliding window and give it a list of
 * filenames with corresponding headers to make filter from.
 * Variant allows users to set a specific filter size before hand
 */
BloomFilterGenerator::BloomFilterGenerator(vector<string> const &filenames,
		unsigned kmerSize, unsigned hashNum, size_t numElements) :
		m_kmerSize(kmerSize), m_hashNum(hashNum), m_expectedEntries(
				numElements), m_filterSize(0), m_totalEntries(0)
{
	//for each file loop over all headers and obtain max number of elements
	for (vector<string>::const_iterator i = filenames.begin();
			i != filenames.end(); ++i)
	{
		WindowedFileParser parser(*i, kmerSize);
		m_fileNamesAndHeaders[*i] = parser.getHeaders();
	}
}

/*
 * Generates a bloom filter outputting it to a filename
 * Returns the redundancy rate of a Bloom Filter generated from a file.
 * Currently only supports fasta files.
 *
 * Outputs to fileName path
 */
size_t BloomFilterGenerator::generate(const string &filename)
{

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	//redundancy metric value
	size_t redundancy = 0;

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	for (boost::unordered_map<string, vector<string> >::iterator i =
			m_fileNamesAndHeaders.begin(); i != m_fileNamesAndHeaders.end();
			++i)
	{
		//let user know that files are being read
		cerr << "Processing File: " << i->first << endl;
		//object to process reads
		WindowedFileParser parser(i->first, m_kmerSize);
		for (vector<string>::iterator j = i->second.begin();
				j != i->second.end(); ++j)
		{
			parser.setLocationByHeader(*j);
			//insert elements into filter
			//read fasta file line by line and split using sliding window
#pragma omp parallel
			while (true) {

				ReadsProcessor proc(m_kmerSize);
				string currentSeq;
#pragma omp critical(currentSeq)
				{
					currentSeq = parser.getNextSeq();
				}

				if (currentSeq == "") {
					break;
				}

				size_t currentLoc = 0;
				while (currentSeq.length() >= currentLoc + m_kmerSize) {
					const unsigned char* currentKmer = proc.prepSeq(currentSeq,
							currentLoc);

					//check to see if string is invalid
					if (currentKmer != NULL) {
						const vector<size_t> &tempHash = multiHash(currentKmer,
								m_hashNum, m_kmerSize);
						if (filter.contains(tempHash)) {
#pragma omp atomic
							++redundancy;
						} else {
							filter.insert(tempHash);
#pragma omp atomic
							++m_totalEntries;
						}
					}
					++currentLoc;
				}
			}
		}
	}
	filter.storeFilter(filename);
	return redundancy;
}

/*
 * Generates a bloom filter outputting it to a filename
 * Input a filename to use as a subtractive filter
 * Returns the redundancy rate of a Bloom Filter generated from a file.
 * Currently only supports fasta files.
 *
 * Outputs to fileName path
 */
//TODO refactor to remove boilerplate-ness to method above
size_t BloomFilterGenerator::generate(const string &filename,
		const string &subtractFilter)
{

	//need the number of hash functions used to be greater than 0
	assert(m_hashNum > 0);

	//need the filter to be greater than the size of the number of expected entries
	assert(m_filterSize > m_expectedEntries);

	//setup bloom filter
	BloomFilter filter(m_filterSize, m_hashNum, m_kmerSize);

	//load other bloom filter info
	string infoFileName = (subtractFilter).substr(0,
			(subtractFilter).length() - 2) + "txt";
	BloomFilterInfo subInfo(infoFileName);

	//load other bloomfilter
	BloomFilter filterSub(subInfo.getCalcuatedFilterSize(),
			subInfo.getHashNum(), subInfo.getKmerSize(), subtractFilter);

	//Check if kmer or subkmers are located in filter
	if (subInfo.getKmerSize() != m_kmerSize) {
		//TODO make compatable with smaller kmer sizes
		cerr << "ERROR: Must use identical size k-mers in subtractive filter"
				<< endl;
	}

	if (subInfo.getKmerSize() > m_kmerSize) {
		cerr
				<< "Error: Subtraction filter's k-mer size is larger than output filter's k-mer size."
				<< endl;
		exit(1);
	}

	//ReadProcessor for subtraction filter
	ReadsProcessor subProc(subInfo.getKmerSize());

	//redundancy metric value
	size_t redundancy = 0;

	size_t kmerRemoved = 0;

	//for each file loop over all headers and obtain seq
	//load input file + make filter
	for (boost::unordered_map<string, vector<string> >::iterator i =
			m_fileNamesAndHeaders.begin(); i != m_fileNamesAndHeaders.end();
			++i)
	{
		//let user know that files are being read
		cerr << "Processing File: " << i->first << endl;
		WindowedFileParser parser(i->first, m_kmerSize);
		for (vector<string>::iterator j = i->second.begin();
				j != i->second.end(); ++j)
		{
			parser.setLocationByHeader(*j);
			//object to process reads
			//insert elements into filter
			//read fasta file line by line and split using sliding window
#pragma omp parallel
			while (true) {
				ReadsProcessor proc(m_kmerSize);
				string currentSeq;
#pragma omp critical(currentSeq)
				{
					currentSeq = parser.getNextSeq();
				}

				if (currentSeq == "") {
					break;
				}

				size_t currentLoc = 0;
				while (currentSeq.length() >= currentLoc + m_kmerSize) {
					const unsigned char* currentKmer = proc.prepSeq(currentSeq,
							currentLoc);
					//check to see if string is invalid
					if (currentKmer != NULL) {
						bool allowKmer = !filterSub.contains(currentKmer);

						if (allowKmer) {
							const vector<size_t> &tempHash = multiHash(
									currentKmer, m_hashNum, m_kmerSize);
							if (filter.contains(tempHash)) {
#pragma omp atomic
								++redundancy;
							} else {
								filter.insert(tempHash);
#pragma omp atomic
								++m_totalEntries;
							}
						} else {
							++kmerRemoved;
						}
					}
					++currentLoc;
				}
			}
		}
	}

	cerr << "Total Number of K-mers not added: " << kmerRemoved << endl;

	filter.storeFilter(filename);
	return redundancy;
}

//setters
void BloomFilterGenerator::setFilterSize(size_t bits)
{
	m_filterSize = bits;
}

//getters

/*
 * Returns the total number of inserted filter entries
 */
size_t BloomFilterGenerator::getTotalEntries() const
{
	return m_totalEntries;
}

/*
 * Returns the maximum possible number of expected filter entries based on inputs
 */
size_t BloomFilterGenerator::getExpectedEntries() const
{
	return m_expectedEntries;
}

//destructor
BloomFilterGenerator::~BloomFilterGenerator()
{
}
