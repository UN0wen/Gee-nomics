#include "provided.h"
#include "Trie.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unordered_map> 
using namespace std;

class GenomeMatcherImpl
{
public:
    GenomeMatcherImpl(int minSearchLength);
    void addGenome(const Genome& genome);
    int minimumSearchLength() const;
    bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
    bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
	int m_minSearchLength;
	vector<Genome> m_genomes;
	struct GenomeLoc
	{
		GenomeLoc(int i, int pos) : index(i), genomePos(pos) {}
		int index;
		int genomePos;
	};
	Trie<GenomeLoc> m_trie;
	vector<DNAMatch> getStrings(vector<GenomeLoc>& matches, int minLength, int maxLength, vector<string>& sequences) const;
	bool checkSnip(string sequence, string fragment, int minimumLength, int& locSnip, bool exactMatchOnly) const;
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength) :
	m_minSearchLength(minSearchLength)
{
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
	m_genomes.push_back(genome);

	string subsequence;
	int lastIndex = genome.length() - m_minSearchLength;

	int pos = m_genomes.size() - 1;

	for (int i = 0; i < lastIndex; i++)
	{
		bool extractSuccess = genome.extract(i, m_minSearchLength, subsequence);
		if (extractSuccess)
			m_trie.insert(subsequence, GenomeLoc(i, pos));
	}
}

int GenomeMatcherImpl::minimumSearchLength() const
{
	return m_minSearchLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
	if ((fragment.length() < minimumLength) || (minimumLength < m_minSearchLength))
		return false;

	// first m_minsearchLength chars of fragment
	string minSearchLengthFrag = fragment.substr(0, m_minSearchLength);
	
	// vector of matches from the trie
	vector<GenomeLoc> trieMatches;

	// Find matches in the trie
	trieMatches = m_trie.find(minSearchLengthFrag, exactMatchOnly);

	vector<string> sequences;
	vector<DNAMatch> DNAMatches;

	// O(H)
	DNAMatches = getStrings(trieMatches, minimumLength, fragment.length(), sequences);

	// Maps a name to its respective closest DNAMatch value
	unordered_map<string, int> nameTable;

	// O(H) time
	for (int i = 0; i < DNAMatches.size(); i++)
	{
		int locSnip = 0;
		// O(F) time
		bool isSnip = checkSnip(sequences[i], fragment, minimumLength, locSnip, exactMatchOnly);
		if (!isSnip)
			continue;
		else if (locSnip != 0)
		// sequence is a snip but only until locSnip
			DNAMatches[i].length = locSnip;
	

		int currentIndex = nameTable[DNAMatches[i].genomeName];
		int x = 0;
		// Uninitialized 
		if (currentIndex == 0)
		{
			// Initializing
			if (DNAMatches[currentIndex].genomeName != DNAMatches[i].genomeName)
				nameTable[DNAMatches[i].genomeName] = i;
			else if (DNAMatches[currentIndex].position > DNAMatches[i].position)
				nameTable[DNAMatches[i].genomeName] = i;
		}
		// New record is closer to beginning than new
		else if (DNAMatches[currentIndex].position > DNAMatches[i].position)
		{
			nameTable[DNAMatches[i].genomeName] = i;
		}
	}

	vector<DNAMatch> results;

	for (auto it : nameTable)
	{
		results.push_back(DNAMatches[it.second]);
	}

	if (!results.empty())
	{
		matches = results;
		return true;
	}
	else
		return false;
}

// Grab strings of maxLength size from the location and genome specified in trieMatches and 
// return them in a vector of DNAMatch in no particular order
// modify sequences so that the index of DNAMatch = string index in sequences
vector<DNAMatch> GenomeMatcherImpl::getStrings(vector<GenomeLoc>& trieMatches, int minLength, int maxLength, vector<string>& sequences) const
{
	vector<DNAMatch> DNAmatches;
	vector<string> sequencestrings;

	for (int i = 0; i < trieMatches.size(); i++)
	{
		DNAMatch currentMatch;
		string sequence;
		GenomeLoc currentGenome = trieMatches[i];
		int genomeLength = m_genomes[currentGenome.genomePos].length();

		// Genome doesn't have enough remaining length
		if (currentGenome.index + minLength >= genomeLength)
			continue;
		// Genome have enough for min but not enough for max, extract everything
		else if (currentGenome.index + maxLength > genomeLength)
		{
			currentMatch.length = genomeLength - 1 - currentGenome.index;
		}
		else
		{
			currentMatch.length = maxLength;
		}
		currentMatch.position = currentGenome.index;
		currentMatch.genomeName = m_genomes[currentGenome.genomePos].name();

		// extract DNA sequence
		m_genomes[currentGenome.genomePos].extract(currentGenome.index, currentMatch.length, sequence);
		
		DNAmatches.push_back(currentMatch);
		sequencestrings.push_back(sequence);
	}

	sequences = sequencestrings;

	return DNAmatches;

}

// Returns true if sequence is a Snip of fragment
bool GenomeMatcherImpl::checkSnip(string sequence, string fragment, int minimumLength, int& locSnip, bool exactMatchOnly) const
{
	// number of mismatches
	int numSnips = 0;

	for (int i = 0; i < sequence.length(); i++)
	{
		if (sequence[i] != fragment[i])
			numSnips++;

		if ((numSnips > 1) || (numSnips > 0 && exactMatchOnly))
		{
			if (i < minimumLength)
				return false;
			else
			{
				locSnip = i;
				return true;
			}
		}
	}

	locSnip = 0;
	return true;
}
bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
	if (fragmentMatchLength < m_minSearchLength)
		return false;

	int numSequences = query.length() / fragmentMatchLength;

	vector<GenomeMatch> localResults;

	unordered_map<string, int> numMatches;

	// for each sequence
	for (int i = 0; i < numSequences; i++)
	{
		// Extract sequence from query genome
		string fragment;
		query.extract(i*fragmentMatchLength, fragmentMatchLength, fragment);

		vector<DNAMatch> matches;
		// Search for sequence across library
		bool genomeFound = findGenomesWithThisDNA(fragment, fragmentMatchLength, exactMatchOnly, matches);
		
		// If a match is found in one or more genomes in the library, then for each
		// such genome, increase the count of matches found thus far for it
		if (genomeFound)
		{
			for (int j = 0; j < matches.size(); j++)
			{
				numMatches[matches[j].genomeName] += 1;
			}
		}
	}

	for (auto it : numMatches)
	{
		double p = (double(it.second) / double(numSequences)) * 100;
		if (p >= matchPercentThreshold)
		{
			GenomeMatch g;
			g.genomeName = it.first;
			g.percentMatch = p;
			localResults.push_back(g);
		}
	}

	if (!localResults.empty())
	{
		results = localResults;
		return true;
	}
	else
	{
		return false;
	}
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
    m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
    delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
    m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
    return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
    return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
    return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
