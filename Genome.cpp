#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
using namespace std;

class GenomeImpl
{
public:
    GenomeImpl(const string& nm, const string& sequence);
    static bool load(istream& genomeSource, vector<Genome>& genomes);
    int length() const;
    string name() const;
    bool extract(int position, int length, string& fragment) const;
private:
	string m_name;
	vector<char> m_sequence;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence) :
	m_name(nm)
{
	for (int i = 0; i < sequence.length(); i++)
	{
		m_sequence.push_back(sequence[i]);
	}
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes) 
{
	char c;
	string name = "";
	string sequence = "";
	char prev = '\0';

	while (genomeSource.get(c))
	{
		if (name == "" && c != '>')
			return false;

		// Empty line
		if (c == '\n' && prev == '\n')
			return false;
		
		if (c == '\n')
			continue;

		if (c == '>')
		{
			if (sequence != "")
			{
				Genome genome(name, sequence);
				genomes.push_back(genome);

				sequence = "";
			}

			getline(genomeSource, name);
			if (name == "")
				return false;
		}
		else 
		{
			char upper_c = toupper(c);
			switch (upper_c)
			{
				case 'A':
				case 'T':
				case 'N':
				case 'G':
				case 'C': sequence += upper_c;
					      break;
				default: return false;
			}
		}
	}
	// Load the last genome
	if (name != "" && sequence != "")
	{
		Genome genome(name, sequence);
		genomes.push_back(genome);
	}

	return true;
}

int GenomeImpl::length() const
{
	return m_sequence.size();
}

string GenomeImpl::name() const
{
	return m_name;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
	if (position + length > m_sequence.size())
		return false;

	string copy;

	for (int i = position; i < position + length; i++)
	{
		copy += m_sequence[i];
	}

	fragment = copy;
	return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
    m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
    delete m_impl;
}

Genome::Genome(const Genome& other)
{
    m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
    GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
    delete m_impl;
    m_impl = newImpl;
    return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes) 
{
    return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
    return m_impl->length();
}

string Genome::name() const
{
    return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
    return m_impl->extract(position, length, fragment);
}
