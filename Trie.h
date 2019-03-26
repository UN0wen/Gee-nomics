#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <vector>

template<typename ValueType>
class Trie
{
public:
    Trie();
    ~Trie();
    void reset();
    void insert(const std::string& key, const ValueType& value);
    std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;

      // C++11 syntax for preventing copying and assignment
    Trie(const Trie&) = delete;
    Trie& operator=(const Trie&) = delete;
private:
	struct Node 
	{
		std::vector<ValueType> data;
		Node* children[26]; // 26 letters of the alphabet
	};

	Node* m_root;

	Node* createNode(); // Initializes an empty node. O(1)
	std::vector<ValueType> findHelper(Node* root, const std::string& key, bool exactMatchOnly) const;
	void deleteNode(Node* root); // Delete all nodes starting from root, O(N)
	bool isLeaf(Node* p) const; // returns true if node has no children. O(1)
	// void addData(Node* p, ValueType data); // adds data to a node p
};

// creates an empty node and set it as root
template<typename ValueType>
Trie<ValueType>::Trie()
{
	m_root = createNode(); 
}

template<typename ValueType>
void Trie<ValueType>::reset()
{
	deleteNode(m_root);
	m_root = createNode();
}

template<typename ValueType>
Trie<ValueType>::~Trie()
{
	deleteNode(m_root);
}

template<typename ValueType>
void Trie<ValueType>::insert(const std::string& key, const ValueType& value)
{
	Node* p = m_root;

	for (int i = 0; i < key.length(); i++) 
	{
		int index = key[i] - 'A'; // index based on ASCII value
		if (p->children[index] == nullptr)
		{
			p->children[index] = createNode();
		}

		p = p->children[index];
	}

	// p points to the last node of the added key
	p->data.push_back(value);;
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string& key, bool exactMatchOnly) const
{
	// check for first char
	char firstChar = key[0];
	int firstCharIndex = firstChar - 'A';
	
	if (m_root->children[firstCharIndex] == nullptr)
		return std::vector<ValueType>(); // Empty vector

	Node* p = m_root->children[firstCharIndex];
	std::string startString = key.substr(1);

	return findHelper(p, startString, exactMatchOnly);
}

// Recursive function to find all values that are 0 or 1 characters off from the key
template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::findHelper(Node* root, const std::string& key, bool exactMatchOnly) const
{

	if (root == nullptr)
		return std::vector<ValueType>();

	// reached the end of the key, return value
	if (key == "")
		return root->data;

	// reached the end of the tree but haven't found the key, return empty vector
	if (isLeaf(root))
		return std::vector<ValueType>();

	std::vector<ValueType> results;

	char firstChar = key[0];
	int firstCharIndex = firstChar - 'A';
	std::string startString = key.substr(1);

	// Explores only the matching key with exactMatchOnly
	if (exactMatchOnly)
	{
		std::vector<ValueType> recResults = findHelper(root->children[firstCharIndex], startString, exactMatchOnly);

		if (!recResults.empty())
			results.insert(results.end(), recResults.begin(), recResults.end());
	}
	// Explores all other subtrees
	else
	{
		std::vector<ValueType> resultsExact = findHelper(root->children[firstCharIndex], startString, exactMatchOnly);

		if (!resultsExact.empty())
			results.insert(results.end(), resultsExact.begin(), resultsExact.end());

		for (int i = 0; i < 26; i++)
		{
			if (i == firstCharIndex)
				continue;
			std::vector<ValueType> recResults = findHelper(root->children[i], startString, true);
			if (!recResults.empty())
				results.insert(results.end(), recResults.begin(), recResults.end());
		}
	}

	return results;
}

template<typename ValueType>
typename Trie<ValueType>::Node* Trie<ValueType>::createNode()
{
	Node* p = new Node;

	for (int i = 0; i < 26; i++)
	{
		p->children[i] = nullptr;
	}

	return p;
}

template<typename ValueType>
bool Trie<ValueType>::isLeaf(Node* p) const
{
	for (int i = 0; i < 26; i++)
	{
		if (p->children[i] != nullptr)
			return false;
	}
	return true;
}

template<typename ValueType>
void Trie<ValueType>::deleteNode(Node* root)
{
	if (root != nullptr)
	{
		// node is a leaf, just delete it and return
		if (isLeaf(root))
		{
			delete root;
			return;
		}

		// recursively delete all nodes
		for (int i = 0; i < 26; i++)
		{
			deleteNode(root->children[i]);
		}

		// delete itself

		delete root;
	}
}

#endif // TRIE_INCLUDED
