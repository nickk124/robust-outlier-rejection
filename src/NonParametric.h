/*
 Robust Chauvenet Rejection (RCR) Official Codebase
 Active Author: Nick C. Konz
 Former Author: Michael Maples
 See license at https://github.com/nickk124/RCR
 */
#include <vector>
class NonParametric
{
public:
	NonParametric();

	virtual void muFunc(std::vector<bool>&, std::vector<int>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
	virtual void muFunc(std::vector<bool>&, std::vector<int>&, std::vector<double>&, std::vector<double>&);
	std::vector<int> indices;

	~NonParametric();
private: 
	
};