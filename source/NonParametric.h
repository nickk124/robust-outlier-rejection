//#pragma once
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

