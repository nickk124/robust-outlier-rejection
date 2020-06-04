//#include "stdafx.h"
#include "NonParametric.h"

NonParametric::NonParametric()
{
}

void NonParametric::muFunc(std::vector<bool> &flags, std::vector<int> &indices, std::vector<double> &w, std::vector<double> &y, std::vector<double> &trueW, std::vector<double> &trueY)
{

	int trueCount = 0, currentIndex;
	std::vector<int> indicesVec;
	std::vector<double> trueWVec, trueYVec;
	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			trueCount += 1;
		}
	}
	trueWVec.resize(trueCount);
	trueYVec.resize(trueCount);
	indicesVec.resize(trueCount);
	currentIndex = 0;
	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			trueWVec[currentIndex] = (w[i]);
			trueYVec[currentIndex] = (y[i]);
			indicesVec[currentIndex] = i;

			currentIndex += 1;
		}
	}
	trueY = trueYVec;
	trueW = trueWVec;
	indices = indicesVec;
}

void NonParametric::muFunc(std::vector<bool> &flags, std::vector<int> &indices, std::vector<double> &y, std::vector<double> &trueY)
{

	int trueCount = 0, currentIndex;
	std::vector<int> indicesVec;
	std::vector<double> trueWVec, trueYVec;
	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			trueCount += 1;
		}
	}
	trueYVec.resize(trueCount);
	indicesVec.resize(trueCount);
	currentIndex = 0;
	for (int i = 0; i < flags.size(); i++)
	{
		if (flags[i])
		{
			trueYVec[currentIndex] = (y[i]);
			indicesVec[currentIndex] = i;

			currentIndex += 1;
		}
	}
	trueY = trueYVec;
	indices = indicesVec;
}

NonParametric::~NonParametric()
{
}
