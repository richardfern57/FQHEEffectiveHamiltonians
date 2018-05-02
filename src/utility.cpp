#include "utility.h"
#include <sstream>
#include <algorithm>
#include <cstring>


//These are some utility functions for dealing with strings and vectors of string/integers.


std::vector<std::string> split(std::string in, char delim)
{
	//Given some input string this splits them by some delimiter.
	std::vector<std::string> out;
	std::stringstream ss(in);
	std::string item;
	while(getline(ss,item,delim))
		out.push_back(item);
	return out;
}


std::vector<int> intsplit(std::string in, char delim){
	//Given some input string we split by a delimiter, converting all parts to integers.
	std::vector<int> out;
	std::stringstream ss(in);
	std::string item;
	while(getline(ss,item,delim))
		out.push_back(atoi(item.c_str()));
	return out;
}


std::vector<int> intsplit(std::string in){
	//Given some input string we split into characters, converting all parts to integers.
	std::vector<int> out;
	for(int i=0; i<in.length(); i++)
		out.push_back(in[i]-'0');
	return out;
}


std::vector<int> unique(std::vector<int> in)
{
	//This produces the vector of unique elements of the input vector.
	std::vector<int> out=in;
	out.erase(std::unique(out.begin(),out.end()),out.end());
	return out;
}


std::string stringefy(std::vector<int> in)
{
	//This turns a vector of integers into a string separated by commas.
	std::stringstream ss;
	int n = in.size();
	if(n>0) ss << in[0];
	for(int i=1; i<n; i++)
		ss << "," << in[i];
	return ss.str();
}


std::string stringefy(std::vector<int> in, std::string delim)
{
	//This turns a vector of integers into a string separated by commas.
	std::stringstream ss;
	int n = in.size();
	if(n>0) ss << in[0];
	for(int i=1; i<n; i++)
		ss << delim << in[i];
	return ss.str();
}


std::string stringefy(std::vector<double> in)
{
	//This turns a vector of integers into a string separated by commas.
	std::stringstream ss;
	int n = in.size();
	if(n>0) ss << in[0];
	for(int i=1; i<n; i++)
		ss << "," << in[i];
	return ss.str();
}


int get_index(std::vector<int> in, int item)
{
	//This returns the index of the first element equal to or less than item.
	int i=0;
	for(int vi:in)
	{
		if(vi<=item) break;
		i++;
	}
	return i;
}

