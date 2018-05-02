#include <vector>
#include <string>


std::vector<std::string> split(std::string in, char delim);


std::vector<int> intsplit(std::string in, char delim);
std::vector<int> intsplit(std::string in);


std::vector<int> unique(std::vector<int> in);


std::string stringefy(std::vector<int> in);
std::string stringefy(std::vector<int> in, std::string delim);
std::string stringefy(std::vector<double> in);


int get_index(std::vector<int> in, int item);
