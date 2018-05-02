#include "Interactions.h"
#include <unordered_map>
#include <vector>
#include <string>


class Jack
{
	protected:
		int dim, LzMax, s;
		std::string file_name;
		Interactions V;
		std::unordered_map<std::string,int> Lookup;

	public:
		std::vector<double> Coefs, VCoefs;
		Jack(int N, std::vector<std::string> krs, std::vector<int> root, std::string filename);
		void Kill();
		void ApplyV(std::vector<double> pseudos);
		void SaveVCs(std::string pot);
		double DotWith(std::vector<double> other);
};


