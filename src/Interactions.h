#include "Binomial.h"
#include "Factorials.h"
#include <vector>

class Interactions
{
	protected:
		Binomial *B;
		Factorials F;
		double ***V;
		double Calculate(int m, int l1, int l2, int l4);

	public:
		Interactions();
		Interactions(std::vector<double> pseudos, int LzMax, int s);
		double Get(int l1, int l2, int l4);
};

