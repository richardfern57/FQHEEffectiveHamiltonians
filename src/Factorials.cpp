#include "Factorials.h"
#include <limits>


//This class calculates and stores binomial coefficients as required and then outputs when asked.


Factorials::Factorials(){}


Factorials::Factorials(int n_max)
{
	//We start by initialising our storage containers. calc stores if we know the value asked
	//for yet or if we must calculate.
	n = n_max;
	bk = new long double[n+1];
	calc = new bool[n+1];
	for(int i=0; i<=n; i++)
		calc[i] = false;
}


long double Factorials::Get(int k)
{
	//If known the coefficient is returned. Otherwise, it is calculated.
	long double out;
	if(k>n) out = std::numeric_limits<long double>::infinity();
	else
	{
		if(calc[k]) {out = bk[k];}
		else {
			out = 1;
			for (int j=2; j<=k; j++) {out = out*j;}
			bk[k] = out;
			calc[k] = true;
		}
	}
	return out;
}

