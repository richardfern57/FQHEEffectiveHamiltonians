#include "Binomial.h"


//This class calculates and stores binomial coefficients as required and then outputs when asked.


Binomial::Binomial(){}


Binomial::Binomial(int n_in)
{
	//We start by initialising our storage containers. calc stores if we know the value asked
	//for yet or if we must calculate.
	n = n_in;
	bk = new double[1+n/2];
	calc = new bool[1+n/2];
	for(int i=0; i<=n/2; i++)
		calc[i] = false;
}


double Binomial::Choose(int k)
{
	//If known the coefficient is returned. Otherwise, it is calculated.
	double out;
	if(k>n) out = 0;
	else
	{
		if(k>n/2) k=n-k;
		if(calc[k]) out = bk[k];
		else if(not calc[k])
		{
			out = 1;
			int ni=n, ki=k;
			while(ki>0)
			{
				out = (out*ni)/ki;
				ni--; ki--;
			}
			bk[k] = out;
			calc[k] = true;
		}
	}
	return out;
}

