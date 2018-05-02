#include "Interactions.h"
#include <math.h>
#include <iostream>

//This class gives the matrix elements for the interaction V = sum Vl4l3l2l1 Adl4 Adl3 Al2 Al1
//arising from pseudopotentials.


double Interactions::Calculate(int m, int l1, int l2, int l4)
{
	//This calculates Vl1l2l3l4 including no counting factors for the m'th pseudopotential.
	double Vllll, B12=0, B43=0, sign=1;
	long double temp=0;
	int L=l1+l2, l3=L-l4;
	if(L<m) Vllll = 0;
	else
	{
		temp = F.Get(m)/pow(2,L+1);
		temp /= sqrt(F.Get(l1));
		temp /= sqrt(F.Get(l2));
		temp /= sqrt(F.Get(l3));
		temp /= sqrt(F.Get(l4));
		temp *= F.Get(L-m);
		Vllll = (double) temp;
	//	Vllll = sqrt(B[L].Choose(l2)*B[L].Choose(l4))/B[L].Choose(m);
	//	Vllll /= pow(2,L+1);
		for(int k=0; k<=m; k++)
		{
			B12 += sign*B[l1].Choose(k)*B[l2].Choose(m-k);
			B43 += sign*B[l4].Choose(k)*B[l3].Choose(m-k);
			sign *= -1;
		}
		Vllll *= B12*B43;
	}
	return Vllll;
}


Interactions::Interactions(){}


Interactions::Interactions(std::vector<double> pseudos, int LzMax, int s)
{
	//This calculates all necessary matrix elements given some maximum cutoff. We start by
	//calculating the Binomial coefficients. We then loop over all l1>=l2, creating the array
	//for V en route. We then loop over l4 such that l2<=l4<=l3<=l1 and calculate the
	//coefficients, including counting factors.
	B = new Binomial[2*LzMax+1];
	F = Factorials(2*LzMax+1);
	for(int i=0; i<=2*LzMax; i++)
		B[i] = Binomial(i);

	int l4Max, p_dim=pseudos.size();
	double fac;
	V = new double**[LzMax+1];
	for(int l1=0; l1<=LzMax; l1++)
	{
		V[l1] = new double*[l1+1];
		for(int l2=0; l2<=l1; l2++)
		{
			l4Max = (l1+l2)/2;
			V[l1][l2] = new double[l4Max+1];
			for(int l4=l2; l4<=l4Max; l4++)
			{
				fac = 4;
				if(l1==l2) fac /= 2;
				if(2*l4==l1+l2) fac /= 2;
				if(l4==l2) fac /= 2;
				V[l1][l2][l4] = 0;
				for(int i=s; i<p_dim; i+=2)
					if(pseudos[i]!=0)
						V[l1][l2][l4] += pseudos[i]*fac*
							Calculate(i,l1,l2,l4);
			}
		}
	}
}


double Interactions::Get(int l1, int l2, int l4)
{
	//Returns the matrix element. Note that not all V's are calculated so use only for
	//l2<=l4<=l3<=l1.
	return V[l1][l2][l4];
}

