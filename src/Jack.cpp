#include "Jack.h"
#include "utility.h"
#include <fstream>
#include <iostream>
#include <cstring>
#include <math.h>
#include <algorithm>
#include <iomanip>
using namespace std;


//This class holds the Jack polynomial. It stores the monomials and coefficients and can apply the
//potential, producing a second set of coefficients.


Jack::Jack(int N, vector<string> krs, vector<int> root, string filename)
{
	//This initialiser reads the Jack from a file. The input file must have lines with format;
	//  coeffient [l1, l2, ..., lfinal]
	//It assumes that any lN=0 is absent and assumes that the very first monomial includes a 0.
	//Furthermore, whilst this routine does not require it, the entire structure of these
	//methods require the partitions to be ordered lexicographically in the file and include
	//every single partition lexicographically inferior to some root with the same total
	//angular momentum and particle number. The monomials are stored in a Lookup table which
	//gives the indices at which that monomial's coefficients are in the vector. We also find
	//the total Hilbert space dimension and LzMax (as found from the first partition).
	this->dim = 0;
	this->LzMax = root.size();
	this->s = 0;
	if(krs[2]=="f") {this->s = 1;}
	this->file_name = filename;

	ifstream file;
	file.open(file_name);
	vector<string> data;
	string line;
	int Ni;
	while(getline(file,line))
	{
		data = split(line,' ');
		Ni = data.size();
		this->Coefs.push_back(stod(data[0]));
		line = "";
		for(int i=1; i<Ni; i++) line += data[i];
		line.erase(line.end()-1);
		line.erase(line.begin());
		if(Ni==N) line += ",0";
		this->Lookup[line] = dim;
		dim++;
	}
	cout << "Read " << file_name << ".\n";
}


void Jack::Kill()
{
	//Clears the data.
	this->Coefs.clear();
	this->VCoefs.clear();
	this->Lookup.clear();
}


void Jack::ApplyV(vector<double> pseudos)
{
	//This applies the potential. We first initialise the interactions and target vector. We
	//then loop over each monomial and apply all possible AdAdAA operators to them. We start
	//by looping over all l1=0,...,LzMax and then all l2<=l1. We then only consider 
	//"squeezing" operators (and their hermitian conjugates) which constitutes all AdAdAA's
	//which keep us in the same Hilbert space in which we began. There's a lot of playing
	//around with vectors of integers at various parts to identify which l's we should be
	//summing over and what the index of the partition is after the operator is applied. A
	//crucial note; because we restrict our operation to our initial Hilbert space only the
	//target vector is not the full application of V. It is instead V*Jack projected back into
	//our original Hilbert space. As such, we can only consider dot products of V*Jack with
	//vectors in the same or smaller Hilbert spaces (for the same N and total angular momentum).
	this->V = Interactions(pseudos, this->LzMax, this->s);
	this->VCoefs = vector<double>(this->dim);

	vector<int> part, part1, u_part, u_part1, aapart, new_part;
	int L, l3, index, new_index;
	double cl1,cl2,cl3,cl4, Vllll;
	int pos1,pos2,pos3,pos4;
	if(this->s) //Fermions
	{
		for(auto mon:Lookup)
		{
			index = mon.second;
			part = intsplit(mon.first,',');
			for(int l1:part)
			{
				pos1 = get_index(part,l1);
				part1 = vector<int>(part.begin()+pos1+1, part.end());
				for(int l2:part1)
				{
					pos2 = get_index(part,l2);
					aapart = part;
					aapart.erase(aapart.begin()+pos2);
					aapart.erase(aapart.begin()+pos1);
					L = l1+l2;
					for(int l4=l2; l4<(L+1)/2; l4++)
					{
						cl4 = count(aapart.begin(),aapart.end(),l4);
						l3 = L-l4;
						cl3 = count(aapart.begin(),aapart.end(),l3);
						if(cl3+cl4==0)
						{
							new_part = aapart;
							pos3 = get_index(new_part,l3);
							new_part.insert(new_part.begin()+pos3,l3);
							pos4 = get_index(new_part,l4);
							new_part.insert(new_part.begin()+pos4,l4);
							new_index = Lookup[stringefy(new_part)];
	
							if((pos1+pos2+pos3+pos4)%2==1)
								Vllll = V.Get(l1,l2,l4);
							else
								Vllll = -V.Get(l1,l2,l4);
							VCoefs[new_index] += Coefs[index]*Vllll;
						}
					}
				}
			}
		}
	}
	else //Bosons
	{
		for(auto mon:Lookup)
		{
			index = mon.second;
			part = intsplit(mon.first,',');
			u_part = unique(part);
			for(int l1:u_part)
			{
				cl1 = count(part.begin(),part.end(),l1);
				pos1 = get_index(part,l1);
				part1 = vector<int>(part.begin()+pos1+1, part.end());
				u_part1 = unique(part1);
				for(int l2:u_part1)
				{
					cl2 = count(part1.begin(),part1.end(),l2);
					pos2 = get_index(part,l2);
					aapart = part;
					aapart.erase(aapart.begin()+pos2);
					aapart.erase(aapart.begin()+pos1);
					L = l1+l2;
					for(int l4=l2; l4<=L/2; l4++)
					{
						l3 = L-l4;
						new_part = aapart;
						pos3 = get_index(new_part,l3);
						new_part.insert(new_part.begin()+pos3,l3);
						cl3 = count(new_part.begin(),new_part.end(),l3);
						pos4 = get_index(new_part,l4);
						new_part.insert(new_part.begin()+pos4,l4);
						cl4 = count(new_part.begin(),new_part.end(),l4);
						new_index = Lookup[stringefy(new_part)];

						Vllll = V.Get(l1,l2,l4)*sqrt(cl1*cl2*cl3*cl4);
						VCoefs[new_index] += Coefs[index]*Vllll;
					}
				}
			}
		}
	}
	cout << "Potential applied.\n";
}


void Jack::SaveVCs(string pot)
{
	//This saves the coefficients of V*Jack to a file specified by the potential's name.
	ofstream file;
	file.open(this->file_name+pot);
	file << setprecision(14);
	for(double c : this->VCoefs) file << c << "\n";
	file.close();
}


double Jack::DotWith(vector<double> other)
{
	//This calculates the dot product of our raw Jack vector with some other vector. It assumes
	//that N and the total angular momentum are the same (and so the Hilbert spaces are the same
	//up to some cutoff) and that the order of elements in the vectors are lexicographic (and
	//so we can match up coefficients by index rather than trawling through lookup tables).
	//Given these assumptions we match coefficients from last element of both vectors (which
	//would be the same monomial) and then we match all the way backwards to whichever cutoff
	//of the two vectors is smallest. The monomials are assumed normalised.
	double out = 0;
	int dim_o = other.size(), dif;
	if(dim>dim_o)
	{
		dif = dim-dim_o;
		for(int i=0; i<dim_o; i++) out += Coefs[i+dif]*other[i];
	} 
	else
	{
		dif = dim_o-dim;
		for(int i=0; i<dim; i++) out += Coefs[i]*other[i+dif];
	}
	return out;
}




