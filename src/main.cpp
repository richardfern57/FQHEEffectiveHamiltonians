#include <iostream>
#include <chrono>
#include <sstream>
#include <math.h>
#include <fstream>

#include "Jack.h"
#include "utility.h"

using namespace std;



void write_VJ(int N, vector<string> krs, vector<int> root, string pot, string filename)
{
	Jack J(N,krs,root, filename);

	vector<double> pseudos;
	if(pot=="local")
		for(int i=0; i<=20; i++) pseudos.push_back(1./pow(5,i+1));
	if(pot=="local2")
		for(int i=0; i<=20; i++) pseudos.push_back(1./pow(10,i+1));
	else if(pot=="exp")
		for(int i=0; i<=20; i++) pseudos.push_back(1./pow(3,i+1));
	else if(pot=="exp2")
		for(int i=0; i<=20; i++) pseudos.push_back(1./pow(2,i+1));
	else if(pot=="exp4")
		for(int i=0; i<=20; i++) pseudos.push_back(1./pow(1.5,i+1));
	else if(pot=="quad")
		for(int i=0; i<=100; i++) pseudos.push_back(double(i));
	else if(pot=="invsq")
	{
		pseudos.push_back(0);
		for(int i=1; i<=100; i++) pseudos.push_back(1./double(i));
	}
	else if(pot=="inv")
	{
		pseudos.push_back(0);
		for(int i=1; i<=100; i++)
		{
			double vm = sqrt(M_PI/2.);
			for(int j=1; j<=i; j++) vm *= (double(j)-0.5)/double(j);
			pseudos.push_back(vm);
		}
	}
	else if(pot.find("pseudo") != string::npos)
	{
		int n=stoi(pot.substr(6));
		for(int i=0; i<n; i++) pseudos.push_back(0);
		pseudos.push_back(1);
	}
	else
		throw invalid_argument("Potential not an implemented option");

	J.ApplyV(pseudos);
	J.SaveVCs(pot);
}




int main(int argc, char* argv[])
{

	if(argc==6)
	{
		int N=atoi(argv[1]);
		vector<string> krs = split(argv[2],',');
		vector<int> root = intsplit(argv[3]);
		string pot = string(argv[4]);
		string filename = string(argv[5]);
		write_VJ(N,krs,root, pot, filename);
	}
	else
		cout << "./Main N krs root pot filename\n";

	//auto t0 = chrono::high_resolution_clock::now();
	//auto t1 = chrono::high_resolution_clock::now();
	//auto dt = chrono::duration<double>(t1-t0);
	//cout << "This took " << dt.count() << " seconds.\n";

	return 0;

}

