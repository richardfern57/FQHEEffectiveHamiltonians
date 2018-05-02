
class Factorials
{
	protected:
		int n;
		long double *bk;
		bool *calc;

	public:
		Factorials();
		Factorials(int n_max);
		long double Get(int k);
};

