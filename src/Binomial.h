
class Binomial
{
	protected:
		int n;
		double *bk;
		bool *calc;

	public:
		Binomial();
		Binomial(int n_in);
		double Choose(int k);
};

