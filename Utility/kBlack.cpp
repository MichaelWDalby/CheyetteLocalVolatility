#include "kBlack.h"
#include "kSolver.h"
#include "kFd1d.h"

class kBlackObj : public kSolverObjective
{
public:

	//	constructor
	kBlackObj(
		double	expiry,
		double	strike,
		double	price,
		double	forward)
		: kSolverObjective(),
		myExpiry(expiry),
		myStrike(strike),
		myPrice(price),
		myForward(forward)
	{}

	//	value
	virtual double	value(double x)
	{
		double res = kBlack::call(myExpiry, myStrike, myForward, x) - myPrice;

		//	done
		return res;
	}

	//	deriv
	virtual double	deriv(double x)
	{
		double res = kBlack::vega(myExpiry, myStrike, myForward, x);

		//	done
		return res;
	}

	//	private parts
private:

	//	expiry
	double	myExpiry;
	double	myStrike;
	double	myPrice;
	double	myForward;

};

// implied vol
double
kBlack::implied(
	double	expiry,
	double	strike,
	double	price,
	double	forward)
{
	//	calc intrinsic
	double intrinc = max(forward - strike,0.0);
	if (price <= intrinc) return 0.0;

	//	objective
	kBlackObj obj(expiry, strike, price, forward);

	//	start guess
	double volatility = 0.1;
	int    numIter = 10;
	double epsilon = (price - intrinc) * kConstants::epsilon();

	//	solve
	kSolver::newtonRapson(obj, volatility, numIter, epsilon, nullptr);

	//	bound
	volatility = max(0.0, volatility);

	//	done
	return volatility;
}

//	fd runner
bool
kBlack::fdRunner(
	const double		s0,
	const double		r,
	const double		mu,
	const double		sigma,
	const double		expiry,
	const double		strike,
	const double		barrier,
	const bool			dig,
	const int			pc,			//	put (-1) call (1)
	const int			ea,			//	european (0), american - price (1), american - excersice (2), american - exercise gride(3) 
	const int			bo,			//  knock out (2) knock out naive (1) non (0)
	const int			smooth,		//	smoothing
	const double		theta,
	const int			wind,
	const double		numStd,
	const int			numT,
	const int			numS,
	const bool			update,
	const int			numPr,
	kVector<kVector<double>>& grid,
	double& res0,
	kVector<double>& s,
	kVector<double>& res,
	string& error)
{
	//	helps
	int h, i, p;

	//	construct s axis
	double t = max(0.0, expiry);
	double std = sigma * sqrt(t);
	double sl = s0 * exp(-numStd * std);
	double su = s0 * exp(numStd * std);
	int    nums = 2 * (max(1, numS) / 2);
	if (numS <= 0 || sl == su)
	{
		nums = 1;
	}
	else
	{
		++nums;
	}
	// double ds = (su - sl) / nums;
	double dx = numStd * std * 2 / (nums-1);
	s.resize(nums);
	for (i = 0; i < nums; ++i)
	{
		s(i) = s0 * exp((i - nums / 2) *dx);
	}

	// Find barrier index
	int bidx;
	if (bo > 0) // Knock out
	{
		for (i = 0; i < nums; i++)
		{
			if (s(i) > barrier) break;
		}
		bidx = i - 1;
	}

	//	construct fd grid
	kFd1d<double> fd;
	fd.init(1, s, false);

	//	set terminal result
	double xl, xu;
	res.resize(nums);
	for (i = 0; i < nums; ++i)
	{
		if (smooth == 0 || i == 0 || i == nums - 1)
		{
			if (dig) res(i) = 0.5 * (kInlines::sign(s(i) - strike) + 1.0);
			else    res(i) = max(0.0, s(i) - strike);

			if (bo > 0 && i<=bidx) res(i) = 0;
		}
		else
		{
			xl = 0.5 * (s(i - 1) + s(i));
			xu = 0.5 * (s(i) + s(i + 1));
			if (dig) res(i) = kFiniteDifference::smoothDigital(xl, xu, strike);
			else	 res(i) = kFiniteDifference::smoothCall(xl, xu, strike);
		}

		if (pc < 0)
		{
			if (dig) res(i) = 1.0 - res(i);
			else    res(i) -= (s(i) - strike);
		}
	}

	//	time steps
	int    numt = max(0, numT);
	double dt = t / max(1, numt);

	// excercise boundary 1d
	kVector<double> boundary;
	if (ea == 2)
	{
		boundary.resize(nums);
		for (i = 0; i < nums; i++) boundary(i) = expiry + 1; // expiry +1 is inf

	}
	// excercis grid
	if (ea == 3)
	{
		grid.resize(numt);
		for (int k = 0; k < numt; ++k) {
			grid[k].resize(nums);
		}
	}

	//	repeat
	int nump = max(1, numPr);
	for (p = 0; p < nump; ++p)
	{
		//	set parameters
		if (bo == 2) {
			for (i = bidx+1; i < nums; ++i)
			{
				fd.r()(i) = r;
				fd.mu()(i) = mu * s(i);
				fd.var()(i) = sigma * sigma * s(i) * s(i);
			}
			for (i = 0; i <= bidx; ++i)
			{
				fd.r()(i) = 0;
				fd.mu()(i) = 0;
				fd.var()(i) = 0;
			}
		}
		else
		{
			for (i = 0; i < nums; ++i)
			{
				fd.r()(i) = r;
				fd.mu()(i) = mu * s(i);
				fd.var()(i) = sigma * sigma * s(i) * s(i);
			}
		}
		//	roll
		fd.res()(0) = res;
		for (h = numt - 1; h >= 0; --h)
		{
			fd.rollBwd(dt, update || h == (numt - 1), theta, wind, fd.res());
			if (ea == 1)
			{
				for (i = 0; i < nums; ++i)
				{
					if (res(i) > fd.res()(0)(i) && fd.res()(0)(i) > 0)
					{
						fd.res()(0)(i) = res(i);
					}
				}
			}
			if (ea == 2)
			{
				for (i = 0; i < nums; ++i)
				{	
					if (res(i) > fd.res()(0)(i) && fd.res()(0)(i) > 0)
					{
						fd.res()(0)(i) = res(i);
						boundary(i) = expiry - (numt-h-1)*dt;
					}
				}
			}
			if (ea == 3)
			{
				for (i = 0; i < nums; ++i)
				{
					if (res(i) > fd.res()(0)(i) && fd.res()(0)(i) > 0)
					{
						fd.res()(0)(i) = res(i);
						grid(h)(i) = 1;
					}
				}
			}
			if (bo == 1) // Knock out, naive
			{
				for (i = 0; i <= bidx; i++)
				{
					fd.res()(0)(i) = 0;
				}
			}
		}
	}

	if (ea == 2)
	{
		fd.res()(0) = boundary;
	}

	//	set result
	res = fd.res()(0);
	res0 = fd.res()(0)(nums / 2);


	//	done
	return true;
}

//	fd Fwd runner
bool
kBlack::fdFwdRunner(
	const double		s0,
	const double		r,
	const double		mu,
	const double		sigma,
	const double		expiry,
	const double		strike,
	const double		barrier,
	const bool			dig,
	const int			pc,			//	put (-1) call (1)
	const int			ea,			//	european (0), american (1)
	const int			bo,			//  knock in (-1) knock out (1) non (0)
	const int			smooth,		//	smoothing
	const double		theta,
	const int			wind,
	const double		numStd,
	const int			numT,
	const int			numS,
	const bool			update,
	const int			numPr,
	double& res0,
	kVector<double>& s,
	kVector<double>& strikes,
	const double timeFactor,
	kVector<kVector<double>>& res,
	string& error)
{
	kVector<kVector<double>> pmat;

	
	//	helps
	int h, i, p;

	//	construct s axis
	double t = max(0.0, expiry);
	double std = sigma * sqrt(t);
	double sl = s0 * exp(-numStd * std);
	double su = s0 * exp(numStd * std);
	int    nums = 2 * (max(1, numS) / 2);
	int	   numk = strikes.size();
	if (numS <= 0 || sl == su)
	{
		nums = 1;
	}
	else
	{
		++nums;
	}
	// double ds = (su - sl) / nums;
	double dx = numStd * std * 2 / (nums-1);
	s.resize(nums);
	for (i = 0; i < nums; ++i)
	{
		s(i) = s0 * exp((i - nums / 2) * dx);
	}

	//	time steps
	int    numt = max(0, numT);
	double dt = t / max(1, numt);
	numt++;

	// Create price matrix
	pmat.resize(numt);
	res.resize(numt);
	for (int k = 0; k < numt; ++k) {
		pmat[k].resize(numk);
		res[k].resize(nums);
	}

	//	set initial result
	for (i = 0; i < nums; ++i) res(0)(i) = 0;
	res(0)(int(nums / 2)) = 1;

	

	//	construct fd grid
	kFd1d<double> fd;
	fd.init(1, s, false);

	//	repeat
	int nump = max(1, numPr);
	for (p = 0; p < nump; ++p)
	{
		//	set parameters
		for (i = 0; i < nums; ++i)
		{
			fd.r()(i) = r;
			fd.mu()(i) = mu * s(i);
			fd.var()(i) = sigma * sigma * s(i) * s(i);
		}

		//	roll
		fd.res()(0) = res(0);
		for (h = 1; h <= numt - 1; ++h)
		{
			fd.rollFwd(dt, update || h == 1, theta, wind, fd.res());
			res(h) = fd.res()(0);
		}
	}
	int j;
	// calc price matrix
	for (h = 0; h < numt; ++h)
	{
		for (i = 0; i < numk; ++i)
		{	
			pmat(h)(i) = 0;
			for (j = nums -1; j>=0; --j)
			{
				if (strikes(i) > s(j)) break;
				pmat(h)(i) += (s(j) - strikes(i)) * res(h)(j);
			}
		}
	}

	//	set result on time grid specified by timeFactor
	res = pmat;

	//	done
	return true;
}