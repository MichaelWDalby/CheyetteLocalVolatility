#include <windows.h>
#include <cmath>

#include "xlcall.h"
#include "framework.h"
#include "xlOper.h"

#include "../Utility/kMatrixAlgebra.h"
#include "../Utility/kBachelier.h"
#include "../Utility/kBlack.h"
#include "../Utility/xlUtils.h"
#include "../Utility/kFd1d.h"

#include "../Utility/kCheyette.h"


//	Wrappers

extern "C" __declspec(dllexport)
double xMultiply2Numbers(double x, double y)
{
	return x * y;
}

extern "C" __declspec(dllexport)
LPXLOPER12 xMatrixMul(LPXLOPER12 A_in, LPXLOPER12 B_in)
{
	FreeAllTempMemory();

	kMatrix<double> A;
	if (!kXlUtils::getMatrix(A_in, A))
		return TempStr12("input 1 is not a matrix");

	kMatrix<double> B;
	if (!kXlUtils::getMatrix(B_in, B))
		return TempStr12("input 2 is not a matrix");

	if (B.rows() != A.cols())
		return TempStr12("input 2 must have number of rows same as input 1 number of cols");

	kMatrix<double> res;
	kMatrixAlgebra::mmult(A, B, res);

	LPXLOPER12 out = TempXLOPER12();
	kXlUtils::setMatrix(res, out);
	return out;
}

extern "C" __declspec(dllexport)
LPXLOPER12
xBachelierCall(
	LPXLOPER12	expiry_,
	LPXLOPER12	strike_,
	LPXLOPER12	forward_,
	LPXLOPER12	volatility_)
{
	FreeAllTempMemory();

	//	helps
	string err;

	//	get expiry
	double expiry;
	if (!kXlUtils::getDbl(expiry_, 0, 0, expiry, &err)) return kXlUtils::setError(err);

	//	get strike
	double strike;
	if (!kXlUtils::getDbl(strike_, 0, 0, strike, &err)) return kXlUtils::setError(err);

	//	get forwad
	double forward;
	if (!kXlUtils::getDbl(forward_, 0, 0, forward, &err)) return kXlUtils::setError(err);

	//	get volatility
	double volatility;
	if (!kXlUtils::getDbl(volatility_, 0, 0, volatility, &err)) return kXlUtils::setError(err);

	//	calc
	double call = kBachelier::call(expiry, strike, forward, volatility);

	//	set output
	LPXLOPER12 out = kXlUtils::getOper(1, 1);
	kXlUtils::setDbl(0, 0, call, out);

	//	done
	return out;
}

extern "C" __declspec(dllexport)
LPXLOPER12
xBachelierImplied(
	LPXLOPER12	expiry_,
	LPXLOPER12	strike_,
	LPXLOPER12	price_,
	LPXLOPER12	forward_)
{
	FreeAllTempMemory();

	//	helps
	string err;

	//	get expiry
	double expiry;
	if (!kXlUtils::getDbl(expiry_, 0, 0, expiry, &err)) return kXlUtils::setError(err);

	//	get strike
	double strike;
	if (!kXlUtils::getDbl(strike_, 0, 0, strike, &err)) return kXlUtils::setError(err);

	//	get volatility
	double price;
	if (!kXlUtils::getDbl(price_, 0, 0, price, &err)) return kXlUtils::setError(err);

	//	get forward
	double forward;
	if (!kXlUtils::getDbl(forward_, 0, 0, forward, &err)) return kXlUtils::setError(err);

	//	calc
	double volatility = kBachelier::implied(expiry, strike, price, forward);

	//	set output
	LPXLOPER12 out = kXlUtils::getOper(1, 1);
	kXlUtils::setDbl(0, 0, volatility, out);

	//	done
	return out;
}

extern "C" __declspec(dllexport)
LPXLOPER12
xBlackCall(
	LPXLOPER12	expiry_,
	LPXLOPER12	strike_,
	LPXLOPER12	forward_,
	LPXLOPER12	volatility_)
{
	FreeAllTempMemory();

	//	helps
	string err;

	//	get expiry
	double expiry;
	if (!kXlUtils::getDbl(expiry_, 0, 0, expiry, &err)) return kXlUtils::setError(err);

	//	get strike
	double strike;
	if (!kXlUtils::getDbl(strike_, 0, 0, strike, &err)) return kXlUtils::setError(err);

	//	get forwad
	double forward;
	if (!kXlUtils::getDbl(forward_, 0, 0, forward, &err)) return kXlUtils::setError(err);

	//	get volatility
	double volatility;
	if (!kXlUtils::getDbl(volatility_, 0, 0, volatility, &err)) return kXlUtils::setError(err);

	//	calc
	double call = kBlack::call(expiry, strike, forward, volatility);

	//	set output
	LPXLOPER12 out = kXlUtils::getOper(1, 1);
	kXlUtils::setDbl(0, 0, call, out);

	//	done
	return out;
}

extern "C" __declspec(dllexport)
LPXLOPER12
xBlackImplied(
	LPXLOPER12	expiry_,
	LPXLOPER12	strike_,
	LPXLOPER12	price_,
	LPXLOPER12	forward_)
{
	FreeAllTempMemory();

	//	helps
	string err;

	//	get expiry
	double expiry;
	if (!kXlUtils::getDbl(expiry_, 0, 0, expiry, &err)) return kXlUtils::setError(err);

	//	get strike
	double strike;
	if (!kXlUtils::getDbl(strike_, 0, 0, strike, &err)) return kXlUtils::setError(err);

	//	get volatility
	double price;
	if (!kXlUtils::getDbl(price_, 0, 0, price, &err)) return kXlUtils::setError(err);

	//	get forward
	double forward;
	if (!kXlUtils::getDbl(forward_, 0, 0, forward, &err)) return kXlUtils::setError(err);

	//	calc
	double volatility = kBlack::implied(expiry, strike, price, forward);

	//	set output
	LPXLOPER12 out = kXlUtils::getOper(1, 1);
	kXlUtils::setDbl(0, 0, volatility, out);

	//	done
	return out;
}

extern "C" __declspec(dllexport)
LPXLOPER12
xTridag(
	LPXLOPER12 A_in,
	LPXLOPER12 b_in)
{
	FreeAllTempMemory();

	kMatrix<double> A;
	if (!kXlUtils::getMatrix(A_in, A))
		return TempStr12("input 1 is not a matrix");

	kVector<double> b;
	if (!kXlUtils::getVector(b_in, b))
		return TempStr12("input 2 is not a vector");

	if (b.size() != A.rows())
		return TempStr12("input 1 must have same number of rows as input 2");

	kVector<double> gam, res;
	kMatrixAlgebra::tridag(A, b, res, gam);

	LPXLOPER12 out = TempXLOPER12();
	kXlUtils::setVector(res, out);
	return out;
}

extern "C" __declspec(dllexport)
LPXLOPER12
xBanmul(
	LPXLOPER12	A_in,
	LPXLOPER12	x_in,
	double		m1_in,
	double		m2_in)
{
	FreeAllTempMemory();

	kMatrix<double> A;
	if (!kXlUtils::getMatrix(A_in, A))
		return TempStr12("input 1 is not a matrix");

	kVector<double> x;
	if (!kXlUtils::getVector(x_in, x))
		return TempStr12("input 2 is not a vector");

	if (x.size() != A.rows())
		return TempStr12("input 1 must have same number of rows as input 2");

	kVector<double> res;
	int m1 = (int)(m1_in + 0.5), m2 = (int)(m2_in + 0.5);
	kMatrixAlgebra::banmul(A, m1, m2, x, res);

	LPXLOPER12 out = TempXLOPER12();
	kXlUtils::setVector(res, out);
	return out;
}

extern "C" __declspec(dllexport)
LPXLOPER12
xFd1d(
	LPXLOPER12 t_in,
	LPXLOPER12 x_in,
	LPXLOPER12 r_in,
	LPXLOPER12 mu_in,
	LPXLOPER12 sigma_in,
	LPXLOPER12 v0_in,
	LPXLOPER12 tech_in)
{
	FreeAllTempMemory();

	//	help
	string err;

	//	get t 
	double t;
	if (!kXlUtils::getDbl(t_in, 0, 0, t, &err)) return kXlUtils::setError(err);

	//	get tech
	kVector<double> tech;
	if (!kXlUtils::getVector(tech_in, tech))
		return kXlUtils::setError("input 1 is not a vector");

	//	standard data
	int    numt = tech.size() > 0 ? (int)std::lround(tech(0)) : 1;
	double theta = tech.size() > 1 ? tech(1) : 0.5;
	int    fb = tech.size() > 2 ? (int)std::lround(tech(2)) : -1;
	int    log = tech.size() > 3 ? (int)std::lround(tech(3)) : 0;
	int    wind = tech.size() > 4 ? (int)std::lround(tech(4)) : 0;

	//	fd grid
	kFd1d<double> fd;
	if (!kXlUtils::getVector(x_in, fd.x()))
		return kXlUtils::setError("input 1 is not a vector");

	int n = fd.x().size();

	//	init fd
	fd.init(1, fd.x(), log > 0);

	if (!kXlUtils::getVector(r_in, fd.r()))
		return kXlUtils::setError("r is not a vector");

	if (n != fd.r().size())
		return kXlUtils::setError("r must have same size as x");

	if (!kXlUtils::getVector(mu_in, fd.mu()))
		return kXlUtils::setError("mu is not a vector");

	if (n != fd.mu().size())
		return kXlUtils::setError("mu must have same size as x");

	if (!kXlUtils::getVector(sigma_in, fd.var()))
		return kXlUtils::setError("sigma is not a vector");

	if (n != fd.var().size())
		return kXlUtils::setError("sigma must have same size as x");

	auto& fd_var = fd.var();
	for (int i = 0; i < fd_var.size(); ++i)
		fd_var(i) *= fd_var(i);

	fd.res().resize(1);
	if (!kXlUtils::getVector(v0_in, fd.res()[0]))
		return kXlUtils::setError("input 2 is not a vector");

	if (n != fd.res()[0].size())
		return kXlUtils::setError("v0 must have same size as x");

	if (theta < 0.0) theta = 0.0;
	if (theta > 1.0) theta = 1.0;

	double dt = t / numt;
	for (int n = 0; n < numt; ++n)
	{
		if (fb <= 0)
		{
			fd.rollBwd(dt, true, theta, wind, fd.res());
		}
		else
		{
			fd.rollFwd(dt, true, theta, wind, fd.res());
		}
	}

	LPXLOPER12 out = TempXLOPER12();
	kXlUtils::setVector(fd.res()[0], out);
	return out;
}

extern "C" __declspec(dllexport)
LPXLOPER12
xBachelierFd(
	LPXLOPER12	params,
	LPXLOPER12	contract,
	LPXLOPER12	gridTech)
{
	FreeAllTempMemory();

	//	help
	string err;
	int numRows, numCols;
	int i, k;

	//	get params
	double s0 = 0.0;
	double r = 0.0;
	double mu = 0.0;
	double sigma = 0.1;
	numRows = (int)getRows(params);
	if (numRows > 0 && !kXlUtils::getDbl(params, 0, 0, s0, &err))		return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getDbl(params, 1, 0, r, &err))		return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(params, 2, 0, mu, &err))		return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getDbl(params, 3, 0, sigma, &err))	return kXlUtils::setError(err);

	//	get contract
	double expiry = 0.0;
	double strike = 0.0;
	int    dig = 0;
	int    pc = 1;
	int	   ea = 0;
	int	   smooth = 0;
	numRows = (int)getRows(contract);
	if (numRows > 0 && !kXlUtils::getDbl(contract, 0, 0, expiry, &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getDbl(contract, 1, 0, strike, &err))	return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getInt(contract, 2, 0, dig, &err))		return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(contract, 3, 0, pc, &err))		return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getInt(contract, 4, 0, ea, &err))		return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getInt(contract, 5, 0, smooth, &err))	return kXlUtils::setError(err);

	//	get grid tech
	double theta = 0.5;
	int	   wind = 0;
	double numStd = 5.0;
	int    numT = 25;
	int    numX = 50;
	int    update = 1;
	int    numPr = 1;
	numRows = (int)getRows(gridTech);
	if (numRows > 0 && !kXlUtils::getDbl(gridTech, 0, 0, theta, &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getInt(gridTech, 1, 0, wind, &err))	return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(gridTech, 2, 0, numStd, &err))return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(gridTech, 3, 0, numT, &err))	return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getInt(gridTech, 4, 0, numX, &err))	return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getInt(gridTech, 5, 0, update, &err))return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getInt(gridTech, 6, 0, numPr, &err))	return kXlUtils::setError(err);

	//	run
	double res0;
	kVector<double> s, res;
	if (!kBachelier::fdRunner(s0, r, mu, sigma, expiry, strike, dig > 0, pc, ea, smooth, theta, wind, numStd, numT, numX, update > 0, numPr, res0, s, res, err)) return kXlUtils::setError(err);

	//	size output
	numRows = 3 + s.size();
	numCols = 2;
	LPXLOPER12 out = kXlUtils::getOper(numRows, numCols);

	//	fill output
	kXlUtils::setStr(0, 0, "res 0", out);
	kXlUtils::setDbl(0, 1, res0, out);
	kXlUtils::setStr(2, 0, "s", out);
	kXlUtils::setStr(2, 1, "res", out);
	for (k = 3, i = 0; i < s.size(); ++i, ++k)
	{
		kXlUtils::setDbl(k, 0, s(i), out);
		kXlUtils::setDbl(k, 1, res(i), out);
	}

	//	done
	return out;
}

extern "C" __declspec(dllexport)
LPXLOPER12
xBlackFd(
	LPXLOPER12	params,
	LPXLOPER12	contract,
	LPXLOPER12	gridTech)
{
	FreeAllTempMemory();

	//	help
	string err;
	int numRows, numCols;
	int i, k;

	//	get params
	double s0 = 0.0;
	double r = 0.0;
	double mu = 0.0;
	double sigma = 0.1;
	numRows = (int)getRows(params);
	if (numRows > 0 && !kXlUtils::getDbl(params, 0, 0, s0, &err))		return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getDbl(params, 1, 0, r, &err))		return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(params, 2, 0, mu, &err))		return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getDbl(params, 3, 0, sigma, &err))	return kXlUtils::setError(err);

	//	get contract
	double expiry = 0.0;
	double strike = 0.0;
	double barrier = 0.0;
	int    dig = 0;
	int    pc = 1;
	int	   ea = 0;
	int	   bo = 0;
	int	   smooth = 0;
	numRows = (int)getRows(contract);
	if (numRows > 0 && !kXlUtils::getDbl(contract, 0, 0, expiry, &err))		return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getDbl(contract, 1, 0, strike, &err))		return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(contract, 2, 0, barrier, &err))	return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(contract, 3, 0, dig, &err))		return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getInt(contract, 4, 0, pc, &err))			return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getInt(contract, 5, 0, ea, &err))			return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getInt(contract, 6, 0, bo, &err))			return kXlUtils::setError(err);
	if (numRows > 7 && !kXlUtils::getInt(contract, 7, 0, smooth, &err))		return kXlUtils::setError(err);

	//	get grid tech
	double theta = 0.5;
	int	   wind = 0;
	double numStd = 5.0;
	int    numT = 25;
	int    numX = 50;
	int    update = 1;
	int    numPr = 1;
	int	   timeFactor = 1;
	numRows = (int)getRows(gridTech);
	if (numRows > 0 && !kXlUtils::getDbl(gridTech, 0, 0, theta, &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getInt(gridTech, 1, 0, wind, &err))	return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(gridTech, 2, 0, numStd, &err))return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(gridTech, 3, 0, numT, &err))	return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getInt(gridTech, 4, 0, numX, &err))	return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getInt(gridTech, 5, 0, update, &err))return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getInt(gridTech, 6, 0, numPr, &err))	return kXlUtils::setError(err);
	if (numRows > 7 && !kXlUtils::getInt(gridTech, 7, 0, timeFactor, &err))	return kXlUtils::setError(err);

	// excersice grid
	kVector<kVector<double>> grid;

	//	run
	double res0;
	kVector<double> s, res;
	if (!kBlack::fdRunner(s0, r, mu, sigma, expiry, strike, barrier, dig > 0, pc, ea, bo, smooth, theta, wind, numStd, numT, numX, update > 0, numPr, grid, res0, s, res, err)) return kXlUtils::setError(err);

	if (ea != 3)
	{
		//	size output
		numRows = 3 + s.size();
		numCols = 2;
		LPXLOPER12 out = kXlUtils::getOper(numRows, numCols);

		//	fill output
		kXlUtils::setStr(0, 0, "res 0", out);
		kXlUtils::setDbl(0, 1, res0, out);
		kXlUtils::setStr(2, 0, "s", out);
		kXlUtils::setStr(2, 1, "res", out);
		for (k = 3, i = 0; i < s.size(); ++i, ++k)
		{
			kXlUtils::setDbl(k, 0, s(i), out);
			kXlUtils::setDbl(k, 1, res(i), out);
		}

		//	done
		return out;
	}
	else
	{
		// time grid for pres
		double dt = expiry / numT;
		int outCols = numT / timeFactor + 1;


		//	size output
		numCols = outCols + 1;
		numRows = s.size() + 1;

		LPXLOPER12 out = kXlUtils::getOper(numRows, numCols);

		//	fill output
		for (i = 1; i < numCols; i++) {
			kXlUtils::setDbl(0, i, ((i - 1) * timeFactor) * dt, out);
		}

		for (k = 1; k < numRows; k++) {
			kXlUtils::setDbl(k, 0, s(numRows - k - 1), out);
			for (i = 1; i < numCols - 1; i++) {
				kXlUtils::setDbl(k, i, grid(int((i - 1) * timeFactor))(numRows - k - 1), out);
			}
			if (s(numRows - k - 1) < strike)
			{
				kXlUtils::setDbl(k, numCols - 1, 0, out);
			}
			else
			{
				kXlUtils::setDbl(k, numCols - 1, 1, out);
			}

		}


		//	done
		return out;
	}

}


extern "C" __declspec(dllexport)
LPXLOPER12
xBlackFwdFd(
	LPXLOPER12	params,
	LPXLOPER12	contract,
	LPXLOPER12	gridTech,
	LPXLOPER12  strikes_in)
{
	FreeAllTempMemory();

	//	help
	string err;
	int numRows, numCols;
	int i, k;

	//	get params
	double s0 = 0.0;
	double r = 0.0;
	double mu = 0.0;
	double sigma = 0.1;
	numRows = (int)getRows(params);
	if (numRows > 0 && !kXlUtils::getDbl(params, 0, 0, s0, &err))		return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getDbl(params, 1, 0, r, &err))		return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(params, 2, 0, mu, &err))		return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getDbl(params, 3, 0, sigma, &err))	return kXlUtils::setError(err);

	//	get contract
	double expiry = 0.0;
	double strike = 0.0;
	double barrier = 0.0;
	int    dig = 0;
	int    pc = 1;
	int	   ea = 0;
	int	   bo = 0;
	int	   smooth = 0;
	numRows = (int)getRows(contract);
	if (numRows > 0 && !kXlUtils::getDbl(contract, 0, 0, expiry, &err))		return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getDbl(contract, 1, 0, strike, &err))		return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(contract, 2, 0, barrier, &err))	return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(contract, 3, 0, dig, &err))		return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getInt(contract, 4, 0, pc, &err))			return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getInt(contract, 5, 0, ea, &err))			return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getInt(contract, 6, 0, bo, &err))			return kXlUtils::setError(err);
	if (numRows > 7 && !kXlUtils::getInt(contract, 7, 0, smooth, &err))		return kXlUtils::setError(err);

	//	get grid tech
	double theta = 0.5;
	int	   wind = 0;
	double numStd = 5.0;
	int    numT = 25;
	int    numX = 50;
	int    update = 1;
	int    numPr = 1;
	int    timeFactor = 1;
	numRows = (int)getRows(gridTech);
	if (numRows > 0 && !kXlUtils::getDbl(gridTech, 0, 0, theta, &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getInt(gridTech, 1, 0, wind, &err))	return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(gridTech, 2, 0, numStd, &err))return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(gridTech, 3, 0, numT, &err))	return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getInt(gridTech, 4, 0, numX, &err))	return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getInt(gridTech, 5, 0, update, &err))return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getInt(gridTech, 6, 0, numPr, &err))	return kXlUtils::setError(err);
	if (numRows > 7 && !kXlUtils::getInt(gridTech, 7, 0, timeFactor, &err))	return kXlUtils::setError(err);


	//	get tech
	kVector<double> strikes;
	if (!kXlUtils::getVector(strikes_in, strikes))
		return kXlUtils::setError("strikes is not a vector");

	//	run
	double res0;
	kVector<double> s;
	kVector<kVector<double>> res;

	if (!kBlack::fdFwdRunner(s0, r, mu, sigma, expiry, strike, barrier, dig > 0, pc, ea, bo, smooth, theta, wind, numStd, numT, numX, update > 0, numPr, res0, s, strikes, timeFactor, res, err)) return kXlUtils::setError(err);

	// time grid for pres
	double dt = expiry / numT;
	int outRows = numT / timeFactor + 1;


	//	size output
	numRows = outRows + 1;
	numCols = res(0).size() + 1;
	LPXLOPER12 out = kXlUtils::getOper(numRows, numCols);

	//	fill output
	for (i = 1; i < numRows; i++) {
		kXlUtils::setDbl(i, 0, ((i - 1) * timeFactor) * dt, out);
	}
	for (k = 1; k < numCols; k++)  kXlUtils::setDbl(0, k, strikes(k - 1), out);

	for (i = 1; i < numRows; i++) {
		for (k = 1; k < numCols; k++) {
			kXlUtils::setDbl(i, k, res(int((i - 1) * timeFactor))(k - 1), out);
		}
	}

	//	done
	return out;
}

extern "C" __declspec(dllexport)
LPXLOPER12
xNumericTest(
	LPXLOPER12	params,
	LPXLOPER12  kIn,
	LPXLOPER12  prices)
{
	FreeAllTempMemory();

	//	help
	string err;
	int numRows;
	int i;

	//	get params
	double s0 = 0.0;
	double r = 0.0;
	double mu = 0.0;
	double sigma = 0.1;

	double dt = 0.0;
	double theta = 0.5;
	int	   wind = 0;
	int equation = 1;

	numRows = (int)getRows(params);
	if (numRows > 0 && !kXlUtils::getDbl(params, 0, 0, s0, &err))		return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getDbl(params, 1, 0, r, &err))		return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(params, 2, 0, mu, &err))		return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getDbl(params, 3, 0, sigma, &err))	return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getDbl(params, 4, 0, dt, &err))	return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getDbl(params, 5, 0, theta, &err))	return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getInt(params, 6, 0, wind, &err))	return kXlUtils::setError(err);
	if (numRows > 7 && !kXlUtils::getInt(params, 7, 0, equation, &err))	return kXlUtils::setError(err);
	//	get contract


	kMatrix<double> strikesM;
	kVector<double> strikes;

	kXlUtils::getMatrix(kIn, strikesM);
	strikes.resize(strikesM.size());

	for (int i = 0; i < strikesM.size(); ++i) strikes(i) = strikesM(0, i);

	kMatrix<double> pricesM;
	kVector<double> outputPrices;
	kVector<double> inputPrices;



	kXlUtils::getMatrix(prices, pricesM);
	int nump = pricesM.cols();

	outputPrices.resize(nump);
	inputPrices.resize(nump);
	for (int i = 0; i < nump; ++i) {
		outputPrices(i) = pricesM(0, i);
		inputPrices(i) = pricesM(1, i);
	}



	kFd1d<double> fd;
	fd.init(1, strikes, false);
	int nums = strikes.size();

	if (nums != outputPrices.size())
	{
		return TempXLOPER12();
	}


	//	set parameters
	for (i = 0; i < nums; ++i)
	{
		fd.r()(i) = r;
		fd.mu()(i) = mu * strikes(i);
		fd.var()(i) = sigma * sigma * strikes(i) * strikes(i);
	}

	kMatrix<double> A1, A2;
	kVector<double>	myVs, myWs;

	fd.calcAx(1, -theta * dt, wind, false, A1);
	fd.calcAx(1, (1 - theta) * dt, wind, false, A2);

	if (equation == 1) {

		if (theta != 0.0)
		{
			myVs = outputPrices;
			kMatrixAlgebra::tridag(A1, myVs, outputPrices, myWs);
		}
		if (theta != 1.0)
		{
			myVs = outputPrices;
			kMatrixAlgebra::banmul(A2, 1, 1, myVs, outputPrices);
		}
	}

	if (equation == -1) {

		if (theta != 1.0)
		{
			myVs = outputPrices;
			kMatrixAlgebra::banmul(A2, 1, 1, myVs, outputPrices);
		}

		if (theta != 0.0)
		{
			myVs = outputPrices;
			kMatrixAlgebra::tridag(A1, myVs, outputPrices, myWs);
		}

	}

	LPXLOPER12 out = TempXLOPER12();

	size_t n = outputPrices.size();
	resize(out, 1, n);
	for (size_t j = 0; j < n; ++j)
		setNum(out, inputPrices((int)j) - outputPrices((int)j), 0, j);
	//	done
	return out;
}

// Multi dim

extern "C" __declspec(dllexport)
LPXLOPER12
xFd2d(
	LPXLOPER12 t_in,
	LPXLOPER12 x_in,
	LPXLOPER12 y_in,
	LPXLOPER12 r_in,
	LPXLOPER12 mux_in,
	LPXLOPER12 muy_in,
	LPXLOPER12 varx_in,
	LPXLOPER12 v0_in,
	LPXLOPER12 gridTech)
{
	FreeAllTempMemory();

	//	help
	int numRows, numCols;
	string err;

	//	get t 
	double t;
	if (!kXlUtils::getDbl(t_in, 0, 0, t, &err)) return kXlUtils::setError(err);

	//	get grid tech
	kVector<double> theta;
	theta.resize(2, 0.5);
	kVector<int>	wind;
	wind.resize(2, 0);
	int	   numt = 100;
	int fb = -1;
	kVector<int> pointDisc(2, 3);
	pointDisc(1) = 5;

	numRows = (int)getRows(gridTech);
	numCols = (int)getCols(gridTech);
	if (numRows > 0 && !kXlUtils::getDbl(gridTech, 0, 0, theta(0), &err))	return kXlUtils::setError(err);
	if (numRows > 0 && !kXlUtils::getDbl(gridTech, 0, 1, theta(1), &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getInt(gridTech, 1, 0, wind(0), &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getInt(gridTech, 1, 1, wind(1), &err))	return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getInt(gridTech, 2, 0, numt, &err))	return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(gridTech, 3, 0, fb, &err))	return kXlUtils::setError(err);



	//	fd grid
	kFd2d<double> fd;
	if (!kXlUtils::getVector(x_in, fd.x()))
		return kXlUtils::setError("x is not a vector");

	int n = fd.x().size();

	if (!kXlUtils::getRowVector(y_in, fd.y()))
		return kXlUtils::setError("y is not a vector");

	int m = fd.y().size();

	//	init fd
	fd.init(1, fd.x(), fd.y(), pointDisc);

	if (!kXlUtils::getMatrix(r_in, fd.r()))
		return kXlUtils::setError("r is not a matrix");

	if (n * m != fd.r().size())
		return kXlUtils::setError("r must have same size as x*y");

	if (!kXlUtils::getMatrix(mux_in, fd.mu()))
		return kXlUtils::setError("mux is not a matrix");

	if (n * m != fd.mu().size())
		return kXlUtils::setError("mux must have same size as x*y");

	if (!kXlUtils::getMatrix(muy_in, fd.nu()))
		return kXlUtils::setError("muy is not a matrix");

	if (n * m != fd.nu().size())
		return kXlUtils::setError("muy must have same size as x*y");

	if (!kXlUtils::getMatrix(varx_in, fd.var()))
		return kXlUtils::setError("var is not a matrix");

	if (n * m != fd.var().size())
		return kXlUtils::setError("var must have same size as x*y");

	fd.res();
	if (!kXlUtils::getMatrix(v0_in, fd.res()[0]))
		return kXlUtils::setError("v0 is not a matrix");

	if (n * m != fd.res()[0].size())
		return kXlUtils::setError("v0 must have same size as x*y");

	for (int i = 0; i < 2; i++)
	{
		if (theta(i) < 0.0) theta(i) = 0.0;
		if (theta(i) > 1.0) theta(i) = 1.0;
	}


	double dt = t / numt;
	for (int n = 0; n < numt; ++n)
	{
		if (fb <= 0)
		{
			fd.rollBwd(dt, theta, wind, fd.res());
		}
		else
		{
			fd.rollFwd(dt, theta, wind, fd.res());
		}
	}

	LPXLOPER12 out = TempXLOPER12();
	kXlUtils::setMatrix(fd.res()[0], out);
	return out;
}


extern "C" __declspec(dllexport)
LPXLOPER12
xCheyetteBwd(
	LPXLOPER12	params,
	LPXLOPER12	contract,
	LPXLOPER12	gridTech,
	LPXLOPER12	timelineIn,
	LPXLOPER12	yieldsIn,
	LPXLOPER12	eta2In,
	LPXLOPER12	eta02In)
{
	FreeAllTempMemory();

	//	help
	string err;
	int numRows, numCols;

	//	get params
	double kappa = 0.0;
	double beta = 0.0;
	double sigma = 0.006;
	double T = 28;
	numRows = (int)getRows(params);
	if (numRows > 0 && !kXlUtils::getDbl(params, 0, 0, kappa, &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getDbl(params, 1, 0, beta, &err))		return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(params, 2, 0, sigma, &err))		return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getDbl(params, 3, 0, T, &err))		return kXlUtils::setError(err);

	//	get contract
	int type = 0;
	int option = 0;
	double strike = 0;
	double maturity = 10;
	double tenor = 10;
	double frequency = 1.0;
	int bermudan = 0;
	double finalUpdateTime = 0;
	numRows = (int)getRows(contract);
	if (numRows > 0 && !kXlUtils::getInt(contract, 0, 0, type, &err))			return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getInt(contract, 1, 0, option, &err))			return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(contract, 2, 0, strike, &err))			return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getDbl(contract, 3, 0, maturity, &err))		return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getDbl(contract, 4, 0, tenor, &err))			return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getDbl(contract, 5, 0, frequency, &err))		return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getInt(contract, 6, 0, bermudan, &err))		return kXlUtils::setError(err);
	if (numRows > 7 && !kXlUtils::getDbl(contract, 7, 0, finalUpdateTime, &err))		return kXlUtils::setError(err);


	//	get grid tech
	kVector<double> theta(2, 0.5);
	kVector<int>	wind(2, 0);
	kVector<int>	pointDisc(2, 3);
	int    n = 100;
	int    m = 10;
	kVector<double> width(2, 1);
	kVector<double> alpha(2, 1);
	int	   numt = 100;
	double dt = 0.3;

	numRows = (int)getRows(gridTech);
	numCols = (int)getCols(gridTech);
	if (numRows > 0 && !kXlUtils::getDbl(gridTech, 0, 0, theta(0), &err))	return kXlUtils::setError(err);
	if (numRows > 0 && !kXlUtils::getDbl(gridTech, 0, 1, theta(1), &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getInt(gridTech, 1, 0, wind(0), &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getInt(gridTech, 1, 1, wind(1), &err))	return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getInt(gridTech, 2, 0, pointDisc(0), &err))	return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getInt(gridTech, 2, 1, pointDisc(1), &err))	return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(gridTech, 3, 0, n, &err))		return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(gridTech, 3, 1, m, &err))		return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getDbl(gridTech, 4, 0, width(0), &err))	return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getDbl(gridTech, 4, 1, width(1), &err))	return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getDbl(gridTech, 5, 0, alpha(0), &err))	return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getDbl(gridTech, 5, 1, alpha(1), &err))	return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getInt(gridTech, 6, 0, numt, &err))	return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getDbl(gridTech, 6, 1, dt, &err))	return kXlUtils::setError(err);

	kVector<double> timeline;
	kXlUtils::getVector(timelineIn, timeline);

	kVector<double> yields;
	kXlUtils::getVector(yieldsIn, yields);

	kMatrix<double> eta2;
	kXlUtils::getMatrix(eta2In, eta2);
	if (eta2.cols() <= 1 || eta2.rows() <= 1) return kXlUtils::setError("No volatilities inserted");

	kVector<double> eta02;
	kXlUtils::getVector(eta02In, eta02);

	// helper
	kFuncHelp<double> helper(timeline, yields);
	helper.type = type;
	helper.option = option > 0;
	helper.strike = strike;
	helper.maturity = maturity;
	helper.tenor = tenor;
	helper.frequency = frequency;
	helper.finalUpdateTime = finalUpdateTime;


	//	run


	kCheyette<double> model;
	kMatrix<double> res;

	LPXLOPER12 out = TempXLOPER12();
	if (numt != 0 && abs(maturity / numt - dt) > 1E-8) return out;

	if (!model.init(kappa, 0, 0, n, m, width, T, sigma, alpha, pointDisc, helper)) return kXlUtils::setError(err);

	model.helper.initEta(eta2, eta02);

	if (!model.bwdRunner(wind, theta, maturity, numt, dt, bermudan > 0, res)) return kXlUtils::setError(err);

	;
	kXlUtils::setMatrix(res, out);

	return out;
}

extern "C" __declspec(dllexport)
LPXLOPER12
xCheyetteFwd(
	LPXLOPER12	params,
	LPXLOPER12	contract,
	LPXLOPER12	gridTech,
	LPXLOPER12	timelineIn,
	LPXLOPER12	yieldsIn,
	LPXLOPER12	eta2In,
	LPXLOPER12	eta02In)
{
	FreeAllTempMemory();

	//	help
	string err;
	int numRows, numCols;

	//	get params
	double kappa = 0.0;
	double beta = 0.0;
	double sigma = 0.006;
	double T = 28;
	numRows = (int)getRows(params);
	if (numRows > 0 && !kXlUtils::getDbl(params, 0, 0, kappa, &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getDbl(params, 1, 0, beta, &err))		return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(params, 2, 0, sigma, &err))		return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getDbl(params, 3, 0, T, &err))		return kXlUtils::setError(err);

	//	get contract
	int type = 0;
	int option = 0;
	double strike = 0;
	double maturity = 10;
	double tenor = 10;
	double frequency = 1.0;
	int bermudan = 0;
	double finalUpdateTime = 0;
	numRows = (int)getRows(contract);
	if (numRows > 0 && !kXlUtils::getInt(contract, 0, 0, type, &err))			return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getInt(contract, 1, 0, option, &err))			return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(contract, 2, 0, strike, &err))			return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getDbl(contract, 3, 0, maturity, &err))		return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getDbl(contract, 4, 0, tenor, &err))			return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getDbl(contract, 5, 0, frequency, &err))		return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getInt(contract, 6, 0, bermudan, &err))		return kXlUtils::setError(err);
	if (numRows > 7 && !kXlUtils::getDbl(contract, 7, 0, finalUpdateTime, &err))		return kXlUtils::setError(err);


	//	get grid tech
	kVector<double> theta(2, 0.5);
	kVector<int>	wind(2, 0);
	kVector<int>	pointDisc(2, 3);
	int    n = 100;
	int    m = 10;
	kVector<double> width(2, 1);
	kVector<double> alpha(2, 1);
	int	   numt = 100;
	double dt = 0.3;

	numRows = (int)getRows(gridTech);
	numCols = (int)getCols(gridTech);
	if (numRows > 0 && !kXlUtils::getDbl(gridTech, 0, 0, theta(0), &err))	return kXlUtils::setError(err);
	if (numRows > 0 && !kXlUtils::getDbl(gridTech, 0, 1, theta(1), &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getInt(gridTech, 1, 0, wind(0), &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getInt(gridTech, 1, 1, wind(1), &err))	return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getInt(gridTech, 2, 0, pointDisc(0), &err))	return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getInt(gridTech, 2, 1, pointDisc(1), &err))	return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(gridTech, 3, 0, n, &err))		return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(gridTech, 3, 1, m, &err))		return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getDbl(gridTech, 4, 0, width(0), &err))	return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getDbl(gridTech, 4, 1, width(1), &err))	return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getDbl(gridTech, 5, 0, alpha(0), &err))	return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getDbl(gridTech, 5, 1, alpha(1), &err))	return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getInt(gridTech, 6, 0, numt, &err))	return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getDbl(gridTech, 6, 1, dt, &err))	return kXlUtils::setError(err);

	kVector<double> timeline;
	kXlUtils::getVector(timelineIn, timeline);

	kVector<double> yields;
	kXlUtils::getVector(yieldsIn, yields);

	kMatrix<double> eta2;
	kXlUtils::getMatrix(eta2In, eta2);
	if (eta2.cols() <= 1 || eta2.rows() <= 1) return kXlUtils::setError("No volatilities inserted");

	kVector<double> eta02;
	kXlUtils::getVector(eta02In, eta02);

	// helper
	kFuncHelp<double> helper(timeline, yields);
	helper.type = type;
	helper.option = option > 0;
	helper.strike = strike;
	helper.maturity = maturity;
	helper.tenor = tenor;
	helper.frequency = frequency;


	//	run


	kCheyette<double> model;
	kMatrix<double> res;

	LPXLOPER12 out = TempXLOPER12();
	if (numt != 0 && abs(maturity / numt - dt) > 1E-8) return out;

	if (!model.init(kappa, 0, 0, n, m, width, T, sigma, alpha, pointDisc, helper)) return kXlUtils::setError(err);

	model.helper.initEta(eta2, eta02);

	// to do fix fwd runner... should maybe price different products at different times
	if (!model.fwdRunner(wind, theta, maturity, numt, dt, res)) return kXlUtils::setError(err);

	kXlUtils::setMatrix(res, out);

	return out;
}


extern "C" __declspec(dllexport)
LPXLOPER12
xCheyetteVol(
	LPXLOPER12	params,
	LPXLOPER12	gridTech,
	LPXLOPER12	timelineIn,
	LPXLOPER12	yieldsIn,
	LPXLOPER12	volParams,
	LPXLOPER12	swptionContract,
	LPXLOPER12	strikesIn,
	LPXLOPER12	impliedVolsIn)
{
	FreeAllTempMemory();

	//	help
	string err;
	int numRows, numCols;

	//	get params
	double kappa = 0.0;
	double beta = 0.0;
	double sigma = 0.006;
	double T = 28;
	numRows = (int)getRows(params);
	if (numRows > 0 && !kXlUtils::getDbl(params, 0, 0, kappa, &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getDbl(params, 1, 0, beta, &err))		return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(params, 2, 0, sigma, &err))		return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getDbl(params, 3, 0, T, &err))		return kXlUtils::setError(err);

	//	get grid tech
	kVector<double> theta(2, 0.5);
	kVector<int>	wind(2, 0);
	kVector<int>	pointDisc(2, 3);
	int    n = 100;
	int    m = 10;
	kVector<double> width(2, 1);
	kVector<double> alpha(2, 1);
	int	   numt = 100;
	double dt = 0.3;

	numRows = (int)getRows(gridTech);
	numCols = (int)getCols(gridTech);
	if (numRows > 0 && !kXlUtils::getDbl(gridTech, 0, 0, theta(0), &err))	return kXlUtils::setError(err);
	if (numRows > 0 && !kXlUtils::getDbl(gridTech, 0, 1, theta(1), &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getInt(gridTech, 1, 0, wind(0), &err))	return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getInt(gridTech, 1, 1, wind(1), &err))	return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getInt(gridTech, 2, 0, pointDisc(0), &err))	return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getInt(gridTech, 2, 1, pointDisc(1), &err))	return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(gridTech, 3, 0, n, &err))		return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(gridTech, 3, 1, m, &err))		return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getDbl(gridTech, 4, 0, width(0), &err))	return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getDbl(gridTech, 4, 1, width(1), &err))	return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getDbl(gridTech, 5, 0, alpha(0), &err))	return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getDbl(gridTech, 5, 1, alpha(1), &err))	return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getInt(gridTech, 6, 0, numt, &err))	return kXlUtils::setError(err);
	if (numRows > 6 && !kXlUtils::getDbl(gridTech, 6, 1, dt, &err))	return kXlUtils::setError(err);

	kVector<double> timeline;
	kXlUtils::getVector(timelineIn, timeline);

	kVector<double> yields;
	kXlUtils::getVector(yieldsIn, yields);

	//	get vol params
	int				method = 1;
	kVector<double> nuBound(2, 0);
	kVector<double>	strikeBound(2, 0);
	int				interpolation = 0;
	int				extrapolation = 0;

	numRows = (int)getRows(volParams);
	numCols = (int)getCols(volParams);
	if (numRows > 0 && !kXlUtils::getInt(volParams, 0, 0, method, &err))		return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getDbl(volParams, 1, 0, nuBound(0), &err))		return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getDbl(volParams, 1, 1, nuBound(1), &err))		return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(volParams, 2, 0, strikeBound(0), &err))	return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(volParams, 2, 1, strikeBound(1), &err))	return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(volParams, 3, 0, interpolation, &err))	return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getInt(volParams, 3, 1, extrapolation, &err))	return kXlUtils::setError(err);

	kVector<double> strikes;
	kXlUtils::getRowVector(strikesIn, strikes);

	int lengthN = (int)getRows(swptionContract);
	kVector<double> maturity(lengthN);
	kVector<double> tenor(lengthN);
	kVector<double> frequency(lengthN);

	for (int i = 0; i < lengthN; ++i)
	{
		kXlUtils::getDbl(swptionContract, i, 0, maturity(i), &err);
		kXlUtils::getDbl(swptionContract, i, 1, tenor(i), &err);
		kXlUtils::getDbl(swptionContract, i, 2, frequency(i), &err);
	}

	kMatrix<double> impliedVols;
	kXlUtils::getMatrix(impliedVolsIn, impliedVols);

	// helper
	kFuncHelp<double> helper(timeline, yields);
	kSwaption<double> swptions = kSwaption<double>(impliedVols, strikes, tenor, frequency, maturity);

	//	run
	kCheyette<double> model;
	kMatrix<double> res;
	if (!model.init(kappa, 0, 0, n, m, width, T, sigma, alpha, pointDisc, helper)) return kXlUtils::setError(err);

	if (method == 1)
	{
		if (!model.initVol(T, numt, swptions, wind, theta, nuBound, strikeBound, interpolation, extrapolation)) return kXlUtils::setError(err);
	}
	else
	{
		if (!model.initVol(sigma, beta, T, numt)) return kXlUtils::setError(err);
	}

	LPXLOPER12 out = TempXLOPER12();

	size_t N = model.helper.myEta2.rows();
	size_t M = model.helper.myEta2.cols();
	resize(out, N + 2, M + 2);
	for (size_t i = 0; i < N; ++i) {
		setNum(out, model.helper.myEta02((int)i), i + 2, M + 1);
		for (size_t j = 0; j < M; ++j) {
			setNum(out, model.helper.myEta2((int)i, (int)j), i + 2, j);
		}
	}

	for (size_t j = 0; j < M; ++j) {
		setNum(out, model.myX((int)j), 0, j);
	}

	return out;
}

extern "C" __declspec(dllexport)
LPXLOPER12
xSwapBachelier(
	LPXLOPER12	contract,
	LPXLOPER12	_timeline,
	LPXLOPER12	_yields
)
{
	FreeAllTempMemory();

	//	help
	string err;
	int numRows, numCols;
	//	get contract
	int type = 1;
	double strike = 0.0;
	double T = 15;
	double tenor = 1.0;
	double frequency = 1.0;
	double sigma = 0.006;

	numRows = (int)getRows(contract);
	if (numRows > 0 && !kXlUtils::getInt(contract, 0, 0, type, &err))				return kXlUtils::setError(err);
	if (numRows > 1 && !kXlUtils::getDbl(contract, 1, 0, strike, &err))		return kXlUtils::setError(err);
	if (numRows > 2 && !kXlUtils::getDbl(contract, 2, 0, T, &err))					return kXlUtils::setError(err);
	if (numRows > 3 && !kXlUtils::getDbl(contract, 3, 0, tenor, &err))			return kXlUtils::setError(err);
	if (numRows > 4 && !kXlUtils::getDbl(contract, 4, 0, frequency, &err))			return kXlUtils::setError(err);
	if (numRows > 5 && !kXlUtils::getDbl(contract, 5, 0, sigma, &err))				return kXlUtils::setError(err);

	kVector<double> timeline;
	kXlUtils::getVector(_timeline, timeline);
	kVector<double> yields;
	kXlUtils::getVector(_yields, yields);

	kFuncHelp<double> helper(timeline, yields);

	helper.type = type;

	kVector<double> x(1, 0.0);
	kVector<double> y(1, 0.0);

	kMatrix<double> floatingLeg, fixedLeg, swpRate;

	helper.swpCalc(0.0, T, 0.0, 0.0, tenor, frequency, x, y, floatingLeg, fixedLeg, swpRate);

	double price = kBachelier::call(T, strike, swpRate(0, 0), sigma);

	if (type < 0)
	{
		price -= swpRate(0, 0) / fixedLeg(0, 0);
	}

	LPXLOPER12 out = TempXLOPER12();

	resize(out, 3, 1);
	setNum(out, fixedLeg(0, 0), 0, 0);
	setNum(out, swpRate(0, 0), 1, 0);
	setNum(out, price * fixedLeg(0, 0), 2, 0);

	return out;
}

extern "C" __declspec(dllexport)
LPXLOPER12
xInterpolation2d(
	LPXLOPER12	settings,
	LPXLOPER12	_xa,
	LPXLOPER12	_ya,
	LPXLOPER12	_za,
	LPXLOPER12	_x,
	LPXLOPER12	_y
)
{
	FreeAllTempMemory();

	//	help
	string err;
	int numRows, numCols;
	//	get contract
	int interpolationType = 1;

	numRows = (int)getRows(settings);
	if (numRows > 0 && !kXlUtils::getInt(settings, 0, 0, interpolationType, &err))	return kXlUtils::setError(err);

	kVector<double> xa;
	kXlUtils::getVector(_xa, xa);

	kVector<double> ya;
	kXlUtils::getRowVector(_ya, ya);

	kMatrix<double> za;
	kXlUtils::getMatrix(_za, za);

	kVector<double> x;
	kXlUtils::getVector(_x, x);

	kVector<double> y;
	kXlUtils::getVector(_y, y);

	if (x.size() != y.size()) return kXlUtils::setError(err);
	int size = x.size();
	int n = xa.size();
	int m = ya.size();

	kVector<double> temp(n);

	kVector<double> res(size, 0.0);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < n; j++)
		{
			kPolation::interpolate(0, 1, interpolationType, y()(i, 1), ya(), za()(j), temp()(j, 1));
		}
		kPolation::interpolate(0, 1, interpolationType, x()(i, 1), xa(), temp(), res()(i, 1));
	}

	LPXLOPER12 out = TempXLOPER12();
	kXlUtils::setVector(res, out);

	return out;
}

//	Registers

extern "C" __declspec(dllexport) int xlAutoOpen(void)
{
	XLOPER12 xDLL;

	Excel12f(xlGetName, &xDLL, 0);

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xMultiply2Numbers"),
		(LPXLOPER12)TempStr12(L"BBB"),
		(LPXLOPER12)TempStr12(L"xMultiply2Numbers"),
		(LPXLOPER12)TempStr12(L"x, y"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Multiplies 2 numbers"),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xMatrixMul"),
		(LPXLOPER12)TempStr12(L"QQQ"),
		(LPXLOPER12)TempStr12(L"xMatrixMul"),
		(LPXLOPER12)TempStr12(L"A, B"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Multiplying 2 matrices."),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xBachelierCall"),
		(LPXLOPER12)TempStr12(L"QQQQQ"),
		(LPXLOPER12)TempStr12(L"xBachelierCall"),
		(LPXLOPER12)TempStr12(L"expiry, strike, forward, volatility"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Price option in Bachelier model"),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xBachelierImplied"),
		(LPXLOPER12)TempStr12(L"QQQQQ"),
		(LPXLOPER12)TempStr12(L"xBachelierImplied"),
		(LPXLOPER12)TempStr12(L"expiry, strike, price, forward"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Compute implied volatility in Bachelier model"),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xBlackCall"),
		(LPXLOPER12)TempStr12(L"QQQQQ"),
		(LPXLOPER12)TempStr12(L"xBlackCall"),
		(LPXLOPER12)TempStr12(L"expiry, strike, forward, volatility"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Price option in Black model"),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xBlackImplied"),
		(LPXLOPER12)TempStr12(L"QQQQQ"),
		(LPXLOPER12)TempStr12(L"xBlackImplied"),
		(LPXLOPER12)TempStr12(L"expiry, strike, price, forward"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Compute implied volatility in Black model"),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xTridag"),
		(LPXLOPER12)TempStr12(L"QQQ"),
		(LPXLOPER12)TempStr12(L"xTridag"),
		(LPXLOPER12)TempStr12(L"A, b"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Solving tri-diagonal system."),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xBanmul"),
		(LPXLOPER12)TempStr12(L"QQQBB"),
		(LPXLOPER12)TempStr12(L"xBanmul"),
		(LPXLOPER12)TempStr12(L"A, x, m1, m2"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Multiplying band-matrix with vector."),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xFd1d"),
		(LPXLOPER12)TempStr12(L"QQQQQQQQ"),
		(LPXLOPER12)TempStr12(L"xFd1d"),
		(LPXLOPER12)TempStr12(L"t, x, r, mu, sigma, v0, tech"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Solve 1d fd."),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xBachelierFd"),
		(LPXLOPER12)TempStr12(L"QQQQ"),
		(LPXLOPER12)TempStr12(L"xBachelierFd"),
		(LPXLOPER12)TempStr12(L"params, contract, gridTech"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Solve fd for Bachelier model."),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xBlackFd"),
		(LPXLOPER12)TempStr12(L"QQQQ"),
		(LPXLOPER12)TempStr12(L"xBlackFd"),
		(LPXLOPER12)TempStr12(L"params, contract, gridTech"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Solve fd for Black model."),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xBlackFwdFd"),
		(LPXLOPER12)TempStr12(L"QQQQQ"),
		(LPXLOPER12)TempStr12(L"xBlackFwdFd"),
		(LPXLOPER12)TempStr12(L"params, contract, gridTech, strikes"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Solve Fwd fd for Black model."),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xNumericTest"),
		(LPXLOPER12)TempStr12(L"QQQQ"),
		(LPXLOPER12)TempStr12(L"xNumericTest"),
		(LPXLOPER12)TempStr12(L"params, strikes, prices"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Numeric test for Dupire"),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xFd2d"),
		(LPXLOPER12)TempStr12(L"QQQQQQQQQQ"),
		(LPXLOPER12)TempStr12(L"xFd2d"),
		(LPXLOPER12)TempStr12(L"t, x, y, r, mux, muy, sigmax, v0, tech"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Solve 2d fd."),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xCheyetteBwd"),
		(LPXLOPER12)TempStr12(L"QQQQQQQQ"),
		(LPXLOPER12)TempStr12(L"xCheyetteBwd"),
		(LPXLOPER12)TempStr12(L"params, contract, gridTech, timeline, yields, swaptions, strikes"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Cheyette backward runner"),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xCheyetteFwd"),
		(LPXLOPER12)TempStr12(L"QQQQQQQQ"),
		(LPXLOPER12)TempStr12(L"xCheyetteFwd"),
		(LPXLOPER12)TempStr12(L"params, contract, gridTech, timeline, yields, swaptions, strikes"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Cheyette forward runner"),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xCheyetteVol"),
		(LPXLOPER12)TempStr12(L"QQQQQQQQQ"),
		(LPXLOPER12)TempStr12(L"xCheyetteVol"),
		(LPXLOPER12)TempStr12(L"params, gridTech, timeline, yields, volParams, maturities, strikes, impliedVols"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Cheyette Volatility Estimator"),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xSwapBachelier"),
		(LPXLOPER12)TempStr12(L"QQQQ"),
		(LPXLOPER12)TempStr12(L"xSwapBachelier"),
		(LPXLOPER12)TempStr12(L"contract, timeline, yields"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Gives the fixed leg, swaprate and price of swaption"),
		(LPXLOPER12)TempStr12(L""));

	Excel12f(xlfRegister, 0, 11, (LPXLOPER12)&xDLL,
		(LPXLOPER12)TempStr12(L"xInterpolation2d"),
		(LPXLOPER12)TempStr12(L"QQQQQQQ"),
		(LPXLOPER12)TempStr12(L"xInterpolation2d"),
		(LPXLOPER12)TempStr12(L"settings, xa, ya, za, x, y"),
		(LPXLOPER12)TempStr12(L"1"),
		(LPXLOPER12)TempStr12(L"myOwnCppFunctions"),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L""),
		(LPXLOPER12)TempStr12(L"Does a two dimensional interpolation, with the given interpolation technique"),
		(LPXLOPER12)TempStr12(L""));

	/* Free the XLL filename */
	Excel12f(xlFree, 0, 1, (LPXLOPER12)&xDLL);

	return 1;
}