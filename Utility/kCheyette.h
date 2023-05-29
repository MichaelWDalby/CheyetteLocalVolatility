#pragma once

//	desc:	cheyette model fd grid
//	auth:	michael worsaae
//	date:	feb 2023	

//	includes
#include "kSpecialFunction.h"
#include "kInlines.h"
#include "kVector.h"
#include "kMatrix.h"
#include "kFd2d.h"
#include "kFuncHelp.h"
#include "kSwaption.h"
#include "kPolation.h"
#include "kColVectorView.h"
#include <cmath>

using std::max;


template <class V>
class kCheyette {
public:
	bool init(
		V				kappa,
		V				x0,
		V				y0,
		int				n,
		int				m,
		kVector<V>&		widthscale,
		V				T,
		V				sigma,
		kVector<V>&		alpha,
		kVector<int>&	pointDisc,
		kFuncHelp<V>&	help
	);

	bool initVol(
		V	etaC,
		V	beta,
		V	T,
		int numt
	);

	bool initVol(
		V				T,
		int				numt,
		kSwaption<V>& swptions,
		kVector<int>&	wind,
		kVector<V>&		theta,
		kVector<V>&		nuBound,
		kVector<V>&		strikeBound,
		int				interpolation,
		int				extrapolation
	);

	bool bwdRunner(
		kVector<int>	wind,
		kVector<V>		theta,
		V				T,
		int				numt,
		V				dt,
		bool			bermudan,
		kMatrix<V>&		res
	);

	bool fwdRunner(
		kVector<int> wind,
		kVector<V> theta,
		V	T,
		int numt,
		V	dt,
		kMatrix<V>& res
	);

	kVector<V>	 myX;
	kFuncHelp<V> helper;
private:
	void setParams(
		int	indext,
		V	t,
		V	dt,
		V	y0t
	);

	
	V					myKappa;
	kFd2d<V>			fd;
	int					myN, myM;
	kVector<V>			myYtil;
	kVector<int>		myXYindex = kVector<int>(2);
	

	
};

template <class V>
bool
kCheyette<V>::init(
	V				kappa,
	V				x0,
	V				y0,
	int				n,
	int				m,
	kVector<V>&		widthscale,
	V				T,
	V				sigma,
	kVector<V>&		alpha,
	kVector<int>&	pointDisc,
	kFuncHelp<V>&	help
)
{
	int i, j;

	myN = int(n / 2) * 2 + 1;
	myM = int(m / 2) * 2 + 1;

	// params
	myKappa = kappa;
	

	// x and y
	double c1, c2, ksi, dksi;

	myX.resize(myN, V(0.0));
	myYtil.resize(myM, V(0.0));

	V qVar = 0.0;
	helper.y0Step(qVar, kappa, 1, T);
	qVar = sigma * sigma * T * T / qVar;
	V gridU = widthscale(0) * sqrt(qVar);

	if (alpha(0) < 1E6 || alpha(1) < 1E6)
	{
		V alphaScaled = alpha(0) * sqrt(qVar);

		c1 = kSpecialFunction::sinhInv(-gridU / alphaScaled);
		c2 = kSpecialFunction::sinhInv(gridU / alphaScaled);

		ksi = 0;
		dksi = 1.0 / (myN - 1.0);

		for (i = 0; i < myN; ++i) {
			myX(i) = alphaScaled * sinh(c2 * ksi + c1 * (1 - ksi));
			ksi += dksi;
		}

		gridU = widthscale(1) * qVar;
		alphaScaled = alpha(1) * qVar;

		c1 = kSpecialFunction::sinhInv(-gridU / alphaScaled);
		c2 = kSpecialFunction::sinhInv(gridU / alphaScaled);

		ksi = 0;
		dksi = 1.0 / (myM - 1.0);

		for (j = 0; j < myM; ++j) {
			myYtil(j) = alphaScaled * sinh(c2 * ksi + c1 * (1 - ksi));
			ksi += dksi;
		}
	}
	else
	{
		ksi = -gridU;
		dksi = 2 * gridU / (myN - 1.0);

		for (i = 0; i < myN; ++i) {
			myX(i) = ksi;
			ksi += dksi;
		}

		gridU = widthscale(1) * qVar;

		ksi = -gridU;
		dksi = 2 * gridU / (myM - 1.0);

		for (i = 0; i < myM; ++i) {
			myYtil(i) = ksi;
			ksi += dksi;
		}
	}

	kPolation::locate(myX(), myN, 0.0, j);
	if (-myX(j) < myX(j + 1))
	{
		myXYindex(0) = j;
		myX(j) = 0.0;
	}
	else
	{
		myXYindex(0) = j + 1;
		myX(j + 1) = 0.0;
	}

	kPolation::locate(myYtil(), myM, 0.0, j);
	if (-myYtil(j) < myYtil(j + 1))
	{
		myXYindex(1) = j;
		myYtil(j) = 0.0;
	}
	else
	{
		myXYindex(1) = j + 1;
		myYtil(j + 1) = 0.0;
	}

	fd.init(1, myX, myYtil, pointDisc);

	helper = help;

	return true;
}


template <class V>
bool
kCheyette<V>::initVol(
	V	etaC,
	V	beta,
	V	T,
	int numt
)
{
	V temp = 0.0;
	helper.y0Step(temp, myKappa, 1, T);
	etaC = sqrt(etaC * etaC * T / temp);
	
	
	int k, i;
	V eta2;
	int size = numt;

	// sizing
	helper.myEta2.resize(size, myN);
	helper.myEta02.resize(size);
	for (k = 0; k < size; ++k)

	// set value
	for (i = 0; i < myN; ++i)
	{
		eta2 = etaC * max(1 + beta * myX(i), 0.0);
		eta2 *= eta2;
		for (k = 0; k < size; ++k)
		{
			helper.myEta2(k, i) = eta2;
		}
	}
	for (k = 0; k < size; ++k)
	{
		helper.myEta02(k) = etaC * etaC;
	}
	
	return true;
}

template <class V>
bool 
kCheyette<V>::initVol(
	V				T,
	int				numt,
	kSwaption<V>&	swptions,
	kVector<int>&	wind,
	kVector<V>&		theta,
	kVector<V>&		nuBound,
	kVector<V>&		strikeBound,
	int				interpolation,
	int				extrapolation
)
{
	// initiation
	int i, j, k, h, l;
	kVector<int> indexStrikeBound(2, 0);
	nuBound(0) *= nuBound(0);
	nuBound(1) *= nuBound(1);
	l = 0;

	V tenor, realisedTenor;
	V maturity;
	V freqInv;
	V missingTenor;

	kMatrix<V> fixedLeg, fixedLegTemp;
	kMatrix<V> floatingLeg, floatingLegTemp;
	kMatrix<V> swpRate(myN, myM);
	V y0t = V(0.0), y0help = V(0.0);
	V discount;
	kVector<V> strikes;
	strikes.resize(myN);


	V dt = min(T, swptions.maturity(swptions.maturity.size() - 1)) / numt;
	V t = 0.0;

	kVector<V> prices;
	kVector<V> inputPrices;
	prices.resize(myN);
	inputPrices.resize(myN);

	int size;

	int dim = 3;
	kVectorView<V>	resVec	= kVectorView<V>();
	kVectorView<V>	eta2 = kVectorView<V>();
	kVectorView<V>	probabilityView = kVectorView<V>();
	kVector<V>		myVs;
	int mmx = (dim - 1) / 2;

	// sizing
	helper.myEta2.resize(numt, myN);
	helper.myEta02.resize(numt);

	// Setting initial boundary value
	for (i = 0; i < myN; ++i)
	{
		for (j = 0; j < myM; ++j)
		{
			fd.res()(0)(i, j) = 0.0;
		}
	}
	fd.res()(0)(myXYindex(0), myXYindex(1)) = 1.0;

	// delta operators
	kMatrix<V> myDx;
	kMatrix<V> myDkk;
	kFiniteDifference::dx(0, myX, myDx);
	

	// A0
	kVector<V> x00 = kVector<V>(1);
	kVector<V> y00 = kVector<V>(1);
	x00(0) = y00(0) = 0.0;

	kMatrix<V> F0, A0, F0Temp, A0Temp;
	kMatrix<V> S0(1, 1);

	// roll
	for (h = 0; h < numt; ++h)
	{		
		// Finding relevant input swaption
		while (!(swptions.maturity(l) + 1E-10 >= t + dt))
		{
			++l;
		}

		maturity = t + dt;

		// Tenor interpolation
		tenor = swptions.tenor(l);
		if (l > 0)
		{
			tenor = ((swptions.maturity(l) - maturity) * swptions.tenor(l - 1)
				+ (maturity - swptions.maturity(l - 1)) * swptions.tenor(l))
				/ (swptions.maturity(l) - swptions.maturity(l - 1));
		}
		
		// No tenor interpolation
		//tenor = swptions.tenor(l);


		freqInv = 1 / swptions.frequency(l);
		realisedTenor = int(tenor * freqInv + 1E-10) * freqInv;
		missingTenor = tenor - realisedTenor;

		// Finding fixed leg, floating leg and swap-rate
		helper.floatingCalc(t, maturity, myKappa, y0t,
			realisedTenor, myX, myYtil, floatingLeg);
		helper.fixedCalc(t, maturity, myKappa, y0t,
			realisedTenor, swptions.frequency(l), myX, myYtil, fixedLeg);

		helper.floatingCalc(0.0, maturity, myKappa, 0.0,
			realisedTenor, x00, y00, F0);
		helper.fixedCalc(0.0, maturity, myKappa, 0.0,
			realisedTenor, swptions.frequency(l), x00, y00, A0);
		
		if (missingTenor > 1E-10)
		{
			helper.floatingCalc(t, maturity + realisedTenor, myKappa, y0t,
				missingTenor, myX, myYtil, floatingLegTemp);
			helper.fixedCalc(t, maturity + realisedTenor, myKappa, y0t,
				missingTenor, 1 / missingTenor, myX, myYtil, fixedLegTemp);
			floatingLeg += floatingLegTemp;
			fixedLeg += fixedLegTemp;
			
			helper.floatingCalc(0.0, maturity + realisedTenor, myKappa, 0.0,
				missingTenor, x00, y00, F0Temp);
			helper.fixedCalc(0.0, maturity + realisedTenor, myKappa, 0.0,
				missingTenor, 1 / missingTenor, x00, y00, A0Temp);
			F0 += F0Temp;
			A0 += A0Temp;
		}

		for (i = 0; i < myN; ++i) {
			for (j = 0; j < myM; ++j)
			{
				swpRate(i, j) = floatingLeg(i, j) / fixedLeg(i, j);
			}
		}

		S0(0, 0) = F0(0, 0) / A0(0, 0);
		
		// Matching strikes to swap-rates
		for (i = 0; i < myN; ++i)
		{
			strikes(i) = swpRate(i, myXYindex(1));
		}

		// Finding observed prices
		swptions.bachSetPrice(l, interpolation, extrapolation, S0(0, 0), maturity, strikes, inputPrices);

		// Finding grid prices
		for (k = 0; k < myN; ++k)
		{
			prices(k) = V(0.0);
			for (i = 0; i < myN; ++i)
			{
				for (j = 0; j < myM; ++j)
				{
					prices(k) += fd.res()(0)(i, j) * fixedLeg(i, j) * max(swpRate(i, j) - strikes(k), 0.0);
				}
			}
			prices(k) /= A0(0, 0);
		}

		// Using Dupire to find Sigma2
		eta2 = helper.myEta2()(h);

		for (i = 0; i < myN; ++i)
		{
			eta2(i) = (2 * (inputPrices(i) - prices(i))) / dt;
			prices(i) = (1 - theta(0)) * prices(i) + theta(0) * inputPrices(i);
		}
		
		kPolation::locate(strikes(), myN, strikeBound(0) + S0(0, 0), indexStrikeBound(0));
		kPolation::locate(strikes(), myN, strikeBound(1) + S0(0, 0), indexStrikeBound(1));

		i = indexStrikeBound(0) + 1 - 1;
		j = indexStrikeBound(1) + 1;

		size = 1 + j - i;

		resVec = strikes()(i, size);
		kFiniteDifference::dxx(resVec, myDkk);


		resVec = prices()(i, size);
		myVs = resVec;
		kMatrixAlgebra::banmul(myDkk, mmx, mmx, myVs, resVec);

		resVec = helper.myEta2()(h);

		for (k = i + 1; k < j; ++k)
		{
			eta2(k) /= prices(k);
		}

		eta2(0) = 0.0;
		for (k = 1; k <= i; ++k)
		{
			eta2(k) = eta2(i + 1);
		}

		eta2(myN - 1) = 0.0;
		for (k = j; k < (myN - 1); ++k)
		{
			eta2(k) = eta2(j - 1);
		}

		

		// finding Sx
		kColVectorView<V> resColVec(swpRate, myXYindex(1));
		myVs = resColVec;
		kMatrixAlgebra::banmul(myDx, mmx, mmx, myVs, resColVec);

		// Computing etas with sigma2/sx2
		for (i = 1; i < myN - 1; ++i)
		{
			eta2(i) /= swpRate(i, myXYindex(1)) * swpRate(i, myXYindex(1));
			eta2(i) = kInlines::bound(nuBound(0), eta2(i), nuBound(1));	
		}


		helper.myEta02(h) = V(0.0);
		discount = V(0.0);
		// computing eta02 under forward measure for our y transformation
		for (i = 0; i < myN; ++i)
		{
			probabilityView = fd.res()(0)()(i);
			for (j = 0; j < myM; ++j)
			{
				helper.myEta02(h) += eta2(i) * probabilityView(j);
				discount += probabilityView(j);
			}
		}
		helper.myEta02(h) /= discount;

		// Setting parameters
		helper.y0Step(y0help, myKappa, helper.myEta02(h), dt);
		y0t = 0.5 * (y0t + y0help);
		setParams(h, t, dt, y0t);
		
		y0t = y0help;
		t += dt;

		// Rolling
		
		fd.rollFwd(dt, theta, wind, fd.res());	
	}

	return true;
}

template <class V>
void
kCheyette<V>::setParams(
	int	indext,
	V	t,
	V	dt,
	V	y0t
)
{
	int i, j;
	V y;
	V T = t + dt;

	V fwdRate = (helper.yield(T) * T - helper.yield(t) * t) / dt;
	V rt, muTemp, eta2diff;
	V eta02t = helper.myEta02(indext);

	// set params
	for (i = 0; i < myN; ++i)
	{
		rt = fwdRate + myX(i);
		muTemp = -myKappa * myX(i);
		eta2diff = helper.myEta2(indext, i) - eta02t;
		for (j = 0; j < myM; ++j) {
			y = myYtil(j) + y0t;
			fd.r()(i, j) = rt;
			fd.mu()(i, j) = muTemp + y;
			fd.var()(i, j) = helper.myEta2(indext, i);
			fd.nu()(i, j) = eta2diff - 2 * myKappa * myYtil(j);
		}
	}
}






template <class V>
bool
kCheyette<V>::bwdRunner(
	kVector<int>	wind,
	kVector<V>		theta,
	V				T,
	int				numt,
	V				dt,
	bool			bermudan,
	kMatrix<V>&		res
)
{
	int h;
	

	V bermudanFinalUpdateTime = helper.finalUpdateTime;
	if (helper.type == 2 || helper.type == -2) {
		helper.finalUpdateTime = min(helper.finalUpdateTime, T);
		
		T += helper.tenor - 1 / helper.frequency;
		numt = T / dt + 0.5;		
	}
	V t = 0.0;
	kMatrix<V> exVal;
	exVal.resize(myN, myM);
	V y0t = V(0.0);
	V y0help;

	for (h = 0; h < numt; ++h)
	{
		helper.y0Step(y0t, myKappa, helper.myEta02(h), dt);
		t += dt;
	}
	y0help = y0t;

	// Setting initial boundary value
	fd.res()(0) = 0.0;
	helper.updateValue(helper.type, T, myKappa, y0t, myX, myYtil, fd.res()(0));
	
	if (helper.option)
	{
		helper.Max(fd.res()(0), 0.0);
	}

	res.resize(1, 1);

	bool willUpdate = false;
	if (helper.type == 2 || helper.type == -2) willUpdate = true;

	V updateTime = helper.nextUpdateTime(T);
	bool updateTimeReached = false;

	// roll
	for (h = numt - 1; h >= 0; --h)
	{
		t -= dt;
		
		helper.y0Step(y0help, myKappa, helper.myEta02(h), -dt);
		y0t = 0.5 * (y0t + y0help);


		setParams(h, t, dt, y0t);
		fd.rollBwd(dt, theta, wind, fd.res());

		y0t = y0help;
		
		updateTimeReached = t - updateTime < 1E-8;

		if (updateTimeReached) {
			updateTime = helper.nextUpdateTime(updateTime);

			if (willUpdate)
			{
				helper.updateValue(helper.type, t, myKappa, y0t, myX, myYtil, fd.res()(0));
			}
			if (bermudan && t + 1E-8 > bermudanFinalUpdateTime) {
				helper.exerciseValue(helper.type, t, myKappa, y0t, myX, myYtil, exVal);
				helper.Max(fd.res()(0), exVal);
			}
		}		
	}
	res(0, 0) = fd.res()(0)(myXYindex(0), myXYindex(1));

	return true;
}

template <class V>
bool
kCheyette<V>::fwdRunner(
	kVector<int>	wind,
	kVector<V>		theta,
	V				T,
	int				numt,
	V				dt,
	kMatrix<V>&		res
)
{
	int i, j, h;
	
	double t = 0.0;
	kMatrix<V> prices;
	prices.resize(myN, myM, 0.0);
	
	V y0t = V(0.0);
	V y0help = V(0.0);

	// Setting initial boundary value
	for (i = 0; i < myN; ++i)
	{
		for (j = 0; j < myM; ++j)
		{
			fd.res()(0)(i, j) = 0.0;
		}
	}
	fd.res()(0)(myXYindex(0), myXYindex(1)) = 1.0;

	
	// roll
	for (h = 0; h < numt; ++h)
	{		
		helper.y0Step(y0help, myKappa, helper.myEta02(h), dt);
		y0t = 0.5 * (y0t + y0help);
		setParams(h, t, dt, y0t);
		fd.rollFwd(dt, theta, wind, fd.res());

		y0t = y0help;
		t += dt;
	}

	res.resize(1, 2, 0);
	
	helper.updateValue(helper.type, T, myKappa, y0t, myX, myYtil, prices);

	if (helper.option)
	{
		helper.Max(prices, 0.0);
	}

	for (i = 0; i < myN; ++i)
	{
		for (j = 0; j < myM; ++j)
		{
			res(0, 0) += fd.res()(0)(i, j) * prices(i, j);
			res(0, 1) += fd.res()(0)(i, j);
		}
	}

	if (helper.type == 42)
	{
		V zcb = res(0, 1);
		res.resize(myN + 1, myM + 1);
		res(0, 0) = y0t;

		for (i = 0; i < myN; ++i)
		{
			res(i + 1, 0) = myX(i);
			for (j = 0; j < myM; ++j) {
				res(i + 1, j + 1) = fd.res()(0)(i, j) / zcb;
			}
		}

		for (j = 0; j < myM; ++j) {
			res(0, j + 1) = myYtil(j);
		}
	}
	

	return true;
}