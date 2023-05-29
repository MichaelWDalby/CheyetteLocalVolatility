#pragma once

//	desc:	helper function for cheyette model
//	auth:	michael worsaae
//	date:	mar 2023	

//	includes
#include "kVector.h"
#include "kMatrix.h"
#include "kPolation.h"
#include <cmath>

using std::max;

template <class V>
class kFuncHelp {
public:
	int type;
	bool option;
	V strike, tenor, frequency, maturity;

	V finalUpdateTime;

	kVector<V> myTimeline, myYields;

	kMatrix<V> myEta2;
	kVector<V> myEta02;

	kFuncHelp() = default;

	kFuncHelp(
		kVector<V>& timeline,
		kVector<V>& yields);

	void initEta(
		kMatrix<V>& _myEta2,
		kVector<V>& _myEta02
	);

	V G(
		V kappa,
		V t,
		V T);

	V bondRecon(
		V bondFactor,
		V valG,
		V x,
		V y
	);


	void fixedCalc(
		V			t,
		V			T,
		V			myKappa,
		V			y0t,
		V			_tenor,
		V			_frequency,
		kVector<V>& myX,
		kVector<V>& myYtil,
		kMatrix<V>& fixedLeg
	);

	void floatingCalc(
		V			t,
		V			T,
		V			myKappa,
		V			y0t,
		V			_tenor,
		kVector<V>& myX,
		kVector<V>& myYtil,
		kMatrix<V>& floatingLeg
	);


	void swpCalc(
		V			t,
		V			T,
		V			myKappa,
		V			y0t,
		V			_tenor,
		V			_frequency,
		kVector<V>& myX,
		kVector<V>& myYtil,
		kMatrix<V>& floatingLeg,
		kMatrix<V>& fixedLeg,
		kMatrix<V>& swpRate
	);

	void exerciseValue(
		int			type,
		V			t,
		V			myKappa,
		V			y0t,
		kVector<V>& myX,
		kVector<V>& myYtil,
		kMatrix<V>& input
	);

	void updateValue(
		int			type,
		V			t,
		V			myKappa,
		V			y0t,
		kVector<V>& myX,
		kVector<V>& myYtil,
		kMatrix<V>& input
	);

	void Max(
		kMatrix<V>& x,
		kMatrix<V>& y
	);

	void Max(
		kMatrix<V>& x,
		V			y
	);

	void y0Step(
		V&	y0t,
		V	myKappa,
		V	etaX02,
		V	dt
	);

	V yield(V t);

	V nextUpdateTime(
		V t);

private:
	kVector<V> yieldtt;
};

template<class V>
kFuncHelp<V>::kFuncHelp(
	kVector<V>& timeline,
	kVector<V>& yields)
{
	myTimeline = timeline;
	myYields = yields;
	yieldtt.resize(myYields.size());
	kPolation::spline(myTimeline, myYields, myTimeline.size(), 1e100, 1e100, yieldtt);
}

template<class V>
void
kFuncHelp<V>::initEta(
	kMatrix<V>& _myEta2,
	kVector<V>& _myEta02
)
{
	myEta2 = _myEta2;
	myEta02 = _myEta02;
	return;
}

template<class V>
V
kFuncHelp<V>::G(
	V kappa,
	V t,
	V T) 
{
	if (kappa == V(0.0))
	{
		return T - t;
	}
	else
	{
		return (1 - exp(-kappa * (T - t))) / kappa;
	}
}

template<class V>
V
kFuncHelp<V>::bondRecon(
	V bondFactor,
	V valG,
	V x,
	V y
)
{
	V res = bondFactor * exp(-valG * x - 0.5 * valG * valG * y);
	return res;
}

template <class V>
void
kFuncHelp<V>::fixedCalc(
	V			t,
	V			T,
	V			myKappa,
	V			y0t,
	V			_tenor,
	V			_frequency,
	kVector<V>& myX,
	kVector<V>& myYtil,
	kMatrix<V>& fixedLeg
)
{
	int i, j, k;
	V y;

	V bondMaturity;
	V valG;
	V bondFactor;


	int myN = myX.size();
	int myM = myYtil.size();
	
	fixedLeg.resize(myN, myM);
	fixedLeg = 0.0;

	V dt = 1.0 / _frequency;
	int payments = (_tenor / dt) + 0.5;
	V yieldt = yield(t);

	bondMaturity = T;
	for (k = 1; k <= payments; ++k)
	{
		bondMaturity += dt;
		valG = G(myKappa, t, bondMaturity);
		bondFactor = exp(yieldt * t - yield(bondMaturity) * bondMaturity);
		for (i = 0; i < myN; ++i)
		{
			for (j = 0; j < myM; ++j)
			{
				y = myYtil(j) + y0t;
				fixedLeg(i, j) += dt * bondRecon(bondFactor, valG, myX(i), y);
			}
		}
	}
	return;
}

template<class V>
void 
kFuncHelp<V>::floatingCalc(
	V			t,
	V			T,
	V			myKappa,
	V			y0t,
	V			_tenor,
	kVector<V>& myX,
	kVector<V>& myYtil,
	kMatrix<V>& floatingLeg
)
{
	int i, j;
	V y;

	V bondMaturity = T + _tenor;
	V valG;
	V bondFactor;

	int  myN = myX.size();
	int  myM = myYtil.size();

	floatingLeg.resize(myN, myM);

	V yieldt = yield(t);
	
	// P(t, Tj)
	valG = G(myKappa, t, T);
	bondFactor = exp(yieldt * t - yield(T) * T);

	for (i = 0; i < myN; ++i) {
		for (j = 0; j < myM; ++j)
		{
			y = myYtil(j) + y0t;
			floatingLeg(i, j) = bondRecon(bondFactor, valG, myX(i), y);
		}
	}

	// -P(t,Tn)
	
	valG = G(myKappa, t, bondMaturity);
	bondFactor = exp(yieldt * t - yield(bondMaturity) * bondMaturity);

	for (i = 0; i < myN; ++i) {
		for (j = 0; j < myM; ++j)
		{
			y = myYtil(j) + y0t;
			floatingLeg(i, j) -= bondRecon(bondFactor, valG, myX(i), y);
		}
	}

	return;
}

// calculates the swaprate and outputs the corresponding fixed leg of the calculation.
template<class V>
void
kFuncHelp<V>::swpCalc(
	V			t,
	V			T,
	V			myKappa,
	V			y0t,
	V			_tenor,
	V			_frequency,
	kVector<V>& myX,
	kVector<V>& myYtil,
	kMatrix<V>& floatingLeg,
	kMatrix<V>& fixedLeg,
	kMatrix<V>& swpRate
	) 
{
	int i, j;

	int myN = myX.size();
	int myM = myYtil.size();

	
	swpRate.resize(myN, myM);

	// Fixed leg
	fixedCalc(t, T, myKappa, y0t, _tenor, _frequency, myX, myYtil, fixedLeg);

	// Floating leg
	floatingCalc(t, T, myKappa, y0t, _tenor, myX, myYtil, floatingLeg);
	

	for (i = 0; i < myN; ++i) {
		for (j = 0; j < myM; ++j)
		{
			swpRate(i, j) = floatingLeg(i, j) / fixedLeg(i, j);
		}
	}

	return;
}


template<class V>
void
kFuncHelp<V>::exerciseValue(
	int			type,
	V			t,
	V			myKappa,
	V			y0t,
	kVector<V>& myX,
	kVector<V>& myYtil,
	kMatrix<V>& input
)
{
	int i, j;
	V y;

	int myN = myX.size();
	int myM = myYtil.size();

	if (type == 0)
	{
		// ZCB
		V valG = G(myKappa, t, maturity);
		V bondFactor = exp(yield(t) * t - yield(maturity) * maturity);
		for (i = 0; i < myN; ++i) {
			for (j = 0; j < myM; ++j)
			{
				y = myYtil(j) + y0t;
				input(i, j) = bondRecon(bondFactor, valG, myX(i), y);
			}

		}
		return;
	}

	if (type == 1 || type == -1)
	{
		// swap
		kMatrix<V> fixedLeg;

		fixedCalc(t, t, myKappa, y0t, tenor, frequency, myX, myYtil, fixedLeg);
		floatingCalc(t, t, myKappa, y0t, tenor, myX, myYtil, input);

		for (i = 0; i < myN; ++i)
		{
			for (j = 0; j < myM; ++j)
			{
				input(i, j) -= strike * fixedLeg(i, j);
			}
		}

		if (type == -1)
		{
			for (i = 0; i < myN; ++i)
			{
				for (j = 0; j < myM; ++j)
				{
					input(i, j) = -input(i, j);
				}
			}
		}
		return;
	}


	if (type == 2 || type == -2)
	{
		input = 0.0;

		return;
	}
}


template<class V>
void
kFuncHelp<V>::updateValue(
	int			type,
	V			t,
	V			myKappa,
	V			y0t,
	kVector<V>& myX,
	kVector<V>& myYtil,
	kMatrix<V>& input
)
{
	int i, j;
	V y;

	int myN = myX.size();
	int myM = myYtil.size();

	if (type == 0)
	{
		// ZCB
		V valG = G(myKappa, t, maturity);
		V bondFactor = exp(yield(t) * t - yield(maturity) * maturity);
		for (i = 0; i < myN; ++i) {
			for (j = 0; j < myM; ++j)
			{
				y = myYtil(j) + y0t;
				input(i, j) = bondRecon(bondFactor, valG, myX(i), y);
			}

		}
		return;
	}

	if (type == 1 || type == -1)
	{
		// swap
		kMatrix<V> fixedLeg;

		fixedCalc(t, t, myKappa, y0t, tenor, frequency, myX, myYtil, fixedLeg);
		floatingCalc(t, t, myKappa, y0t, tenor, myX, myYtil, input);

		for (i = 0; i < myN; ++i)
		{
			for (j = 0; j < myM; ++j)
			{
				input(i, j) -= strike * fixedLeg(i, j);
			}
		}

		if (type == -1)
		{
			for (i = 0; i < myN; ++i)
			{
				for (j = 0; j < myM; ++j)
				{
					input(i, j) = -input(i, j);
				}
			}
		}
		return;
	}


	if (type == 2 || type == -2)
	{
		V bondMaturity = t + 1 / frequency;
		V valG = G(myKappa, t, bondMaturity);
		V bondFactor = exp(yield(t) * t - yield(bondMaturity) * bondMaturity);
		
		if (type == 2)
		{
			for (i = 0; i < myN; ++i) {
				for (j = 0; j < myM; ++j)
				{
					y = myYtil(j) + y0t;
					input(i, j) += 1 - (1 + strike / frequency) * bondRecon(bondFactor, valG, myX(i), y);
				}
			}
		}
		if (type == -2)
		{
			for (i = 0; i < myN; ++i) {
				for (j = 0; j < myM; ++j)
				{
					y = myYtil(j) + y0t;
					input(i, j) -= 1 - (1 + strike / frequency) * bondRecon(bondFactor, valG, myX(i), y);
				}
			}
		}

		return;
	}
}

// max between all elements returned in x
template<class V>
void 
kFuncHelp<V>::Max(
	kMatrix<V>& x,
	kMatrix<V>& y
)
{
	int i, j;
	int myN = x.rows();
	int myM = y.cols();

	for (i = 0; i < myN; ++i)
	{
		for (j = 0; j < myM; ++j)
		{
			x(i, j) = max(x(i, j), y(i, j));
		}
	}
}

template<class V>
void
kFuncHelp<V>::Max(
	kMatrix<V>& x,
	V			y
) 
{
	int i, j;
	
	int myN = x.rows();
	int myM = x.cols();

	for (i = 0; i < myN; ++i)
	{
		for (j = 0; j < myM; ++j)
		{
			x(i, j) = max(x(i, j), y);
		}
	}
}


template<class V>
void 
kFuncHelp<V>::y0Step(
	V&	y0t,
	V	myKappa,
	V	etaX02,
	V	dt
) 
{
	if (myKappa == 0.0) {
		y0t += etaX02 * (dt);
	}
	else{
		V factor = exp(-2 * myKappa * dt);
		V addedV = etaX02 * ((1 - factor) / (2 * myKappa));
		y0t = factor * y0t + addedV;
	}
}

template<class V>
V
kFuncHelp<V>::yield(
	V t)
{
	int i;
	V res;
	int n = myYields.size() - 1;
	kPolation::locate(myTimeline(), myTimeline.size(), t, i);

	if (i == -1) {
		return myYields(0);
	}
	if (i == n) {
		return myYields(n);
	}
	
	kPolation::splint(myTimeline, myYields, yieldtt, i, t, res);
	return res;
}

template<class V>
V
kFuncHelp<V>::nextUpdateTime(
	V t)
{
	t -= 1 / frequency;
	if (t < finalUpdateTime - 1E-8)
	{
		return -1;
	}

	return t;
}