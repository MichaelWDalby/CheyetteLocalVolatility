#pragma once

#include "kMatrix.h"
#include "kVector.h"
#include "kSpecialFunction.h"
#include "kPolation.h"
#include <cmath>
#include "kBachelier.h"

template<class V>
class kSwaption {
public:	
	kMatrix<V> impliedVol;
	kVector<V> relativeStrike, tenor, frequency, maturity;

	kSwaption(	kMatrix<V>& _impliedVol,
				kVector<V>& _relativeStrike,
				kVector<V>& _tenor,
				kVector<V>& _frequency,
				kVector<V>& _maturity);

	void bachSetPrice(
		int				index,
		int				interpolation,
		int				extrapolation,
		V				s0,
		V				_maturity,
		kVector<V>&		_strikes,
		kVector<V>&		res
	);
private:
	kMatrix<V> impliedVoltt;
};

template<class V>
kSwaption<V>::kSwaption(kMatrix<V>& _impliedVol, 
						kVector<V>& _relativeStrike, 
						kVector<V>& _tenor, 
						kVector<V>& _frequency, 
						kVector<V>& _maturity)
	: impliedVol(_impliedVol), relativeStrike(_relativeStrike), tenor(_tenor), frequency(_frequency), maturity(_maturity)
{
	int n = impliedVol.rows();
	int m = impliedVol.cols();
	
	impliedVoltt.resize(n, m);
	for (int i = 0; i < m; ++i)
	{
		kColVectorView<V> impliedVolCol(impliedVol, i);
		kColVectorView<V> impliedVolttCol(impliedVoltt, i);
		kPolation::spline(maturity, impliedVolCol, n, 1e100, 1e100, impliedVolttCol);
	}
}

template<class V>
void
kSwaption<V>::bachSetPrice(
	int				index,
	int				interpolation,
	int				extrapolation,
	V				s0,
	V				_maturity,
	kVector<V>&		_strikes,
	kVector<V>&		res
)
{
	int i, j;
	int n = _strikes.size();
	int m = relativeStrike.size();
	
	kVector<V> impliedVolV(m);
	j = 0;
	if (_maturity > maturity(index) || _maturity < maturity(0))
	{
		_maturity = maturity(index);
	}

	for (i = 0; i < m; ++i)
	{
		kColVectorView<V> impliedVolCol(impliedVol, i);
		kColVectorView<V> impliedVolttCol(impliedVoltt, i);
		kPolation::splint(maturity, impliedVolCol, impliedVolttCol, j, _maturity, impliedVolV(i));
	}

	res.resize(n);

	kVector<V> sigma;
	kVector<V> strikes;

	strikes.resize(m);
	for (i = 0; i < m; ++i)
	{
		strikes(i) = s0 + relativeStrike(i);
	}

	// extra- and interpolation
	kPolation::polate(interpolation, extrapolation, _strikes, strikes, impliedVolV(), sigma);
	
	for (i = 0; i < n; ++i) 
	{
		res(i) = kBachelier::call(_maturity, _strikes(i), s0, sigma(i));
	}

	return;
}