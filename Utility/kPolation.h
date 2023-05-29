#pragma once

#include "kVector.h"
#include "kMatrix.h"
#include "kColVectorView.h"

namespace kPolation 
{
	// Locates the index for which interval the given x lies between in xx
	template<class V>
	void locate(
		kVectorView<V> xx,
		int n,
		V x,
		int& j
	)
	{
		int ju, jm, jl;

		jl = 0;
		ju = n + 1;

		while (ju-jl > 1)
		{
			// Bit shifting, fast way of dividing by 2
			jm = (ju + jl) >> 1;
			if (x >= xx[jm - 1])	jl = jm;
			else					ju = jm;
		}

		if (x == xx[0])				j = 0;
		else if (x == xx[n - 1])	j = n - 2;
		else						j = jl - 1;
	}

	// first part of finding spline values. Finding all the second derivatives
	template<class V>
	void spline(
		kVectorView<V>	x,
		kVectorView<V>	y,
		int				n,
		V				yp1,
		V				ypn,
		kVectorView<V>  y2
	)
	{
		int i, k;
		V p, qn, sig, un;
		kVector<V> u(n);

		if (yp1 > 0.99e30)
		{
			y2[0] = u[0] = V(0.0);
		}
		else 
		{
			y2[0] = V(-0.5);
			u[0] = (V(3.0) / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
		}
		for (i = 1; i < n - 1; ++i)
		{
			sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
			p = sig * y2[i - 1] + V(2.0);
			y2[i] = (sig - V(1.0)) / p;
			u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
			u[i] = (V(6.0) * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
		}
		if (ypn > 0.99e30)
		{
			qn = un = V(0.0);
		}
		else
		{
			qn = V(0.5);
			un = (V(3.0) / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
		}
		y2[n-1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + V(1.0));

		for (k = n - 2; k >= 0; --k)
			y2[k] = y2[k] * y2[k + 1] + u[k];
	}

	// first part of finding spline values. Finding all the second derivatives
	template<class V>
	void spline(
		kVector<V>&		x,
		kVector<V>&		y,
		int				n,
		V				yp1,
		V				ypn,
		kVector<V>&		y2
	)
	{
		spline(x(), y(), n, yp1, ypn, y2());
	}

	// first part of finding spline values. Finding all the second derivatives
	template<class V>
	void spline(
		kVector<V>&		x,
		kVectorView<V>	y,
		int				n,
		V				yp1,
		V				ypn,
		kVector<V>&		y2
	)
	{
		spline(x, y, n, yp1, ypn, y2());
	}

	// first part of finding spline values. Finding all the second derivatives
	template<class V>
	void spline(
		kVector<V>& x,
		kColVectorView<V>	y,
		int				n,
		V				yp1,
		V				ypn,
		kColVectorView<V>  y2
	)
	{
		int i, k;
		V p, qn, sig, un;
		kVector<V> u(n);

		if (yp1 > 0.99e30)
		{
			y2[0] = u[0] = V(0.0);
		}
		else
		{
			y2[0] = V(-0.5);
			u[0] = (V(3.0) / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
		}
		for (i = 1; i < n - 1; ++i)
		{
			sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
			p = sig * y2[i - 1] + V(2.0);
			y2[i] = (sig - V(1.0)) / p;
			u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
			u[i] = (V(6.0) * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
		}
		if (ypn > 0.99e30)
		{
			qn = un = V(0.0);
		}
		else
		{
			qn = V(0.5);
			un = (V(3.0) / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
		}
		y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + V(1.0));

		for (k = n - 2; k >= 0; --k)
			y2[k] = y2[k] * y2[k + 1] + u[k];
	}

	// computes spline value in given point
	template<class V>
	void splint(
		kVectorView<V>	xa,
		kVectorView<V>  ya,
		kVectorView<V>	y2a,
		int&			klo,
		V				x,
		V&				y
	) 
	{
		int khi;
		V h, b, a;

		while (xa[klo + 1] < x)
		{
			++klo;
		}
		khi = klo + 1;		
		
		h = xa[khi] - xa[klo];
		a = (xa[khi] - x) / h;
		b = (x - xa[klo]) / h;
		y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / V(6.0);
	}

	// computes spline value in given point
	template<class V>
	void splint(
		kVector<V>&			xa,
		kColVectorView<V>	ya,
		kColVectorView<V>	y2a,
		int&				klo,
		V					x,
		V&					y
	)
	{
		int khi;
		V h, b, a;

		while (xa[klo + 1] < x)
		{
			++klo;
		}
		khi = klo + 1;

		h = xa[khi] - xa[klo];
		a = (xa[khi] - x) / h;
		b = (x - xa[klo]) / h;
		y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / V(6.0);
	}

	// computes spline value in given point
	template<class V>
	void splint(
		kVector<V>&		xa,
		kVectorView<V>  ya,
		kVector<V>&		y2a,
		int&			klo,
		V				x,
		V&				y
	)
	{
		splint(xa, ya, y2a(), klo, x, y);
	}

	// computes spline value in given point
	template<class V>
	void splint(
		kVector<V>& xa,
		kVector<V>& ya,
		kVector<V>& y2a,
		int&		klo,
		V			x,
		V&			y
	)
	{
		splint(xa(), ya(), y2a(), klo, x, y);
	}

	// extrapolates with the given method
	template<class V>
	void extrapolate(
		int il,
		int iu,
		bool lower,
		int method,
		kVector<V>& trueX,
		kVector<V>& inX,
		kVectorView<V> inY,
		kVector<V>& res
	)
	{
		int i;
		if (method == 0)
		{
			V deriv;
			if (lower)
			{
				deriv = (inY(1) - inY(0)) / (inX(1) - inX(0));
				for (i = il; i < iu; ++i)
				{
					res(i) = inY(0) - deriv * (inX(0) - trueX(i));
				}
			}
			else
			{
				int N = inX.size() - 1;
				deriv = (inY(N) - inY(N - 1)) / (inX(N) - inX(N - 1));
				for (i = il; i < iu; ++i)
				{
					res(i) = inY(N) + deriv * (trueX(i) - inX(N));
				}
			}
			return;
		}


	}

	// interpolates with the given method
	template<class V>
	void interpolate(
		int il,
		int iu,
		int method,
		kVectorView<V> trueX,
		kVectorView<V> inX,
		kVectorView<V> inY,
		kVectorView<V> res
	)
	{
		int i, j;

		if (method == 0)
		{
			V deriv;
			j = 0;
			for (i = il; i < iu; ++i)
			{
				while (j < inX.size())
				{
					if (inX(j + 1) >= trueX(i)) {
						break;
					}
					++j;
				}
				deriv = (inY(j + 1) - inY(j)) / (inX(j + 1) - inX(j));
				res(i) = inY(j) + deriv * (trueX(i) - inX(j));
			}
			return;
		}

		if (method == 1)
		{
			int n = inX.size();
			V y;
			kVector<V> y2;
			y2.resize(n);
			int xl;
			locate(inX, n, trueX[il], xl);
			spline(inX, inY, n, 1.0e100, 1.0e100, y2());
			for (i = il; i < iu; ++i)
			{
				splint(inX, inY, y2(), xl, trueX[i], res[i]);
			}
			return;
		}
		
	}

	// Finds the right places to do inter- and extrapolation and returns res as a mix between the two methods
	template<class V>
	void polate(
		int			interpolation,
		int			extrapolation,
		kVector<V>& trueX,
		kVector<V>& inX,
		kVectorView<V> inY,
		kVector<V>& res
	)
	{

		int N = trueX.size();
		int M = inX.size();

		if (res.size() != N)
		{
			res.resize(N);
		}

		// extrapolate bounds
		int jl, ju;

		locate(trueX(), N, inX[0], jl);
		locate(trueX(), N, inX[M - 1], ju);

		extrapolate(0, jl + 1, true, extrapolation, trueX, inX, inY, res);
		extrapolate(ju + 1, N, false, extrapolation, trueX, inX, inY, res);

		interpolate(jl + 1, ju + 1, interpolation, trueX(), inX(), inY, res());

		return;
	}


	// Probability interpolation
	template<class V>
	void probInterpolation(
		int i,
		int j,
		kVector<V>& x,
		kVector<V>& y,
		kMatrix<V>& res
	)
	{
		V p1 = x[i + 1] / (x[i + 1] - x[i]);
		V p2 = y[j + 1] / (y[j + 1] - y[j]);

		res(i, j) = p1 * p2;
		res(i, j + 1) = p1 * (1 - p2);
		res(i + 1, j) = (1 - p1) * p2;
		res(i + 1, j + 1) = (1 - p1) * (1 - p2);
	}

	// Bilinear interpolation
	template<class V>
	V bilinearInterpolation(
		int i,
		int j,
		kVector<V>& x,
		kVector<V>& y,
		kMatrix<V>& z
	)
	{
		V res = V(0.0);
		
		V p1 = x[i + 1] / (x[i + 1] - x[i]);
		V p2 = y[j + 1] / (y[j + 1] - y[j]);

		res += p1 * p2 * z(i, j);
		res += p1 * (1 - p2)* z(i, j + 1);
		res += (1 - p1) * p2 * z(i + 1, j);
		res += (1 - p1) * (1 - p2) * z(i + 1, j + 1);

		return res;
	}
}
