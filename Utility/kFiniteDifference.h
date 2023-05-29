#pragma once

//	desc:	various tools for finite difference

//	includes
#include "kVector.h"
#include "kMatrix.h"
#include "kInlines.h"

//	class declaration
class kFiniteDifference
{
public:
	//	1st order diff operator
	template <class V>
	static void	dx5(
		int					wind,
		const kVector<V>&	x,
		kMatrix<V>&			out)
	{
		out.resize(x.size(), 5);
		if (!x.size()) return;
		int n = x.size();
		
		// helpers
		V a, b, c, d, a2, b2, c2, d2, denominator;
		int i, j;

		// Lower boundary
		i = 0;

		if (wind >= 0) {
			WindStep5(i, true, x, out);
		}

		if (wind < 0) {
			for (j = 0; j < 5; ++j) out(i, j) = V(0.0);			
		}

		// next step
		++i;

		if (wind == 0) {
			a = x(i - 1) - x(i);
			b = x(i + 1) - x(i);
			c = x(i + 2) - x(i);			

			a2 = a * a;
			b2 = b * b;
			c2 = c * c;

			denominator = a * b * c * (a - b) * (a - c) * (b - c);

			out(i, 0) = V(0.0);
			out(i, 1) = b2 * c2 * (b - c) / denominator;
			out(i, 3) = a2 * c2 * (c - a) / denominator; 
			out(i, 4) = a2 * b2 * (a - b) / denominator;
			out(i, 2) = -(out(i, 0) + out(i, 1) + out(i, 3) + out(i, 4));
		}
		
		if (wind > 0) { 
			WindStep5(i, true, x, out); 
		}
		
		if (wind < 0)
		{
			a = x(i - 1) - x(i);

			denominator = a;

			out(i, 0) = V(0.0);
			out(i, 1) = 1 / denominator;
			out(i, 3) = V(0.0);
			out(i, 4) = V(0.0);
			out(i, 2) = -(out(i, 0) + out(i, 1) + out(i, 3) + out(i, 4));
		}
		

		// Middle
		++i;
		for (i; i < n - 2; ++i)
		{
			if (wind == 0) {
				a = x(i - 2) - x(i);
				b = x(i - 1) - x(i);
				c = x(i + 1) - x(i);
				d = x(i + 2) - x(i);
				a2 = a * a;
				b2 = b * b;
				c2 = c * c;
				d2 = d * d;

				denominator = a * b * c * d * (a - b) * (a - c) * (b - c) * (a - d) * (b - d) * (c - d);

				out(i, 0) = (b2 * c2 * d2 * (b2 * (-c + d) + c2 * (b - d) + d2 * (-b + c))) / denominator;
				out(i, 1) = (a2 * c2 * d2 * (a2 * (c - d) + c2 * (-a + d) + d2 * (a - c))) / denominator;
				out(i, 3) = (a2 * b2 * d2 * (a2 * (-b + d) + b2 * (a - d) + d2 * (-a + b))) / denominator;
				out(i, 4) = (a2 * b2 * c2 * (a2 * (b - c) + b2 * (-a + c) + c2 * (a - b))) / denominator;
				out(i, 2) = -(out(i, 0) + out(i, 1) + out(i, 3) + out(i, 4));
			}
			else { WindStep5(i, wind > 0, x, out); }
		}

		// Upper boundary
		if (wind == 0) {
			a = x(i - 2) - x(i);
			b = x(i - 1) - x(i);
			c = x(i + 1) - x(i);
			a2 = a * a;
			b2 = b * b;
			c2 = c * c;

			denominator = a * b * c * (a - b) * (a - c) * (b - c);

			out(i, 0) = b2 * c2 * (b - c) / denominator;
			out(i, 1) = a2 * c2 * (c - a) / denominator;
			out(i, 3) = a2 * b2 * (a - b) / denominator;
			out(i, 4) = V(0.0);
			out(i, 2) = -(out(i, 0) + out(i, 1) + out(i, 3) + out(i, 4));
		}
		
		
		if (wind > 0){
			a = x(i + 1) - x(i);

			denominator = a;

			out(i, 0) = V(0.0);
			out(i, 1) = V(0.0);
			out(i, 3) = 1 / denominator;
			out(i, 4) = V(0.0);
			out(i, 2) = -(out(i, 0) + out(i, 1) + out(i, 3) + out(i, 4));
		}
		
		if (wind < 0){
			WindStep5(i, false, x, out);
		}

		// last step
		++i;
		if (wind <= 0) {
			WindStep5(i, false, x, out);
		}

		if (wind > 0) {
			for (j = 0; j < 5; ++j) out(i, j) = V(0.0);
		}

		//	done
		return;
	}


	//	1st order diff operator
	template <class V>
	static void	dx(
		int					wind,
		const kVector<V>&	x,
		kMatrix<V>&			out)
	{
		out.resize(x.size(), 3);
		if(!x.size()) return;
		int n = x.size() - 1;
		if(wind>=0)
		{
			V dxu   =  x(1)-x(0);
			out(0,0) =  V(0.0);
			out(0,1) = -1.0/dxu;
			out(0,2) =  1.0/dxu;
		}
		else
		{
			out(0,0) = V(0.0);
			out(0,1) = V(0.0);
			out(0,2) = V(0.0);
		}

		for(int i = 1;i<n;++i)
		{
			V dxl = x(i)-x(i-1);
			V dxu = x(i+1)-x(i);
			if(wind<0)
			{
				out(i,0) = -1.0/dxl;
				out(i,1) =  1.0/dxl;
				out(i,2) =  V(0.0);
			}
			else if(wind==0)
			{
				out(i,0) = -dxu/dxl/(dxl+dxu);
				out(i,1) = (dxu/dxl-dxl/dxu)/(dxl+dxu);
				out(i,2) =  dxl/dxu/(dxl+dxu);
			}
			else
			{
				out(i,0) =  V(0.0);
				out(i,1) = -1.0/dxu;
				out(i,2) =  1.0/dxu;
			}
		}

		if(wind<=0)
		{
			V dxl   =  x(n)-x(n-1);
			out(n,0) = -1.0/dxl;
			out(n,1) =  1.0/dxl;
			out(n,2) =  V(0.0);
		}
		else
		{
			out(n,0) = V(0.0);
			out(n,1) = V(0.0);
			out(n,2) = V(0.0);
		}

		//	done
		return;
	}

	//	2nd order diff operator
	template <class V>
	static void	dxx(
		const kVectorView<V> x,
		kMatrix<V>& out)
	{
		out.resize(x.size(), 3);
		if (!x.size()) return;

		int n = x.size() - 1;

		out(0, 0) = V(0.0);
		out(0, 1) = V(0.0);
		out(0, 2) = V(0.0);

		for (int i = 1; i < n; ++i)
		{
			V dxl = x(i) - x(i - 1);
			V dxu = x(i + 1) - x(i);
			out(i, 0) = 2.0 / (dxl * (dxl + dxu));
			out(i, 1) = -(2.0 / dxl + 2.0 / dxu) / (dxl + dxu);
			out(i, 2) = 2.0 / (dxu * (dxl + dxu));
		}

		out(n, 0) = V(0.0);
		out(n, 1) = V(0.0);
		out(n, 2) = V(0.0);

		//	done
		return;
	}

	//	2nd order diff operator
	template <class V>
	static void	dxx(
		const kVector<V>&	x,
		kMatrix<V>&			out)
	{
		dxx(x(), out);
		return;
	}

	//	2nd order diff operator
	template <class V>
	static void	dxx5(
		const kVectorView<V> x,
		kMatrix<V>& out)
	{
		int i;
		V a, b, c, d, a2, b2, c2, d2, denominator;
		
		out.resize(x.size(), 5);
		if (!x.size()) return;

		int n = x.size() - 1;

		out(0, 0) = V(0.0);
		out(0, 1) = V(0.0);
		out(0, 2) = V(0.0);
		out(0, 3) = V(0.0);
		out(0, 4) = V(0.0);

		i = 1;

		a = x(i - 1) - x(i);
		b = x(i + 1) - x(i);
		c = x(i + 2) - x(i);
		a2 = a * a;
		b2 = b * b;
		c2 = c * c;

		denominator = a * b * c * (a - b) * (a - c) * (b - c);

		out(i, 0) = V(0.0);
		out(i, 1) = -2 * (b * c * (b2 - c2)) / denominator;
		out(i, 3) = -2 * (a * c * (-a2 + c2)) / denominator;
		out(i, 4) = -2 * (a * b * (a2 - b2)) / denominator;
		out(i, 2) = -(out(i, 0) + out(i, 1) + out(i, 3) + out(i, 4));

		for (i = 2; i < n - 1; ++i)
		{
			a = x(i - 2) - x(i);
			b = x(i - 1) - x(i);
			c = x(i + 1) - x(i);
			d = x(i + 2) - x(i);
			a2 = a * a;
			b2 = b * b;
			c2 = c * c;
			d2 = d * d;

			denominator = a * b * c * d * (a - b) * (a - c) * (b - c) * (a - d) * (b - d) * (c - d);

			out(i, 0) = -2 * b * c * d * (b * b2 * (-c2 + d2) + c * c2 * (b2 - d2) + d * d2 * (-b2 + c2)) / denominator;
			out(i, 1) = -2 * a * c * d * (a * a2 * (c2 - d2) + c * c2 * (-a2 + d2) + d * d2 * (a2 - c2)) / denominator;
			out(i, 3) = -2 * a * b * d * (a * a2 * (-b2 + d2) + b * b2 * (a2 - d2) + d * d2 * (-a2 + b2)) / denominator;
			out(i, 4) = -2 * a * b * c * (a * a2 * (b2 - c2) + b * b2 * (-a2 + c2) + c * c2 * (a2 - b2)) / denominator;
			out(i, 2) = -(out(i, 0) + out(i, 1) + out(i, 3) + out(i, 4));
		}

		i = n - 1;

		a = x(i - 2) - x(i);
		b = x(i - 1) - x(i);
		c = x(i + 1) - x(i);
		a2 = a * a;
		b2 = b * b;
		c2 = c * c;

		denominator = a * b * c * (a - b) * (a - c) * (b - c);

		
		out(i, 0) = -2 * (b * c * (b2 - c2)) / denominator;
		out(i, 1) = -2 * (a * c * (-a2 + c2)) / denominator;
		out(i, 3) = -2 * (a * b * (a2 - b2)) / denominator;
		out(i, 0) = V(0.0);
		out(i, 2) = -(out(i, 0) + out(i, 1) + out(i, 3) + out(i, 4));


		out(n, 0) = V(0.0);
		out(n, 1) = V(0.0);
		out(n, 2) = V(0.0);
		out(n, 3) = V(0.0);
		out(n, 4) = V(0.0);

		//	done
		return;
	}

	//	2nd order diff operator
	template <class V>
	static void	dxx5(
		const kVector<V>& x,
		kMatrix<V>& out)
	{
		dxx5(x(), out);
		return;
	}

	//	smooth call payoff function
	template <class V>
	static V	smoothCall(
		const V&	xl,
		const V&	xu,
		const V&	strike)
	{
		V res;
		if(xu<=strike)
		{
			res = 0.0;
		}
		else if(strike<=xl)
		{
			res = 0.5 * (xl + xu) - strike;
		}
		else
		{
			res = 0.5 * kInlines::sqr(xu - strike) / (xu - xl);
		}

		//	done
		return res;
	}

	//	smooth digital payoff function
	template <class V>
	static V	smoothDigital(
		const V& xl,
		const V& xu,
		const V& strike)
	{
		V res;
		if (xu <= strike)
		{
			res = 0.0;
		}
		else if (strike <= xl)
		{
			res = 1.0;
		}
		else
		{
			res = (xu-strike)/(xu-xl);
		}

		//	done
		return res;
	}

private:
	// One step when winding
	template <class V>
	static void WindStep5(
		int i,
		bool forward,
		const kVector<V>& x,
		kMatrix<V>& out
	) 
	{
		V a, b, a2, b2, denominator;

		if (forward)
		{
			a = x(2 + i) - x(i);
			b = x(1 + i) - x(i);
		}
		else
		{
			a = x(-2 + i) - x(i);
			b = x(-1 + i) - x(i);
		}
		a2 = a * a;
		b2 = b * b;

		denominator = a * b * (a - b);

		if (forward)
		{
			out(i, 0) = V(0.0);
			out(i, 1) = V(0.0);
			out(i, 3) = a2 / denominator;
			out(i, 4) = -b2 / denominator;
			out(i, 2) = -(out(i, 0) + out(i, 1) + out(i, 3) + out(i, 4));
		}
		else
		{
			out(i, 0) = -b2 / denominator;
			out(i, 1) = a2 / denominator;
			out(i, 3) = V(0.0);
			out(i, 4) = V(0.0);
			out(i, 2) = -(out(i, 0) + out(i, 1) + out(i, 3) + out(i, 4));
		}
	}
};
