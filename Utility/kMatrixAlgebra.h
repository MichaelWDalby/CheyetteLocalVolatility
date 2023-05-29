#pragma once

//	includes
#include "kMatrix.h"
#include "kInlines.h"

//	class
namespace kMatrixAlgebra
{
	//	Decomposition of matrix a to upper and lower. Upper returned in a and lower in al
	template <class V>
	void bandec(
		kMatrix<V>&		a, 
		int				m1, 
		int				m2, 
		kMatrix<V>&		al, 
		kVector<int>&	indx, 
		int&			d)
	{
		int n = a.rows();
		al.resize(n, a.cols());
		indx.resize(n);
		
		int i, j, k, l;
		V dum;
		
		int mm = m1 + m2 + 1;
		l = m1;
		for (i = 0; i < m1; ++i)
		{
			for (j = m1 - i; j < mm; ++j) a(i, j - l) = a(i, j);
			--l;
			for (j = mm - l - 1; j < mm; ++j) a(i, j) = V(0.0);
		}
		d = 1;
		l = m1;
		for (k = 0; k < n; ++k) {
			dum = a(k, 0);
			i = k;
			if (l < n) ++l;
			for (j = k + 1; j < l; ++j)
			{
				if (kInlines::abs(a(j, 0)) > kInlines::abs(dum)) {
					dum = a(j, 0);
					i = j;
				}
			}
			indx(k) = i;
			
			if (dum == 0.0) {
				a(k, 0) = 1.0E-21;
			}
			

			if (i != k) {
				d = -d;
				for (j = 0; j < mm; ++j) kInlines::swap(a(k, j), a(i, j));
			}
			for (i = k + 1; i < l; ++i) {
				dum = a(i, 0) / a(k, 0);
				al(k, i - k) = dum;
				for (j = 1; j < mm; ++j) a(i, j - 1) = a(i, j) - dum * a(k, j);
				a(i, mm - 1) = V(0.0);
			}
		}

	}

	// Solving matrix equation on compact form by using the given LU decomposition
	template <class V>
	void banbks(
		kMatrix<V>&		a,
		int				m1,
		int				m2,
		kMatrix<V>&		al,
		kVector<int>&	indx,
		kVectorView<V>&	b)
	{
		int n = a.rows();

		int i, k, l;
		int mm;
		V dum;

		mm = m1 + m2 + 1;
		l = m1;
		for (k = 0; k < n; ++k) {
			i = indx(k);
			if (i != k) kInlines::swap(b(k), b(i));
			if (l < n) ++l;
			for (i = k + 1; i < l; ++i) b(i) -= al(k, i - k) * b(k);
		}
		l = 1;
		for (i = n - 1; i >= 0; --i) {
			dum = b(i);
			for (k = 1; k < l; k++) dum -= a(i, k) * b(k + i);
			b(i) = dum / a(i, 0);
			if (l < mm) ++l;
		}
	}
	
	// Solving matrix equation on compact form by using the given LU decomposition
	template <class V>
	void banbks(
		kMatrix<V>&		a, 
		int				m1, 
		int				m2, 
		kMatrix<V>&		al, 
		kVector<int>&	indx, 
		kVector<V>&		b) 
	{
		banbks(a, m1, m2, al, indx, b());
	}

	// Solving matrix equation on compact form by using the given LU decomposition
	template <class V>
	void banbks(
		kMatrix<V>&			a,
		int					m1,
		int					m2,
		kMatrix<V>&			al,
		kVector<int>&		indx,
		kColVectorView<V>&	b)
	{
		int n = a.rows();

		int i, k, l;
		int mm;
		V dum;

		mm = m1 + m2 + 1;
		l = m1;
		for (k = 0; k < n; ++k) {
			i = indx(k);
			if (i != k) kInlines::swap(b(k), b(i));
			if (l < n) ++l;
			for (i = k + 1; i < l; ++i) b(i) -= al(k, i - k) * b(k);
		}
		l = 1;
		for (i = n - 1; i >= 0; --i) {
			dum = b(i);
			for (k = 1; k < l; k++) dum -= a(i, k) * b(k + i);
			b(i) = dum / a(i, 0);
			if (l < mm) ++l;
		}
	}
	
	//	mat mult res = A * B
	template <class U, class V, class W>
	void mmult(
		const kMatrix<U>& a,
		const kMatrix<V>& b,
		kMatrix<W>& ab)
	{
		//	dims
		int m1 = a.rows();
		int m2 = a.cols();
		int m3 = b.cols();

		ab.resize(m1, m3, 0.0);

		//	calc
		for (int i = 0; i < m1; ++i)
		{
			for (int k = 0; k < m2; ++k)
			{
				for (int j = 0; j < m3; ++j)
					ab(i, j) += a(i, k) * b(k, j);
			}
		}

		//	done
		return;
	}

	//	tridag: solves A u = r when A is tridag
	template <class V>
	void	tridag(
		const kMatrixView<V>	A,		//	n x 3
		const kVectorView<V>	r,
		kVectorView<V>			u,
		kVectorView<V>			gam)
	{
		//	helps
		V bet;
		int j;

		//	dim
		int n = A.rows();

		//	go
		u(0) = r(0) / (bet = A(0, 1));
		for (j = 1; j < n; ++j)
		{
			gam(j) = A(j - 1, 2) / bet;
			bet = A(j, 1) - A(j, 0) * gam(j);
			u(j) = (r(j) - A(j, 0) * u(j - 1)) / bet;
		}
		for (j = n - 2; j >= 0; --j)
		{
			u(j) -= gam(j + 1) * u(j + 1);
		}

		//	done
		return;
	}


	//	tridag: solves A u = r when A is tridag
	template <class V>
	void	tridag(
		const kMatrix<V>		A,		//	n x 3
		const kVector<V>		r,
		kVectorView<V>			u,
		kVector<V>				gam)
	{
		tridag(A(), r(), u, gam());

		//	done
		return;
	}


	//	tridag: solves A u = r when A is tridag
	template <class V>
	void	tridag(
		const kMatrix<V>& A,		//	n x 3
		const kVector<V>& r,
		kVector<V>& u,
		kVector<V>& gam)
	{
		//	dim
		int n = A.rows();

		//	check dim
		if (u.size() < n) u.resize(n);
		if (gam.size() < n) gam.resize(n);

		tridag(A(), r(), u(), gam());
		//	done
		return;
	}

	//	tridag: solves A u = r when A is tridag
	template <class V>
	void	tridag(
		const kMatrixView<V>	A,		//	n x 3
		const kVectorView<V>	r,
		kColVectorView<V>			u,
		kVectorView<V>		gam)
	{
		//	helps
		V bet;
		int j;

		//	dim
		int n = A.rows();

		//	go
		u(0) = r(0) / (bet = A(0, 1));
		for (j = 1; j < n; ++j)
		{
			gam(j) = A(j - 1, 2) / bet;
			bet = A(j, 1) - A(j, 0) * gam(j);
			u(j) = (r(j) - A(j, 0) * u(j - 1)) / bet;
		}
		for (j = n - 2; j >= 0; --j)
		{
			u(j) -= gam(j + 1) * u(j + 1);
		}

		//	done
		return;
	}

	//	tridag: solves A u = r when A is tridag
	template <class V>
	void	tridag(
		const kMatrix<V>		A,		//	n x 3
		const kVector<V>		r,
		kColVectorView<V>		u,
		kVector<V>				gam)
	{
		tridag(A(), r(), u, gam());

		//	done
		return;
	}

	//	band diagonal matrix vector multiplication
	template <class V>
	void banmul(
		const kMatrixView<V> A,
		int					 m1,
		int					 m2,
		const kVectorView<V> b,
		kVectorView<V>		 x)
	{
		int n = A.rows() - 1;
		V xi;
		for (int i = 0; i <= n; ++i)
		{
			int jl = max<int>(0, i - m1);
			int ju = min<int>(i + m2, n);
			xi = 0.0;
			for (int j = jl; j <= ju; ++j)
			{
				int k = j - i + m1;
				xi += A(i, k) * b(j);
			}
			x(i) = xi;
		}

		//	done
		return;
	}

	//	band diagonal matrix vector multiplication
	template <class V>
	void banmul(
		const kMatrixView<V> A,
		int					 m1,
		int					 m2,
		const kVectorView<V> b,
		kColVectorView<V>	 x)
	{
		int n = A.rows() - 1;
		V xi;
		for (int i = 0; i <= n; ++i)
		{
			int jl = max<int>(0, i - m1);
			int ju = min<int>(i + m2, n);
			xi = 0.0;
			for (int j = jl; j <= ju; ++j)
			{
				int k = j - i + m1;
				xi += A(i, k) * b(j);
			}
			x(i) = xi;
		}

		//	done
		return;
	}

	//	band diagonal matrix vector multiplication
	template <class V>
	void banmul(
		const kMatrix<V>		A,		
		int						m1,
		int						m2,
		const kVector<V>		b,
		kVectorView<V>			x)
	{
		banmul(A(), m1, m2, b(), x);
		return;
	}

	//	band diagonal matrix vector multiplication
	template <class V>
	void banmul(
		const kMatrix<V>&	A,
		int					m1,
		int					m2,
		const kVector<V>&	b,
		kColVectorView<V>&	x)
	{
		banmul(A(), m1, m2, b(), x);
		return;
	}

	//	band diagonal matrix vector multiplication
	template <class V>
	void banmul(
		const kMatrix<V>&	A,		
		int					m1,
		int					m2,
		const kVector<V>&	b,
		kVector<V>&			x)
	{
		int n = A.rows();
		if(x.size() < n) x.resize(n);

		banmul(A(), m1, m2, b(), x());

		//	done
		return;
	}


	//	transposing a banddiagonal matrix
	template <class V>
	void	transpose(
		const int			mm,
		kMatrix<V>& A)
	{
		int n = A.rows() - 1;
		int i, j, k, l;
		int jl;
		for (i = 0; i <= n; ++i)
		{
			jl = max(0, i - mm);
			for (j = jl; j < i; ++j)
			{
				k = j - i + mm;
				l = i - j + mm;
				kInlines::swap(A(i, k), A(j, l));
			}
		}

		//	done
		return;
	}


	// Transposes makes a transposed copy of a matrix. (Not for compact form)
	template<class V>
	void
	transposeCopy(
		kMatrix<V>& M,
		kMatrix<V>& Mt
		) 
	{
		int n = M.rows();
		int m = M.cols();
		Mt.resize(m, n);

		for (int i = 0; i < n; ++i){
			for (int j = 0; j < m; ++j) {
				Mt(j, i) = M(i, j);
			}
		}
		return;
	}

	// sum product of two matrices assumed to have same dimensions
	template<class V>
	V
	sumProduct(
			kMatrix<V>& A,
			kMatrix<V>& B
		)
	{
		int n = A.rows();
		int m = B.cols();

		V res = V(0.0);

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				res += A(i, j) * B(i, j);
			}
		}
	}
}


