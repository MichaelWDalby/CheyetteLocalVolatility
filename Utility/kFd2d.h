#pragma once
//	includes
#include "kFiniteDifference.h"
#include "kMatrixAlgebra.h"
#include "kFd1d.h"

//	class declaration
template <class V>
class kFd2d
{
public:
	//	init 
	void	init(
		int						numV,
		const kVector<V>&		x,
		const kVector<V>&		y,
		kVector<int>&			pointDisc);

	const kMatrix<V>&					r()		const { return myR; }
	const kMatrix<V>&					mu()	const { return myMu; }
	const kMatrix<V>&					var()	const { return myVar; }
	const kVector<V>&					x()		const { return myX; }

	kMatrix<V>&							r()		{ return myR; }
	kMatrix<V>&							mu()	{ return myMu; }
	kMatrix<V>&							var()	{ return myVar; }
	kVector<V>&							x()		{ return myX; }

	const kMatrix<V>&					nu()	const { return myNu; }
	const kVector<V>&					y()		const { return myY; }
	const kVector<kMatrix<V>>&			res()	const { return myRes; }

	kMatrix<V>&							nu()	{ return myNu; }
	kVector<V>&							y()		{ return myY; }
	kVector<kMatrix<V>>&				res()	{ return myRes; }

	//	operator
	void	calcAx(
		V						one,
		V						dtTheta,
		int						wind,
		bool					tr,
		int						yIndex,
		kMatrix<V>&				A) const;

	//	operator
	void	calcAy(
		V						one,
		V						dtTheta,
		int						wind,
		bool					tr,
		int						xIndex,
		kMatrix<V>&				A) const;

	void	rollBwd(
		V						dt,
		kVector<V>				theta,
		kVector<int>			wind,
		kVector<kMatrix<V>>&	res);

	//	roll fwd
	void	rollFwd(
		V						dt,
		kVector<V>				theta,
		kVector<int>			wind,
		kVector<kMatrix<V>>&	res);


private:
	
	//	r, mu, var, nu, x, y
	kVector<V>			myX, myY; 
	kMatrix<V>			myR, myVar, myMu, myNu;		

	//	diff operators
	kMatrix<V>	myDxd, myDxu, myDx, myDxx;
	kMatrix<V>	myDyd, myDyu, myDy;

	//	operator matrix
	kMatrix<V>	myAxe, myAxi, myAye, myAyi;

	//	helper
	kVector<V>				myVs, myWns, myWms;

	//	matrix of results
	kVector<kMatrix<V>> myRes;

	// dimensions
	int myN, myM = 0;

	void Inverse2d(
		const kMatrix<V>&	A,
		int					m1,
		int					m2,
		kVectorView<V>&		input,
		kVector<V>&			myWs
	);

	void Inverse2d(
		const kMatrix<V>& A,
		int					m1,
		int					m2,
		kColVectorView<V>& input,
		kVector<V>& myWs
	);
};


// Makes a linear solve with matrix and vector. Assumes x dim = 0 and tridag while y dim = 1 and five point
template <class V>
void
kFd2d<V>::Inverse2d(
	const kMatrix<V>&	A,
	int					m1,
	int					m2,
	kVectorView<V>&		input,
	kVector<V>&			myWs
) {	
	if (A.cols() == 3)
	{
		myVs = input;
		kMatrixAlgebra::tridag(A, myVs, input, myWs);
		return;
	}

	if (A.cols() == 5)
	{
		kMatrix<V> au = kMatrix<V>(A);
		kMatrix<V> al;
		kVector<int> indx;
		int d;

		kMatrixAlgebra::bandec(au, m1, m2, al, indx, d);
		kMatrixAlgebra::banbks(au, m1, m2, al, indx, input);
		return;
	}

	throw std::invalid_argument("Input is neither a 3 or 5 point discretisation");
}

// Makes a linear solve with matrix and vector. Assumes x dim = 0 and tridag while y dim = 1 and five point
template <class V>
void
kFd2d<V>::Inverse2d(
	const kMatrix<V>&	A,
	int					m1,
	int					m2,
	kColVectorView<V>&	input,
	kVector<V>&			myWs
) {
	if (A.cols() == 3)
	{
		myVs = input;
		kMatrixAlgebra::tridag(A, myVs, input, myWs);
		return;
	}

	if (A.cols() == 5)
	{
		kMatrix<V> au = kMatrix<V>(A);
		kMatrix<V> al;
		kVector<int> indx;
		int d;

		kMatrixAlgebra::bandec(au, m1, m2, al, indx, d);
		kMatrixAlgebra::banbks(au, m1, m2, al, indx, input);
		return;
	}

	throw std::invalid_argument("Input is neither a 3 or 5 point discretisation");
}

// init
template <class V>
void
kFd2d<V>::init(
	int						numV,
	const kVector<V>&		x,
	const kVector<V>&		y,
	kVector<int>&			pointDisc
	)
{
	myX = x;
	myY = y;
	int n = myX.size();
	int m = myY.size();
	myN = n;
	myM = m;

	// Result sizing
	myRes.resize(numV);
	for (int k = 0; k < numV; ++k) {
		myRes[k].resize(n, m);
	}

	//	resize params
	myR.resize(n, m, 0.0);
	myVar.resize(n, m, 0.0);
	myMu.resize(n, m, 0.0);

	myNu.resize(n, m, 0.0);


	if (pointDisc(0) < 5)
	{
		kFiniteDifference::dx(-1, myX, myDxd);
		kFiniteDifference::dx(0, myX, myDx);
		kFiniteDifference::dx(1, myX, myDxu);
		kFiniteDifference::dxx(myX, myDxx);
	}
	else
	{
		kFiniteDifference::dx5(-1, myX, myDxd);
		kFiniteDifference::dx5(0, myX, myDx);
		kFiniteDifference::dx5(1, myX, myDxu);
		kFiniteDifference::dxx5(myX, myDxx);
	}

	if (pointDisc(1) < 5)
	{
		kFiniteDifference::dx(-1, myY, myDyd);
		kFiniteDifference::dx(0, myY, myDy);
		kFiniteDifference::dx(1, myY, myDyu);
	}
	else
	{
		kFiniteDifference::dx5(-1, myY, myDyd);
		kFiniteDifference::dx5(0, myY, myDy);
		kFiniteDifference::dx5(1, myY, myDyu);
	}
	

	if (myX.empty()) return;
	if (myY.empty()) return;

	int numCx = myDx.cols();
	int numCy = myDy.cols();

	myAxe.resize(n, numCx);
	myAxi.resize(n, numCx);

	myAye.resize(m, numCy);
	myAyi.resize(m, numCy);
	
	myWns.resize(n);
	myWms.resize(m);
	

	//	done
	return;
}


//	construct operator
template <class V>
void
kFd2d<V>::calcAx(
	V						one,
	V						dtTheta,
	int						wind,
	bool					tr,
	int						yIndex,
	kMatrix<V>&				A) const
{
	//	helps
	int i, j;
	
	//	dims
	int n = myN;
	int numCols = myDx.cols();
	int m = myM;
	int mm = numCols / 2;

	//	resize
	A.resize(n, numCols);		

	//	wind
	const kMatrix<V>*	Dx = 0;
	if (wind < 0)		Dx = &myDxd;
	else if (wind == 0)	Dx = &myDx;
	else if (wind == 1)	Dx = &myDxu;

	//	fill
	for (i = 0; i < n; ++i)
	{
		if (wind > 1)
		{
			Dx = myMu(i, yIndex) < 0.0 ? &myDxd : &myDxu;
		}
		for (j = 0; j < numCols; ++j)
		{
			A(i, j) = dtTheta * (myMu(i, yIndex) * (*Dx)(i, j) + 0.5 * myVar(i, yIndex) * myDxx(i, j));
		}
		A(i, mm) += one - dtTheta * myR(i, yIndex);
	}

	if (tr) kMatrixAlgebra::transpose(mm, A);

	//	done
	return;
}

template <class V>
void
kFd2d<V>::calcAy(
	V					one,
	V					dtTheta,
	int					wind,
	bool				tr,
	int					xIndex,
	kMatrix<V>&			A) const
{
	//	dims
	int n = myN;
	int numCols = myDy.cols();
	int m = myM;
	int mm = numCols / 2;

	//	helps
	int i, j;


	//	resize
	A.resize(m ,numCols);

	//	wind
	const kMatrix<V>*	Dy = 0;
	if (wind < 0)		Dy = &myDyd;
	else if (wind == 0)	Dy = &myDy;
	else if (wind == 1)	Dy = &myDyu;

	//	fill
	for (i = 0; i < m; ++i)
	{
		if (wind > 1)
		{
			Dy = myNu(xIndex)(i) < 0.0 ? &myDyd : &myDyu;
		}
		for (j = 0; j < numCols; ++j)
		{
			A(i, j) = dtTheta * (myNu(xIndex)(i) * (*Dy)(i, j));
		}
		A(i, mm) += one;
	}

	if (tr) kMatrixAlgebra::transpose(mm, A);
	//	done
	return;
}

//	roll bwd
template <class V>
void
kFd2d<V>::rollBwd(
	V								dt,
	kVector<V>						theta,
	kVector<int>					wind,
	kVector<kMatrix<V>>&			res)
{
	//	helps
	int i, j, k;

	//	dims
	int n = myN;
	int m = myM;
	int mmx = myDx.cols() / 2;
	int mmy = myDy.cols() / 2;
	int numV = (int)res.size();
	
	kVectorView<V> resVec = kVectorView<V>();

	kMatrix<V> resX ,resY;
	
	for (k = 0; k < numV; ++k)
	{
		kColVectorView<V> resColVec(res[k], 0);	
		resY = resX = res[k];
		
		
		// Finding U
		// y step
		for (i = 0; i < n; ++i) {
			calcAy(1.0, dt, wind(1), false, i, myAye);
			resVec = res[k]()(i);
			myVs = resVec;
			kMatrixAlgebra::banmul(myAye, mmy, mmy, myVs, resVec);
		}

		// x step
		if (theta(0) != 1.0)
		{
			kColVectorView<V> resColVec(resX, 0);
		
			for (j = 0; j < m; ++j) {
				calcAx(0.0, dt * (1 - theta(0)), wind(0), false, j, myAxe);
				resColVec.setCol(j);
				myVs = resColVec;
				kMatrixAlgebra::banmul(myAxe, mmx, mmx, myVs, resColVec);
			}

			// combine
			res[k] += resX;
		}

		// implicit x step giving U
		if (theta(0) != 0.0)
		{
			kColVectorView<V> resColVec(res[k], 0);
			for (j = 0; j < m; ++j) {
				calcAx(1.0, -dt * theta(0), wind(0), false, j, myAxi);
				resColVec.setCol(j);
				Inverse2d(myAxi, mmx, mmx, resColVec, myWns);
			}
		}

		// Solving for V
		// explicit
		if (theta(1) != 0.0)
		{
			for (i = 0; i < n; ++i) {
				calcAy(0.0, - dt * theta(1), wind(1), false, i, myAye);
				resVec = resY()(i);
				myVs = resVec;
				kMatrixAlgebra::banmul(myAye, mmy, mmy, myVs, resVec);
			}

			// combine
			res[k] += resY;
		}

		// implicit
		if (theta(1) != 0.0)
		{
			for (i = 0; i < n; ++i) {
				calcAy(1.0, -dt * theta(1), wind(1), false, i, myAyi);
				resVec = res[k]()(i);
				Inverse2d(myAyi, mmy, mmy, resVec, myWms);
			}
		}
	}
	
	return;
}


//	roll fwd
template <class V>
void
kFd2d<V>::rollFwd(
	V								dt,
	kVector<V>						theta,
	kVector<int>					wind,
	kVector<kMatrix<V>>&			res)
{
	//	helps
	int i, j, k;

	//	dims
	int n = myN;
	int m = myM;
	int mmx = myDx.cols() / 2;
	int mmy = myDy.cols() / 2;
	int numV = (int)res.size();

	kVectorView<V> resVec = kVectorView<V>();

	kMatrix<V> resX, resY;

	for (k = 0; k < numV; ++k)
	{
		kColVectorView<V> resColVec(res[k], 0);

		// implicit
		if (theta(1) != 0.0)
		{
			for (i = 0; i < n; ++i) {
				calcAy(1.0, -dt * theta(1), wind(1), true, i, myAyi);
				resVec = res[k]()(i);
				Inverse2d(myAyi, mmy, mmy, resVec, myWms);
			}
		}

		resY = res[k];

		// Through U
		// implicit 
		if (theta(0) != 0.0)
		{
			kColVectorView<V> resColVec(res[k], 0);
			for (j = 0; j < m; ++j) {
				calcAx(1.0, -dt * theta(0), wind(0), true, j, myAxi);
				resColVec.setCol(j);
				Inverse2d(myAxi, mmx, mmx, resColVec, myWns);
			}
		}

		resX = res[k];

		// y step
		for (i = 0; i < n; ++i) {
			calcAy(1.0, dt, wind(1), true, i, myAye);
			resVec = res[k]()(i);
			myVs = resVec;
			kMatrixAlgebra::banmul(myAye, mmy, mmy, myVs, resVec);
		}

		// x step
		if (theta(0) != 1.0)
		{
			kColVectorView<V> resColVec(resX, 0);

			for (j = 0; j < m; ++j) {
				calcAx(0.0, dt * (1 - theta(0)), wind(0), true, j, myAxe);
				resColVec.setCol(j);
				myVs = resColVec;
				kMatrixAlgebra::banmul(myAxe, mmx, mmx, myVs, resColVec);
			}

			// combine
			res[k] += resX;
		}

		

		// Combining for V
		// explicit
		if (theta(1) != 0.0)
		{
			for (i = 0; i < n; ++i) {
				calcAy(0.0, -dt * theta(1), wind(1), true, i, myAye);
				resVec = resY()(i);
				myVs = resVec;
				kMatrixAlgebra::banmul(myAye, mmy, mmy, myVs, resVec);
			}

			// combine
			res[k] += resY;
		}
	}

	return;
}



