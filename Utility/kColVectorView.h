#pragma once
#include "kMatrix.h"

template<typename T>
class kColVectorView
{
public:
	//	declarations
	using value_type = T;

	//	trivi c'tors
	kColVectorView()noexcept = default;
	kColVectorView(const kColVectorView&)noexcept = default;
	kColVectorView(kColVectorView&&) noexcept = default;
	~kColVectorView()noexcept = default;

	//	c'tors will make a view on the rhs (i.e. updating values will update the rhs)
	kColVectorView(kMatrix<T>& _matrix, int _col) : matrix(_matrix), col(_col) {}
	

	//	element access
	const T& operator[](int i) const
	{
#ifdef _DEBUG
		if (i < 0 || i >= size()) throw std::runtime_error("kVectorView subscript out of range");
#endif
		return matrix(i, col);
	}

	T& operator[](int i)
	{
#ifdef _DEBUG
		if (i < 0 || i >= size()) throw std::runtime_error("kVectorView subscript out of range");
#endif
		return matrix(i, col);
	}

	const T& operator()(int i)	const { return operator[](i); }
	T& operator()(int i) { return operator[](i); }

	//	funcs
	void setMatrix(kMatrix<T>& _matrix) {
		matrix = _matrix;
	}
	void setCol(int _col) {
		col = _col;
	}

	int			size()				const { return (int)matrix.rows(); }
	bool		empty()				const { return size() == 0; }

private:
	kMatrix<T>& matrix;
	int			col;
};
