// X = done
// * = to discuss
// E = exercise for students
// N = don't do
#ifndef _PMATRIX_
#define _PMATRIX_
// TODO mostrare i seguenti metodi nel lab del 12/11
// pmatrix<double,3> M;
// M = {1.,2.,3.,4.,5.,6.,7.,8.,9.};
//
// [E] constructor with initializer_list
// [E] overloading of () operator (m(i,j))
// [E] overloading of << operator to print a matrix
// [E] matrix addition and subtraction
//
// =================================================
//
// mostrare/far implementare nel lab del 12/11
// [E] += and -= operators
// [E] matrix divided by scalar (also /=)
// [E] matrix times scalar (also *=)
// [E] scalar times matrix <--- with friend function
// [E] matrix times vector
// [E] vector times matrix
// [E] matrix times matrix (also *=)
// [E] transpose
// [E] operator == to pmatrix class
//
#include "./pvector.hpp"
template <typename ntype, int NT>
class pmatrixq
{
	ntype m[NT][NT];

public:
	~pmatrixq()
	{
	}
	pmatrixq()
	{
		for (int i = 0; i < NT; i++)
			for (int j = 0; j < NT; j++)
				m[i][j] = 0;
	}

	pmatrixq(std::initializer_list<ntype> list)
	{
		int cc = 0;
		// WRITE THIS
		for (auto el : list)
		{
			m[cc / NT][cc % NT] = el;
			cc++;
		}
	}

	// ============================================================
	// overloading of () operator to get/set (i,j)-th element
	// m(i,j)=2.0;
	// std::cout << m(i, j);
	ntype &operator()(int i, int j)
	{
		// write your code here
		return m[i][j];
	}

	ntype operator()(int i, int j) const
	{
		// write your code here
		return m[i][j];
	}
	// =============================================================
	//      std::cout << M;
	friend std::ostream &operator<<(std::ostream &os, const pmatrixq &m)
	{
		// write your code here
		int i, j;
		os << "{";
		for (i = 0; i < NT; i++)
		{
			os << "{";
			for (j = 0; j < NT; j++)
			{
				os << m(i, j);
				if (j < NT - 1)
					os << ",";
			}
			os << "}";
			if (i < NT - 1)
				os << ",";
		}
		os << "}";
		return os;
	}
	// addition
	// note that omitting the const modifier in the argument, m2 can not be a temporary object like pmatrixq<double,3>()
	pmatrixq operator+(const pmatrixq &m2) const
	{
		// mt = m1+m2
		pmatrixq<ntype, NT> mt;
		for (auto i = 0; i < NT; i++)
		{
			for (auto j = 0; j < NT; j++)
			{
				mt(i, j) = m[i][j] + m2.m[i][j];
			}
		}
		return mt;
	}
	// subtraction
	// note that omitting the const modifier in the argument, m2 can not be a temporary object like pmatrixq<double,3>()
	pmatrixq operator-(const pmatrixq &m2) const
	{
		// mt = m1+m2
		pmatrixq mt;
		for (auto i = 0; i < NT; i++)
		{
			for (auto j = 0; j < NT; j++)
			{
				mt(i, j) = m[i][j] - m2.m[i][j];
			}
		}
		return mt;
	}
	pmatrixq<ntype, NT> &operator+=(const pmatrixq<ntype, NT> &m2)
	{
		return (*this = *this + m2);
	}

	pmatrixq<ntype, NT> &operator-=(const pmatrixq<ntype, NT> &m2)
	{
		return (*this = *this - m2);
	}

	pmatrixq<ntype, NT> &operator/=(const ntype &s)
	{
		for (auto i = 0; i < NT; i++)
		{
			for (auto j = 0; j < NT; j++)
			{
				m[i][j] /= s;
			}
		}
		return (*this);
	}

	pmatrixq<ntype, NT> &operator*=(const ntype &s)
	{
		for (auto i = 0; i < NT; i++)
		{
			for (auto j = 0; j < NT; j++)
			{
				m[i][j] *= s;
			}
		}
		return (*this);
	}
	//
	// matrix times scalar
	pmatrixq<ntype, NT> operator*(const ntype &s) const
	{
		pmatrixq<ntype, NT> mt;
		for (auto i = 0; i < NT; i++)
		{
			for (auto j = 0; j < NT; j++)
			{
				mt(i, j) = m[i][j];
				mt(i, j) *= s;
			}
		}
		return mt;
	}

	// scalar times matrix
	friend pmatrixq<ntype, NT> operator*(const ntype &s, const pmatrixq<ntype, NT> &m2)
	{
		pmatrixq<ntype, NT> mt;
		for (auto i = 0; i < NT; i++)
		{
			for (auto j = 0; j < NT; j++)
			{
				mt(i, j) = m2.m[i][j] * s;
			}
		}
		return mt;
	}

	// matrix times vector
	// note that omitting the const modifier in the argument, v2 can not be a temporary object like
	// pvector<double,3>(1,2,3)
	pvector<ntype, NT> operator*(const pvector<ntype, NT> &v2) const
	{
		// mt = (*this)*v2
		pvector<ntype, NT> vt;
		for (auto i = 0; i < NT; i++)
		{
			vt(i) = 0.0;
			for (auto j = 0; j < NT; j++)
			{
				vt(i) += m[i][j] * v2(j);
			}
		}
		return vt;
	}

	// vector times matrix
	// note that omitting the const modifier in the argument, v1 and m2 can not be temporary objects
	friend pvector<ntype, NT> operator*(const pvector<ntype, NT> &v1, const pmatrixq &m2)
	{
		// vt = m1*v2
		pvector<ntype, NT> vt;
		for (auto i = 0; i < NT; i++)
		{
			vt(i) = 0.0;
			for (auto j = 0; j < NT; j++)
			{
				vt(i) += v1(j) * m2.m[j][i];
			}
		}
		return vt;
	}

	// matrix times matrix
	// m1 è *this
	// m1_{il}*m2_{lj}
	// note that omitting the const modifier in the argument, m2 can not be a temporary object like pmatrixq<double,3>()
	pmatrixq<ntype, NT> operator*(const pmatrixq<ntype, NT> &m2) const
	{
		// mt = m1*m2
		pmatrixq<ntype, NT> mt;
		for (auto i = 0; i < NT; i++)
		{
			for (auto j = 0; j < NT; j++)
			{
				mt(i, j) = 0.0;
				for (auto l = 0; l < NT; l++)
				{
					mt(i, j) += (*this).m[i][l] * m2.m[l][j];
				}
			}
		}
		return mt;
	}

	pmatrixq<ntype, NT> transpose(void)
	{
		pmatrixq<ntype, NT> mt;
		for (auto i = 0; i < NT; i++)
		{
			for (auto j = 0; j < NT; j++)
			{
				mt(i, j) = m[j][i];
			}
		}
		return mt;
	}

	bool operator==(const pmatrixq<ntype, NT> &m2) const
	{
		for (auto i = 0; i < NT; i++)
		{
			for (auto j = 0; j < NT; j++)
			{
				if (m2.m[i][j] != m[i][j])
				{
					return 0;
				}
			}
		}
		return 1;
	}
};
#if 0
// =============================================================
  template<typename ntype, int NT>
  std::ostream& operator<<(std::ostream& os, const pmatrixq<ntype,NT>& m)
    {
      // write your code here
      int i, j;
      os << "{";
      for (i=0; i < NT; i++)
        {
          os << "{";
          for (j=0; j < NT; j++)
            {
              os << m(i,j);
              if (j < NT-1)
                os << ",";
            }
          os << "}";
          if (i < NT-1)
            os << ",";
        }
      os << "}";
      return os;
    }
#endif
#endif
