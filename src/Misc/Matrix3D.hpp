// -*- Mode : c++; tab-width: 3; c-tab-always-indent: t; indent-tabs-mode: nil; c-basic-offset: 3 -*-
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 18-Dec-02 at 11:08:22
//  LAST-MOD: 22-Jan-03 at 17:31:19 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file Matrix3D.hpp

3*3 Matric class, row-oriented, 0-based [i][j] indexing.

*/
// DESCRIP-END.
//



#ifndef Matrix3D_hpp
#define Matrix3D_hpp

#ifdef USE_STD_INCLUDES
#include <iostream>
#include <strstream>
#else
#include <iostream.h>
#include <strstream.h>
#endif

#ifdef USE_C_PREFIX_INCLUDES
#include <cassert>
#include <cstdlib>
#else
#include <assert.h>
#include <stdlib.h>
#endif

#include "Mesquite.hpp"
#include "Vector3D.hpp"

namespace Mesquite
{

  /*! \class Matrix3D
      \brief 3*3 Matric class, row-oriented, 0-based [i][j] indexing.*/
  class Matrix3D 
  {
  protected:
    double v_[9];                  
    double* row_[3];

    // internal helper function to create the array
    // of row pointers

    void initialize()
    {
      row_[0] = v_;
      row_[1] = v_+3;
      row_[2] = v_+6;
    }
   
    void copy(const double*  v)
    { memcpy(v_, v, 9*sizeof(double)); }

    void set(const double& val)
    {
      v_[0]=val;  v_[1]=val;  v_[2]=val;
      v_[3]=val;  v_[4]=val;  v_[5]=val;
      v_[6]=val;  v_[7]=val;  v_[8]=val;
    }

    void set_values(const char *s)
    {
      //        initialize(M,N);
      std::istrstream ins(s);
      ins>>v_[0];  ins>>v_[1];  ins>>v_[2]; 
      ins>>v_[3];  ins>>v_[4];  ins>>v_[5]; 
      ins>>v_[6];  ins>>v_[7];  ins>>v_[8]; 
    }
    
  public:

    operator double**(){ return  row_; }

    size_t size() const { return 9; }

    // constructors
    Matrix3D() // : row_[0](v_), row_[1](v_+3), row_[2](v_+6) {}
    { initialize(); }
    
    Matrix3D(const Matrix3D &A)
    {
      initialize();
      copy(A.v_);
    }

    Matrix3D(const double& value)
    {
      initialize();
      set(value);
    }

    Matrix3D(const double* v)
    {
      initialize();
      copy(v);
    }

    Matrix3D(const char *s)
    {
      initialize();
      set_values(s);
    }

    // destructor
    ~Matrix3D() { }

    // assignments
    Matrix3D& operator=(const Matrix3D &A)
    {
      if (v_ == A.v_)
        return *this;
      copy(A.v_);
      return *this;
    }
        
    Matrix3D& operator=(const double& scalar)
    { 
      set(scalar); 
      return *this;
    }

    Matrix3D& operator=(const char* s)
    { 
      set_values(s); 
      return *this;
    }

    // Matrix Operators
    friend bool operator==(const Matrix3D &lhs, const Matrix3D &rhs);
    friend bool operator!=(const Matrix3D &lhs, const Matrix3D &rhs);
    friend Matrix3D operator+(const Matrix3D &A, const Matrix3D &B);
    friend Matrix3D operator-(const Matrix3D &A, const Matrix3D &B);
    friend Matrix3D operator*(const Matrix3D &A, const Matrix3D &B);
    friend Matrix3D mult_element(const Matrix3D &A, const Matrix3D &B);
    friend Matrix3D transpose(const Matrix3D &A);
    friend int matmult(Matrix3D& C, const Matrix3D  &A, const Matrix3D &B);
    friend Vector3D operator*(const Matrix3D  &A, const Vector3D &x);
    friend Vector3D operator*(const Vector3D &x, const Matrix3D  &A);
    Matrix3D& operator+=(const Matrix3D &rhs);
    Matrix3D plus_transpose(const Matrix3D &B);

    size_t num_rows() const { return 3; }
    size_t num_cols() const { return 3; }

    //! returns a pointer to a row.
    inline double* operator[](size_t i)
    {
#ifdef MSQ_DEBUG
      assert(0<=i);
      assert(i < 3) ;
#endif
      return row_[i];
    }

    //! returns a pointer to a row.
    inline const double* operator[](size_t i) const
    {
#ifdef MSQ_DEBUG
      assert(0<=i);
      assert(i < 3) ;
#endif
      return row_[i];
    }

  };


  /* ***********  I/O  **************/

  inline std::ostream& operator<<(std::ostream &s, const Matrix3D &A)
  {
    for (size_t i=0; i<3; ++i)
      {
        for (size_t j=0; j<3; ++j)
          s << A[i][j] << " ";
        s << "\n";
      }
    return s;
  }

  inline std::istream& operator>>(std::istream &s, Matrix3D &A)
  {
    for (size_t i=0; i<3; i++)
      for (size_t j=0; j<3; j++)
        {
          s >>  A[i][j];
        }
    return s;
  }

  // *********** matrix operators *******************

  // comparison functions
  inline bool operator==(const Matrix3D &lhs, const Matrix3D &rhs)
  {
    return (memcmp(lhs.v_, rhs.v_, 9*sizeof(double)) == 0);
  }
  inline bool operator!=(const Matrix3D &lhs, const Matrix3D &rhs)
  {
    return (memcmp(lhs.v_, rhs.v_, 9*sizeof(double)) != 0);
  }

  inline Matrix3D operator+(const Matrix3D &A, 
                            const Matrix3D &B)
  {
    Matrix3D tmp;
    size_t i;
    for (i=0; i<3; ++i) {
      tmp[i][0] = A[i][0] + B[i][0];
      tmp[i][1] = A[i][1] + B[i][1];
      tmp[i][2] = A[i][2] + B[i][2];
    }
    return tmp;
  }

  inline Matrix3D operator-(const Matrix3D &A, 
                            const Matrix3D &B)
  {
    Matrix3D tmp;
    size_t i;
    for (i=0; i<3; ++i) {
      tmp[i][0] = A[i][0] - B[i][0];
      tmp[i][1] = A[i][1] - B[i][1];
      tmp[i][2] = A[i][2] - B[i][2];
    }
    return tmp;
  }

  inline Matrix3D mult_element(const Matrix3D &A, 
                               const Matrix3D &B)
  {
    Matrix3D tmp;
    size_t i;
    for (i=0; i<3; ++i) {
      tmp[i][0] = A[i][0] * B[i][0];
      tmp[i][1] = A[i][1] * B[i][1];
      tmp[i][2] = A[i][2] * B[i][2];
    }
    return tmp;
  }

  inline Matrix3D transpose(const Matrix3D &A)
  {
    Matrix3D S;
    size_t i;
    for (i=0; i<3; ++i) {
        S[size_t(0)][i] = A[i][0];
        S[size_t(1)][i] = A[i][1];
        S[size_t(2)][i] = A[i][2];
    }
    return S;
  }

  inline Matrix3D& Matrix3D::operator+=(const Matrix3D &rhs)
  {
      v_[0] += rhs.v_[0]; v_[1] += rhs.v_[1]; v_[2] += rhs.v_[2];
      v_[3] += rhs.v_[3]; v_[4] += rhs.v_[4]; v_[5] += rhs.v_[5];
      v_[6] += rhs.v_[6]; v_[7] += rhs.v_[7]; v_[8] += rhs.v_[8];

      return *this;
  }

  inline Matrix3D Matrix3D::plus_transpose(const Matrix3D &B)
  {
    Matrix3D tmp;
    size_t i;
    for (i=0; i<3; ++i) {
      tmp[i][0] = v_[3*i+0] + B[0][i];
      tmp[i][1] = v_[3*i+1] + B[1][i];
      tmp[i][2] = v_[3*i+2] + B[2][i];
    }
    return tmp;
  }

  inline Matrix3D operator*(const Matrix3D  &A, 
                            const Matrix3D &B)
  {
    Matrix3D tmp;
    double sum;
    for (size_t i=0; i<3; ++i)
      for (size_t k=0; k<3; ++k)
        {
          sum = 0;
          for (size_t j=0; j<3; j++)
            sum = sum +  A[i][j] * B[j][k];
          tmp[i][k] = sum; 
        }
    return tmp;
  }

  inline int matmult(Matrix3D& C, const Matrix3D  &A, 
                     const Matrix3D &B)
  {
    double sum;
    const double* row_i;
    const double* col_k;
    for (size_t i=0; i<3; ++i)
      for (size_t k=0; k<3; ++k)
        {
          row_i  = &(A[i][0]);
          col_k  = &(B[0][k]);
          sum = 0;
          for (size_t j=0; j<3; ++j)
            {
              sum  += *row_i * *col_k;
              row_i++;
              col_k += 3;
            }
          C[i][k] = sum; 
        }
    return 0;
  }

  /*! \brief Computes \f$ A v \f$ . */
  inline Vector3D operator*(const Matrix3D  &A, const Vector3D &x)
  {
    Vector3D tmp;
    double sum;
    for (size_t i=0; i<3; ++i)
      {
        sum = 0;
        const double* rowi = A[i];
        for (size_t j=0; j<3; ++j)
          sum = sum +  rowi[j] * x[j];
        tmp[i] = sum; 
      }
    return tmp;
  }

  /*! \brief Computes \f$ v^T A \f$ .
    
      This function implicitly considers the transpose of vector x times
      the matrix A and it is implicit that the returned vector must be
      transposed. */
  inline Vector3D operator*(const Vector3D &x, const Matrix3D  &A)
  {
    Vector3D res(0., 0., 0.);
    for (size_t i=0; i<3; ++i)
      {
        const double* rowi = A[i];
        for (size_t j=0; j<3; ++j)
          res[j] += rowi[j] * x[i];
      }
    return res;
  }

} // namespace Mesquite

#endif // Matrix3D_hpp
