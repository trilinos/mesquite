#ifndef MESQUITE_VECTOR3D_HPP
#define MESQUITE_VECTOR3D_HPP

#include <iostream>

#include "Mesquite.hpp"
#ifdef USE_C_PREFIX_INCLUDES
#include <cassert>
#else
#include <assert.h>
#endif

namespace Mesquite
{
  class Vector3D
  {
  public:
    // Constructors
    Vector3D();
    Vector3D(const double& x, const double& y, const double& z);
    Vector3D(const double xyz[3]);
    Vector3D(const Vector3D& to_copy);
    
    // Functions to get the coordinates
    double x() const;
    double y() const;
    double z() const;
    void get_coordinates(double& x, double& y, double& z) const;
    void get_coordinates(double xyz[3]) const;
    const double& operator[](size_t index) const; // 0-based
    
    // Functions to set the coordinates.
    void x(const double x);
    void y(const double y);
    void z(const double z);
    void set(const double x, const double y, const double z);
    void set(const double xyz[3]);
    void set(const Vector3D& to_copy);
    // Subscripts on non-consts both get and set coords
    double& operator[](size_t index); // 0-based
    Vector3D& operator=(const Vector3D &to_copy);
    
    // Functions to modify existing coordinates
    Vector3D operator-() const;  //- unary negation.
    Vector3D& operator*=(const double scalar);
    Vector3D& operator/=(const double scalar);
    Vector3D& operator*=(const Vector3D &rhs); //- cross product
    Vector3D& operator+=(const Vector3D &rhs);
    Vector3D& operator-=(const Vector3D &rhs);
    
    // Binary operators (like a+b).
    friend Vector3D operator+(const Vector3D &lhs,
                              const Vector3D &rhs);
    friend Vector3D operator-(const Vector3D &lhs,
                              const Vector3D &rhs);
    friend Vector3D operator*(const Vector3D &lhs,
                              const double scalar); //- lhs * scalar
    friend Vector3D operator*(const double scalar,
                              const Vector3D &rhs); //- scalar * rhs
    friend Vector3D operator/(const Vector3D &lhs,
                              const double scalar); //- lhs / scalar
    friend double operator%(const Vector3D &v1,
                            const Vector3D &v2); //- dot product
    friend Vector3D operator*(const Vector3D &v1, 
                              const Vector3D &v2); //- cross product

    // Comparison functions
    friend int operator==(const Vector3D &lhs, const Vector3D &rhs);
    friend int operator!=(const Vector3D &lhs, const Vector3D &rhs);
    static double distance_between(const Vector3D& p1,
                                   const Vector3D& p2);
    int within_tolerance_box(const Vector3D &compare_to,
                             double tolerance) const;
    //- Compare two Vector3Ds to see if they are spatially equal.  
    // Return TRUE if difference in x, y, and z are all within tolerance.
    // Essentially checks to see if 'this' lies within a box centered on
    // 'compare_to' with sides of length ('tolerance' * 2).

    // Length functions
    inline double length_squared() const;
    inline double length() const;
    inline void set_length(const double new_length);
    inline void normalize();
    
    // Utility functions.  All angle functions work in radians.
    static double interior_angle(const Vector3D &a,
                                 const Vector3D &b);
    //- Interior angle: acos((a%b)/(|a||b|))
    static Vector3D interpolate(const double param, const Vector3D &p1,
                                const Vector3D &p2);
    //- Interpolate between two points. Returns (1-param)*v1 + param*v2.

  protected:
    double mCoords[3];
  };
  
  // Constructors
  inline Vector3D::Vector3D() 
  {
    mCoords[0] = 0;
    mCoords[1] = 0;
    mCoords[2] = 0;
  }
  inline Vector3D::Vector3D(const double &x,
                            const double &y,
                            const double &z) 
  {
    mCoords[0] = x;
    mCoords[1] = y;
    mCoords[2] = z;
  }
  inline Vector3D::Vector3D(const double xyz[3]) 
  { memcpy(mCoords, xyz, 3*sizeof(double)); }
  inline Vector3D::Vector3D(const Vector3D& to_copy) 
  { memcpy(mCoords, to_copy.mCoords, 3*sizeof(double)); }
  
  // Functions to get coordinates
  inline double Vector3D::x() const
  { return mCoords[0]; }
  inline double Vector3D::y() const
  { return mCoords[1]; }
  inline double Vector3D::z() const
  { return mCoords[2]; }
  inline void Vector3D::get_coordinates(double &x, double &y, double &z) const
  {
    x = mCoords[0]; 
    y = mCoords[1]; 
    z = mCoords[2];
  }
  inline void Vector3D::get_coordinates(double xyz[3]) const
  { memcpy(xyz, mCoords, 3*sizeof(double)); }
  inline const double& Vector3D::operator[](size_t index) const
  {
#ifdef MSQ_DEBUG
    assert(0<=(int)index);
    assert(index<=2);     
#endif
    return mCoords[index]; }
  
  // Functions to set coordinates
  inline void Vector3D::x( const double x )
  { mCoords[0] = x; }
  inline void Vector3D::y( const double y )
  { mCoords[1] = y; }
  inline void Vector3D::z( const double z )
  { mCoords[2] = z; }
  inline void Vector3D::set(const double x,
                            const double y,
                            const double z)
  {
    mCoords[0] = x;
    mCoords[1] = y;
    mCoords[2] = z;
  }
  inline void Vector3D::set(const double xyz[3])
  { memcpy(mCoords, xyz, 3*sizeof(double)); }
  inline void Vector3D::set(const Vector3D& to_copy)
  { memcpy(mCoords, to_copy.mCoords, 3*sizeof(double)); }
  inline double& Vector3D::operator[](size_t index)
  {
#ifdef MSQ_DEBUG
    assert(0<=(int)index);
    assert(index<=2);     
#endif
    return mCoords[index]; }
  inline Vector3D& Vector3D::operator=(const Vector3D &to_copy)  
  {
    memcpy(mCoords, to_copy.mCoords, 3*sizeof(double));
    return *this;
  }
  
  // Functions that modify existing coordinates
  inline Vector3D Vector3D::operator-() const
  {
    return Vector3D(-mCoords[0], -mCoords[1], -mCoords[2]);
  }
  inline Vector3D& Vector3D::operator*=(const double scalar)
  {
    mCoords[0] *= scalar;
    mCoords[1] *= scalar;
    mCoords[2] *= scalar;
    return *this;
  }
  inline Vector3D& Vector3D::operator/=(const double scalar)
  {
    assert (scalar != 0);
    mCoords[0] /= scalar;
    mCoords[1] /= scalar;
    mCoords[2] /= scalar;
    return *this;
  }
  inline Vector3D& Vector3D::operator*=(const Vector3D &rhs)
  {
    double new_coords[3] = 
      {mCoords[1]*rhs.mCoords[2] - mCoords[2]*rhs.mCoords[1],
       mCoords[2]*rhs.mCoords[0] - mCoords[0]*rhs.mCoords[2],
       mCoords[0]*rhs.mCoords[1] - mCoords[1]*rhs.mCoords[0]
      };
    memcpy(mCoords, new_coords, 3*sizeof(double));
    return *this;
  }
  inline Vector3D& Vector3D::operator+=(const Vector3D &rhs)
  {
    mCoords[0] += rhs.mCoords[0];
    mCoords[1] += rhs.mCoords[1];
    mCoords[2] += rhs.mCoords[2];
    return *this;
  }
  inline Vector3D& Vector3D::operator-=(const Vector3D &rhs)
  {
    mCoords[0] -= rhs.mCoords[0];
    mCoords[1] -= rhs.mCoords[1];
    mCoords[2] -= rhs.mCoords[2];
    return *this;
  }

  // Binary operators
  inline Vector3D operator+(const Vector3D &lhs,
                            const Vector3D &rhs)
  {
    return Vector3D(lhs.x() + rhs.x(),
                    lhs.y() + rhs.y(),
                    lhs.z() + rhs.z());
  }
  inline Vector3D operator-(const Vector3D &lhs,
                            const Vector3D &rhs)
  {
    return Vector3D(lhs.x() - rhs.x(),
                    lhs.y() - rhs.y(),
                    lhs.z() - rhs.z());
  }
  inline Vector3D operator*(const Vector3D &lhs,
                            const double scalar)
  {
    return Vector3D(lhs.x() * scalar,
                    lhs.y() * scalar,
                    lhs.z() * scalar);
  }
  inline Vector3D operator*(const double scalar,
                            const Vector3D &rhs)
  {
    return Vector3D(rhs.x() * scalar,
                    rhs.y() * scalar,
                    rhs.z() * scalar);
  }
  inline Vector3D operator/(const Vector3D &lhs,
                            const double scalar)
  {
    assert (scalar != 0);
    return Vector3D(lhs.x() / scalar,
                    lhs.y() / scalar,
                    lhs.z() / scalar);
  }
  inline double operator%(const Vector3D &lhs,
                          const Vector3D &rhs) // Dot Product
  {
    return( lhs.mCoords[0] * rhs.mCoords[0] +
            lhs.mCoords[1] * rhs.mCoords[1] +
            lhs.mCoords[2] * rhs.mCoords[2] );
  }
  inline Vector3D operator*(const Vector3D &lhs, 
                            const Vector3D &rhs) // Cross Product
  {
    return Vector3D(lhs.mCoords[1]*rhs.mCoords[2]-lhs.mCoords[2]*rhs.mCoords[1],
                    lhs.mCoords[2]*rhs.mCoords[0]-lhs.mCoords[0]*rhs.mCoords[2],
                    lhs.mCoords[0]*rhs.mCoords[1]-lhs.mCoords[1]*rhs.mCoords[0]);
  }
  
  // output operator
  inline std::ostream& operator<<(std::ostream &s,
                                  const Mesquite::Vector3D &v)
  {
    return s << v[0] << ' ' << v[1] << ' ' << v[2] << ' ' << std::endl;
  }
  
  // comparison functions
  inline int operator==(const Vector3D &lhs, const Vector3D &rhs)
  {
    return (memcmp(lhs.mCoords, rhs.mCoords, 3*sizeof(double)) == 0);
  }
  inline int operator!=(const Vector3D &lhs, const Vector3D &rhs)
  {
    return (memcmp(lhs.mCoords, rhs.mCoords, 3*sizeof(double)) != 0);
  }
  inline double Vector3D::distance_between(const Vector3D &p1,
                                           const Vector3D &p2)
  {
    Vector3D v = p2 - p1;
    return v.length();
  }
  inline int Vector3D::within_tolerance_box(const Vector3D &compare_to,
                                            double tolerance) const
  {
    return ((fabs(this->mCoords[0] - compare_to.mCoords[0]) < tolerance) &&
            (fabs(this->mCoords[1] - compare_to.mCoords[1]) < tolerance) &&
            (fabs(this->mCoords[2] - compare_to.mCoords[2]) < tolerance));
  }
  
  // Length functions
  inline double Vector3D::length_squared() const
  {
    return (mCoords[0]*mCoords[0] +
            mCoords[1]*mCoords[1] +
            mCoords[2]*mCoords[2]);
  }
  inline double Vector3D::length() const
  {
    return sqrt(mCoords[0]*mCoords[0] +
                mCoords[1]*mCoords[1] +
                mCoords[2]*mCoords[2]);
  }
  inline void Vector3D::set_length(const double new_length)
  {
    // Make sure at least one coordinate is non-zero so
    // we've got a non-zero length.
    assert(mCoords[0] != 0.0 ||
           mCoords[1] != 0.0 ||
           mCoords[2] != 0.0);
    
    double factor = new_length / length();
    *this *= factor;
  }
  inline void Vector3D::normalize()
  { set_length(1.0); }

  // Utility functions.
  inline double Vector3D::interior_angle(const Vector3D &lhs,
                                         const Vector3D &rhs)
  {
    double len1 = lhs.length();
    double len2 = rhs.length();
    
    assert(len1 > 0);
    assert(len2 > 0);
    
    double angle_cos = (lhs % rhs)/(len1 * len2);
    
    // Adjust the cosine if slightly out of range
    if ((angle_cos > 1.0) && (angle_cos < 1.0001))
      {
        angle_cos = 1.0;
      }
    else if (angle_cos < -1.0 && angle_cos > -1.0001)
      {
        angle_cos = -1.0;
      }
    assert(angle_cos >= -1.0 && angle_cos <= 1.0);
    
    return acos(angle_cos);
  }
  inline Vector3D Vector3D::interpolate(const double param,
                                        const Vector3D &p1,
                                        const Vector3D &p2)
  {
    return (1-param)*p1 + param*p2;
  }

} // namespace Mesquite

#endif
