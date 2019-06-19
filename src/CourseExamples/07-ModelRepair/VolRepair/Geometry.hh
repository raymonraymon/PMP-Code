//=============================================================================
//                                                                            
//   Example code for the full-day course
//
//   M. Botsch, M. Pauly, L. Kobbelt, P. Alliez, B. Levy,
//   "Geometric Modeling Based on Polygonal Meshes"
//   held at SIGGRAPH 2007, San Diego, and Eurographics 2008, Crete.
//
//   Copyright (C) 2007 by  Computer Graphics Laboratory, ETH Zurich, 
//                      and Computer Graphics Group,      RWTH Aachen
//
//                                                                            
//-----------------------------------------------------------------------------
//                                                                            
//                                License                                     
//                                                                            
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License
//   as published by the Free Software Foundation; either version 2
//   of the License, or (at your option) any later version.
//   
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//   
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 51 Franklin Street, Fifth Floor, 
//   Boston, MA  02110-1301, USA.
//                                                                            
//=============================================================================
//=============================================================================
#ifndef GEOMETRY_HH
#define GEOMETRY_HH
//=============================================================================

#include <OpenMesh/Core/Math/VectorT.hh>

//=============================================================================

//=============================================================================
//
// This file provides the self-explaining classes
//
//  Triangle : representation of a triangle
//  Plane    : representation of a plane
//
//=============================================================================

class Triangle
{
  
public:
  
  typedef OpenMesh::Vec3d Vec3d;
  
  Triangle() {}
  
  Triangle( const Vec3d & _point0,
	    const Vec3d & _point1,
	    const Vec3d & _point2 )
  {
    point_[0] = _point0;
    point_[1] = _point1;
    point_[2] = _point2;
  }

  const Vec3d & point( int i ) const
  {
    return point_[ i ];
  }

  double intersect( const Vec3d & _origin,
		    const Vec3d & _direction ) const
  {
    //
    // Note: This code was adapted from
    //
    //    Tomas Möller, Ben Trumbore
    //    Fast, Minimum Storage Ray-Triangle Intersection
    //    Journal of Graphics Tools, 2(1), 21-28, 1997
    //
   
    const double eps = 0.000001;

    Vec3d edge1( point(1) - point(0) );
    Vec3d edge2( point(2) - point(0) );

    Vec3d pvec = cross( _direction, edge2 );
    double det =  dot( edge1, pvec );

    if ( det > -eps && det < eps )
      return std::numeric_limits<double>::infinity();
    
    double inv_det = 1.0f / det;
    
    Vec3d tvec( _origin - point( 0 ) );

    double u = dot( tvec, pvec ) * inv_det;
    if ( u < 0.0f || u > 1.0f )
      return std::numeric_limits<double>::infinity();

    Vec3d qvec = cross( tvec, edge1 );

    double v = dot( _direction, qvec ) * inv_det;
    if ( v < 0.0f || u + v > 1.0f )
      return std::numeric_limits<double>::infinity();

    return dot( edge2, qvec ) * inv_det;
  }


private:

  Vec3d  point_[ 3 ];


};


//=============================================================================


class Plane
{
public:

  typedef OpenMesh::Vec3d Vec3d;
  
  enum Classification {
    FRONT    = 0,
    INCIDENT = 1,
    BACK     = 2,
  };

  Plane() : normal_( Vec3d( 0, 0, 0 ) ), distance_( 0 ) {}

  Plane( const Triangle & _triangle )
  {
    Vec3d e0( _triangle.point(1) - _triangle.point(0) );
    Vec3d e1( _triangle.point(2) - _triangle.point(0) );

    normal_ = cross( e0, e1 );

    double l = normal_.norm();

    if ( fabs(l) < 0.0001 )
    {
      normal_ = Vec3d( 0, 0, 0 );
      distance_ = 0;
    }
    else
    {
      normal_ /= l;
      distance_ = ( _triangle.point( 0 ) | normal_ );
    }
  }
  
  Plane( const Vec3d & _normal,
	 double         _distance ) :
    normal_( _normal ), distance_( _distance )
  {}
  
  ~Plane() {}


  const Vec3d & normal() const { return normal_; }

  double distance() const { return distance_; }

  double intersect( const Vec3d & _origin,
		    const Vec3d & _direction ) const
  {
    return ( ( distance_ - dot( _origin, normal_ ) ) /
	     dot( _direction, normal_ ) );
  }


  Classification classify( const Vec3d & _point ) const
  {
    double v = dot( _point, normal_ ) - distance_;
    if ( v > 0 )
      return FRONT;
    if ( v < 0 )
      return BACK;
    return INCIDENT;
  }
  
  Classification classify( const Triangle & _triangle ) const
  {
    Classification c0 = classify( _triangle.point( 0 ) );
    Classification c1 = classify( _triangle.point( 1 ) );
    Classification c2 = classify( _triangle.point( 2 ) );
    
    if ( c0 == FRONT && c1 == FRONT && c2 == FRONT )
      return FRONT;

    if ( c0 == BACK && c1 == BACK && c2 == BACK )
      return BACK;

    return INCIDENT;
  }
 
private:

  Vec3d  normal_;
  double distance_;

};


//=============================================================================
#endif
//=============================================================================
