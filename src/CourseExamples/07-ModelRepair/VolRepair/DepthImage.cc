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
#include "DepthImage.hh"
#include <math.h>
//=============================================================================


//=============================================================================
//
// Initialize a depth image with a scene
//
// The depth image plane is orthogonal to the given normal and
// intersects the scene center. Axis 0 and 1 span the image plane,
// axis 2 is the normal. The image itself is scaled and translated
// such that it encompasses the scene radius.
//
//=============================================================================

DepthImage::DepthImage( const Vec3d & _normal,
			const Vec3d & _scene_center,
			double        _scene_radius,
			int _width, int _height )
{
  axis_[0] = cross( _normal, Vec3d(_normal[1]+1, _normal[2]+2, _normal[0]+3 ) );
  axis_[0].normalize();
  
  axis_[1] = cross( _normal, axis_[0] );
  axis_[1].normalize();
  
  axis_[2] = cross( axis_[0], axis_[1] );
  axis_[2] = _normal;
  axis_[2].normalize();
  
  width_  =  _width;
  height_ = _height;
  
  scene_center_ = _scene_center;
  scene_radius_ = _scene_radius;
  
  front_.resize( width_ * height_ );
  back_ .resize( width_ * height_ );
  
  for ( int y = 0; y < height_; ++y )
    for ( int x = 0; x < width_; ++x )
    {
      front( x, y ) = - 2 * _scene_radius;
      back( x, y )  =   2 * _scene_radius;
    }

}



//=============================================================================
//
// Take a depth image of a scene
//
// Let p be the center of a pixel of the image. We shoot a ray
// orthogonal to the image plane through p and record the first and
// last intersection with the scene.
//
//=============================================================================

void
DepthImage::record( const Tree & _tree )
{
  for ( int y = 0; y < height_; ++y )
    for ( int x = 0; x < width_; ++x )
    {
      // Normalized coordinates of pixel (x,y)
      
      double nx = ( x + 0.5 ) / width_;
      double ny = ( y + 0.5 ) / height_;
      
      
      // Local coordinates in image plane
      
      double lx = ( 2.0 * nx - 1.0 ) * scene_radius_;
      double ly = ( 2.0 * ny - 1.0 ) * scene_radius_;
      
      
      // Instead of a ray, we use a line segment which is orthogonal
      // to the image plane and intersects the whole scene.
      
      Vec3d segmentA = scene_center_
		     + lx * axis_[0]
		     + ly * axis_[1]
		     + 2 * scene_radius_ * axis_[2];
      Vec3d segmentB = scene_center_
		     + lx * axis_[0]
		     + ly * axis_[1]
		     - 2 * scene_radius_ * axis_[2];
      
      // Test for intersection of the segment with the scene ...

      double ft = _tree.intersect( segmentA, segmentB, true );
      double bt = _tree.intersect( segmentB, segmentA, true );
      
      // ... and if yes ...
      
      if ( ft != std::numeric_limits< double >::infinity() &&
	   bt != std::numeric_limits< double >::infinity() )
      {
	// ... compute the intersecttion points ...
	
	Vec3d front_is = ( 1 - ft ) * segmentA + ft * segmentB;
	Vec3d back_is  = ( 1 - bt ) * segmentB + bt * segmentA;
	
	// ... and store them in the image.

	front( x, y ) = dot( front_is - scene_center_, axis_[2] );
	back ( x, y ) = dot( back_is  - scene_center_, axis_[2] );
      }
      
    }
}


//=============================================================================
//
// Classify a point
//
// Classify a point as solid, if it lies between the front and back
// depth sample.
//
//=============================================================================

DepthImage::PointClass
DepthImage::classify( const Vec3d & _point ) const
{
  // Coordinates in image plane
  
  Vec3d local( dot( _point - scene_center_, axis_[0] ),
	       dot( _point - scene_center_, axis_[1] ),
	       dot( _point - scene_center_, axis_[2] ) );
  
  
  // Normalized image coordinates
  
  double nx = ( local[0] / scene_radius_ + 1.0 ) / 2.0;
  double ny = ( local[1] / scene_radius_ + 1.0 ) / 2.0;
  
  
  // Float pixel coordinates
  
  double fx = nx * width_;
  double fy = ny * height_;
  
  // Bilinear interpolation
  
  double ax = fx - floor( fx );
  double ay = fy - floor( fy );
  
  int x0 = static_cast<int>( floor( fx ) );
  int x1 = static_cast<int>( ceil( fx ) );
  int y0 = static_cast<int>( floor( fy ) );
  int y1 = static_cast<int>( ceil( fy ) );
  
  double f00 = front(x0,y0);
  double f01 = front(x0,y1);
  double f10 = front(x1,y0);
  double f11 = front(x1,y1);

  if ( f00 < -scene_radius_ && f01 < -scene_radius_ &&
       f10 < -scene_radius_ && f11 < -scene_radius_ )
    return EXTERIOR;

  if ( f00 < -scene_radius_ || f01 < -scene_radius_ ||
       f10 < -scene_radius_ || f11 < -scene_radius_ )
    return UNKNOWN;

  double b00 = back(x0,y0);
  double b01 = back(x0,y1);
  double b10 = back(x1,y0);
  double b11 = back(x1,y1);

  if ( b00 > scene_radius_ || b01 > scene_radius_ ||
       b10 > scene_radius_ || b11 > scene_radius_ )
    return UNKNOWN;

  double f = (1-ax)*(1-ay)*f00
	   + (1-ax)*   ay *f01
	   +    ax *(1-ay)*f10
	   +    ax *   ay *f11;

  double b = (1-ax)*(1-ay)*b00
	   + (1-ax)*   ay *b01
	   +    ax *(1-ay)*b10
	   +    ax *   ay *b11;

  if ( f > local[2] && local[2] > b )
    return INTERIOR;
  
  return EXTERIOR;
}


//=============================================================================
