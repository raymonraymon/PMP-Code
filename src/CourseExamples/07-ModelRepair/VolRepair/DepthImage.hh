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
#ifndef DEPTH_IMAGE_HH
#define DEPTH_IMAGE_HH
//=============================================================================

#include <vector>
#include <limits>

#include <OpenMesh/Core/Math/VectorT.hh>

#include "Tree.hh"

//=============================================================================


class DepthImage
{
public:

  typedef OpenMesh::Vec3d Vec3d;

  // Initialize a depth image with a scene
  DepthImage( const Vec3d & _normal,
	      const Vec3d & _scene_center,
	      double        _scene_radius,
	      int _width, int _height );

  // Take a depth image of a scene
  void record( const Tree & _tree );


  enum PointClass {
    INTERIOR,
    EXTERIOR,
    UNKNOWN,
  };

  // Classify a point
  PointClass classify( const Vec3d & _point ) const;


private:

  // Return the front depth sample of a pixel

  double & front( int _i0, int _i1 )
  { return front_[ width_ * _i1 + _i0 ]; }

  double front( int _i0, int _i1 ) const
  { return front_[ width_ * _i1 + _i0 ]; }


  // Return the back depth sample of a pixel

  double & back( int _i0, int _i1 )
  { return back_[ width_ * _i1 + _i0 ]; }

  double back( int _i0, int _i1 ) const
  { return back_[ width_ * _i1 + _i0 ]; }


  // Center and radius of the bounding sphere of the scene

  Vec3d  scene_center_;
  double scene_radius_;

  // Axes of the image coordinate system:
  // - the origin of the system is the scene center
  // - the image plane is spanned by axis 0 and 1 (x- and y-coordinate)
  // - axis 2 gives the projection direction

  Vec3d  axis_[3];

  // Dimensions of the image (in pixels)

  int width_;
  int height_;

  // Image data for the first and last depth sample along a ray.

  std::vector< double > front_;
  std::vector< double > back_;

};


//=============================================================================
#endif
//=============================================================================
