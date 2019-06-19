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
#ifndef TREE_HH
#define TREE_HH
//=============================================================================

#include <vector>
#include <limits>

#include <OpenMesh/Core/Math/VectorT.hh>

#include "Geometry.hh"

//=============================================================================


class Tree
{
public:

  typedef OpenMesh::Vec3d Vec3d;

  typedef std::vector< Triangle >             TriVec;
  typedef std::vector< Triangle >::iterator   TriVecIter;
  typedef std::vector< Triangle >::const_iterator CTriVecIter;
  typedef std::pair< TriVecIter, TriVecIter > TriVecIterPair;
  typedef std::pair< Vec3d, Vec3d >           BoundingBox;

  struct Node {

    Node() :
      front( 0 ), incident( 0 ), back( 0 ) {}

    Plane plane;
    
    Node * front;
    Node * incident;
    Node * back;
    
    TriVecIter begin;
    TriVecIter end;
  };
  
  
  Tree();
  ~Tree();

  // Return the bounding box
  Vec3d bbmin() const { return bbox_.first; }
  Vec3d bbmax() const { return bbox_.second; }

  // Add a triangle to the tree
  void push_back( const Triangle & _triangle );

  // Build up a kd-tree
  void build_kdtree();

  // Check for intersection with a line segment
  double intersect( const Vec3d & _p0,
		    const Vec3d & _p1, int _level = 0 ) const;

private:

  // Delete a subtree
  void delete_node( Node * _node );

  // Helper
  int maxarg( const Vec3d & _v ) const
  {
    return ( _v[0]>_v[1] ? ( _v[0]>_v[2] ? 0 : 2 ) :
                           ( _v[1]>_v[2] ? 1 : 2 ) );
  }
  
  
  // The bounding box of the scene
  BoundingBox   bbox_;

  // Root node
  Node        * root_;

  // All triangles that are stored in the tree
  TriVec        triangle_;

  // Helper members for the segment/tree intersection routines
  mutable Vec3d origin_;
  mutable Vec3d direction_;
  
  // The three axes for splitting nodes
  Vec3d axis_[3];


  // Main functionality for segment/tree intersection
  double intersect( Node * _node,
		    double _t0, double _t1, int _level = 0 ) const;

  // Main functionality for building a kd-tree
  void build_kdtree( Node     * _node,
		     TriVecIter   _begin,
		     TriVecIter   _end );
  
  
  // Partition triangles according to their position w.r.t. a plane
  TriVecIterPair partition( const Plane & _plane,
			    TriVecIter    _begin,
			    TriVecIter    _end );

  // Compute the bounding box for a set of triangles
  BoundingBox bounding_box( TriVecIter _begin,
			    TriVecIter _end );

};



//=============================================================================
#endif
//=============================================================================
