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
#include "Tree.hh"
//=============================================================================
template<typename T>
T Tree_min(T& a, T& b)
{
	if (a<b)
	{
		return a;
	} 
	else
	{
		return b;
	}
	
}

Tree::Tree()
{
  root_ = 0;
  
  bbox_.first = Vec3d( std::numeric_limits<double>::max(),
		       std::numeric_limits<double>::max(),
		       std::numeric_limits<double>::max() );
  bbox_.second = - bbox_.first;
  
  axis_[0] = Vec3d( 1, 0, 0 );
  axis_[1] = Vec3d( 0, 1, 0 );
  axis_[2] = Vec3d( 0, 0, 1 );
}


//=============================================================================


Tree::~Tree()
{
  delete_node( root_ );
}



//=============================================================================
//
// Delete a subtree
//
//=============================================================================

void
Tree::delete_node( Node * _node )
{
  if ( _node )
  {
    if ( _node->front )
      delete_node( _node->front );
    if ( _node->incident )
	delete_node( _node->incident );
    if ( _node->back )
      delete_node( _node->back );
    delete _node;
  }
}




//=============================================================================
//
// Add a triangle to the tree
//
// Note: After adding all triangles one has to call build_kdtree() to
// actually build the tree.
//
//=============================================================================

void
Tree::push_back( const Triangle & _triangle )
{
  triangle_.push_back( _triangle );
  
  bbox_.first.minimize( _triangle.point( 0 ) );
  bbox_.first.minimize( _triangle.point( 1 ) );
  bbox_.first.minimize( _triangle.point( 2 ) );
  
  bbox_.second.maximize( _triangle.point( 0 ) );
  bbox_.second.maximize( _triangle.point( 1 ) );
  bbox_.second.maximize( _triangle.point( 2 ) );
}


//=============================================================================
//
// Build up a kd-tree
//
//=============================================================================


void
Tree::build_kdtree()
{
  root_ = new Node();
  build_kdtree( root_, triangle_.begin(), triangle_.end() );
}




//=============================================================================
//
// Check for intersection with a line segment
//
// This function checks wether the line segment (_p0,_p1) intersects a
// triangle in the tree. If so it returns the barycentric coordinate t
// of the intersection point q that is nearest to _p0, i.e. q =
// (1-t)*_p0+t*_p1. If there is no intersection this function returns
// infinity.
//
//=============================================================================

double
Tree::intersect( const Vec3d & _p0,
		 const Vec3d & _p1, int _level ) const
{
  // The origin and direction of the ray are globally stored. A line
  // segment (a,b) is then given as two scalar parameters ta, tb such that
  //   a = origin_ + ta * direction_ and
  //   b = origin_ + tb * direction_
  // All further processing is done with ta and tb instead of a and b.

  origin_    = _p0;
  direction_ = _p1 - _p0;
  
  return intersect( root_, 0, 1 );
}




//=============================================================================
//
// Main functionality for segment/tree intersection
//
//=============================================================================

double
Tree::intersect( Node * _node,
		 double _t0, double _t1, int _level ) const
{

#ifdef BRUTE_FORCE
  double mindist = std::numeric_limits<double>::infinity();
  int idx = -1;
  for ( CTriVecIter i = triangle_.begin(); i != triangle_.end(); ++i )
  {
    double dist = i->intersect( origin_, direction_ );
    if ( dist < mindist )
    {
      mindist = dist;
      idx = i - triangle_.begin();
    }
  }
    
  return mindist;
#endif


  // If _node is a leaf, we have to test all triangles in this leaf
  // for intersections with the line segment.

  if ( ! _node->front && 
       ! _node->incident && 
       ! _node->back )
  {
    double mindist = std::numeric_limits<double>::infinity();

    for ( TriVecIter i = _node->begin; i != _node->end; ++i )
    {
      double dist = i->intersect( origin_, direction_ );
      if ( dist < mindist )
	mindist = dist;
    }
    
    return mindist;
  }
  

  // If _node is not a leaf, we classify the endpoints of the line
  // segment with respect to the plane that is tsored in _node. If
  // both endpoints lie in the same halfspace, we do a simple
  // recursion. If not, we visit the corresponding spaces in the
  // correct order.

  Plane::Classification c0 =
    _node->plane.classify( origin_ + _t0 * direction_ );
  Plane::Classification c1 =
    _node->plane.classify( origin_ + _t1 * direction_ );
  

  if ( c0 != Plane::BACK && c1 != Plane::BACK )
  {
    double dist1 = std::numeric_limits< double >::infinity();
    double dist2 = std::numeric_limits< double >::infinity();
    
    if ( _node->front )
      dist1 = intersect( _node->front, _t0, _t1, _level+1 );
    
    if ( _node->incident )
      dist2 = intersect( _node->incident, _t0, _t1, _level+1 );
    
    return Tree_min( dist1, dist2 );
  }
  
  if ( c0 != Plane::FRONT && c1 != Plane::FRONT )
  {
    double dist1 = std::numeric_limits< double >::infinity();
    double dist2 = std::numeric_limits< double >::infinity();
    
    if ( _node->back )
      dist1 = intersect( _node->back, _t0, _t1, _level+1 );
    
    if ( _node->incident )
      dist2 = intersect( _node->incident, _t0, _t1, _level+1 );
    
    return Tree_min( dist1, dist2 );
  }
  
  
  double t = _node->plane.intersect( origin_, direction_ );
  
  
  if ( c0 == Plane::FRONT )
  {
    double dist1 = std::numeric_limits< double >::infinity();
    double dist2 = std::numeric_limits< double >::infinity();
    double dist3 = std::numeric_limits< double >::infinity();
    
    if ( _node->front )
      dist1 = intersect( _node->front, _t0, t, _level+1);

    if ( _node->incident )
      dist2 = intersect( _node->incident, _t0, _t1, _level+1 );
    
    if ( dist1 != std::numeric_limits< double >::infinity() )
      return Tree_min( dist1, dist2 );
    
    if ( _node->back )
      dist3 = intersect( _node->back, t, _t1, _level+1 );

    return Tree_min( dist2, dist3 );
  }
  
  
  if ( c0 == Plane::BACK )
  {
    double dist1 = std::numeric_limits< double >::infinity();
    double dist2 = std::numeric_limits< double >::infinity();
    double dist3 = std::numeric_limits< double >::infinity();
    
    if ( _node->back )
      dist1 = intersect( _node->back, _t0, t, _level+1 );

    if ( _node->incident )
      dist2 = intersect( _node->incident, _t0, _t1, _level+1 );
    
    if ( dist1 != std::numeric_limits< double >::infinity() )
      return Tree_min( dist1, dist2 );
    
    if ( _node->front )
      dist3 = intersect( _node->front, t, _t1, _level+1 );
    
    return Tree_min( dist2, dist3 );
  }
  
  return std::numeric_limits< double >::infinity();
}




//=============================================================================
//
// Build up the tree.
//
// We use a very simple heuristic to build up a kd-tree, namely: split
// along the longest edge of the bounding box of the triangles. More
// sophisticated heuristics are given in e.g.
//
//=============================================================================

void
Tree::build_kdtree( Node     * _node,
		    TriVecIter   _begin,
		    TriVecIter   _end )
{
  // If there are less than a given threshold triangles, we do not
  // bother to split the node.
  
  if ( _end - _begin < 10 )
  {
    _node->begin = _begin;
    _node->end   = _end;
    return;
  }
  
  // Choose splitting plane: For simplicity, we compute the bounding
  // box of the triangles and split it along the longest axis.
  
  BoundingBox bb = bounding_box( _begin, _end );
  
  int dir = maxarg( bb.second - bb.first );
  
  _node->plane = Plane( axis_[ dir ],
			0.5 * ( bb.first[dir] + bb.second[dir] ) );
  
  // Partition the triangles according to that plane and distribute
  // them to the child nodes.
  
  TriVecIterPair p = partition( _node->plane, _begin, _end );
  
  // In case that all triangles have been distributed into different
  // partitions, we can stop.
  
  if ( p.first == _begin && p.second == _end )
  {
    _node->begin = _begin;
    _node->end   = _end;
    return;
  }
  
  // FRONT triangles
  
  if ( _begin != p.first )
  {
    _node->front = new Node();
    build_kdtree( _node->front, _begin, p.first );
  }
  
  
  // INCIDENT triangles
  
  if ( p.first != p.second )
  {
    _node->incident = new Node();
    build_kdtree( _node->incident, p.first, p.second );
  }
  
  
  // BACK triangles
  if ( p.second != _end )
  {
    _node->back = new Node();
    build_kdtree( _node->back, p.second, _end );
  }
  
}


  
//=============================================================================
//
// Partition triangles according to their position w.r.t. a plane
//
// Partitions all triangles in the range [_begin,_end) according
// to _plane and returns a pair (i1,i2) of iterators such that all
// triangles in
//  [_begin, i1) are FRONT triangles
//  [i1,i2)      are INCIDENT triangles
//  [i2,end)     are BACK triangles
//
//=============================================================================

Tree::TriVecIterPair
Tree::partition( const Plane & _plane,
		 TriVecIter    _begin,
		 TriVecIter    _end )
{
  
  TriVecIter it1( _begin );
  
  for ( TriVecIter i = it1; i != _end; ++i )
    if ( _plane.classify( *i ) == Plane::FRONT )
    {
      std::swap( *it1, *i );
      ++it1;
    }
  
  TriVecIter it2( it1 );
  for ( TriVecIter i = it2; i != _end; ++i )
    if ( _plane.classify( *i ) == Plane::INCIDENT )
    {
      std::swap( *it2, *i );
      ++it2;
    }
  
  return TriVecIterPair( it1, it2 );
}



//=============================================================================
//
// Compute the bounding box for a set of triangles
//
//=============================================================================

Tree::BoundingBox
Tree::bounding_box( TriVecIter _begin,
		    TriVecIter _end )
{
  Vec3d bbmin( std::numeric_limits<double>::max(),
	       std::numeric_limits<double>::max(),
	       std::numeric_limits<double>::max() );
  Vec3d bbmax = - bbmin;
  
  for ( TriVecIter i = _begin; i != _end; ++i )
  {
    bbmin.minimize( i->point( 0 ) );
    bbmin.minimize( i->point( 1 ) );
    bbmin.minimize( i->point( 2 ) );
    
    bbmax.maximize( i->point( 0 ) );
    bbmax.maximize( i->point( 1 ) );
    bbmax.maximize( i->point( 2 ) );
  }
  
  return BoundingBox( bbmin, bbmax );
}



//=============================================================================

