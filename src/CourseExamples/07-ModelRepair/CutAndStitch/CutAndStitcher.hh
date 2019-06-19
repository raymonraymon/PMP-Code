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
#ifndef CUTANDSTITCH_HH
#define CUTANDSTITCH_HH
//=============================================================================

#include <map>
#include <vector>
#include <float.h>
#include "OpenMeshUtils.hh"
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>

//=============================================================================


template< class TheMesh >
class CutAndStitcher
{
public:

  typedef TheMesh Mesh;

  import_om_abbreviations( typename Mesh );


  //
  // An unordered pair of indices
  //

  class UnorderedPair {
  public:
    
    UnorderedPair( int _v0, int _v1 ) : v0_( _v0 ), v1_( _v1 )
    {
      if ( _v1 < _v0 ) std::swap( v0_, v1_ );
    }

    bool operator<( const UnorderedPair & _rhs ) const {
      return ( v0() < _rhs.v0() ||
	       ( v0() == _rhs.v0() &&
		 v1() <  _rhs.v1() ) );		 
    }

    int v0() const { return v0_; }
    int v1() const { return v1_; }

  private:

    int v0_, v1_;

  };


  //
  // Data associated with each edge:
  // - The number of incident faces
  // - Up to two incident faces (if there are more than two, then this edge
  //   is a complex edge anyway)
  //

  class EdgeData {
  public:
    
    EdgeData() : count_( 0 )
    {}

    EdgeData( FH _fh ) : count_( 1 ), face0_( _fh )
    {}

    void add_face( FH _fh )
    {
      ++count_;
      face1_ = _fh;
    }

    int count() const { return count_; }

    FH face0() const { return face0_; }
    FH face1() const { return face1_; }


  private:

    int count_;
    FH face0_;
    FH face1_;    
  };



  CutAndStitcher( Mesh & _mesh ) : mesh_( _mesh )
  {
    mesh_.request_vertex_status();
    mesh_.request_edge_status();
    mesh_.request_face_status();
    mesh_.request_halfedge_status();
  }

  ~CutAndStitcher()
  {
    mesh_.release_vertex_status();
    mesh_.release_edge_status();
    mesh_.release_face_status();
    mesh_.release_halfedge_status();
  }


  void cut_and_stitch()
  {
    std::cerr << "Cutting and stitching ... ";

    // Retrive all faces of the mesh as a list of points (3 successive
    // points make up a face)

    std::vector< Point > p;

    for ( FI fi = mesh_.faces_begin(); fi != mesh_.faces_end(); ++fi )
      for ( FVI fvi = mesh_.fv_iter( fi ); fvi; ++fvi )
	p.push_back( mesh_.point( fvi ) );



    // Re-insert faces into mesh and identify vertices at the same
    // spatial position
    //
    // For identifying vertices we use a simple map
    //
    //      vertex position -> ID
    //
    // Before we add a new point, we look up its ID in the map. If the
    // vertex already exists we do not need to insert it anymore, but
    // can directly use its index.

    typedef std::map< Point, int > Point2ID;
    Point2ID id_from_pos;

    mesh_.clear();

    mesh_.add_property( vertex_id );

    int next_id = 0;

    std::vector< VH > face( 3 );

    for ( int i = 0; i < (int) p.size(); i += 3 )
    {
      for ( int j = 0; j < 3; ++j )
      {
	typename Point2ID::iterator it;
	it = id_from_pos.find( p[i+j] );
	if ( it == id_from_pos.end() )
	{
	  it = id_from_pos.insert
	     ( typename Point2ID::value_type( p[i+j], next_id ) ).first;
	  ++next_id;
	}
	
	face[j] = mesh_.add_vertex( p[i+j] );
	mesh_.property( vertex_id, face[j] ) = it->second;
      }

      
      mesh_.add_face( face );
    }


    // Colect all edges
    //
    // Similar to the vertex case, we use a map
    //
    //   unordered pair (of vertex indices) -> edge data
    //
    // to associate with each edge its corresponding data. If an edge
    // exists more than twice, it a complex edge and will be handled
    // later.

    typedef std::map< UnorderedPair, EdgeData > Pair2Edge;

    Pair2Edge pair2edge;

    for ( FI fi = mesh_.faces_begin(); fi != mesh_.faces_end(); ++fi )
    {
      FVI fvi = mesh_.fv_iter( fi );

      int id0 = mesh_.property( vertex_id, fvi ); ++fvi;
      int id1 = mesh_.property( vertex_id, fvi ); ++fvi;
      int id2 = mesh_.property( vertex_id, fvi ); ++fvi;

      if ( pair2edge.find( UnorderedPair( id0, id1 ) ) == pair2edge.end() )
	pair2edge[ UnorderedPair(id0,id1) ] = EdgeData( fi );
      else
	pair2edge[ UnorderedPair(id0,id1) ].add_face( fi );


      if ( pair2edge.find( UnorderedPair( id1, id2 ) ) == pair2edge.end() )
	pair2edge[ UnorderedPair(id1,id2) ] = EdgeData( fi );
      else
	pair2edge[ UnorderedPair(id1,id2) ].add_face( fi );

      if ( pair2edge.find( UnorderedPair( id2, id0 ) ) == pair2edge.end() )
	pair2edge[ UnorderedPair(id2,id0) ] = EdgeData( fi );
      else
	pair2edge[ UnorderedPair(id2,id0) ].add_face( fi );
    }


    //
    // Edges that do have exactly two incident faces, are manifold and
    // can be identified.
    //

    for ( typename Pair2Edge::iterator i = pair2edge.begin();
	  i != pair2edge.end(); ++i )
      if ( i->second.count() == 2 )
	identify( i->first.v0(), i->first.v1(),
		  i->second.face0(), i->second.face1() );



    //
    // Pinching
    //
    // We search for configurations like the following:
    //
    //   | hole 
    //   |     /
    //   a    /
    //   |   c
    //   |  /  mesh
    //   | /
    //   |/
    //   b						
    //
    // where vertex a and c have the same spatial position. We then
    // pinch the mesh along the edges (a,b) and (b,c) to get
    //
    //   |hole/
    //   |   /
    //   |  /
    //   | /
    //   |/     
    //   ac    
    //   ||  
    //   ||   mesh
    //   || 
    //   ||
    //   b						
    //


    mesh_.garbage_collection();

    // Collect all boundary halfedges
    
    std::vector< HH > boundary_halfedge;
  
    for ( HI hi = mesh_.halfedges_begin(); hi != mesh_.halfedges_end(); ++hi )
      if ( mesh_.is_boundary( hi ) )
	boundary_halfedge.push_back( hi );



    // Pinch boundaries

    for ( int i = 0; i < (int) boundary_halfedge.size(); ++i )
      if ( mesh_.is_boundary( boundary_halfedge[ i ] ) )
      {
	HH bh = boundary_halfedge[ i ];
	do {

	  HH nh = mesh_.next_halfedge_handle( bh );
	  HH ph = mesh_.prev_halfedge_handle( bh );
	  
	  int id0 = mesh_.property( vertex_id, mesh_.from_vertex_handle( ph ) );
	  int id1 = mesh_.property( vertex_id, mesh_.from_vertex_handle( bh ) );
	  int id2 = mesh_.property( vertex_id, mesh_.to_vertex_handle( bh ) );
	  
	  if ( id0 == id2 )
	    identify( id0, id1,
		      mesh_.face_handle( mesh_.opposite_halfedge_handle( ph ) ),
		      mesh_.face_handle( mesh_.opposite_halfedge_handle( bh ) ) );

	  bh = nh;
	} while ( mesh_.is_boundary( bh ) && bh != boundary_halfedge[i] );
      }
    
    std::cerr << "ok\n";

    //
    // Clean up
    //

    mesh_.remove_property( vertex_id );
    mesh_.garbage_collection();
  }


  // identify
  //
  // Stitch two faces along an edge. The edge is given by the vertex
  // IDs of its endpoints, the two faces are given as face
  // handles. The problem is that OpenMesh provides no explicit
  // operation for doing this. As we do not want to mess with
  // OpenMesh's internal data structure, we use the following trick:
  // We bride the gap between the two faces by a triangulated
  // quad. Then we decimate away this quad.
  //
  // Suppose we want to identify two edges: 
  //
  //     \        /
  //      *      *
  //      |      |  
  //      |      |  
  //   A  |      |  B
  //      |      |
  //      |      |
  //      |      |
  //      *      *
  //     /        \_
  //
  // First the edges are connected by a triagulated quad:
  //
  //     \        /
  //      *------*
  //      |     /|  
  //      |    / |  
  //   A  |   /  |  B
  //      |  /   |
  //      | /    |
  //      |/     |
  //      *------*
  //     /        \_
  //
  // Then the quad is decimated away
  //
  //     \  /
  //      **
  //      ||  
  //      ||  
  //   A  ||  B
  //      ||
  //      ||
  //      ||
  //      **
  //     /  \_
  //

  bool identify( int _id0, int _id1, FH _f0, FH _f1 )
  {
    // get edge handle of f0

    HH h0;
    for ( FHI fhi = mesh_.fh_iter( _f0 ); fhi; ++fhi )
    {
      int fid = mesh_.property( vertex_id, mesh_.from_vertex_handle( fhi ) );
      int tid = mesh_.property( vertex_id, mesh_.to_vertex_handle( fhi ) );

      if ( ( fid == _id0 && tid == _id1 ) ||
	   ( fid == _id1 && tid == _id0 ) )
	h0 = mesh_.opposite_halfedge_handle( fhi );
    }

    if ( ! h0.is_valid() || ! mesh_.is_boundary( h0 ) )
      return false;

    // get edge handle of f1

    HH h1;
    for ( FHI fhi = mesh_.fh_iter( _f1 ); fhi; ++fhi )
    {
      int fid = mesh_.property( vertex_id, mesh_.from_vertex_handle( fhi ) );
      int tid = mesh_.property( vertex_id, mesh_.to_vertex_handle( fhi ) );

      if ( ( fid == _id0 && tid == _id1 ) ||
	   ( fid == _id1 && tid == _id0 ) )
	h1 = mesh_.opposite_halfedge_handle( fhi );
    }

    if ( ! h1.is_valid() || ! mesh_.is_boundary( h1 ) )
      return false;


    // get adjacent vertices

    VH v0a = mesh_.from_vertex_handle( h0 );
    VH v0b = mesh_.to_vertex_handle( h0 );

    VH v1a = mesh_.from_vertex_handle( h1 );
    VH v1b = mesh_.to_vertex_handle( h1 );

    // Do the faces have compatible orientations?

    if ( mesh_.property( vertex_id, v0a ) !=
	 mesh_.property( vertex_id, v1b ) )
      return false;

    if ( mesh_.property( vertex_id, v0b ) !=
	 mesh_.property( vertex_id, v1a ) )
      return false;

    // Add the trianglulated quad

    HH bh0;
    FH fh0;

    if ( v0a != v1b )
    {
      fh0 = mesh_.add_face( v1b, v0a, v0b );

      if ( ! fh0.is_valid() )
	return false;

      for ( FHI fhi = mesh_.fh_iter( fh0 ); fhi; ++fhi )
	if ( mesh_.from_vertex_handle( fhi ) == v1b )
	  bh0 = fhi;
    }


    HH bh1;
    FH fh1;

    if ( v0b != v1a )
    {
      fh1 = mesh_.add_face( v0b, v1a, v1b );

      if ( ! fh1.is_valid() )
      {
	if ( fh0.is_valid() ) mesh_.delete_face( fh0 );
	return false;
      }

      for ( FHI fhi = mesh_.fh_iter( fh1 ); fhi; ++fhi )
	if ( mesh_.from_vertex_handle( fhi ) == v0b )
	  bh1 = fhi;
    }



    // Now it gets ugly because OpenMesh does not allow to collapse
    // edges if a hole is closed by this collapse. Thus we possibly
    // need to insert additional triangles, to close this hole
    // ourselves.

    if ( bh0.is_valid() )
    {
      HH h0 = mesh_.next_halfedge_handle( mesh_.opposite_halfedge_handle( bh0 ) );
      HH h1 = mesh_.next_halfedge_handle( h0 );
      HH h2 = mesh_.next_halfedge_handle( h1 );

      if ( h2 == mesh_.opposite_halfedge_handle( bh0 ) )
	mesh_.add_face( mesh_.to_vertex_handle( h0 ),
			mesh_.to_vertex_handle( h1 ),
			mesh_.to_vertex_handle( h2 ) );
    }

    if ( bh1.is_valid() )
    {
      HH h0 = mesh_.next_halfedge_handle( mesh_.opposite_halfedge_handle( bh1 ) );
      HH h1 = mesh_.next_halfedge_handle( h0 );
      HH h2 = mesh_.next_halfedge_handle( h1 );

      if ( h2 == mesh_.opposite_halfedge_handle( bh1 ) )
	mesh_.add_face( mesh_.to_vertex_handle( h0 ),
			mesh_.to_vertex_handle( h1 ),
			mesh_.to_vertex_handle( h2 ) );
    }



    if ( bh0.is_valid() && ! mesh_.is_collapse_ok( bh0 ) )
    {
      if ( fh0.is_valid() ) mesh_.delete_face( fh0 );
      if ( fh1.is_valid() ) mesh_.delete_face( fh1 );
      return false;
    }

    if ( bh1.is_valid() && ! mesh_.is_collapse_ok( bh1 ) )
    {
      if ( fh0.is_valid() ) mesh_.delete_face( fh0 );
      if ( fh1.is_valid() ) mesh_.delete_face( fh1 );
      return false;
    }

    if ( bh0.is_valid() )
      mesh_.collapse( bh0 );
    if ( bh1.is_valid() )
      mesh_.collapse( bh1 );

    return true;
  }

private:

  Mesh & mesh_;  

  OpenMesh::VPropHandleT< int > vertex_id;

};

//=============================================================================
#endif // CUTANDSTITCH_HH defined
//=============================================================================

