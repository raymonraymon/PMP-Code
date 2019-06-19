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
//
// This program provides a simplified version of
//
// F.S. Nooruddin, G. Turk
// Simplification and Repair of Polygonal Models Using Volumetric Techniques
// IEEE Trans. on Visualization and Computer Graphics, 9(2), 2003,191-205
//
//=============================================================================

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/IO/BinaryHelper.hh>
#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>

#include <IsoEx/Grids/ScalarGridT.hh>
#include <IsoEx/Extractors/MarchingCubesT.hh>

#include <limits>

#include <iostream>

#include "Geometry.hh"
#include "Tree.hh"
#include "DepthImage.hh"

//=============================================================================

using namespace IsoEx;
using namespace OpenMesh;
using namespace OpenMesh::IO;

//=============================================================================

typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::DefaultTraits> MyMesh;

//=============================================================================

int main( int argc, char ** argv )
{
  
  if ( argc != 7 )
  {
    std::cerr << "Usage:\n"
    << "\n"
    << "  volrepair <volres> <imgres> <ssf> <smooth> <input> <output>\n"
    << "\n"
    << "where\n"
    << "  <volres> : resolution of the volume (voxels)\n"
    << "  <imgres> : resolution of the images (pixels)\n"
    << "  <ssf>    : supersampling factor\n"
    << "  <smooth> : number of smoothing steps\n"
    << "  <input>  : input model\n"
    << "  <output> : output volume (has to be post-processed\n"
    << "             with marching cubes)\n";
    
    std::cerr << "\n"
    << "Example:\n"
    << "  volrepair 100 300 3 1 octahedron.off result.off\n";
    
    exit( EXIT_FAILURE );
  }
  
  int volres = atoi( argv[1] );
  int imgres = atoi( argv[2] );
  int ssf    = atoi( argv[3] );
  int smooth = atoi( argv[4] );
  
  
  //---------------------------------------------------------------------------
  //
  // Read a mesh
  //
  //---------------------------------------------------------------------------
  
  std::cerr << "Reading the mesh ... ";
  
  MyMesh mesh;
  if ( ! read_mesh( mesh, argv[5] ) )
  {
    std::cerr << "error: Could not open input\n";
    exit( EXIT_FAILURE );
  }
  
  std::cerr << "ok\n";
  
  
  
  //---------------------------------------------------------------------------
  //
  // Build up the tree
  //
  //---------------------------------------------------------------------------
  
  
  std::cerr << "Building up the kd-tree ... ";
  
  Tree tree;
  
  for ( MyMesh::FaceIter fi = mesh.faces_begin();
        fi != mesh.faces_end(); ++fi )
  {
    MyMesh::FaceVertexIter fvi = mesh.fv_iter( fi );
    
    Vec3d p0( mesh.point( fvi ) ); ++fvi;
    Vec3d p1( mesh.point( fvi ) ); ++fvi;
    Vec3d p2( mesh.point( fvi ) ); ++fvi;
    
    tree.push_back( Triangle( p0, p1, p2 ) );
  }
  
  tree.build_kdtree();
  
  std::cerr << "ok\n";
  
  
  
  
  //---------------------------------------------------------------------------
  //
  // Determine scene geometry
  //
  //---------------------------------------------------------------------------
  
  // Make the bounding box square and slightly enlarge it to encompass
  // the whole model.
  
  Vec3d bbmin = tree.bbmin();
  Vec3d bbmax = tree.bbmax();
  
  Vec3d scene_center = 0.5 * ( bbmin + bbmax );
  
  Vec3d diag = bbmax - bbmin;
  
  double r = std::max( diag[0], std::max( diag[1], diag[2] ) ) / 2;
  
  bbmin = scene_center - 1.1 * r * Vec3d( 1, 1, 1 );
  bbmax = scene_center + 1.1 * r * Vec3d( 1, 1, 1 );
  diag  = bbmax - bbmin;
  
  double scene_radius = diag.norm() / 2.0;
  
  
  //---------------------------------------------------------------------------
  //
  // Record depth images
  //
  //---------------------------------------------------------------------------
  
  std::cerr << "Recording depth images ... ";
  
  Vec3d projection_dir[] = {
    Vec3d( 1, 1, 1 ),  Vec3d( 1, 1,-1 ),  Vec3d( 1,-1, 1 ),
    Vec3d(-1, 1, 1 ),
    Vec3d( 0, 1,-1 ),  Vec3d(-1, 0, 1 ),  Vec3d( 1,-1, 0 ),
    Vec3d( 0, 1, 1 ),  Vec3d( 1, 0, 1 ),  Vec3d( 1, 1, 0 ),
    Vec3d( 1, 0, 0 ),  Vec3d( 0, 1, 0 ),  Vec3d( 0, 0, 1 ),
  };
  
  int images = 13;
  
  std::vector< DepthImage > depth_image;
  
  for ( int i = 0; i < images; ++i )
  {
    std::cerr << i << " ";
    DepthImage image( projection_dir[i],
                      scene_center,
                      scene_radius,
                      imgres, imgres );
    image.record( tree );
    
    depth_image.push_back( image );
  }
  
  std::cerr << "ok\n";
  
  
  //---------------------------------------------------------------------------
  //
  // Convert the input into a volumetric representation
  //
  //---------------------------------------------------------------------------
  
  
  std::cerr << "Convert to volume ... ";
  
  
  typedef ScalarGridT< float > ScalarGrid;
  
  ScalarGrid grid( (Vec3f) bbmin,
                   Vec3f( float(bbmax[0]-bbmin[0]), 0.0f, 0.0f ),
                   Vec3f( 0.0f, float(bbmax[1]-bbmin[1]), 0.0f ),
                   Vec3f( 0.0f, 0.0f, float(bbmax[2]-bbmin[2]) ),
                   volres, volres, volres );
  
  // Assume all voxels to be exterior (=0)  
  for ( int i2 = 0; i2 < volres; ++i2 )
    for ( int i1 = 0; i1 < volres; ++i1 )
      for ( int i0 = 0; i0 < volres; ++i0 )
        grid( i0, i1, i2 ) = 0;
  
  double sinv  = 1.0 / ssf;
  double sinv2 = 0.5 / ssf;
  
  for ( int i2 = 0; i2 < volres; ++i2 )
  {
    std::cerr << i2 << " ";
    for ( int i1 = 0; i1 < volres; ++i1 )
      for ( int i0 = 0; i0 < volres; ++i0 )
      {
        Vec3d corner0( bbmin[0] + i0 / double( volres ) * diag[0],
                       bbmin[1] + i1 / double( volres ) * diag[1],
                       bbmin[2] + i2 / double( volres ) * diag[2] );
        Vec3d corner1( bbmin[0] + (i0+1) / double( volres ) * diag[0],
                       bbmin[1] + (i1+1) / double( volres ) * diag[1],
                       bbmin[2] + (i2+1) / double( volres ) * diag[2] );
        
        // Supersampling
        
        float votes_for_interior = 0;
        
        
        for ( double f0 = sinv2; f0 < 1.0; f0 += sinv )
          for ( double f1 = sinv2; f1 < 1.0; f1 += sinv )
            for ( double f2 = sinv2; f2 < 1.0; f2 += sinv )
            {
              Vec3d sample( (1-f0)*corner0[0] + f0*corner1[0],
                            (1-f1)*corner0[1] + f1*corner1[1],
                            (1-f2)*corner0[2] + f2*corner1[2] );
              
              int i = 0;
              for ( ; i < images; ++i )
                if ( depth_image[i].classify( sample ) == DepthImage::EXTERIOR )
                  break;
              
              if ( i == images )
                votes_for_interior += 1;
            }
              
              grid( i0, i1, i2 ) += votes_for_interior;
      }
  }  
    
    
    float f = 1.0f / ( ssf * ssf * ssf );
    for ( int i2 = 0; i2 < volres; ++i2 )
      for ( int i1 = 0; i1 < volres; ++i1 )
        for ( int i0 = 0; i0 < volres; ++i0 )
          grid( i0, i1, i2 ) = 0.5f - f * grid( i0, i1, i2 );
    
    std::cerr << "ok\n";
    
    
    //---------------------------------------------------------------------------
    //
    // Smooth the volume
    //
    //---------------------------------------------------------------------------
    
    
    std::cerr << "Smooth the volume ... ";
    
    ScalarGrid grid2( (Vec3f) bbmin,
                      Vec3f( float(bbmax[0]-bbmin[0]), 0, 0 ),
                      Vec3f( 0, float(bbmax[1]-bbmin[1]), 0 ),
                      Vec3f( 0, 0, float(bbmax[2]-bbmin[2]) ),
                      volres, volres, volres );
    
    // Gaussian smoothing, width = 2, sigma = 0.7
    
    float mask0 = 0.7748f / ( 2.0f * ( 0.7748f + 0.1f ) );
    float mask1 = 0.1f / ( 2.0f * ( 0.7748f + 0.1f ) );
    
    // smooth only once
    
    for ( int i = 0; i < smooth; ++i )
    {
      std::cerr << (i+1) << " ";
      
      // Smooth along x-axis
      
      for ( int i2 = 0; i2 < volres; ++i2 )
        for ( int i1 = 0; i1 < volres; ++i1 )
          for ( int i0 = 0; i0 < volres; ++i0 )
          {
            grid2( i0, i1, i2 )
            = mask1 * grid( ( i0 + volres  - 2 ) % volres, i1, i2 )
            + mask0 * grid( ( i0 + volres  - 1 ) % volres, i1, i2 )
            + mask0 * grid( ( i0               ) % volres, i1, i2 )
            + mask1 * grid( ( i0           + 1 ) % volres, i1, i2 );
          }
            
            for ( int i2 = 0; i2 < volres; ++i2 )
              for ( int i1 = 0; i1 < volres; ++i1 )
                for ( int i0 = 0; i0 < volres; ++i0 )
                  grid( i0, i1, i2 ) = grid2( i0, i1, i2 );
      
      // Smooth along y-axis
      
      for ( int i2 = 0; i2 < volres; ++i2 )
        for ( int i1 = 0; i1 < volres; ++i1 )
          for ( int i0 = 0; i0 < volres; ++i0 )
          {
            grid2( i0, i1, i2 )
            = mask1 * grid( i0, ( i1 + volres  - 2 ) % volres, i2 )
            + mask0 * grid( i0, ( i1 + volres  - 1 ) % volres, i2 )
            + mask0 * grid( i0, ( i1               ) % volres, i2 )
            + mask1 * grid( i0, ( i1           + 1 ) % volres, i2 );
          }
            
            for ( int i2 = 0; i2 < volres; ++i2 )
              for ( int i1 = 0; i1 < volres; ++i1 )
                for ( int i0 = 0; i0 < volres; ++i0 )
                  grid( i0, i1, i2 ) = grid2( i0, i1, i2 );
      
      // Smooth along z-axis
      
      for ( int i2 = 0; i2 < volres; ++i2 )
        for ( int i1 = 0; i1 < volres; ++i1 )
          for ( int i0 = 0; i0 < volres; ++i0 )
          {
            grid2( i0, i1, i2 )
            = mask1 * grid( i0, i1, ( i2 + volres  - 2 ) % volres )
            + mask0 * grid( i0, i1, ( i2 + volres  - 1 ) % volres )
            + mask0 * grid( i0, i1, ( i2               ) % volres )
            + mask1 * grid( i0, i1, ( i2           + 1 ) % volres );
          }
            
            for ( int i2 = 0; i2 < volres; ++i2 )
              for ( int i1 = 0; i1 < volres; ++i1 )
                for ( int i0 = 0; i0 < volres; ++i0 )
                  grid( i0, i1, i2 ) = grid2( i0, i1, i2 );
    }
    
    std::cerr << "ok\n";
    
    
    
    //---------------------------------------------------------------------------
    //
    // Apply marching cubes and write result
    //
    //---------------------------------------------------------------------------
    
    
    
    std::cerr << "Extract result ... ";
    
    MyMesh result_mesh;
    marching_cubes( grid, result_mesh ); 
    
    
    if ( ! write_mesh( result_mesh, argv[6]) )
    {
      std::cerr << "Error: Could not open output file\n";
      exit( EXIT_FAILURE );
    }
    
    
    std::cerr << "ok\n";
    
    exit( EXIT_SUCCESS );
    
}


//=============================================================================
