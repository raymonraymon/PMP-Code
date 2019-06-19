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

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>

#include "CutAndStitcher.hh"

//=============================================================================

typedef OpenMesh::TriMesh_ArrayKernelT<OpenMesh::DefaultTraits> MyMesh;

//=============================================================================


int main( int argc, char ** argv )
{
  if ( argc != 3 )
  {
    std::cerr << "\n"
	      << "Usage :\n"
	      << "\n"
	      << "  " << argv[0] << " <infile> <outfile>\n"
	      << "\n"
	      << "Read a mesh, identify vertices at identical positions,\n"
	      << "perform cutting and stitching and write result."
	      << std::endl;
    exit( EXIT_FAILURE );
  }



  // read a mesh

  std::cerr << "Reading mesh ... ";

  MyMesh mesh;
  if ( ! OpenMesh::IO::read_mesh( mesh, argv[1] ) )
  {
    std::cerr << "Error: Could not read mesh\n";
    exit( EXIT_FAILURE );
  }

  std::cerr << "ok\n";


  // stitch the mesh

  CutAndStitcher< MyMesh > cas( mesh );
  cas.cut_and_stitch();


  // write result

  std::cerr << "Saving result ... ";

  if ( ! OpenMesh::IO::write_mesh( mesh, argv[2] ) )
  {
    std::cerr << "Error: Could not write mesh\n";
    exit( EXIT_FAILURE );
  }

  std::cerr << "ok\n";


  exit( EXIT_SUCCESS );
}

//=============================================================================
