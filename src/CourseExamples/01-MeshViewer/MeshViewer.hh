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
//  CLASS MeshViewerWidget
//
//=============================================================================


#ifndef MESH_VIEWER_WIDGET_HH
#define MESH_VIEWER_WIDGET_HH


//== INCLUDES =================================================================


#include "GlutExaminer.hh"
#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>


//== CLASS DEFINITION =========================================================


class MeshViewer : public GlutExaminer
{
public:
   
  /// default constructor
  MeshViewer(const char* _title, int _width, int _height);

  /// open mesh
  virtual bool open_mesh(const char* _filename);

  /// update buffer with face indices
  void update_face_indices();

  /// draw the scene
  virtual void draw(const std::string& _draw_mode);
  

protected:

  typedef OpenMesh::TriMesh_ArrayKernelT<>  Mesh;
  

protected:

  Mesh                       mesh_;
  std::vector<unsigned int>  indices_;
};


//=============================================================================
#endif // MESH_VIEWER_WIDGET_HH defined
//=============================================================================

