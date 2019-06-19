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
//  CLASS DeformationViewer
//
//=============================================================================


#ifndef DEFORMATION_VIEWER_HH
#define DEFORMATION_VIEWER_HH


//== INCLUDES =================================================================


#include "../04-Fairing/FairingViewer.hh"



//== CLASS DEFINITION =========================================================

	      

class DeformationViewer : public FairingViewer
{
public:
   
  // default constructor
  DeformationViewer(const char* _title, int _width, int _height);

  // open mesh
  virtual bool open_mesh(const char* _filename);



private:

	virtual void keyboard(int key, int x, int y);
  virtual void draw(const std::string& _draw_mode);
  virtual void motion(int x, int y);
  virtual void mouse(int button, int state, int x, int y);

	enum Mode { MOVE, PICK, DRAG } mode_;
	void set_mode( Mode _mode );

	void glText(int x, int y, const std::string& _text);

	bool pick(int x, int y, Vec3f& _p);
	
  Mesh::Point& orig_point(Mesh::VertexHandle _vh) 
  { return mesh_.property(orig_point_, _vh); }

	void precompute_basis();
	void deform_mesh();


private:

	OpenMesh::VPropHandleT<Mesh::Point>  orig_point_;

	std::vector< std::vector< Mesh::Scalar > >  basis_;
	bool  compute_basis_;
	
	std::vector<Mesh::VHandle>  handle_vertices_;
	std::vector<Vec3f>          orig_handles_, moved_handles_;
	int                         active_handle_;
};


//=============================================================================
#endif // CURVVIEWERWIDGET_HH defined
//=============================================================================

