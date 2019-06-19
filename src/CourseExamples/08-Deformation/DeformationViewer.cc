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
//  CLASS DeformationViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include "DeformationViewer.hh"
#include "../04-Fairing/TaucsSolver.hh"
#include <OpenMesh/Tools/Utils/Timer.hh>
#include <vector>
#include <float.h>



//== IMPLEMENTATION ========================================================== 


DeformationViewer::
DeformationViewer(const char* _title, int _width, int _height)
: FairingViewer(_title, _width, _height)
{ 
  mesh_.add_property(orig_point_);
	compute_basis_ = true;
  init();
	set_mode(MOVE);
}


//-----------------------------------------------------------------------------


void
DeformationViewer::
keyboard(int key, int x, int y)
{
  switch (key)
  {
    case ' ':
    {
      std::cout << "Deforming..." << std::flush;
			deform_mesh();
      glutPostRedisplay();
      std::cout << "done\n";
      break;
    }
			
		case 'm':
		{
			set_mode(MOVE);
      glutPostRedisplay();
			break;
		}
			
		case 'p':
		{
			set_mode(PICK);
      glutPostRedisplay();
			break;
		}
			
		case 'd':
		{
			set_mode(DRAG);
      glutPostRedisplay();
			break;
		}
			
		case 'r': // reset
		{
			set_mode(MOVE);

			Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
			for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
				mesh_.set_point(v_it, orig_point(v_it));

			mesh_.update_normals();

			moved_handles_ = orig_handles_;
			
      glutPostRedisplay();
			break;
		}
			
    default:
    {
      MeshViewer::keyboard(key, x, y);
      break;
    }
  }
}


//-----------------------------------------------------------------------------


void
DeformationViewer::
set_mode( Mode _mode )
{
	switch(mode_ = _mode)
	{
		case MOVE:
			glutSetCursor( GLUT_CURSOR_LEFT_ARROW );
			break;

		case PICK:
			glutSetCursor( GLUT_CURSOR_CROSSHAIR );
			break;

		case DRAG:
			glutSetCursor( GLUT_CURSOR_INFO );
			break;
	}	
}


//-----------------------------------------------------------------------------


void
DeformationViewer::
glText(int x, int y, const std::string& _text)
{
  GLint viewport[4];
  glGetIntegerv(GL_VIEWPORT, viewport);
	
  // set raster pos
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0.0, (GLfloat) viewport[2], 0.0, (GLfloat) viewport[3]);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glRasterPos2i(x, y);
	
  // draw characters
  std::string::const_iterator s_it(_text.begin()), s_end(_text.end());
  for (; s_it!=s_end; ++s_it)
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *s_it);
	
  // restore matrices
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}


//-----------------------------------------------------------------------------


bool
DeformationViewer::
open_mesh(const char* _filename)
{
  // load mesh
  if (MeshViewer::open_mesh(_filename))
  {
    // store original vertex positions
		Mesh::VertexIter v_it, v_end(mesh_.vertices_end());
		for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
			orig_point(v_it) = mesh_.point(v_it);

		// compute cotang weights  (method of CurvatureViewer)
		calc_weights();
		
		// need to (re-)compute basis
		compute_basis_ = true;
		
    return true;
  }
  return false;
}


//-----------------------------------------------------------------------------


void 
DeformationViewer::
draw(const std::string& _draw_mode)
{
	// draw mesh
  CurvatureViewer::draw(_draw_mode);
	
	
	// draw spheres
	GLfloat mat_sphere[4] = {1.0, 0.0, 0.0, 1.0};
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_sphere);
	
	for (unsigned int i=0; i<moved_handles_.size(); ++i)
	{
		glPushMatrix();
		glTranslatef( moved_handles_[i][0],
									moved_handles_[i][1],
									moved_handles_[i][2] );
		glutSolidSphere(0.05*radius_, 20, 20);
		glPopMatrix();
	}


	GLfloat mat_mesh[4] = {0.4, 0.4, 0.4, 1.0};
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_mesh);

	
	// draw text
	glDisable(GL_LIGHTING);
	glColor3f(1,1,1);
	switch (mode_)
	{
		case MOVE:
			glText(10, 10, "Move");
			break;

		case PICK:
			glText(10, 10, "Pick");
			break;

		case DRAG:
			glText(10, 10, "Drag");
			break;
	}
}


//-----------------------------------------------------------------------------


bool
DeformationViewer::
pick(int x, int y, Vec3f& _p)
{
	GLdouble  modelview[16], projection[16];
	GLint     viewport[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);
	
	
	// read depth buffer value at (x, y_new)
	float  z;
	int    y_new = viewport[3] - y; // in OpenGL y is zero at the 'bottom'
	glReadPixels(x, y_new, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z );
	
	
	// reverse projection to get 3D point
	double pos[3];
	gluUnProject(x, y_new, z, 
							 modelview, 
							 projection, 
							 viewport, 
							 &pos[0], &pos[1], &pos[2]);
	
	
	if (z != 1.0f)
	{
		_p = Vec3f(pos[0], pos[1], pos[2]);
		return true;
	}
	
	return false;
}


//-----------------------------------------------------------------------------


void 
DeformationViewer::
mouse(int button, int state, int x, int y)
{
	switch (mode_)
	{
		// move the mesh
		case MOVE:
		{
			MeshViewer::mouse(button, state, x, y);
			break;
		}
		
			
		// create a new handle point 
		case PICK:
		{
		  if (state == GLUT_DOWN)
			{	
				Vec3f p;
				if (pick(x, y, p))
				{
					// find closest vertex
					Mesh::VertexIter  v_it, v_end(mesh_.vertices_end());
					Mesh::VHandle     vh;
					Mesh::Scalar      d, dmin(FLT_MAX);
					for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
					{
						d = (mesh_.point(v_it) - p).sqrnorm();
						if (d < dmin) { dmin = d; vh = v_it; }
					}

					// lock vertex to later recognize it as handle
					mesh_.status(vh).set_locked(true);

					// store handle and remember its index
					handle_vertices_.push_back(vh);
					idx(vh) = handle_vertices_.size()-1;

					// snap to this vertex & store sphere centers
					p = mesh_.point(vh);
					
					// create handle sphere
					orig_handles_.push_back(p);
					moved_handles_.push_back(p);
					active_handle_ = moved_handles_.size()-1;

					
					// need to (re-)compute basis
					compute_basis_ = true;
				}
				
				glutPostRedisplay();
			}

			break;
		}
			

			
		case DRAG:
		{
			if (state == GLUT_DOWN)
			{	
				if (!moved_handles_.empty())
				{
					Vec3f p;
					if (pick(x, y, p))
					{
						// find closest handle
						int    imin(-1);
						float  d, dmin(FLT_MAX);
						for (unsigned int i=0; i<moved_handles_.size(); ++i)
						{
							d = (moved_handles_[i] - p).sqrnorm();
							if (d<dmin) { dmin=d; imin=i; }
						}
						
						// mark this one as active handle
						active_handle_ = imin;
					}
				}
			}
			break;
		}
	}
}


//-----------------------------------------------------------------------------


void 
DeformationViewer::
motion(int x, int y)
{
	switch (mode_)
	{
		case MOVE:
		{
			MeshViewer::motion(x, y);
			break;
		}

			
		case DRAG:
		{
			GLdouble  modelview[16], projection[16];
			GLint     viewport[4];
			glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
			glGetDoublev(GL_PROJECTION_MATRIX, projection);
			glGetIntegerv(GL_VIEWPORT, viewport);
			double y_new = viewport[3] - y; // in OpenGL y is zero at the 'bottom'


			GLdouble ox, oy, oz, wx, wy, wz;

			gluProject( moved_handles_[active_handle_][0],
									moved_handles_[active_handle_][1],
									moved_handles_[active_handle_][2],
									modelview,
									projection,
									viewport,
									&wx, &wy, &wz );
			
			gluUnProject( x, y_new, wz,
										modelview, 
										projection, 
										viewport, 
										&ox, &oy, &oz );
			
			moved_handles_[active_handle_] = Mesh::Point(ox, oy, oz);


			// compute deformation on the fly...
			deform_mesh();

			glutPostRedisplay();
			break;
		}


	  case PICK:
	    break;
	}
}


//-----------------------------------------------------------------------------


void
DeformationViewer::
precompute_basis()
{
	std::cout << "Initialize deformation..." << std::flush;
	OpenMesh::Utils::Timer t; t.start();

	
	Mesh::VertexIter  v_it, v_end(mesh_.vertices_end());
	
	
	// collect free vertices, assign indices to them
	// handle vertices have been collected when picking them
	std::vector<Mesh::VertexHandle>  free_vertices;
  for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
	{
		if (!mesh_.status(v_it).locked())
		{
			free_vertices.push_back(v_it);
			idx(v_it) = free_vertices.size()-1;
		}
	}

	
	
	int nV = mesh_.n_vertices();
	int nH = handle_vertices_.size();
	int nF = free_vertices.size();
	
	assert(nH == (int)orig_handles_.size());
	

	// allocate basis functions, initialize them to zero
	// one for each vertex with #handles entries
	basis_.clear();
	basis_.resize(nV);
	for (int i=0; i<nV; ++i)
		basis_[i].resize(nH, 0.0);
	

	// allocate right hand sides (one for each handle)
	std::vector< std::vector< double > >  b(nH);
	for (int h=0; h<nH; ++h)
		b[h].resize(nF, 0.0);
	
	
	// setup bilaplacian matrix row by row (use function from FairingViewer)
	// setup all right-hand sides for the different basis functions
	TaucsSolver  solver;
	std::map<VertexHandle,double>  row;
  std::map<VertexHandle,double>::const_iterator r_it, r_end;
	
  for (int i=0; i<nF; ++i)
  {
		row.clear();
		setup_matrix_row(free_vertices[i], 2, 1.0, row);
		solver.begin_row();
		
		for (r_it=row.begin(), r_end=row.end(); r_it!=r_end; ++r_it)
			if (!mesh_.status(r_it->first).locked())
				solver.add_value(idx(r_it->first), r_it->second);
			else
				b[idx(r_it->first)][i] -= r_it->second;
				
		solver.end_row();
	}
		
		
	// factorize matrix
	if (!solver.factorize())
		exit(1);
	
	
	
	// precompute basis functions, solve system for each rhs
	std::vector<double>  x(nF);
	for (int h=0; h<nH; ++h)
	{
		if (!solver.solve(b[h], x))
			exit(1);
		
		// store basis function h for all free vertices
	  for (int i=0; i<nF; ++i)
			basis_[ free_vertices[i].idx() ][h] = x[i];
		
		// store basis function h for handle vertices
		basis_[ handle_vertices_[h].idx() ][h] = 1.0;
		
		std::cout << '.' << std::flush;
	}
	
	
	compute_basis_ = false;


	t.stop();
	std::cout << "done (" << t.seconds() << "s)\n";
}


//-----------------------------------------------------------------------------


void 
DeformationViewer::
deform_mesh()
{
	// do we have to compute the deformation basis functions?
	if (compute_basis_) precompute_basis();


	Mesh::VertexIter  v_it, v_end(mesh_.vertices_end());
	Mesh::Point       p;
	int               h, nH=orig_handles_.size();

	
	// collect handle displacements
	std::vector<Vec3f> displacements(nH);
	for (h=0; h<nH; ++h)
		displacements[h] = moved_handles_[h] - orig_handles_[h];

	
  for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
	{
		const std::vector<Mesh::Scalar>& basis = basis_[v_it.handle().idx()];
		
		p = orig_point(v_it);
		for (h=0; h<nH; ++h)
			p += displacements[h] * basis[h];
		
		mesh_.set_point(v_it, p);
	}
	
	
	mesh_.update_normals();
}


//=============================================================================
