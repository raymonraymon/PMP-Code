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
//  CLASS ParameterizationViewer - IMPLEMENTATION
//
//=============================================================================


//== INCLUDES =================================================================

#include "ParameterizationViewer.hh"
#include "../04-Fairing/TaucsSolver.hh"


//== IMPLEMENTATION ========================================================== 


ParameterizationViewer::
ParameterizationViewer(const char* _title, int _width, int _height)
  : FairingViewer(_title, _width, _height)
{ 
  mesh_.request_vertex_texcoords2D();

  set_draw_mode( add_draw_mode("Textured") );
}


//-----------------------------------------------------------------------------


bool
ParameterizationViewer::
open_mesh(const char* _filename)
{
  if (FairingViewer::open_mesh(_filename))
  {
    parameterize();
    return true;
  }
  return false;
}


//-----------------------------------------------------------------------------


void
ParameterizationViewer::
keyboard(int key, int x, int y)
{
  switch (key)
  {
    case ' ':
    {
      parameterize();
      glutPostRedisplay();
      break;
    }

    default:
    {
      CurvatureViewer::keyboard(key, x, y);
      break;
    }
  }
}


//-----------------------------------------------------------------------------


void
ParameterizationViewer::
parameterize()
{
  std::cout << "Harmonic parameterization: " << std::flush;

  Mesh::VertexIter      v_it, v_end(mesh_.vertices_end());
  Mesh::VertexHandle    vh;
  Mesh::HalfedgeHandle  hh;

  std::vector<Mesh::VertexHandle>  loop;



  // find 1st boundary vertex
  for (v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
    if (mesh_.is_boundary(v_it.handle()))
      break;


  // boundary found ?
  if (v_it == v_end)
  {
    std::cerr << "No boundary found\n";
    return;
  }


  // collect boundary loop
  vh = v_it.handle();
  hh = mesh_.halfedge_handle(vh);
  do 
  { 
    loop.push_back(mesh_.to_vertex_handle(hh));
    hh = mesh_.next_halfedge_handle(hh);
  }
  while (hh != mesh_.halfedge_handle(vh));



  // map boundary loop to unit circle in texture domain
  unsigned int i, n = loop.size();
  Mesh::Scalar  angle, l, length;

  for (i=0, length=0.0; i<n; ++i)
    length += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();

  for (i=0, l=0.0; i<n; ++i)
  {
    angle = l/length*2.0*M_PI;
    mesh_.set_texcoord2D(loop[i],
			 Vec2f(0.5*cos(angle)+0.5, 0.5*sin(angle)+0.5));
    l += (mesh_.point(loop[i]) - mesh_.point(loop[(i+1)%n])).norm();
  }


  
  // setup free vertices
  std::vector<Mesh::VertexHandle>   free_vertices;
  for (i=0, v_it=mesh_.vertices_begin(); v_it!=v_end; ++v_it)
  {
    mesh_.status(v_it).set_locked(mesh_.is_boundary(v_it));

    if (!mesh_.status(v_it).locked())
    {
      idx(v_it) = i++;
      free_vertices.push_back(v_it);
    }
  }



  TaucsSolver          solver;
  unsigned int         N = free_vertices.size();
  std::vector<Vec2d>   b(N), x(N);

  std::map<VertexHandle,double>  row;
  std::map<VertexHandle,double>::const_iterator r_it, r_end;



  // setup matrix and rhs
  for (i=0; i<N; ++i)
  {
    row.clear();
    setup_matrix_row(free_vertices[i], 1, -1.0, row);

    solver.begin_row();

    b[i] = Vec2d(0,0);

    for (r_it=row.begin(), r_end=row.end(); r_it!=r_end; ++r_it)
      if (mesh_.status(r_it->first).locked())
	b[i] -= r_it->second * (Vec2d) mesh_.texcoord2D(r_it->first);
      else
	solver.add_value(idx(r_it->first), r_it->second);

    solver.end_row();
  }


  // factorize & solve
  if (!solver.factorize() || !solver.vec_solve(b, x))
  {
    std::cerr << "Solving failed\n";
    exit(1);
  }


  // copy solution
  for (i=0; i<N; ++i)
    mesh_.set_texcoord2D(free_vertices[i], (Mesh::TexCoord2D)x[i]);


  std::cout << "done\n";
}


//-----------------------------------------------------------------------------


void 
ParameterizationViewer::
draw(const std::string& _draw_mode)
{
  if (_draw_mode == "Textured")
  {
    glEnable( GL_TEXTURE_2D ); // use texture from CurvatureViewer
    glEnable(GL_LIGHTING);
    glShadeModel(GL_SMOOTH);

    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    GL::glVertexPointer(mesh_.points());
    GL::glNormalPointer(mesh_.vertex_normals());
    GL::glTexCoordPointer(mesh_.texcoords2D());

    glDrawElements(GL_TRIANGLES, indices_.size(), GL_UNSIGNED_INT, &indices_[0]);

    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);

    glDisable( GL_TEXTURE_2D );
  }



  else FairingViewer::draw(_draw_mode);
}


//=============================================================================
