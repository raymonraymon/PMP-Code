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
//== INCLUDES =================================================================

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include "QuadricT.hh"
#include <iostream>
#include <set>
#include <float.h>


//== IMPLEMENTATION ===========================================================


typedef OpenMesh::TriMesh_ArrayKernelT<>        Mesh;

Mesh                                            mesh;
OpenMesh::VPropHandleT<Quadricd>                vquadric;
OpenMesh::VPropHandleT<float>                   vprio;
OpenMesh::VPropHandleT<Mesh::HalfedgeHandle>    vtarget;


//-----------------------------------------------------------------------------


void  init();
bool  is_collapse_legal(Mesh::HalfedgeHandle _hh);
float priority(Mesh::HalfedgeHandle _heh);
void  decimate(unsigned int _n_vertices);
void  enqueue_vertex(Mesh::VertexHandle vh);


//-----------------------------------------------------------------------------


// access quadric of vertex _vh
Quadricd& quadric(Mesh::VertexHandle _vh) 
{ return mesh.property(vquadric, _vh); }

// access priority of vertex _vh
float& priority(Mesh::VertexHandle _vh)
{ return mesh.property(vprio, _vh); }

// access target halfedge of vertex _vh
Mesh::HalfedgeHandle& target(Mesh::VertexHandle _vh)
{ return mesh.property(vtarget, _vh); }


//-----------------------------------------------------------------------------


// compare functor for priority queue
struct VertexCmp
{
  bool operator()(Mesh::VertexHandle _v0, Mesh::VertexHandle _v1) const
  {
    // std::set needs UNIQUE keys -> handle equal priorities
    return ((priority(_v0) == priority(_v1)) ? 
	    (_v0.idx() < _v1.idx()) :
	    (priority(_v0) < priority(_v1)));
  }
};


std::set<Mesh::VertexHandle, VertexCmp>  queue;
  


//-----------------------------------------------------------------------------


int main(int argc, char **argv)
{
  if (argc < 4) 
  {
    std::cerr << "Usage: \n" 
	      << argv[0] << " <percentage>  <input_mesh>  <output_mesh>\n\n";
    exit(1);
  }


  // add required properties
  mesh.request_vertex_status();
  mesh.request_edge_status();
  mesh.request_face_status();
  mesh.request_face_normals();
  mesh.add_property(vquadric);
  mesh.add_property(vprio);
  mesh.add_property(vtarget);


  // read mesh
  OpenMesh::IO::read_mesh(mesh, argv[2]);
  std::cout << "#vertices: " << mesh.n_vertices() << std::endl;

  // compute normals & quadrics
  init();

  // decimate
  decimate((int)(atof(argv[1])*mesh.n_vertices()));
  std::cout << "#vertices: " << mesh.n_vertices() << std::endl;

  // write mesh
  OpenMesh::IO::write_mesh(mesh, argv[3]);
}


//-----------------------------------------------------------------------------


void init()
{
  // compute face normals
  mesh.update_face_normals();


  Mesh::VertexIter  v_it, v_end = mesh.vertices_end();
  Mesh::Point       n;
  double              a, b, c, d;

  for (v_it=mesh.vertices_begin(); v_it != v_end; ++v_it)
  {
    priority(v_it) = -1.0;
    quadric(v_it).clear();


    // calc vertex quadrics from incident triangles
    for (Mesh::VFIter vf_it(mesh, v_it.handle()); vf_it; ++vf_it)
    {
      // plane equation
      n = mesh.normal(vf_it);
      a = n[0];
      b = n[1];
      c = n[2];
      d = -(n | mesh.point(v_it));

      // plane -> quadric
      quadric(v_it) += Quadricd(a, b, c, d);
    }
  }
}


//-----------------------------------------------------------------------------


bool
is_collapse_legal(Mesh::HalfedgeHandle _hh)
{
  // collect vertices
  Mesh::VertexHandle v0, v1;
  v0 = mesh.from_vertex_handle(_hh);
  v1 = mesh.to_vertex_handle(_hh);


  // collect faces
  Mesh::FaceHandle fl = mesh.face_handle(_hh);
  Mesh::FaceHandle fr = mesh.face_handle(mesh.opposite_halfedge_handle(_hh));


  // backup point positions
  Mesh::Point p0 = mesh.point(v0);
  Mesh::Point p1 = mesh.point(v1);


  // topological test
  if (!mesh.is_collapse_ok(_hh))
    return false;


  // test boundary
  if (mesh.is_boundary(v0) && !mesh.is_boundary(v1))
    return false;

  

  // test normal flipping:
  //   if normal vector of a (non-degenerate) triangle changes by 
  //   more than pi degrees, return false.
  
  // simulate collapse
  mesh.set_point(v0, p1);

  // check for flipping normals
  Mesh::Scalar c(1.0);
  const float min_cos = (float) cos(0.25*M_PI);

  for (Mesh::ConstVertexFaceIter vf_it(mesh.cvf_iter(v0)); vf_it; ++vf_it) 
  {
    if (vf_it.handle() != fl && vf_it.handle() != fr)
    {
      Mesh::Point n0 = mesh.normal(vf_it);
      Mesh::Point n1 = mesh.calc_face_normal(vf_it);
      
      if ((c=(n0|n1)) < min_cos)
	break;
    }
  }

  // undo simulation changes
  mesh.set_point(v0, p0);

  if (c < 0.0) return false;

  
  // collapse passed all tests -> ok
  return true;
}


//-----------------------------------------------------------------------------


float priority(Mesh::HalfedgeHandle _heh)
{
  // return priority: the smaller the better
  // use quadrics to estimate approximation error

  Mesh::VertexHandle v0(mesh.from_vertex_handle(_heh));
  Mesh::VertexHandle v1(mesh.to_vertex_handle(_heh));

  Quadricd q = quadric(v0);  
  q += quadric(v1);

  return (float) q(mesh.point(v1));
}


//-----------------------------------------------------------------------------


void enqueue_vertex(Mesh::VertexHandle _vh)
{
  float                 prio, min_prio(FLT_MAX);
  Mesh::HalfedgeHandle  min_hh;


  // find best out-going halfedge
  for (Mesh::VOHIter vh_it(mesh, _vh); vh_it; ++vh_it)
  {
    if (is_collapse_legal(vh_it))
    {
      prio = priority(vh_it);
      if (prio != -1.0 && prio < min_prio)
      {
	min_prio = prio;
	min_hh   = vh_it.handle();
      }
    }
  }
    
    
  // update queue
  if (priority(_vh) != -1.0) 
  {
    queue.erase(_vh);
    priority(_vh) = -1.0;
  }

  if (min_hh.is_valid()) 
  {
    priority(_vh) = min_prio;
    target(_vh)   = min_hh;
    queue.insert(_vh);
  }
}


//-----------------------------------------------------------------------------


void  decimate(unsigned int _n_vertices)
{
  unsigned int nv(mesh.n_vertices());

  Mesh::HalfedgeHandle hh;
  Mesh::VertexHandle   to, from;
  Mesh::VVIter         vv_it;

  std::vector<Mesh::VertexHandle>            one_ring;
  std::vector<Mesh::VertexHandle>::iterator  or_it, or_end;



  // build priority queue
  Mesh::VertexIter  v_it  = mesh.vertices_begin(), 
                      v_end = mesh.vertices_end();
  
  queue.clear();
  for (; v_it!=v_end; ++v_it)
    enqueue_vertex(v_it.handle());



  while (nv > _n_vertices && !queue.empty())
  {
    // Decimate using priority queue:
    //   1) take 1st element of queue
    //   2) collapse this halfedge
    //   3) update queue


    // get 1st element
    Mesh::VertexHandle vh = *queue.begin();
    queue.erase(queue.begin());
    
    hh   = target(vh);
    to   = mesh.to_vertex_handle(hh);
    from = mesh.from_vertex_handle(hh);


    // store one-ring
    one_ring.clear();
    for (vv_it=Mesh::VVIter(mesh, from); vv_it; ++vv_it)
      one_ring.push_back(vv_it.handle());


    // perform collapse
    if (is_collapse_legal(hh))
    {
      mesh.collapse(hh);
      quadric(to) += quadric(from);
      --nv;
      if (nv % 100 == 0) std::cerr << nv << "\r";
    }


    // update queue
    for (or_it=one_ring.begin(), or_end=one_ring.end(); or_it!=or_end; ++or_it)
      enqueue_vertex(*or_it);
  }



  // clean up
  queue.clear();

  // now, delete the items marked to be deleted
  mesh.garbage_collection();
}


//-----------------------------------------------------------------------------
