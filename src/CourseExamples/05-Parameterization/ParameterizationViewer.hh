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
//  CLASS ParameterizationViewer
//
//=============================================================================


#ifndef PARAM_VIEWER_HH
#define PARAM_VIEWER_HH


//== INCLUDES =================================================================


#include "../04-Fairing/FairingViewer.hh"



//== CLASS DEFINITION =========================================================

	      

class ParameterizationViewer : public FairingViewer
{
public:
   
  /// default constructor
  ParameterizationViewer(const char* _title, int _width, int _height);

  /// open mesh
  virtual bool open_mesh(const char* _filename);

  /// discrete harmonic parameterizaton
  void parameterize();



private:

  virtual void keyboard(int key, int x, int y);
  virtual void draw(const std::string& _draw_mode);

};


//=============================================================================
#endif // PARAM_VIEWER_HH defined
//=============================================================================

