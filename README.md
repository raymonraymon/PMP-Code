# PMP-Code

===============================================================================

 This package contains example code for the full-day course

 M. Botsch, M. Pauly, L. Kobbelt, P. Alliez, B. Levy,
 "Geometric Modeling Based on Polygonal Meshes"

 held at SIGGRAPH 2007, San Diego, and Eurographics 2008, Crete.

===============================================================================



0. Files & Directories
----------------------

ReadMe.txt    This file.
src           Source code for the mesh data structure OpenMesh, the 
              isosurface extraction IsoEx, and the mesh course.
lib           The TAUCS library for a sparse direct Cholesky solver
acgmake       The compile tool acgmake (Linux)
MSVC          Project files for MS Visual Studio 2005 (Windows)
Xcode         Project files for Apple's Xcode 3.0 (Mac)



1. License
----------

The code examples are released under the GNU GPL, see src/CodeExamples/LICENSE.

OpenMesh is released under the GNU LGPL, see src/OpenMesh/LICENSE and 
http://www.openmesh.org

IsoEx is released under the GNU LGPL, see src/IsoEx/LICENSE and 
http://www-i8.informatik.rwth-aachen.de/software/software.html

TAUCS is released under the GNU LGPL, see lib/taucs-source/doc/taucs.pdf and
http://www.tau.ac.il/~stoledo/taucs/



2. Compiling
------------

Under Linux, please use the tool acgmake. Simply go to "cd src" 
and call "../acgmake/bin/acgmake". Platform-specific
configurations can be done in acgmake/configs/config.Linux. More
details on acgmake can be found in acgmake/docu/html/index.html.

For Windows the directory MSVC contains a solution file for Microsoft
Visual Studio .NET 2005.

The directory Xcode contains project files for Apple's Xcode 3.0 on MacOS.



3. Example Applications
-----------------------

The code examples include the following applications:

1. A simple GLUT-based triangle mesh viewer.
2. Visualization of curvatures and reflection lines for quality analysis.
3. Iterative Laplacian smoothing for noise removal.
4. Mesh fairing based on minimizing the surface thin-plate energy.
5. Discrete harmonic parameterization of disk-shaped surfaces.
6. Iterative mesh decimation based on the quadric error metric.
7. Three different algorithms for mesh repair.
8. Interactive mesh deformation.



4. 3D Models
------------

For the following triangle meshes are included in the package. Please
refer to their sources to check the respective terms of useage.

bunny.off: 
Courtesy of the Stanford 3D Scanning Repository. See
http://graphics.stanford.edu/data/3Dscanrep/

scanned_face.off: 
Courtesy Leif Kobbelt, RWTH Aachen.

max.off: 
Part of the Max Planck model, courtesy Leif Kobbelt and
Max-Planck-Institute for Computer Science, Saarbruecken.

fan.off:
Courtesy Hughes Hoppe, Microsoft Research USA.

teapot.off: 
The famous Utah teapot. See e.g. http://www.sjbaker.org/teapot/



--
Mario Botsch (botsch@inf.ethz.ch)
April 2008
