Raytrace Monte Carlo Code for SEE/modified modified for Xe+ on W or Ar+ on W
To compile go into source and type 
g++ -std=c++11 -w -o raytrace raytrace.cxx
must have gcc version 7.2.0, other versions may work.

execute raytrace to view list of tags
./raytrace

RTMC outputs:
  outputs.dat file contains each ray shot daughter yield and subsequent daughters grandaughter yield etc.
  in_out.dat file that contains each elements number of deposited and sputtered specie atom. 
  angle_dis.dat contains angle of incidence of every ray
  penetration.dat contains energy and z-position of rays that do not generate yield.
  outangles.dat contains outgoing energy and angle of only escaped rays wrt to z-axis. 
  ray_data.dat is optional under the tag print, shows every start and end position of rays. 
  distribution.dat is optional under the tag print, shows the energy and angle distribtuion of initial and outgoing rays. 
Contains optional tags 
STL file input name example:
     filename=foam.stl
Rays shot example: 
     Nsteps=10000
Initial Ray Energy:
     ParentE=100
Randomize beam orientation:
     random=yes
Cutoff Energy:
     cutoff=8.68
bounding volume hierarchy:
     bvh=yes
distribution of outgoing angles example:
     distribution=cosine
     distribution=table.txt
if distribution tab is used for table please ensure the name of file is species_"filename"
for example if you have a several tables to use, please name it to its initial species. For example if you have one table for Xe on W and another for W on W then write the tables as  Xe_AngDis.txt and W_AngDis.txt. So for the distribution tag write in distribution=AngDis.txt. if the file appears elsewhere write the path distribution=/home/tables/AngDis.txt 
species tag up to 4:
     species1=Xe
     species2=W
print of additional files for visual use or otherwise
     print=yes

--------------------------------------------------------------------
Remesh code located in remeshing directory:
to compile while in directory:
 g++ -std=c++11 -o remesh remesh.cxx
remesh.cxx contains evolutionary algorithm, a tri-tri self-intersections, deletion, and island removal.
remesh.cxx optional tags
STL file input name:
     filename=foam.stl 
Total number of rays shot:
     Nsteps=1000000
a scaling factor for the evolution algorithm
     scale=10
laplacian smoothing:
     smooth=yes
self-intersect deletion:
     del=yes
island removal (input is an element connected to atleast this many other elements):
     island=100
--------------------------------------------------------------------

Example folder contains two examples of raytrace code
Directories simple and for_hsge_use
simple is a 1 core process, execute ./run1.sh to see an example run. 
for_hsge_use is a 5 core parallel process, only the raytracing monte carlo part is parallel, remeshing is executed after. 
this should be executed in schedular cluster like hoffman2 in UCLA, this is just an example of a quick way to run multiple RTMC,
combine results at the end, and iterate. execute ./wrapperfunc.sh to see an example run
