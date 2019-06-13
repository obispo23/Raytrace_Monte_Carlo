//Read STL file, store x, y, and z for each vertex of a triangle mesh.
//Raytracing Monte Carlo for rectangular mesh
//Andrew Alvarado
//UCLA Department of Materials Science
//Marian Research group
//started July 24, 2017
//current edited dat February 5, 2019
#include<fstream>
#include<string>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<random>
#include<ctime>
#include"vec_bvh.h"
#include"source_f.h"
using namespace std;

/*
vec_bvh.h holds functions. Create vector class for simplification and function traingularintersection based off of Moller-Trumbore algorithm
*/

//#define eps 0.000001

//---------------------------------------------------------------

int main(int argc, char **argv)
{
  clock_t t;
  t = clock();
  //input and output declaration
  ofstream outfile, outfile1, outfile2, outfile3, outfile4, outfile5;
  ifstream infile;
  //Place holders for chucking excess from stl files
  string temp1, temp2, temp3, temp4;
  int tn = 0;
  int lines = 0;
  char buffer[32];

  string input, first, last, filename = "mesh.stl", distribution;
  double ParentE = 100;
  int Nsteps = 10000; // 1E4  # of beams shot
  bool random = false, print = false, bvh = false, substrate=false;
  double cutoff = 1.00;
  int boxes = 1; //BVH
  string species[4] = {"Xe", "W", "Ar", "Mo"};
  if(argc > 1){
    for(int i = 1; i < argc; i++){
      input = argv[i];
      grabber(input, first, last, '=');
      if(first == "ParentE")
	ParentE = stod(last);
      if(first == "Nsteps")
	Nsteps = stod(last);
      if(first == "filename")
	filename = last;
      if(first == "random")
	random = true;
      if(first == "cutoff")
	cutoff = stod(last);
      if(first == "print")
	print = true;
      if(first == "bvh"){
	bvh = true;
	boxes = stoi(last);
      }
      if(first == "species1")
	species[0] = last;
      if(first == "species2")
	species[1] = last;
      if(first == "species3")
	species[2] = last;
      if(first == "species4")
	species[3] = last;
      if(first == "distribution")
	distribution = last;
      if(first == "substrate")
	substrate = true;
    }
    
  }
  else{
    cout << "list of arguments: ParentE, Nsteps, filename, random, cutoff, print, bvh, species1, species2, species3, species4, distribution, substrate" << endl;
      return 0;
  }
  //string backdirect = "../../";
  //filename = backdirect+filename;
  snprintf(buffer, sizeof(char)*32, filename.c_str()); //flatmesh.stl");
  cout << " file used: " << buffer << endl;
  //----------------------------mesh input--------------------------------------
   infile.open(buffer); //<----STL file here
   if(infile.good() == false){
     cout << buffer << " not found " << endl;
     return 1;
   }
   //count number of triangles and lines in file
  while(!infile.eof()){
    infile >> temp1;
    if(temp1 == "facet"){
      tn++;
    }
    getline(infile,temp2);
    lines++;
  }
  infile.close();
  cout << tn << " triangles in " << lines << " lines" << endl;
  infile.open(buffer);

  //Data storing declaration --pointers
  ray rays;
  elements *face;
  face = new elements[tn];
  for(int i =0; i < tn; i++)
    face[i] = {};
  int n=0; //counter for facets

  double xmax = 0, ymax = 0,  zmax = 0; //solid object size dimensions
  double xmin = 0, ymin = 0,  zmin = 0; //for BC:width, length, and height
  //vector class store
  while(n < tn){
    infile >> temp1;
    if(temp1 == "facet"){
      infile >> temp2;
      infile >> face[n].norm.x >> face[n].norm.y >> face[n].norm.z;
      if(fabs(face[n].norm.x) < eps) face[n].norm.x = 0; //making sure zeros are zeros.
      if(fabs(face[n].norm.y) < eps) face[n].norm.y = 0;
      if(fabs(face[n].norm.z) < eps) face[n].norm.z = 0;
    }
    if(temp1 == "loop"){
	for(int i = 0; i < 3; i++){
	  infile >> temp4 >> face[n].vert[i].x >> face[n].vert[i].y >> face[n].vert[i].z;
	  getline(infile,temp1);
	  //---------------find length of mesh object------------------------
	  if(face[n].vert[i].x > xmax)
	    xmax = face[n].vert[i].x; 
	  if(face[n].vert[i].x < xmin)
	    xmin = face[n].vert[i].x;
	  
	  if(face[n].vert[i].y > ymax)
	    ymax = face[n].vert[i].y; 
	  if(face[n].vert[i].y < ymin)
	    ymin = face[n].vert[i].y;
	  
	  if(face[n].vert[i].z > zmax)
	    zmax = face[n].vert[i].z;
	  if(face[n].vert[i].z < zmin)
	    zmin = face[n].vert[i].z;
	}
      n++;
    }
  }
  infile.close();
  if(boxes != 8)
    boxes = 8;
  vector< vector<int> > bvh_id;
  Vector bvmin[boxes];
  Vector bvmax[boxes];
  if(bvh == true){    
    double length = abs(xmax - xmin);
    double width  = abs(ymax - ymin);
    double height = abs(zmax - zmin);
    //split into even and smaller rectangles
    //calculate centers of each triangle sort into the bounding volumes
    bvmin[0] = Vector(0.0*length+xmin, 0.0*width+ymin, 0.0*height+zmin);
    bvmin[1] = Vector(0.5*length+xmin, 0.0*width+ymin, 0.0*height+zmin);
    bvmin[2] = Vector(0.0*length+xmin, 0.5*width+ymin, 0.0*height+zmin);
    bvmin[3] = Vector(0.5*length+xmin, 0.5*width+ymin, 0.0*height+zmin);
    bvmin[4] = Vector(0.0*length+xmin, 0.0*width+ymin, 0.5*height+zmin);
    bvmin[5] = Vector(0.5*length+xmin, 0.0*width+ymin, 0.5*height+zmin);
    bvmin[6] = Vector(0.0*length+xmin, 0.5*width+ymin, 0.5*height+zmin);
    bvmin[7] = Vector(0.5*length+xmin, 0.5*width+ymin, 0.5*height+zmin);

    bvmax[0] = Vector(0.5*length+xmin, 0.5*width+ymin, 0.5*height+zmin);
    bvmax[1] = Vector(1.0*length+xmin, 0.5*width+ymin, 0.5*height+zmin);
    bvmax[2] = Vector(0.5*length+xmin, 1.0*width+ymin, 0.5*height+zmin);
    bvmax[3] = Vector(1.0*length+xmin, 1.0*width+ymin, 0.5*height+zmin);
    bvmax[4] = Vector(0.5*length+xmin, 0.5*width+ymin, 1.0*height+zmin);
    bvmax[5] = Vector(1.0*length+xmin, 0.5*width+ymin, 1.0*height+zmin);
    bvmax[6] = Vector(0.5*length+xmin, 1.0*width+ymin, 1.0*height+zmin);
    bvmax[7] = Vector(1.0*length+xmin, 1.0*width+ymin, 1.0*height+zmin);
    
    for(int i = 0; i < boxes; i++){
      vector<int> row;
      bvh_id.push_back(row);
    }
    for(int i = 0; i < tn; i++){
      bool used = false;
      for(int j = 0; j < boxes; j++){
	if((face[i].vert[0].x >= bvmin[j].x) &&
	   (face[i].vert[0].y >= bvmin[j].y) &&
	   (face[i].vert[0].z >= bvmin[j].z) &&
	   (face[i].vert[0].x <= bvmax[j].x) &&
	   (face[i].vert[0].y <= bvmax[j].y) &&
	   (face[i].vert[0].z <= bvmax[j].z)){
	  bvh_id[j].push_back(i);
	  used = true;
	  continue;
	}
	if((face[i].vert[1].x >= bvmin[j].x) &&
	   (face[i].vert[1].y >= bvmin[j].y) &&
	   (face[i].vert[1].z >= bvmin[j].z) &&
	   (face[i].vert[1].x <= bvmax[j].x) &&
	   (face[i].vert[1].y <= bvmax[j].y) &&
	   (face[i].vert[1].z <= bvmax[j].z)){
	  bvh_id[j].push_back(i);
	  used = true;
	  continue;
	}
	if((face[i].vert[2].x >= bvmin[j].x) &&
	   (face[i].vert[2].y >= bvmin[j].y) &&
	   (face[i].vert[2].z >= bvmin[j].z) &&
	   (face[i].vert[2].x <= bvmax[j].x) &&
	   (face[i].vert[2].y <= bvmax[j].y) &&
	   (face[i].vert[2].z <= bvmax[j].z)){
	  used = true;
	  bvh_id[j].push_back(i);
	  continue;
	}
      }
    }
    
  }

  cout << "File read: triangle data stored and sorted";
  if(bvh == true)
    cout << " | Bounding Volume Hierarchy added. ";
  cout << endl;
  cout << "min point ( " << xmin << ", " << ymin << ", " << zmin << " )" << endl;
  cout << "max point ( " << xmax << ", " << ymax << ", " << zmax << " )" << endl; //data stored from mesh

  //-------------------Electron Beam------------------------
  //BC
  double u = 0; // a possbile intersection P = O + u* D
  double max[3] = {xmax+eps, ymax+eps, zmax+eps};
  double min[3] = {xmin-eps, ymin-eps, zmin-eps};
  

  uniform_real_distribution<double> unifrndx(xmin, xmax);//-beam randomizer dimensions
  uniform_real_distribution<double> unifrndy(ymin, ymax);//-
  uniform_real_distribution<double> unifrndz(zmin, zmax);//-----
  
  int see = 0; //count of secondary particle emitted
  int thru = 0; // daughters that made it through
  int init = 0; // initial rays that missed;
  int err = 0;
  double deg2rad = M_PI/180.0;
  double MAXYIELD = 10;
  int trapped = 0; //error for rays bouncing from wall to wall.
  double angle = 0, angle1, angle2, psi;
  double r1, r2; // random numbers
  double MAXDIST = 5;//check file for scale of mesh (normalize first to avoid issues)
  cout << "cutoff energy = " << cutoff << " eV" << endl;
  outfile.open("outputs.dat"); //tracks all rays created (for debugging)
  if(print == true){
    outfile1.open("ray_data.dat");  //for visualization. The start and end points of rays
    outfile2.open("distribution.dat"); //Energy and angle distributions of initial and outgoing rays
  }
  outfile3.open("angle_dis.dat"); // angle of incidence
  outfile4.open("penetration.dat"); // rays that do not generate more rays (zero yield or < Ecut energy).
  outfile5.open("outangles.dat"); //outgoing energy and angle of only escaped rays wrt z-axis
  //-------------------MC part---------------
  for(int steps = 0; steps < Nsteps; steps++){
    double Energy = ParentE, Energy1, Energy2 = 0; //eV
    int gencount = 0;
    int found = 0;
    int found1 = 0;
    int hold = 0;
    int count = 0;
    double numdaught = 0;
    double randx = unifrndx(generator);
    double randy = unifrndy(generator);
    double randz = unifrndz(generator);
    //-------------------------randomize ray origin and direction----------------
    Vector POT; //Point on Triangle, the intersection
    rays.Origin = Vector(randx, randy, eps+zmax); //initially just above the mesh in random location
    rays.Direction = Vector(0.0, 0.0, -1.0); //0 degree incidence w.r.t z-axis
    rays.spec = species[0];
    if(random == true){
      angle = asin(sqrt(unifrnd(generator)));
      psi = 2*M_PI*unifrnd(generator);
      rays.Direction = (cos(angle)*Vector(0.0, 0.0, -1.0)+sin(angle)*cos(psi)*Vector(0.0, 1.0, 0.0)+sin(angle)*sin(psi)*Vector(1.0, 0.0, 0.0)).UnitVector();
    }
    if(print == true)
      outfile2 << "\t" << Energy << "\t" << angle/deg2rad << "\t"; //initial outgoing energy and outgoing angle
    //-------------------------Initial intersection search----------------------
    if(bvh == true)
      found1 = ray_box_intersect(rays, face, POT, bvmax, bvmin, boxes, bvh_id, MAXDIST, hold, found);
    else
      search(tn, rays, face, found, found1, hold, POT, MAXDIST);
    trapped = 0;
    while(found1 != 1){//initial beam did not make an intersection begin B.C.
      trapped++; //check if beam is stuck somewhere
      if(trapped > 20){
	rays.Direction.z = rays.Direction.z + -0.1;
	if(trapped > 30){
	  cout << "error 1" << endl;
	  err++;
	  break;
	}
      }
      if(print == true)
	outfile1 << fixed << setprecision(5) << rays.Origin.x << "\t" << rays.Origin.y << "\t" << rays.Origin.z << "\t";
      bc(rays, u, max, min, found);
      if(print == true)      
	outfile1 << setprecision(5) << POT.x << "\t" << POT.y << "\t" << POT.z  << "\t" << "0 \t" << gencount << endl; //output initial rays that went through for visual purposes
      if(u == 0.0){
	cout << "error 2" << endl;
	err++;
	break;
      }
      POT = rays.Origin + u * rays.Direction;
      //B.C, ray is bounced back
      bounce(POT, rays, found);
      if(found == 7){//hits bottom on mesh
	init++;
	if(substrate == true){
	  found1 = 1;
	  hold = -1;
	  found = 1;
	}
	break;
      }
      //check if there is a ray-mesh collision - use M-T algorithm
      if(hold != -1)
	if(bvh == true)
	  found1 = ray_box_intersect(rays, face, POT, bvmax, bvmin, boxes, bvh_id, MAXDIST, hold, found);
	else
	  search(tn, rays, face, found, found1, hold, POT, MAXDIST);
    }
    if(found1 == 1){ //M-T Made an intersection
      Vector Normal;
      //------------------------angle of incidence-----------------------------
      if(hold == -1)
	Normal = Vector(0.0,0.0,1.0);
      else
	Normal = face[hold].norm; 
      angle = acos(DotProduct(-rays.Direction, Normal)); //angle of incidence
      if(angle >= M_PI/2){
	angle = angle - M_PI/2;
	if(angle >= M_PI/2){
	  angle = angle - M_PI/2;
	}
      }
      outfile3 << rays.spec << "\t" << Energy << "\t" <<  angle/deg2rad << endl;
      //-------------------------Initialize list--------------------------------
      node *root = newNode(Energy, angle, rays, gencount);
      node *look;
      look = root;
      //------------------------write out ------------------------------------
      /*if(hold == -1){
	outfile << "Hit bottom of mesh \t";
	outfile << setprecision(6) << POT.x << "\t" << POT.y << "\t" << POT.z << "\t";
      }
      else{
	outfile << hold << "\t"<< setprecision(6) << face[hold].vert[0].x << "\t" <<  face[hold].vert[0].y << "\t" <<  face[hold].vert[0].z << "\t"; 
	outfile << setprecision(6) <<  face[hold].vert[1].x << "\t" <<  face[hold].vert[1].y << "\t" <<  face[hold].vert[1].z << "\t"; 
	outfile << setprecision(6) <<  face[hold].vert[2].x << "\t" <<  face[hold].vert[2].y << "\t" <<  face[hold].vert[2].z << "\t"; 
	outfile << setprecision(6) << Normal.x << "\t" << Normal.y << "\t" << Normal.z << "\t";
	outfile << setprecision(6) << POT.x << "\t" << POT.y << "\t" << POT.z << "\t";
	}*/
      //-----------------------ray data --------------------------------------
      if(print == true){
	outfile1 << rays.Origin.x << "\t" << rays.Origin.y << "\t" << rays.Origin.z << "\t";
	outfile1 << POT.x << "\t" << POT.y << "\t" << POT.z  << "\t";
      }
      if(rays.Origin.x == POT.x && rays.Origin.y == POT.y && rays.Origin.z == POT.z){
	cout << "exiting... status 0" << endl;
	break;
      }
      //-----------------------Daughter generation---------------------------
      //draw from tables given energy and incident angle.
      //new energy calculated from old rays energy and angle of incidence
      
      source_f_yield(look->Energy, angle, numdaught, look->rays.spec, species[1]); //electron yield
      if(numdaught > MAXYIELD)
	numdaught = MAXYIELD;
      source_f_energy(look->Energy, angle, Energy, look->rays.spec, species[1]); // Emitted particle energy
      
      r1 = unifrnd(generator);
      if(Energy > cutoff){
	gencount = look->Gen;
	gencount++;
	rays.Origin = POT;
	double daughtholder = numdaught, daughtint, daughtfrac;
	Energy2 = 0;
	daughtfrac = modf(daughtholder, &daughtint);
	if(daughtfrac > r1)
	  daughtint += 1;
	count = count + int(daughtint);
	if(look->rays.spec == species[0] && hold != -1){
	  face[hold].in += 0;
	  face[hold].out += int(daughtint);
	}
	if(look->rays.spec == species[1] && hold != -1){
	  face[hold].in += 1;
	  face[hold].out += int(daughtint);
	}
	if(print == true)
	  outfile1 << daughtint << "\t" << look->Gen <<  endl; //daughter generated and generation	    
	outfile << int(daughtint) << "\t";
	if(daughtint == 0.0) // zero yield
	  outfile4 << look->Energy << "\t" << POT.z << endl;  //cout << "no daughters generated" << endl;
	double Eholders[int(daughtint)];
	for(int i = 0; i < int(daughtint); i++){
	  r2 = unifrnd(generator);	    
	  Eholders[i] = r2;	
	  Energy2 = Eholders[i] + Energy2;
	  if(int(daughtint) == 1)
	      Energy2 = 1.0;  
	}
	for(int i = 0; i < int(daughtint); i++){
	  r2 = Eholders[i] / Energy2;	    
	  angle1 = angle_distribution(distribution, rays, face, look->Energy, angle, min, max, Normal, hold);
	  Energy1 = Energy*r2;
	  if(int(daughtint) == 1)
	    Energy1 = Energy;
	  insertKey(look, Energy1, angle1, rays, gencount);
	}
      }
      else{//energy was lower than cutoff no generation made
	outfile4 << Energy << "\t" << POT.z << endl; // no energy to produce more
	if(print == true)	
	  outfile1 << "0 \t" << look->Gen << endl;
      }
      //---------------while loop for generations after initial----------------
      while(look->next != NULL){ //while the next in line of the chain is not empty
	look = look->next;
	Energy = look->Energy;
	found = 0; //set to not found
	found1 = 0; //set to not found
	hold = 0;
	rays.Origin = look->rays.Origin;
	rays.Direction = look->rays.Direction;
	rays.spec = species[1];
	//----------------------next node calculation M-T--------------------------
	if(bvh == true)
	  found1 = ray_box_intersect(rays, face, POT, bvmax, bvmin, boxes, bvh_id, MAXDIST, hold, found);
	else
	  search(tn, rays, face, found, found1, hold, POT, MAXDIST);	
	trapped = 0;
	while(found1 != 1){//an intersection is not found... ray must go somewhere...
	  look->rays = rays;
	  trapped++;
	  if(trapped > 20){
	    if(rays.Direction.z >= eps)
	      rays.Direction.z = rays.Direction.z + 0.1;
	    else
	      rays.Direction.z = rays.Direction.z -0.1;
	    //cout << "error 3" << endl;
	    if(trapped > 30){
	      cout << "error 1.5" << endl;
	      err++;
	      break;
	    }
	  }
	  if(look->rays.Direction.z >= eps){ //<--------------change for direction, for beam starting above mesh do > 0 | below mesh z < 0	    
	    if(print == true)	    
	      outfile1 << setprecision(5) << look->rays.Origin.x << "\t" << look->rays.Origin.y << "\t" << look->rays.Origin.z  << "\t";
	    //---------------chance of escapes-------------------------------   see++; 
	    bc(rays, u, max, min, found);
	    if(found == 2){
	      see++;
	      POT = look->rays.Origin + u * look->rays.Direction;
	      double anglehold;
	      anglehold = acos(DotProduct(look->rays.Direction, Vector(0.0, 0.0, 1.0)));
	      outfile5 << look->Energy << "\t" << (anglehold)/deg2rad << endl;
	      if(print == true)	      
		outfile1 << POT.x << "\t" << POT.y << "\t" << POT.z  << "\t" << "0 \t" << look->Gen << endl; //output daughters that escaped
	      break;
	    }
	    POT = look->rays.Origin + u * look->rays.Direction;
	    if(print == true)	    
	      outfile1 << setprecision(5) << POT.x << "\t" << POT.y << "\t" << POT.z  << "\t" << "0 \t" << look->Gen << endl; //for visual purpose only
	    if(POT.x == look->rays.Origin.x && POT.y == look->rays.Origin.y && POT.z == look->rays.Origin.z){//escaped because it was at max value.
	      see++;
	      double anglehold;
	      anglehold = acos(DotProduct(look->rays.Direction, Vector(0.0, 0.0, 1.0)));
	      outfile5 << look->Energy << "\t" << (look->Angle)/deg2rad << endl;
	      break;
	    }
	    //B.C
	    bounce(POT, rays, found);
	    if(found == 7){//hits bottom on mesh
	      thru++;
	      if(substrate == true){
		found1 = 1;
		hold = -1;
		found = 1;
	      }
	      break;
	    }
	    //M-T Algorithm after correction from boundary condition
	    if(hold != -1)
	      if(bvh == true)
		found1 = ray_box_intersect(rays, face, POT, bvmax, bvmin, boxes, bvh_id, MAXDIST, hold, found);
	      else
		search(tn, rays, face, found, found1, hold, POT, MAXDIST);
	  }
	  else{ //direction is down so intersection must occur begin B.C. // thru++;
	    if(print == true)	    
	      outfile1 << setprecision(5) << look->rays.Origin.x << "\t" << look->rays.Origin.y << "\t" << look->rays.Origin.z  << "\t";
	    bc(rays, u, max, min, found);
	    POT = look->rays.Origin + u * look->rays.Direction;
	    if(print == true)	    
	      outfile1 << setprecision(5) << POT.x << "\t" << POT.y << "\t" << POT.z  << "\t" << "0 \t" << look->Gen << endl; //output daughters that went through
	    if(POT.x == look->rays.Origin.x && POT.y == look->rays.Origin.y && POT.z == look->rays.Origin.z){ //errors?
	      cout << "error 4" << endl;
	      err++;
	      break;
	    }
	    if(u == 0){ //causes POT = rays.Origin
	      err++;
	      cout << "error 5" << endl;
	      break;
	    }
	    //B.C.
	    bounce(POT, rays, found);
	    if(found == 7){ //hit the bottom
	      thru++;
	      if(substrate == true){
		found1 = 1;
		hold = -1;
		found = 1;
	      }
	      break;
	    }
	    //M-T Algorithm
	    if(hold != -1)
	      if(bvh == true)
		found1 = ray_box_intersect(rays, face, POT, bvmax, bvmin, boxes, bvh_id, MAXDIST, hold, found);
	      else
		search(tn, rays, face, found, found1, hold, POT, MAXDIST);	    
	    
	  }
	  if(POT.z >= max[2]+eps){//somehow out of bounds, check B.C if so
	    err++;
	    cout << "error 6 up" << endl;
	    break;
	  }
	  if(POT.z <= min[2]-eps){ //somehow out of bounds, check B.C if so
	    err++;
	    cout << "error 6 down" << endl;
	    break;
	  }
	}
	// if intersection made, generation begins
	if(found1 == 1){
	  Vector Normal;
	  if(hold == -1)
	    Normal = Vector(0.0, 0.0, 1.0);
	  else
	    Normal = face[hold].norm;
	  angle = acos(DotProduct(-rays.Direction, Normal)); //angle of incidence
	  if(angle >= M_PI/2){
	    angle = angle - M_PI/2;
	    if(angle >= M_PI/2){	
	      angle = angle - M_PI/2;
	    }
	  }
	  outfile3 << rays.spec << " " << Energy << " " <<  angle/deg2rad << endl;
	  //----------------write out------------------------
	  /*if(hold == -1){
	    outfile << "Hit bottom of mesh \t";
	    outfile << setprecision(6) << POT.x << "\t" << POT.y << "\t" << POT.z << "\t";
	  }
	  else{
	    outfile << hold << "\t" << setprecision(6) <<  face[hold].vert[0].x << "\t" <<  face[hold].vert[0].y << "\t" <<  face[hold].vert[0].z << "\t";
	    outfile << setprecision(6) << face[hold].vert[1].x << "\t" << face[hold].vert[1].y << "\t" <<  face[hold].vert[1].z << "\t";
	    outfile << setprecision(6) << face[hold].vert[2].x << "\t" << face[hold].vert[2].y << "\t" << face[hold].vert[2].z << "\t";
	    outfile << setprecision(6) << Normal.x << "\t" << Normal.y << "\t" << Normal.z << "\t";
	    outfile << setprecision(6) << POT.x << "\t" << POT.y << "\t" << POT.z << "\t";
	    }*/
	  //-----------------ray data-----------------------
	  if(print == true){	  
	    outfile1 << setprecision(5) << look->rays.Origin.x << "\t" << look->rays.Origin.y << "\t" << look->rays.Origin.z  << "\t";
	    outfile1 << setprecision(5) << POT.x << "\t" << POT.y << "\t" << POT.z  << "\t";
	  }	  
	  if(look->rays.Origin.x == POT.x && look->rays.Origin.y == POT.y && look->rays.Origin.z == POT.z)
	    break;
	  //if incoming ray is species[0] use first set of equations
	  //if incoming ray is species[1] use second set of equations

	  source_f_yield(look->Energy, angle, numdaught, look->rays.spec, species[1]);
	  if(numdaught > MAXYIELD)
	    numdaught = MAXYIELD;
	  source_f_energy(look->Energy, angle, Energy, look->rays.spec, species[1]);	  

	  r1 = unifrnd(generator);
	  if(Energy > ParentE){
	    cout << "Daughter energy is greater than Parent!" << steps << "\t" << count << endl;
	    break;
	  }
	  //----------------grandaughter+ generation----------------------
	  if(Energy > cutoff){
	    gencount = look->Gen;
	    gencount++;
	    rays.Origin = POT;
	    double daughtholder = numdaught;
	    double daughtint;
	    double daughtfrac;
	    Energy2 = 0;
	    daughtfrac = modf(daughtholder, &daughtint);
	    if(daughtfrac > r1 )
	      daughtint += 1;
	    count = count + int(daughtint);
	    if(look->rays.spec == species[0] && hold != -1){
	      face[hold].in += 0;
	      face[hold].out += int(daughtint);
	    }
	    if(look->rays.spec == species[1] && hold != -1){
	      face[hold].in += 1;
	      face[hold].out += int(daughtint);
	    }
	    if(print == true)
	      outfile1 << daughtint << "\t" << look->Gen <<  endl;	    
	    outfile << int(daughtint) << "\t";  
	    if(daughtint == 0.0) //zero yields
	      outfile4 << look->Energy << "  " << POT.z << endl;  //cout << "no daughters generated"  
	    double Eholders[int(daughtint)];
	    for(int i = 0; i < int(daughtint); i++){
	      r2 = unifrnd(generator);
	      Eholders[i] = r2;
	      Energy2 = Eholders[i] + Energy2;
	      if(int(daughtint) == 1)
		Energy2 = 1.0;
	    }
	    for(int i = 0; i < int(daughtint); i++){
	      r2 = Eholders[i] / Energy2;
	      Energy1 = Energy*r2;
	      angle1 = angle_distribution(distribution, rays, face, look->Energy, angle, min, max, Normal, hold);
	      if(int(daughtint) == 1)
		Energy1 = Energy;     
	      //outfile5 << Energy1 << "\t" << angle1/deg2rad << endl;
	      insertKey(look, Energy1, angle1, rays, gencount);   
	    }
	  }
	  else{//no daughters generated not enough energy
	    outfile4 << look->Energy << "\t" <<  POT.z << endl;
	    if(print == true)	    
	      outfile1 << "0 \t" << look->Gen << endl;
	  }
	}
      }//close daughter generation while loop
      outfile << endl;
      clear(root);   
    }// close initial found intersection if statement
    //--------- there must be an intersection!--------------
  }
  outfile.close();
  if(print == true){
    outfile1.close();
    outfile2.close();
  }
  outfile3.close();
  outfile4.close();
  outfile.open("in_out.dat");
  for(int i = 0 ; i < tn; i++)
    outfile << i << "\t" << face[i].in << "\t" << face[i].out << endl;
  outfile.close();
  cout << Nsteps << " # of rays shot at target." << endl;
  cout << see << " # of secondary electrons emitted." << endl;
  cout << see/double(Nsteps) << " yield" << endl;
  cout << init << " initial rays that hit bottom of material" << endl;
  cout << thru << " # of emitted particles that hit bottom of material | " << err << endl;
  cout << "done in ";
  t = clock() - t;
  cout << float(t)/CLOCKS_PER_SEC << " seconds" << endl;

  delete [] face;
  return 0;

}

