#include<fstream>
#include<iostream>
#include<vector>
#include<cmath>
#include<string>
#include"func.h"
using namespace std;

int main( int argc, char **argv){
  ofstream outfile;
  ifstream infile;
  string filename, temp1,temp2, first, last;
  int n = 0;
  int bin = 25;
  int ybin = 1;
  char buffer[30];
  vector<elements> face;
  double max = 0.5, min = -0.5;
  if(argc > 1)
    for(int i = 1; i < argc; i++){
      grabber(argv[i], first, last, '=');
      if(first == "filename")
	 filename = last;
      if(first == "bin")
	bin = stoi(last);
      if(first == "ybin")
	ybin = stoi(last);
      if(first == "max")
	max = stod(last);
      if(first == "min")
	min = stod(last);
    }
  snprintf(buffer, sizeof(char)*32, filename.c_str());
  infile.open(buffer);
  if(infile.good() == false){
    cout << buffer << " not found" << endl;
    return 1;
  }
  double xmax=0, ymax=0, zmax=0;
  double xmin=0, ymin=0, zmin=0;
  while(!infile.eof()){
    infile >> temp1;
    if( temp1 == "facet"){
      infile >> temp2;
      face.push_back(elements());
      infile >> face[n].norm.x >> face[n].norm.y >> face[n].norm.z;
      if(fabs(face[n].norm.x) < eps) face[n].norm.x = 0; //making sure zeros are zer0
      if(fabs(face[n].norm.y) < eps) face[n].norm.y = 0;
      if(fabs(face[n].norm.z) < eps) face[n].norm.z = 0;
    }
    if(temp1 == "loop"){
      for(int i = 0; i < 3; i++){
	infile >> temp2 >> face[n].vert[i].x >> face[n].vert[i].y >> face[n].vert[i].z;
	getline(infile,temp1);
	face[n].center.x += face[n].vert[i].x;
	face[n].center.y += face[n].vert[i].y;
	face[n].center.z += face[n].vert[i].z; 
	if(face[n].vert[i].x > xmax)	  xmax = face[n].vert[i].x;
	if(face[n].vert[i].x < xmin)	  xmin = face[n].vert[i].x;
	
	if(face[n].vert[i].y > ymax)	  ymax = face[n].vert[i].y;
	if(face[n].vert[i].y < ymin)	  ymin = face[n].vert[i].y;
	
	if(face[n].vert[i].z > zmax)	  zmax = face[n].vert[i].z;
	if(face[n].vert[i].z < zmin)	  zmin = face[n].vert[i].z;
      }
      face[n].center.x /= 3.0;
      face[n].center.y /= 3.0;
      face[n].center.z /= 3.0;
      n++; 
    }
  }
  cout << "min point ( " << xmin << ", " << ymin << ", " << zmin << " )" << endl;
  cout << "max point ( " << xmax << ", " << ymax << ", " << zmax << " )" << endl;
  cout << n << endl;
  infile.close();
  Vector Origin, Origin2;
  Vector Direction, Direction2;
  Vector POT;
  double inpoint;
  int found = 0;
  double lengthz = zmax - zmin - eps;
  double lengthy = ymax - ymin - eps;
  double lengthx = xmax - xmin - eps;
  double binstep = (max-min)/bin;
  double bins[bin+1];
  double value[bin];
  for(int i =0; i < bin; i++)
    value[i] = 0;
  for(int i = 0; i <= bin; i++){
    bins[i] = max-binstep*i;
  }

  int ceps = 0;
  for(int x = 0; x < bin; x++){
    vector<Vector> POTS;
    for(int y = 0; y < ybin; y++){
      Origin = Vector(2*xmin-eps, ymax-eps-y*lengthy/ybin, bins[x]-0.5*binstep);
      //Origin = Vector(0, 0, bins[x]-binstep*0.5);
      Direction = Vector(1.0, 0.0, 0.0).UnitVector();
      int i = 0;
      while(i < n){
	found = 0;
	if(POTS.size() %2 == 0)
	  intersect_triangle(Origin, Direction, face[i].vert[0], face[i].vert[1], face[i].vert[2], inpoint, found);
	else
	  intersect_triangle(Origin, Direction, face[i].vert[0], face[i].vert[2], face[i].vert[1], inpoint, found);
	if(found == 1){
	  POT = Origin + inpoint*Direction;
	  POTS.push_back(POT);
	    i = 0;
	    Origin = POT;
	  }
	i++;
      }
    }
    ceps = 0;  
    for(int i = 0; i < POTS.size(); i=i+2){
      Vector subtr = POTS[i+1] - POTS[i];
      if(subtr.Magnitude() < 0.01)
	ceps++;
      else
	value[x] += subtr.Magnitude();
    }
    cout << "ceps " << ceps << endl;
    value[x] = value[x]/(0.5*(double(POTS.size())-eps));
    cout << x << " " << POTS.size() << " " << value[x] <<  endl;
  }
  
    outfile.open("ligsize.dat");
  for(int i = 0; i < bin; i++)
    outfile << bins[i] << "\t" << value[i] << endl; 
  
return 0;
}
