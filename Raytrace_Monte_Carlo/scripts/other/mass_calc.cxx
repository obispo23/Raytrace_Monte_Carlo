#include<fstream>
#include<string>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<random>
#include<ctime>
#include"vec.h"

using namespace std;

//#define eps 0.00001

int main(int argc, char **argv)
{
  clock_t t;
  t = clock();
  ofstream outfile, outfile1, outfile2, outfile3, outfile4, outfile5;
  ifstream infile;
  string temp1, temp2, temp3, temp4;
  int tn = 0;
  int lines = 0;
  char buffer[32];
  string input, first, last;
  string filename = "flatmesh.stl";
  if(argc > 1){
    for(int i = 1; i < argc; i++){
      input = argv[i];
      grabber(input, first, last, '=');
      if(first == "filename")
        filename = last;
    }
  }
  snprintf(buffer, sizeof(char)*32, filename.c_str());
  cout << " file used: " << buffer << endl;
 //----------------------------mesh input--------------------------------------               
  infile.open(buffer);
  if(infile.good() == false){
    cout << buffer << " not found " << endl;
    return 1;
  }
  //count number of triangles and lines in file                                                
  while(!infile.eof()){
    infile >> temp1 >> in >> out;
    if(temp1 == "facet"){
      tn++;
    }
    getline(infile,temp2);
    lines++;
  }
  infile.close();
  cout << tn << " triangles in " << lines << " lines" << endl;
  infile.open(buffer);
  elements *face;
  face = new elements[tn];
  elements *original;
  original = new elements[tn];
  for(int i =0; i < tn; i++){
    face[i] = {};
    original[i] = {};
    for(int j = 3; j < 3; j++)
      original[i].used[j] = false;
  }
  int n=0; //counter for facets                 
  double xmax = 0, ymax = 0,  zmax = 0; //solid object size dimensions
  double xmin = 0, ymin = 0,  zmin = 0; //for cwall:width, length, and height
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
  cout << "File read: triangle data stored and sorted" << endl;
  cout << "min point:" << xmin << "\t" << ymin << "\t" << zmin << endl;
  cout << "max point:" << xmax << "\t" << ymax << "\t" << zmax << endl; //data stored from mesh
  double average = 0;

  for(int i = 0 ; i < tn; i++){
    double angle;
    Vector AB = face[i].vert[0]-face[i].vert[1];
    Vector AC = face[i].vert[0]-face[i].vert[2];
    double dotprod = DotProduct(AB,AC);
    angle = acos(dotprod / (abs(AB.Magnitude())*abs(AC.Magnitude())));
    average = average + 0.5 * (abs(AB.Magnitude())*abs(AC.Magnitude()))*sin(angle);
  }
  average = average/tn;
  cout << "Triangle Surface Area Average: " << average << endl;
  double length  = zmax - zmin;
  double volume = 0.0;
  
  for(int i = 0; i < tn; i++){
    double v0 = -(face[i].vert[2].x * face[i].vert[1].y * face[i].vert[0].z);
    double v1 = (face[i].vert[1].x * face[i].vert[2].y * face[i].vert[0].z);
    double v2 = (face[i].vert[2].x * face[i].vert[0].y * face[i].vert[1].z);
    double v3 = -(face[i].vert[0].x * face[i].vert[2].y * face[i].vert[1].z);
    double v4 = -(face[i].vert[1].x * face[i].vert[0].y * face[i].vert[2].z);
    double v5 = (face[i].vert[0].x * face[i].vert[1].y * face[i].vert[2].z);
    volume = volume + 1.0/6.0*(v0 + v1 + v2 + v3 + v4 + v5); //signed volume
  }
  double volumetot = (xmax-xmin)*(ymax-ymin)*(zmax-zmin);
  cout << "Mesh Volume: " << volume << ". Total Volume:" << volumetot << ". Volume Fraction" << volume/volumetot << endl;

  infile.open("in_out.dat");
  for(int i = 0; i < tn; i++)
    infile >> temp1
    return 0;
}
