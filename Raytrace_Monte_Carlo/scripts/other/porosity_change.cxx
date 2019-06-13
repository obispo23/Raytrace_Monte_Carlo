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
  double scale = 0.0001;
  if(argc > 1){
    for(int i = 1; i < argc; i++){
      input = argv[i];
      grabber(input, first, last, '=');
      if(first == "filename")
        filename = last;
      if(first == "scale")
	scale = stod(last);
    }
  }
  else{
    cout << "arguements: filename and scale" << endl;
    return 0;
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
  cout << "Triangle Average: " << average << endl;
  double length  = zmax - zmin;
  double volume = 0.0;
  /*  for(int p = 0; p < 100; p++)
    for(int i = 0; i < tn; i++)
      if((face[i].vert[0].z > zmin+length*p) && (face[i].vert[0].z < zmin+length*(p+1))){
	double angle;
	Vector AB = face[i].vert[0]-face[i].vert[1];
	Vector AC = face[i].vert[0]-face[i].vert[2];
	double dotprod = DotProduct(AB,AC);
	angle = acos(dotprod / (abs(AB.Magnitude())*abs(AC.Magnitude())));
	average = average + 0.5 * (abs(AB.Magnitude())*abs(AC.Magnitude()))*sin(angle);
	}*/
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
  cout << volume << " " << volumetot << " " << volume/volumetot << endl;
 

  //while(input != 'q')
  
  for(int i = 0; i < tn; i++)
    for(int j = 0; j < 3; j++)
      original[i].vert[j] = face[i].vert[j];
  for(int i = 0 ; i < tn; i++){
    for(int j = 0; j < 3; j++){
      face[i].vert[j] = face[i].vert[j] +scale*face[i].norm.UnitVector();
    }
  }
  for(int i = 0; i < tn; i++)
    for(int p = 0; p < 3; p++)
      if(original[i].used[p] == false){
	vector<int> vlist;
	vector<int> vindex;
	vlist.push_back(i);
	vindex.push_back(p);
	original[i].used[p] = true;
	for(int j = 0; j < tn; j++)
	  for(int q = 0; q < 3; q++)
	    if(original[j].used[q] == false)
	      if((original[j].vert[q].x == original[i].vert[p].x) && (original[j].vert[q].y == original[i].vert[p].y) && (original[j].vert[q].z == original[i].vert[p].z)){
		original[j].used[q] = true;
		vlist.push_back(j);
		vindex.push_back(q);
	      }
	Vector nvert = {};
	for(int k = 0; k < vlist.size(); k++)
	  nvert = nvert + face[vlist[k]].vert[vindex[k]];
	nvert = 1.0/double(vlist.size())*nvert;
	for(int k = 0; k < vlist.size(); k++){
	  face[vlist[k]].vert[vindex[k]] = nvert;
	}
      }
  for(int i = 0; i < tn; i++)
    for(int j = 0; j < 3; j++)
      if(face[i].used[j] == true){
        Vector line1 = face[i].vert[1] - face[i].vert[0];
        Vector line2 = face[i].vert[2] - face[i].vert[0];
        face[i].norm = CrossProduct(line2, line1).UnitVector();
      }
  double nvol = 0.0;
  for(int i = 0; i < tn; i++){
    double v0 = -(face[i].vert[2].x * face[i].vert[1].y * face[i].vert[0].z);
    double v1 = (face[i].vert[1].x * face[i].vert[2].y * face[i].vert[0].z);
    double v2 = (face[i].vert[2].x * face[i].vert[0].y * face[i].vert[1].z);
    double v3 = -(face[i].vert[0].x * face[i].vert[2].y * face[i].vert[1].z);
    double v4 = -(face[i].vert[1].x * face[i].vert[0].y * face[i].vert[2].z);
    double v5 = (face[i].vert[0].x * face[i].vert[1].y * face[i].vert[2].z);
    nvol = nvol + 1.0/6.0*(v0 + v1 + v2 + v3 + v4 + v5); //signed volume
  }
  cout << "changed to : " << nvol/volumetot << endl;
  cout << " Writing output to newmesh.stl " << endl;
  //  snprintf(buffer, sizeof(char)*32, filename.c_str());
  outfile.open("newmesh.stl");
  outfile << "#output by Andrew Alvarado's code remesh" << endl;
  outfile << fixed << setprecision(6);
  for(int i = 0; i < tn; i++){
    int found = face[i].in;
    if(found == 0){
      outfile << setprecision(6) << "facet normal \t" << face[i].norm.x << "\t" <<  face[i].norm.y << "\t" << face[i].norm.z << endl;
      outfile << "outer loop" << endl;
      for(int j = 0; j < 3; j++)
        outfile << setprecision(6) << "vertex \t" << face[i].vert[j].x << "\t" << face[i].vert[j].y << "\t" << face[i].vert[j].z << endl;
      outfile << "endloop" << endl << "endfacet" << endl;
    }
  }
  outfile.close();
  delete [] face;
  delete [] original;
  return 0;
}
