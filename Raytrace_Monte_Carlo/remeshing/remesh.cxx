#include<fstream>
#include<string>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<random>
#include<ctime>
#include"vec_remesh.h"
#include"tri_tri_inter.h"
//#include"opttritri.h"
using namespace std;

#define eps 0.000001
#define eps1 0.001
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
  string input, first, last;
  string filename = "flatmesh.stl";
  double atom_volume=9.53E3/6.02E23; //(mm)^3
  double flux=1E15; //(#/mm^2/s)
  double Nsteps=100;
  bool smooth = false;
  bool del = false;
  bool island = false;
  bool remesh = false;
  bool experimental = false;
  int islecon = 0;
  double scale = 1;
  if(argc > 1){
    for(int i = 1; i < argc; i++){
      input = argv[i];
      grabber(input, first, last, '=');
      if(first == "filename")
        filename = last;
      if(first == "Nsteps")
	Nsteps = stod(last);
      if(first == "atom_volume")
	atom_volume = stod(last);
      if(first == "flux")
	flux = stod(last); 
      if(first == "scale")
	scale = stod(last);
      if(first == "smooth"){
	smooth = true;
	if(last == "no")
	  smooth = false;
      }
      if(first == "del"){
	del = true;
	if(last == "no")
	  del = false;
      }
      if(first == "island"){
	island = true;
	if(last == "no")
	  island = false;
	else
	  islecon = stoi(last);
      }
      if(first == "remesh"){
	remesh = true;
	if(last == "no")
	  remesh = false;
      }
      if(first == "experimental"){
	experimental = true;
	if(last == "no")
	  experimental = false;
      }
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
  for(int i = 0; i < tn; i++)
    for(int j = 0; j < 3; j++)
      original[i].vert[j] = face[i].vert[j];
  cout << "File read: triangle data stored and sorted" << endl;
  cout << "min point:" << xmin << "\t" << ymin << "\t" << zmin << endl;
  cout << "max point:" << xmax << "\t" << ymax << "\t" << zmax << endl; //data stored from mesh
  infile.open("in_out.dat");
  double Area = abs(xmax-xmin)*abs(ymax-ymin);
  double timestep = Nsteps/(flux*Area)*scale;
  cout << " Timestep: " << timestep << endl;
  int id;
  int in;
  int out;
  int intot = 0;
  int outtot = 0;
  //  double deltah = 0;
  Vector line1;
  Vector line2;
  double angle;
  double a_i;
  for(int i = 0; i < tn; i++){
    infile >> id >> in >> out;
    intot += in;
    outtot += out;
    getline(infile, temp1);
    if((in - out) != 0 ){
      line1 = face[id].vert[1]-face[id].vert[0];
      line2 = face[id].vert[2]-face[id].vert[0];
      angle = acos(DotProduct(line1,line2)/abs(line1.Magnitude()*line2.Magnitude()));
      a_i = 0.5*line1.Magnitude()*line2.Magnitude()*sin(angle);
      face[id].deltah=(double(in - out)/abs(a_i))*(atom_volume)*scale;
    }
    else
      face[id].deltah = 0;
    for(int j = 0; j < 3; j++)
      face[id].vert[j] = face[id].vert[j] + face[id].deltah*face[id].norm.UnitVector();
  }
  infile.close();
  cout << "total deposit: " << intot << endl << "total sputter: " << outtot << endl;
  cout << "total mass change: " << (outtot-intot) << endl;
  if(smooth == true){
    cout << "smoothing..." << endl;
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
  }
  
  for(int i = 0; i < tn; i++)
    for(int j = 0; j < 3; j++)
      if(face[i].used[j] == true){
	line1 = face[i].vert[1] - face[i].vert[0];
	line2 = face[i].vert[2] - face[i].vert[0];
	face[i].norm = CrossProduct(line2, line1).UnitVector();
      }

  for(int i = 0; i < tn; i++){
    face[i].in = 0; //if 1 delete
    face[i].out = 0; // if 1 connected to faces > islecon
  }


  if(del == true){
    cout << "deleting... ";
    int intcount = 0;
    for(int i = 0; i < tn; i++)
      for(int j = 0; j < tn; j++)
	if(i != j){
	  int found = 0;
	  float V0[3] = {float(face[i].vert[0].x), float(face[i].vert[0].y), float(face[i].vert[0].z)};
	  float V1[3] = {float(face[i].vert[1].x), float(face[i].vert[1].y), float(face[i].vert[1].z)};
	  float V2[3] = {float(face[i].vert[2].x), float(face[i].vert[2].y), float(face[i].vert[2].z)};
	  float U0[3] = {float(face[j].vert[0].x), float(face[j].vert[0].y), float(face[j].vert[0].z)};
	  float U1[3] = {float(face[j].vert[1].x), float(face[j].vert[1].y), float(face[j].vert[1].z)};
	  float U2[3] = {float(face[j].vert[2].x), float(face[j].vert[2].y), float(face[j].vert[2].z)};
	  //found = NoDivTriTriIsect(V0, V1, V2, U0, U1, U2);
	  int coplanar=0;
	  float source[3];
	  float target[3];
	  found =  tri_tri_intersection_test_3d(V0, V1, V2, U0, U1, U2, &coplanar, source, target);
	  //found = tri_tri_overlap_test_3d(V0, V1, V2, U0, U1, U2);
	  if(found == 1){
	    for(int p = 0; p < 3; p++)
	      for(int q = 0; q < 3; q++)
                if((face[i].vert[p].x < face[j].vert[q].x + eps1 && face[i].vert[p].x > face[j].vert[q].x - eps1) &&
                   (face[i].vert[p].y < face[j].vert[q].y + eps1 && face[i].vert[p].y > face[j].vert[q].y - eps1) &&
                   (face[i].vert[p].z < face[j].vert[q].z + eps1 && face[i].vert[p].z > face[j].vert[q].z - eps1))
	    	  found = 0;
	    if(coplanar == 1)
	      found = 0;
	  }
	  if(found == 1){
	    intcount++;
	    face[i].in = 1;
	    face[j].in = 1;
	  }
	}
    cout << intcount << " intersections " << endl;
  }

  vector<elements> facetoadd; 
  int fta_index = 0;
  if(remesh == true){
    cout << "Remeshing the structure... ";
    int intcount = 0;
    vector<Vector> seg_chain;
    int ichain = 0; 
    for(int i = 0; i < tn; i++)
      for(int j = 0; j < tn; j++)
	if(i != j){
	  int found = 0;
	  float V0[3] = {float(face[i].vert[0].x), float(face[i].vert[0].y), float(face[i].vert[0].z)};
	  float V1[3] = {float(face[i].vert[1].x), float(face[i].vert[1].y), float(face[i].vert[1].z)};
	  float V2[3] = {float(face[i].vert[2].x), float(face[i].vert[2].y), float(face[i].vert[2].z)};
	  float U0[3] = {float(face[j].vert[0].x), float(face[j].vert[0].y), float(face[j].vert[0].z)};
	  float U1[3] = {float(face[j].vert[1].x), float(face[j].vert[1].y), float(face[j].vert[1].z)};
	  float U2[3] = {float(face[j].vert[2].x), float(face[j].vert[2].y), float(face[j].vert[2].z)};
	  int coplanar=0;
	  float source[3];
	  float target[3];
	  found =  tri_tri_intersection_test_3d(V0, V1, V2, U0, U1, U2, &coplanar, source, target);
	  if(found == 1){
	    for(int p = 0; p < 3; p++)
	      for(int q = 0; q < 3; q++)
		if((face[i].vert[p].x == face[j].vert[q].x) && (face[i].vert[p].y == face[j].vert[q].y) && (face[i].vert[p].z == face[j].vert[q].z))
		  found = 0;
	    if(coplanar == 1)
	      found = 0;
	  }
	  if(found == 1){
	    intcount++;
	    face[i].in = 1;
	    face[j].in = 1;
	    bool inchain = false;
	    if(ichain > 2){
	      for(int k = 0; k < ichain-1; k++){
		if(((source[0] == seg_chain[k].x && source[1] == seg_chain[k].y && source[2] == seg_chain[k].z) &&
		    (target[0] == seg_chain[k+1].x && target[1] == seg_chain[k+1].y && target[2] == seg_chain[k+1].z)) ||
		   ((source[0] == seg_chain[k+1].x && source[1] == seg_chain[k+1].y && source[2] == seg_chain[k+1].z) &&
		    (target[0] == seg_chain[k].x && target[1] == seg_chain[k].y && target[2] == seg_chain[k].z)))
		  inchain = true;
		k++;
	      }
	      if(inchain == false){
		seg_chain.push_back(Vector());
		seg_chain[ichain].x = source[0];
		seg_chain[ichain].y = source[1];
		seg_chain[ichain].z = source[2];
		ichain++;
		seg_chain.push_back(Vector());
		seg_chain[ichain].x = target[0];
		seg_chain[ichain].y = target[1];
		seg_chain[ichain].z = target[2];
		ichain++;
	      }
	    }
	    else{
	      seg_chain.push_back(Vector());
	      seg_chain[ichain].x = source[0];
	      seg_chain[ichain].y = source[1];
	      seg_chain[ichain].z = source[2];
	      ichain++;
	      seg_chain.push_back(Vector());
	      seg_chain[ichain].x = target[0];
	      seg_chain[ichain].y = target[1];
	      seg_chain[ichain].z = target[2];
	      ichain++;
	    }
	  }
	}
    outfile.open("seglines.txt");
    for(int i = 0; i < seg_chain.size(); i++)
    outfile << seg_chain[i].x << " " << seg_chain[i].y << " " << seg_chain[i].z << endl;
    outfile.close();
    cout << intcount << " intersections " << endl;
    //cout << "intersection chain length " << ichain << endl;
    //------------------------island removal-----------------------
    if(island == true){
      cout << "Removing isolated faces with a connectivity of less than " << islecon << "..." << endl;
      for(int i = 0; i < tn; i++)
	if(face[i].in != 1 && face[i].out != 1){
	  vector<int> connectivity;
	  connectivity.push_back(i);
	  face[i].out = 1;
	  for(int a = 0; a < connectivity.size(); a++){
	    for(int j = 0; j < tn; j++){
	      if(face[j].in != 1 && face[j].out != 1){
		int found = 0;
		for(int p = 0; p < 3; p++)
		  for(int q = 0; q < 3; q++)
		    if((face[connectivity[a]].vert[p].x == face[j].vert[q].x) && (face[connectivity[a]].vert[p].y == face[j].vert[q].y) && (face[connectivity[a]].vert[p].z == face[j].vert[q].z))
		      found++;  
		if(found == 2){
		  connectivity.push_back(j);
		  face[j].out = 1;
		}
	      }
	    }
	  }
	  if(connectivity.size() < islecon){
	    for(int k = 0; k < connectivity.size(); k++){
	      face[connectivity[k]].out = 1;
	      face[connectivity[k]].in = 2;
	    }
	  }
	}
    }
    //------------------------end island removal-----------------------
    //experimental functions
    if(experimental == true){    
      //------------------------find only those with 1 or 2 shared borders---------------------------------
      vector<int> incomplete_index;
      vector< vector<int> > incomplete_edges;
      for(int i = 0; i < tn; i++)
	if(face[i].in == 0){
	  vector<Vector> edgehold;
	  int edges = 0;
	  for(int j = 0; j < tn; j++)
	    if(face[j].in == 0 && i != j){
	      int edgehold[3] = {-1, -1, -1};
	      for(int p = 0; p < 3; p++)
		if((face[i].vert[0].x == face[j].vert[p].x) &&
		   (face[i].vert[0].y == face[j].vert[p].y) &&
		   (face[i].vert[0].z == face[j].vert[p].z)){
		  for(int q = 0; q < 3; q++)
		    if((face[i].vert[1].x == face[j].vert[q].x) &&
		       (face[i].vert[1].y == face[j].vert[q].y) &&
		       (face[i].vert[1].z == face[j].vert[q].z)){
		      edgehold[0]++; //vert0 and vert1 edge
		    }
		}
	      for(int p = 0; p < 3; p++)
		if((face[i].vert[1].x == face[j].vert[p].x) &&
		   (face[i].vert[1].y == face[j].vert[p].y) &&
		   (face[i].vert[1].z == face[j].vert[p].z)){
		  for(int q = 0; q < 3; q++)
		    if((face[i].vert[2].x == face[j].vert[q].x) &&
		       (face[i].vert[2].y == face[j].vert[q].y) &&
		       (face[i].vert[2].z == face[j].vert[q].z)){
		      edgehold[1]++; //vert1 and vert2 edge
		    }
		}
	      for(int p = 0; p < 3; p++)
		if((face[i].vert[2].x == face[j].vert[p].x) &&
		   (face[i].vert[2].y == face[j].vert[p].y) &&
		   (face[i].vert[2].z == face[j].vert[p].z)){
		  for(int q = 0; q < 3; q++)
		    if((face[i].vert[0].x == face[j].vert[q].x) &&
		       (face[i].vert[0].y == face[j].vert[q].y) &&
		       (face[i].vert[0].z == face[j].vert[q].z)){
		      edgehold[2]++; //vert2 and vert0 edge
		    }
		}
	    }

	  for(int j = 0; j < tn; j++)
	    if(face[j].in == 0 && i != j){
	      int found = 0;	
	      for(int p = 0; p < 3; p++)
		for(int q = 0; q < 3; q++)
		  if((face[i].vert[p].x == face[j].vert[q].x) &&
		     (face[i].vert[p].y == face[j].vert[q].y) &&
		     (face[i].vert[p].z == face[j].vert[q].z)){
		    found++;
		    //	    edgehold.push_back(face[i].vert[p]);
		  }
	      if(found == 2)
		edges++;
	    }
	  if(edges < 3){
	    incomplete_index.push_back(i);
	    //incomplete_found.push_back(edges);
	  }
	}
      
      
      cout << " Rebuilding... " << endl;
      //incomplete triangle chain && intersect line-segment chain done...
      //Build gap-filling triangles
      //idea1 start from incomplete triangles
      int holdj = 0;
      int holdk = 0;
      for(int i = 0; i < incomplete_index.size(); i++){
	double smallestdist = 1000;
	for(int j = 0; j < seg_chain.size(); j++)
	  for(int k =0; k < 3; k++){
	    Vector distance = seg_chain[j] - face[incomplete_index[i]].vert[k];
	    if(distance.Magnitude() < smallestdist){
	      smallestdist = distance.Magnitude();
	      holdj = j;
	      holdk = k;
	    }
	  }
	Vector vechold[3];
	int found = 1;
	while(found == 1)
	  for(int k = 0; k < 3; k++)
	    if(k != holdk){
	      vechold[0] = face[incomplete_index[i]].vert[holdk];
	      vechold[1] = face[incomplete_index[i]].vert[k];
	      vechold[2] = seg_chain[holdj+1];
	      float V0[3] = {float(vechold[0].x), float(vechold[0].y), float(vechold[0].z)};
	      float V1[3] = {float(vechold[1].x), float(vechold[1].y), float(vechold[1].z)};
	      float V2[3] = {float(vechold[2].x), float(vechold[2].y), float(vechold[2].z)};
	      float U0[3] = {float(face[incomplete_index[i]].vert[0].x), float(face[incomplete_index[i]].vert[0].y), float(face[incomplete_index[i]].vert[0].z)};
	      float U1[3] = {float(face[incomplete_index[i]].vert[1].x), float(face[incomplete_index[i]].vert[1].y), float(face[incomplete_index[i]].vert[1].z)};
	      float U2[3] = {float(face[incomplete_index[i]].vert[2].x), float(face[incomplete_index[i]].vert[2].y), float(face[incomplete_index[i]].vert[2].z)};
	      int coplanar;
	      float source[3];
	      float target[3];
	      found =  tri_tri_intersection_test_3d(V0, V1, V2, U0, U1, U2, &coplanar, source, target);
	      if(found == 1){
		for(int p = 0; p < 3; p++)
		  for(int q = 0; q < 3; q++)
		    if((vechold[p].x == face[incomplete_index[i]].vert[q].x) &&
		       (vechold[p].y == face[incomplete_index[i]].vert[q].y) &&
		       (vechold[p].z == face[incomplete_index[i]].vert[q].z))
		      found++;
		if(coplanar == 1)
		  found = 0;
	      }
	    }
	//add new element
	facetoadd.push_back(elements());      
	facetoadd[fta_index].vert[0] = vechold[0];
	facetoadd[fta_index].vert[1] = vechold[1];
	facetoadd[fta_index].vert[2] = vechold[2];
	line1 = facetoadd[fta_index].vert[1] - facetoadd[fta_index].vert[0];
	line2 = facetoadd[fta_index].vert[2] - facetoadd[fta_index].vert[0];
	facetoadd[fta_index].norm = CrossProduct(line2, line1).UnitVector();
	fta_index++;
      }
    
	//idea2 start from line-segment chain
	/*
      for(int i = 0; i < seg_chain.size(); i++){
	double smallestdist = 1000;
	for(int j = 0; j < incomplete_index.size(); j++)
	  if(face[incomplete_index[j]].in == 0)
	    for(int k = 0; k < 3; k++){
	      Vector distance = seg_chain[i] - face[incomplete_index[j]].vert[k];
	      if(distance.Magnitude() < smallestdist){
		smallestdist = distance.Magnitude();
		holdj = incomplete_index[j];
		holdk = k;
	      }
	    }
	bool intersects = false;
	Vector vechold[3];
	vechold[0] = face[holdj].vert[holdk];
	vechold[1] = seg_chain[i];
	vechold[2] = seg_chain[i+1];
	if(facetoadd.size() > 0){
	  for(int j = 0; j < facetoadd.size(); j++){
	    int found = 0;
	    float V0[3] = {float(vechold[0].x), float(vechold[0].y), float(vechold[0].z)};
	    float V1[3] = {float(vechold[1].x), float(vechold[1].y), float(vechold[1].z)};
	    float V2[3] = {float(vechold[2].x), float(vechold[2].y), float(vechold[2].z)};
	    float U0[3] = {float(facetoadd[j].vert[0].x), float(facetoadd[j].vert[0].y), float(facetoadd[j].vert[0].z)};
	    float U1[3] = {float(facetoadd[j].vert[1].x), float(facetoadd[j].vert[1].y), float(facetoadd[j].vert[1].z)};
	    float U2[3] = {float(facetoadd[j].vert[2].x), float(facetoadd[j].vert[2].y), float(facetoadd[j].vert[2].z)};
	    int coplanar;
	    float source[3];
	    float target[3];
	    found =  tri_tri_intersection_test_3d(V0, V1, V2, U0, U1, U2, &coplanar, source, target);
	    if(found == 1){
	      for(int p = 0; p < 3; p++)
		for(int q = 0; q < 3; q++)
		  if((vechold[p].x == facetoadd[j].vert[q].x) &&
		     (vechold[p].y == facetoadd[j].vert[q].y) &&
		     (vechold[p].z == facetoadd[j].vert[q].z))
		    found = 0;
	      if(coplanar == 1)
		found = 0;
	    }
	    if(found == 1){
	      intersects = true;
	      break;
	    }
	  }
	  if(intersects == false){
	     facetoadd.push_back(elements());      
	     facetoadd[fta_index].vert[0] = vechold[0];
	     facetoadd[fta_index].vert[1] = vechold[1];
	     facetoadd[fta_index].vert[2] = vechold[2];
	     line1 = facetoadd[fta_index].vert[1] - facetoadd[fta_index].vert[0];
	     line2 = facetoadd[fta_index].vert[2] - facetoadd[fta_index].vert[0];
	     facetoadd[fta_index].norm = CrossProduct(line2, line1).UnitVector();
	     fta_index++;
	     }
	     else{
	     //facetoadd.push_back(elements());
	     //facetoadd[fta_index].vert[0] =  
	     }
	  }
	  else{
	  facetoadd.push_back(elements());      
	  facetoadd[fta_index].vert[0] = vechold[0];
	  facetoadd[fta_index].vert[1] = vechold[1];
	  facetoadd[fta_index].vert[2] = vechold[2];
	  line1 = facetoadd[fta_index].vert[1] - facetoadd[fta_index].vert[0];
	  line2 = facetoadd[fta_index].vert[2] - facetoadd[fta_index].vert[0];
	  facetoadd[fta_index].norm = CrossProduct(line2, line1).UnitVector();
	  fta_index++;
	  }
	
	i++;
	}*/
  
    }
  }
  //------------------------island removal-----------------------
  if((island == true) && (remesh == false)){
    cout << "Removing isolated faces with a connectivity of less than " << islecon << "..." << endl;
    for(int i = 0; i < tn; i++)
      if(face[i].in != 1 && face[i].out != 1){
	vector<int> connectivity;
	connectivity.push_back(i);
	face[i].out = 1;
	for(int a = 0; a < connectivity.size(); a++){
	  for(int j = 0; j < tn; j++){
	    if(face[j].in != 1 && face[j].out != 1){
	      int found = 0;
	      for(int p = 0; p < 3; p++)
		for(int q = 0; q < 3; q++)
		  if((face[connectivity[a]].vert[p].x == face[j].vert[q].x) &&
		     (face[connectivity[a]].vert[p].y == face[j].vert[q].y) &&
		     (face[connectivity[a]].vert[p].z == face[j].vert[q].z))
		    found++;  
	      if(found == 2){
		connectivity.push_back(j);
		face[j].out = 1;
	      }
	    }
	  }
	}
	if(connectivity.size() < islecon){
	  for(int k = 0; k < connectivity.size(); k++){
	    face[connectivity[k]].out = 1;
	    face[connectivity[k]].in = 2;
	  }
	}
      }
  }
  //------------------------end island removal-----------------------
  for(int i = 0; i < tn; i++)
    if(face[i].in == 0){
      for(int j = 0; j < 3; j++){
	if(face[i].vert[j].x < xmin)
	  face[i].vert[j].x = xmin;
      	if(face[i].vert[j].y < ymin)
	  face[i].vert[j].y = ymin;
	if(face[i].vert[j].z < zmin)
	  face[i].vert[j].z = zmin;
	if(face[i].vert[j].x > xmax)
	  face[i].vert[j].x = xmax;
	if(face[i].vert[j].y > ymax)
	  face[i].vert[j].y = ymax;
	if(face[i].vert[j].z > zmax)
	  face[i].vert[j].z = zmax;
      }
    }

  cout << " Writing output to newmesh.stl " << endl;
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
  for(int i = 0; i < fta_index; i++){
    outfile << setprecision(6) << "facet normal \t" << facetoadd[i].norm.x << "\t" <<  facetoadd[i].norm.y << "\t" << facetoadd[i].norm.z << endl;
    outfile << "outer loop" << endl;
    for(int j = 0; j < 3; j++)
      outfile << setprecision(6) << "vertex \t" << facetoadd[i].vert[j].x << "\t" << facetoadd[i].vert[j].y << "\t" << facetoadd[i].vert[j].z << endl;
    outfile << "endloop" << endl << "endfacet" << endl;
  }
  outfile.close();
  t = clock() - t;
  cout << float(t)/CLOCKS_PER_SEC << " seconds" << endl;

  delete [] face;
  delete [] original;
  return 0;
}


