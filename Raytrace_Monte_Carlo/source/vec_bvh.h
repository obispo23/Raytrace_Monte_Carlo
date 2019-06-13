//Vector class and Functions stored here. 
//Andrew Alvarado
//Feb 4 2019
//Professor Marian Research group

#include<fstream>
#include<string>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<random>
#include<ctime>
using namespace std;
#define eps 0.000001
class Vector
{
    public:
  double x;
  double y;
  double z;
  
        // Default constructor: create a vector whose 
        // x, y, z components are all zero.
  Vector()
    : x(0.0)
    , y(0.0)
    , z(0.0)
  {
  }
  
  // This constructor initializes a vector 
  // to any desired component values.
  Vector(double _x, double _y, double _z)
    : x(_x)
    , y(_y)
    , z(_z)
  {
  }
  
  // Returns the square of the magnitude of this vector.
  // This is more efficient than computing the magnitude itself,
  // and is just as good for comparing two vectors to see which
  // is longer or shorter.
  const double MagnitudeSquared() const
  {
    return (x*x) + (y*y) + (z*z);
  }
  
  const double Magnitude() const
  {
    return sqrt(MagnitudeSquared());
  }
  
  const Vector UnitVector() const
  {
    const double mag = Magnitude();
    return Vector(x/mag, y/mag, z/mag);
  }
  
  Vector& operator *= (const double factor)
  {
    x *= factor;
    y *= factor;
    z *= factor;
    return *this;
  }
  
  Vector& operator += (const Vector& other)
  {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }
};


//------------------------------------------------------------------------

inline Vector operator + (const Vector &a, const Vector &b)
{
  return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vector operator - (const Vector &a, const Vector &b)
{
  return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Vector operator - (const Vector& a)
{
  return Vector(-a.x, -a.y, -a.z);
}

inline double DotProduct (const Vector& a, const Vector& b) 
{
  return (a.x*b.x) + (a.y*b.y) + (a.z*b.z);
}

inline Vector CrossProduct (const Vector& a, const Vector& b)
{
  return Vector(
		(a.y * b.z) - (a.z * b.y), 
		(a.z * b.x) - (a.x * b.z), 
		(a.x * b.y) - (a.y * b.x));
}

inline Vector operator * (double s, const Vector& v)
{
  return Vector(s*v.x, s*v.y, s*v.z);
}

inline Vector operator / (const Vector& v, double s)
{
  return Vector(v.x/s, v.y/s, v.z/s);
}
//randomizers

default_random_engine generator(random_device{}()); //stochastic random engine
uniform_real_distribution<double> unifrnd(0, 1); //-----------
uniform_real_distribution<double> unifrndEnergy(100, 1000);//-

//--------------------------M-T algorithm-----------------                                                                                             
int intersect_triangle(Vector orig, Vector dir, Vector vert0, Vector vert1, Vector vert2, double& inpoint, int& found)
{
  Vector edge1, edge2, tvec, pvec, qvec;
  double det, inv_det;
  double t, u, v;
  found = 0;
  //inpoint = 0;
  //find vector for two edges sharing vert0
  edge1 = vert1 - vert0;
  edge2 = vert2 - vert0;
  //Start calculating determinant - also used to calculate U parameter 
  pvec = CrossProduct(dir, edge2);
  // if determinant is near zero, ray lies in plane of triangle
  det =  DotProduct(edge1, pvec);
  //#ifdef TEST_CULL // define test_cull if culling is desired
  //if(det < eps)
  //  return 0;
  //#else
  if(det > -eps && det < eps)
    return 0;
  //#endif
  //calculate distance from vert0 to ray origin
  tvec = orig - vert0;
  //calculate U paramter and test bounds
  u = DotProduct(tvec, pvec);
  if(u < 0.0 || u > det)
    return 0;
  //prepare to test V parameter   
  qvec = CrossProduct(tvec, edge1);
  //calculate V parameter and test bounds
  v = DotProduct(dir, qvec);
  if(v < 0.0 || u + v > det)
    return 0;
  // calculate t, scale parameters, ray intersects triangle  
  t = DotProduct(edge2, qvec);
  inv_det = 1.0 / det;
  t *= inv_det;
  u *= inv_det;
  v *= inv_det;
  if(t > eps){
    t = DotProduct(edge2, qvec) * inv_det;
    inpoint = t;
    //cout << "Great shot kid that was one in a million" << endl; 
    found = 1;
    return 1;
  }
  else{
    return 0;
  }
}
//---------------------------------------------------  
struct elements{
  Vector norm;
  Vector vert[3] = {};
  int in = 0;
  int out = 0;
  bool used[3] = {};
  double deltah;
};
struct ray{
  Vector Origin;
  Vector Direction;
  string spec;
 };

struct node{
  double Energy;
  double Angle;
  ray rays;
  double Gen;
  node* next;
};
struct node *newNode(double Energy, double Angle, ray rays, double Gen){
  node *trail = new node;
  trail->Energy = Energy;
  trail->Angle = Angle;
  trail->rays = rays;
  trail->Gen = Gen;
  trail->next = NULL;
  return trail;
}
struct node* insertKey(struct node *trail, double Energy, double Angle, ray rays, double Gen){
  if(trail == NULL)
    return(newNode(Energy, Angle, rays, Gen));
  else{
    trail->next = insertKey(trail->next, Energy, Angle, rays, Gen);
  }
  return trail;
}

void clear(struct node *trail){
  if(trail != NULL){
    clear(trail->next);
    delete trail;
  }
}

//Cuboid intersect find: boundary condition solve
//---For rays that miss the foam!
int bc(ray rays, double &u, double min[], double max[], int &found){
  Vector invdir = Vector(1.0/rays.Direction.x, 1.0/rays.Direction.y, 1.0/rays.Direction.z);
  Vector find;
  double ux, uy, uz;
  double tmin, tmax;
  bool escape = false;
  if(rays.Direction.x >= 0){
    tmin = (min[0] - rays.Origin.x) * invdir.x;
    tmax = (max[0] - rays.Origin.x) * invdir.x;
  }
  else{
    tmin = (max[0] - rays.Origin.x) * invdir.x;
    tmax = (min[0] - rays.Origin.x) * invdir.x;
  }
  if(tmax <= 0)
    ux = tmin;
  else
    ux = tmax;
  
  if(rays.Direction.y >= 0){
    tmin = (min[1] - rays.Origin.y) * invdir.y;
    tmax = (max[1] - rays.Origin.y) * invdir.y;
  }
  else{
    tmin = (max[1] - rays.Origin.y) * invdir.y;
    tmax = (min[1] - rays.Origin.y) * invdir.y;
  }
  if(tmax <= 0)
    uy = tmin;
  else
    uy = tmax;
  
  if(rays.Direction.z >= 0){
    tmin = (min[2] - rays.Origin.z) * invdir.z;
    tmax = (max[2] - rays.Origin.z) * invdir.z;
    escape = true;
  }
  else{
    tmin = (max[2] - rays.Origin.z) * invdir.z;
    tmax = (min[2] - rays.Origin.z) * invdir.z;
  }
  if((tmin > 0) && (tmax > 0)){//for initial rays
    if(tmin < tmax)
      uz = tmax;
    else
      uz = tmin;
  }
  else{
    if(tmax <= 0)
      uz = tmin;
    else 
      uz = tmax;
  }
  if((ux <= uy) && (ux <= uz)){ //ux is smallest distance
    u = ux-eps;
    found = 3;
  }
  if((uy <= ux) && (uy <= uz)){ //uy is smallest distance
    u = uy-eps;
    found = 5;
  }
  if((uz <= ux) && (uz <= uy)){ //uz is smallest distance
    u = uz-eps;
    found = 7;
    if(escape == true)
      found = 2;    
  }
  //if(found == 0)
  
  
  return found;
}
int bounce(Vector POT, ray &rays, int found){
  if(found == 0 || found == 1 || found == 2 || found ==7)
    return 0;
  rays.Origin.x = POT.x;
  rays.Origin.y = POT.y;
  rays.Origin.z = POT.z;
  
  if(found == 3 || found == 4) //bouncing off left or right face
    rays.Direction.x *= -1;
  if(found == 5 || found == 6) // bouncing off front or back face
    rays.Direction.y *= -1;
    //rays.Direction.z *= -1;
  return 0;   
}
void grabber(string input, string &first, string &last, char mid){
  string hold;
  hold = input;
  int index = 0;
  bool midthere = false;
  for(int i = 0 ; i < input.length(); i++){
    if(hold[i] == mid){
      index = i;
      midthere = true;
    }
  }
  if(midthere == false)
    first = input;
  else{
    last = input.substr(index+1);
    first = input.substr(0, index);
  }
}
int ray_box_intersect(ray rays, elements face[], Vector &POT, Vector bvmax[], Vector bvmin[], int boxes, vector< vector<int> > bvh_id, double MAXDIST, int &hold, int &found){
  Vector invdir = Vector(1.0/rays.Direction.x, 1.0/rays.Direction.y, 1.0/rays.Direction.z);
  //int hold = -1;
  double distmax = MAXDIST;
  double dist = 1000;
  found = 0;
  int found1 = 0;
  double tmin, tmax, tymin, tymax, tzmin, tzmax;
  for(int i = 0; i < boxes; i++){
    if(rays.Direction.x >= 0){
      tmin = (bvmin[i].x - rays.Origin.x) * invdir.x;
      tmax = (bvmax[i].x - rays.Origin.x) * invdir.x;
    }
    else{
      tmin = (bvmax[i].x - rays.Origin.x) * invdir.x;
      tmax = (bvmin[i].x - rays.Origin.x) * invdir.x;
    }
    if(rays.Direction.y >= 0){
      tymin = (bvmin[i].y - rays.Origin.y) * invdir.y;
      tymax = (bvmax[i].y - rays.Origin.y) * invdir.y;
    }
    else{
      tymin = (bvmax[i].y - rays.Origin.y) * invdir.y;
      tymax = (bvmin[i].y - rays.Origin.y) * invdir.y;
    }
    if((tmin > tymax) || (tymin > tmax))
      continue;
    if(tymin > tmin)
      tmin = tymin;
    if(tymax < tmax)
      tmax = tymax;
    if(rays.Direction.z >= 0){
      tzmin = (bvmin[i].z - rays.Origin.z) * invdir.z;
      tzmax = (bvmax[i].z - rays.Origin.z) * invdir.z;
    }
    else{
      tzmin = (bvmax[i].z - rays.Origin.z) * invdir.z;
      tzmax = (bvmin[i].z - rays.Origin.z) * invdir.z;
    }
    if((tmin > tzmax) || (tzmin > tmax))
      continue;
    if(tzmin > tmin)
      tmin = tzmin;
    if(tzmax < tmax)
      tmax = tzmax;
    found = 0;
    for(int j = 0; j < bvh_id[i].size(); j++){
      intersect_triangle(rays.Origin, rays.Direction, face[bvh_id[i][j]].vert[0], face[bvh_id[i][j]].vert[1], face[bvh_id[i][j]].vert[2], dist, found);
      if(found == 1)
	if(dist < distmax){
	  distmax = dist;
	  hold = bvh_id[i][j];
	  POT = rays.Origin + dist*rays.Direction;
	  found1 = found;
	}
    }
  }
  found = found1;
  return found1;
}
void search(int tn, ray rays, elements face[], int &found, int &found1, int &hold, Vector &POT, double MAXDIST){
  double  dist = 1000;
  double  distmax = MAXDIST;
  found = 0;
  found1 = 0;
  for(int i = 0; i < tn; i++){
    intersect_triangle(rays.Origin, rays.Direction, face[i].vert[0], face[i].vert[1], face[i].vert[2], dist, found);
    if(found == 1){
      if(dist < distmax){
	distmax = dist;
	hold = i;
	POT = rays.Origin + dist*rays.Direction;
	found1 = found;
      }
    }
  }
}
//angle distribution for outgoing angles, can choose cosine or a txt file with species#_AngDis.txt (see example of table)
double angle_distribution(string distribution, ray &rays, elements face[], double Energy, double angle, double min[], double max[], Vector Normal, int hold){
  double angle1;
  double sin1, cos1, a1, b1, c1;
  Vector line1, line2;
  ifstream read;
  ofstream out;
  double psi;
  if(distribution == "cosine"){
  here:
    if(hold == -1)
      line1 = (Vector(min[0]+eps, min[1]+eps, min[2]+eps) - Vector(max[0]-eps, max[1]-eps, min[2]+eps)).UnitVector();
    else
      line1 = (face[hold].vert[1] - face[hold].vert[0]).UnitVector();
    line2 = CrossProduct(Normal,line1).UnitVector();
    sin1 = sqrt(unifrnd(generator));
    cos1 = sqrt(1.0 - sin1*sin1); 
    psi = 2.0*M_PI*unifrnd(generator);
    a1 = sin1*cos(psi);
    b1 = sin1*sin(psi);
    c1 = cos1;
    rays.Direction = ((a1*line1) + (b1*line2) + (c1*Normal)).UnitVector();
    angle1 = acos(cos1);
  }
  else if(distribution != "cosine"){
    double energy_table, angle_table;
    string temp, first, last;
    int rowcount = 0, rowhold = 0;
    double randnum = unifrnd(generator);
    int index = 0;
    for(int i = 0 ; i < distribution.length(); i++){
      if(distribution[i] == '/'){
	index = i;
      }
    }
    if(index > 0){
      last = distribution.substr(index+1);
      first = distribution.substr(0, index+1);
      distribution = first+rays.spec+"_"+last;
    }
    else{
      distribution = rays.spec+"_"+distribution;
    }
    read.open(distribution.c_str());
    if(read.good() == false){
      cout << distribution << "file not found" << endl;
      goto here;
    }
        
    double c_E = 1000.0;
    double c_A = 1000.0;
    double interval[20];
    double binsize[20];
    double size = 0;
    double frequency;
    for(int i = 0; i < 20; i++){
      read >> interval[i];
      interval[i] = (interval[i] - 90.0)*M_PI/180.0;
    }
    getline(read, temp);
    getline(read, temp);
    while(!read.eof()){
      read >> energy_table >> angle_table;
      getline(read, temp);
      rowcount++;
      if((abs(Energy - energy_table) <= c_E) && (abs(angle - angle_table*M_PI/180.0) <= c_A)){
	c_E = abs(Energy - energy_table);
	c_A = abs(angle - angle_table*M_PI/180.0);
	rowhold = rowcount;
      }
    }
    read.close();
    read.open(distribution.c_str());
    getline(read, temp);
    getline(read, temp);
    for(int i = 0; i < rowhold; i++)
      getline(read, temp);

    read >> energy_table >> angle_table;
    //if((abs(Energy - energy_table) == c_E) && (abs(angle - angle_table*M_PI/180.0) == c_A))
      for(int i = 0; i < 20; i++){
	read >> frequency;
	size = size+frequency;
	binsize[i] = size;
      }
    read.close();
    int binhold = 0;
    for(int i  = 0; i < 20; i++)
      if(randnum <= binsize[i]){
	binhold = i;
	break;
      }

    if(hold == -1)
      line1 = (Vector(min[0]+eps, min[1]+eps, min[2]+eps) - Vector(max[0]-eps, max[1]-eps, min[2]+eps)).UnitVector();
    else
      line1 = (face[hold].vert[1] - face[hold].vert[0]).UnitVector();
    line2 = CrossProduct(Normal,line1).UnitVector();
    angle1 = interval[binhold]*randnum/binsize[binhold];
    
    sin1 = sin(angle1);
    cos1 = sqrt(1.0 - sin1*sin1); 
    psi = 2.0*M_PI*unifrnd(generator);
    a1 = sin1*cos(psi);
    b1 = sin1*sin(psi);
    c1 = cos1;
    rays.Direction = ((a1*line1) + (b1*line2) + (c1*Normal)).UnitVector();
    //    out.open("ain_aout_300.dat", ios::app);
    //if(Energy == 300.0)
    //  if((angle*180.0/M_PI > 40.0) && (angle*180.0/M_PI < 50.0))
    //	out << angle*180.0/M_PI << " " << angle1*180.0/M_PI << endl;
    //out.close();
  }
  
  return angle1;
}
