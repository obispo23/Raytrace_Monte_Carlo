//Vector class and Functions stored here. 
//Andrew Alvarado
//Sept 25 2018
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

//--------------------------M-T algorithm-----------------                                                                                             
int intersect_triangle(Vector orig, Vector dir, Vector vert0, Vector vert1, Vector vert2, double& inpoint, int& found)
{
  Vector edge1, edge2, tvec, pvec, qvec;
  double det, inv_det;
  double t, u, v;
  found = 0;
  inpoint = 0;
  //find vector for two edges sharing vert0
  edge1 = vert1 - vert0;
  edge2 = vert2 - vert0;
  //Start calculating determinant - also used to calculate U parameter 
  pvec = CrossProduct(dir, edge2);
  // if determinant is near zero, ray lies in plane of triangle
  det =  DotProduct(edge1, pvec);
  ///calculate distance from vert0 to ray origin
  tvec = orig - vert0;
  inv_det = 1.0 / det;
  if(det > eps){
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
  }
  else if(det < -eps){
    //calculate U paramter and test bounds
    u = DotProduct(tvec, pvec);
    if(u > 0.0 || u < det)
      return 0;
    //prepare to test V parameter
    qvec = CrossProduct(tvec, edge1);
    //calculate V parameter and test bounds
    v = DotProduct(dir, qvec);
    if(v > 0.0 || u + v < det)
      return 0;
  }
  else{
    return 0;
  }
  // calculate t, scale parameters, ray intersects triangle
  t = DotProduct(edge2, qvec)*inv_det;
  u *= inv_det;
  v *= inv_det;
  return 1; 
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
double cwall(ray rays, double& u, double a, double b, double c, int& found){
  u = 0.0;
  found = 0;
  Vector find; 
  if(rays.Direction.x <= eps){
    u = (-a-rays.Origin.x)/rays.Direction.x;  //left face       
    if(u > eps){
      find = rays.Origin + u * rays.Direction;
      if((fabs(find.x) <= a+eps) && (fabs(find.y) <= b+eps) && (fabs(find.z) <= c+eps)){
	found = 3;
	return 0;
      }
    }
  }
  if(rays.Direction.x >= eps){
    u = (a-rays.Origin.x)/rays.Direction.x;   //right face   
    if(u > eps){
      find = rays.Origin + u * rays.Direction;
      if((fabs(find.x) <= a+eps) && (fabs(find.y) <= b+eps) && (fabs(find.z) <= c+eps)){
	found = 4;
	return 0;
      }
    }
  }
  if(rays.Direction.y <= eps){
    u = (-b-rays.Origin.y)/rays.Direction.y;  //front face    
    if(u > eps){
      find = rays.Origin + u * rays.Direction;
      if((fabs(find.x) <= a+eps) && (fabs(find.y) <= b+eps) && (fabs(find.z) <= c+eps)){
	found = 5;
	return 0;
      }
    }
  }
  if(rays.Direction.y >= eps){
    u = (b-rays.Origin.y)/rays.Direction.y;   //back face
    if(u > eps){
      find = rays.Origin + u * rays.Direction;
      if((fabs(find.x) <= a+eps) && (fabs(find.y) <= b+eps) && (fabs(find.z) <= c+eps)){
	found = 6;
	return 0;
      }
    }
  }
  if(rays.Direction.z >= eps){
    u = (c-rays.Origin.z)/rays.Direction.z;   //top face                      
    if(u > eps){
      find = rays.Origin + u * rays.Direction;
      if((fabs(find.x) <= a+eps) && (fabs(find.y) <= b+eps) && (fabs(find.z) <= c+eps)){
	found = 2;
	return 0;
      }
    }
  }
  if(rays.Direction.z <= eps){
    u = (-c-rays.Origin.z)/rays.Direction.z;  //bottom face         
    if(u > eps){
      find = rays.Origin + u * rays.Direction;
      if((fabs(find.x) <= a+eps) && (fabs(find.y) <= b+eps) && (fabs(find.z) <= c+eps)){
	found = 7;
	return 0;
      }
    }
  }
  //                                                                                                
  if((rays.Origin.z >= c + eps) && (rays.Direction.z <= eps)){ //initial beams only
    if(rays.Direction.x <= eps){
      u = (-a-rays.Origin.x)/rays.Direction.x;  //left face 
      if(u > eps){
        find = rays.Origin + u * rays.Direction;
	if((fabs(find.x) <= a+eps) && (fabs(find.y) <= b+eps)){
          found = 3;
          return 0;
        }
      }
    }
    if(rays.Direction.x >= eps){
      u = (a-rays.Origin.x)/rays.Direction.x;   //right face 
      if(u > eps){
        find = rays.Origin + u * rays.Direction;
        if((fabs(find.x) <= a+eps) && (fabs(find.y) <= b+eps)){
          found = 4;
          return 0;
        }
      }
    }
    if(rays.Direction.y <= eps){
      u = (-b-rays.Origin.y)/rays.Direction.y;  //front face 
      if(u > eps){
        find = rays.Origin + u * rays.Direction;
        if((fabs(find.x) <= a+eps) && (fabs(find.y) <= b+eps)){
          found = 5;
          return 0;
        }
      }
    }
    if(rays.Direction.y >= eps){
      u = (b-rays.Origin.y)/rays.Direction.y;   //back face
      if(u > eps){
        find = rays.Origin + u * rays.Direction;
        if((fabs(find.x) <= a+eps) && (fabs(find.y) <= b+eps)){
          found = 6;
          return 0;
        }
      }
    }
    if(rays.Direction.z <= eps){
      u = (-c-rays.Origin.z)/rays.Direction.z;  //bottom face
      if(u > eps){
        find = rays.Origin + u * rays.Direction;
        if((fabs(find.x) <= a+eps) && (fabs(find.y) <= b+eps)){
          found = 7;
          return 0;
        }
      }
    }
  }
  if(found ==0)
    u = 0;
  return 0;
}
double ofind(Vector POT, ray &rays, int found){
  if(found == 0)
    return 0;
  if(found == 1)
    return 0;
  if(found == 2)
    return 0;
  if(found == 3){ // left face = right face 
    rays.Origin.x = -POT.x;
    rays.Origin.y = POT.y;
    rays.Origin.z = POT.z;
  }
  if(found == 4){ // right face = left face
    rays.Origin.x = -POT.x;
    rays.Origin.y = POT.y;
    rays.Origin.z = POT.z;
  }
  if(found == 5){ // front face = back face
    rays.Origin.x = POT.x;
    rays.Origin.y = -POT.y;
    rays.Origin.z = POT.z;
  }
  if(found == 6){ // back face = front face
    rays.Origin.x = POT.x;
    rays.Origin.y = -POT.y;
    rays.Origin.z = POT.z;
  }
  if(found == 7){ // bottom face = top face
    return 0;
  }
  return 0;
}

void grabber(string input, string &first, string &last, char mid){
  string hold;
  hold = input;
  int index = 0;
  for(int i = 0 ; i < input.length(); i++){
    if(hold[i] == mid)
      index = i;
  }
  last = input.substr(index+1);
  first = input.substr(0, index);

}

void search(int tn, ray rays, elements face[], int &found, int &found1, int &hold, Vector &POT, double MAXDIST){
  double  dist = 0;
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

