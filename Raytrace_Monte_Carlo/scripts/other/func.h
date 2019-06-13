//Vector class and Functions stored here. 
//Andrew Alvarado
//May 1 2019
//Professor Marian Research group

#include<fstream>
#include<string>
#include<iostream>
#include<iomanip>
#include<cmath>
#include<vector>
using namespace std;
#define eps 0.00001
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
  //inpoint = 0;
  //find vector for two edges sharing vert0
  edge1 = vert1 - vert0;
  edge2 = vert2 - vert0;
  //Start calculating determinant - also used to calculate U parameter 
  pvec = CrossProduct(dir, edge2);
  // if determinant is near zero, ray lies in plane of triangle
  det =  DotProduct(edge1, pvec);
  //#ifdef TEST_CULL // define test_cull if culling is desired
  //*if(det < eps)
  //  return 0;*/
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
  Vector center = {};
};
struct ray{
  Vector Origin;
  Vector Direction;
 };

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

        

