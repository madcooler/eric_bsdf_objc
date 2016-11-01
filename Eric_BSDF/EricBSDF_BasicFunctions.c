//
//  MicrosurfaceHeight.h
//  imp
//
//  Created by Chi Wang on 30/10/16.
//
//

#include "EricBSDF_BasicFunctions.h"

double generateRandomNumber()
{
    return 0;
}

double eric_erf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
 
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);
 
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
 
    return sign*y;
}

float eric_erfinv(float x)
{
    float w, p;
    w = - logf((1.0f-x)*(1.0f+x));
    if ( w < 5.000000f ) {
        w = w - 2.500000f;
        p = 2.81022636e-08f;
        p = 3.43273939e-07f + p*w;
        p = -3.5233877e-06f + p*w;
        p = -4.39150654e-06f + p*w;
        p = 0.00021858087f + p*w;
        p = -0.00125372503f + p*w;
        p = -0.00417768164f + p*w;
        p = 0.246640727f + p*w;
        p = 1.50140941f + p*w;
    }
    else {
        w = sqrtf(w) - 3.000000f;
        p = -0.000200214257f;
        p = 0.000100950558f + p*w;
        p = 0.00134934322f + p*w;
        p = -0.00367342844f + p*w;
        p = 0.00573950773f + p*w;
        p = -0.0076224613f + p*w;
        p = 0.00943887047f + p*w;
        p = 1.00167406f + p*w;
        p = 2.83297682f + p*w;
    }
    return p*x;
}

double  eric_abgam (double x)
{
  double  gam[10],
          temp;

  gam[0] = 1./ 12.;
  gam[1] = 1./ 30.;
  gam[2] = 53./ 210.;
  gam[3] = 195./ 371.;
  gam[4] = 22999./ 22737.;
  gam[5] = 29944523./ 19733142.;
  gam[6] = 109535241009./ 48264275462.;
  temp = 0.5*log (2*M_PI) - x + (x - 0.5)*log (x)
    + gam[0]/(x + gam[1]/(x + gam[2]/(x + gam[3]/(x + gam[4] /
	  (x + gam[5]/(x + gam[6]/x))))));

  return temp;
}

double  eric_gamma (double x)
{
  double  result;
  result = exp (eric_abgam (x + 5))/(x*(x + 1)*(x + 2)*(x + 3)*(x + 4));
  return result;
}

double  eric_beta (double m, double n)
{
  return (eric_gamma (m)* eric_gamma (n) / eric_gamma (m + n));
}

double sign(double x)
{
    if(x>0)
    return 1;
    else
    return -1;
}

int IsFiniteNumber(double x)
{
    if( x < INFINITY)
        return 1;
    else return 0;
}



 // build orthonormal basis (Building an Orthonormal Basis from a 3D Unit Vector Without Normalization, [Frisvad2012])
void buildOrthonormalBasis(struct vec3 * omega_1, struct vec3 * omega_2, const struct vec3 * omega_3)
{
	if(omega_3->z < -0.9999999f)
	{
	   *omega_1 = vec3 ( 0.0f , -1.0f , 0.0f );
	   *omega_2 = vec3 ( -1.0f , 0.0f , 0.0f );
	} else {
	   const float a = 1.0f /(1.0f + omega_3->z );
	   const float b = -omega_3->x*omega_3->y*a ;
	   *omega_1 = vec3 (1.0f - omega_3->x*omega_3->x*a , b , -omega_3->x );
	   *omega_2 = vec3 (b , 1.0f - omega_3->y*omega_3->y*a , -omega_3->y );
	}
}


// constructor of vec2
struct vec2 vec2(double x, double y)
{
    struct vec2 temp;
    temp.x = x;
    temp.y = y;
    return temp;
}

// constructor of vec3
struct vec3 vec3(double x, double y, double z)
{
    struct vec3 temp;
    temp.x = x;
    temp.y = y;
    temp.z = z;
    return temp;
}

float length(struct vec3 v)
{
    return sqrtf(v.x*v.x+v.y*v.y+v.z*v.z);
}

struct vec3 normalize(struct vec3 v)
{
    double l = length(v);
    return vec3(v.x/l,v.y/l,v.z/l);
}

struct vec3 vec3_add (struct vec3 v1,struct vec3 v2)
{
    return vec3(v1.x+v2.x,v1.y+v2.y,v1.z+v2.z);
}

struct vec3 vec3_mul (struct vec3 v,double l)
{
    return vec3(v.x*l,v.y*l,v.z*l);
}


struct vec3 negtive(struct vec3 v)
{
    return vec3(-v.x,-v.y,-v.z);
}

struct vec3 Vec3D_to_vec3(const Vec3D * v)
{
    double x = VEC3D_I(*v, 0);
    double y = VEC3D_I(*v, 1);
    double z = VEC3D_I(*v, 2);
    return vec3(x,y,z);
}

Vec3D vec3_to_Vec3D(struct vec3 v)
{
    Vec3D v3d = VEC3D(v.x, v.y, v.z);
    return v3d;
}

struct vec2 Vec2D_to_vec2(const Vec2D * v)
{
    double x = VEC2D_I(*v, 0);
    double y = VEC2D_I(*v, 1);
    return vec2(x,y);
}

Vec2D vec2_to_Vec2D(struct vec2 v)
{
    Vec2D v2d = VEC2D(v.x, v.y);
    return v2d;
}




