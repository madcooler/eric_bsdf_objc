//
//  MicrosurfaceHeight.h
//  imp
//
//  Created by Chi Wang on 30/10/16.
//
//

#ifndef EricBSDF_BasicFunctions_h
#define EricBSDF_BasicFunctions_h

#include "Vec3D.h"
#include "Vec2D.h"


#include <math.h>
//#include <iostream.h>
//#include <algorithm.h>
//#include <random.h>

#ifndef M_PI
#define M_PI			3.14159265358979323846f	/* pi */
#endif

#define SQRT_M_PI		1.77245385090551602729f /* sqrt(pi) */
#define SQRT_2			1.41421356237309504880f /* sqrt(2) */
#define INV_SQRT_M_PI	0.56418958354775628694f /* 1/sqrt(pi) */
#define INV_2_SQRT_M_PI	0.28209479177387814347f /* 0.5/sqrt(pi) */
#define INV_SQRT_2_M_PI 0.3989422804014326779f /* 1/sqrt(2*pi) */
#define INV_SQRT_2		0.7071067811865475244f /* 1/sqrt(2) */



double eric_erf(double x);
float  eric_erfinv(float x);


/*
 * A method to compute the gamma() function.
 *
 */  

double  eric_abgam (double x);

double  eric_gamma (double x);
double  eric_beta (double m, double n);
int     IsFiniteNumber(double x);
double  sign(double x);



// vec3 and vec2

struct vec2
{
    double x;
    double y;
};

// constructor of vec2
struct vec2 vec2(double x, double y);


struct vec3
{
    double x;
    double y;
    double z;
};

// constructor of vec3
struct vec3 vec3(double x, double y, double z);

float length(struct vec3 v);

struct vec3 normalize(struct vec3 v);

struct vec3 vec3_add (struct vec3 v1,struct vec3 v2);

struct vec3 negtive(struct vec3 v);

struct vec3 Vec3D_to_vec3(const Vec3D * v);

Vec3D vec3_to_Vec3D(struct vec3 v);

struct vec2 Vec2D_to_vec2(const Vec2D * v);

Vec2D vec2_to_Vec2D(struct vec2 v);

void buildOrthonormalBasis(struct vec3 * omega_1, struct vec3 * omega_2, const struct vec3 * omega_3);

#endif /* EricBSDF_BasicFunctions_h */
