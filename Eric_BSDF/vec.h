//
//  Header.h
//  imp
//
//  Created by Administrator on 30/10/16.
//
//

#include <math.h>


#ifndef VEC_H
#define VEC_H

#define vec2(_x,_y,_z) initVec2(_x,_y);
#define vec3(_x,_y,_z) initVec3(_x,_y,_z);

struct vec2
{
    float x;
    float y;
}

vec2 initVec2(float x, float y)
{
    vec2 temp;
    temp.x = x;
    temp.y = y;
    return temp;
}


struct vec3
{
    float x;
    float y;
    float z;
}

vec3 initVec3(float x, float y, float z)
{
    vec3 temp;
    temp.x = x;
    temp.y = y;
    temo.z = z;
    return temp;
}

float length(vec3 v)
{
    return sqrtf(x*x+y*y+z*z);
}

vec3 normalize(vec3 v)
{
    float l = length(v);
    
    return vec3(v.x/l,v.y/l,v.z/l);
}




#endif /* VEC_H */
