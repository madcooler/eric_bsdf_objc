//
//  basic.h
//  eric_original
//
//  Created by Administrator on 07/11/2016.
//  Copyright © 2016 Chi Wang. All rights reserved.
//

#ifndef basic_h
#define basic_h


#include <random>
#include "Vector.h"
#include <iostream>
#include "MicrosurfaceScattering.h"

using namespace std;

//static std::random_device rd;
//static std::mt19937 gen(rd);
//static std::uniform_real_distribution<> dis(0, 1);


inline double getUniformRandomNumber ( double min , double max )
{
	std::random_device rd;   // non-deterministic generator
	std::mt19937 gen ( rd () );
	std::uniform_real_distribution<> distr ( min , max );
	return distr ( gen );
}

float generateRandomNumber()
{
    return getUniformRandomNumber(0,1);
}

vec3 generateRandomDirection()
{
    const float theta = M_PI * generateRandomNumber();
    const float phi = 2.0f * M_PI * generateRandomNumber();
    const vec3 w(cosf(phi)*sinf(theta), sinf(phi)*sinf(theta), cosf(theta));
    cout << "w␣=␣(" << w.x << ",␣" << w.y << ",␣" << w.z << ")" << endl;
    cout << endl;
    return w;
}

vec3 generateRandomDirectionUp()
{
    const float theta = 0.5f * M_PI * generateRandomNumber();
    const float phi = 2.0f * M_PI * generateRandomNumber();
    const vec3 w(cosf(phi)*sinf(theta), sinf(phi)*sinf(theta), cosf(theta));
    cout << "w␣=␣(" << w.x << ",␣" << w.y << ",␣" << w.z << ")" << endl;
    cout << endl;
    return w;
    
}

Microsurface * generateRandomMicrosurface()
{
    // height
    const bool height_uniform = generateRandomNumber() > 0.5f;
    
    // roughness
    const float alpha_x = generateRandomNumber() + 0.1f;
    const float alpha_y = generateRandomNumber() + 0.1f;
    // distribution (Beckmann or GGX)
    const bool Beckmann = generateRandomNumber() > 0.5f;
    
    // material
    
    Microsurface * m;
    
    const int material = (int)floor(generateRandomNumber()*3.0);
    switch(material)
    {
        case 0: m = new MicrosurfaceDiffuse(height_uniform, Beckmann, alpha_x, alpha_y);
        cout << "Material:␣Diffuse" << endl;
            break;

        case 1: m = new MicrosurfaceConductor(height_uniform, Beckmann, alpha_x, alpha_y);
        cout << "Material:␣Conductor" << endl;
            break;
            
        default: m = new MicrosurfaceDielectric(height_uniform, Beckmann, alpha_x, alpha_y);
        
        cout << "Material:␣Dielectric" << endl;
        break;
        
        
    }
    
    cout<<"height␣distribution:␣" << ((height_uniform)?"Uniform":"Gaussian") << endl; cout<<"slope␣distribution:␣" << ((Beckmann)?"Beckmann":"GGX") << endl;
        cout<<"alpha_x␣=␣" << m->m_microsurfaceslope->m_alpha_x << endl;
        cout<<"alpha_y␣=␣" << m->m_microsurfaceslope->m_alpha_y << endl;

    
    return m;
}

#endif /* basic_h */
