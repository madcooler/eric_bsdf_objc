//
//  test_single_scattering.h
//  eric_original
//
//  Created by Administrator on 07/11/2016.
//  Copyright © 2016 Chi Wang. All rights reserved.
//

#ifndef test_single_scattering_h
#define test_single_scattering_h

void test_single_scattering()
{
    const vec3 wi = vec3(0,0,1);
    //generateRandomDirectionUp();
    const vec3 wo = generateRandomDirectionUp();
    Microsurface * m = generateRandomMicrosurface();
    
    const int N = 100000;
    // eval truncated random walk (loop because eval is stochastic)
    double V = 0;
    for(int n=0 ; n<N ; ++n)
    {
        V += m->eval(wi, wo, 1) / (double)N;
    }
    
    // eval single (loop because single of diffuse is also stochastic)
    double V_single = 0;
    for(int n=0 ; n<N ; ++n)
    {
        V_single += m->evalSingleScattering(wi, wo) / (double)N;

    }

    cout << "random␣walk␣cut␣after␣1st␣bounce␣=␣\t" << V << endl;
    cout << "single␣scattering␣=␣\t" << V_single << endl;

    delete m;
}

#endif /* test_single_scattering_h */
