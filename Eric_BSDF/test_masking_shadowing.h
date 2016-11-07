//
//  test_masking_shadowing.h
//  eric_original
//
//  Created by Administrator on 07/11/2016.
//  Copyright © 2016 Chi Wang. All rights reserved.
//

#ifndef test_masking_shadowing_h
#define test_masking_shadowing_h

void test_masking_shadowing()
{
    const vec3 wi = generateRandomDirectionUp();
    const vec3 wo = generateRandomDirectionUp();
    Microsurface * m = generateRandomMicrosurface();
    
    // analytic
    float G2 = 1.0f / (1.0f + m->m_microsurfaceslope->Lambda(wi) + m->m_microsurfaceslope->Lambda(wo));
    // sampling
    const int N = 100000; double V = 0;
    for(int n=0 ; n<N ; ++n)
    {
        // generate random height
        const float U1 = generateRandomNumber();
        const float h = m->m_microsurfaceheight->invC1(U1);
        // next heights
        const float U2 = generateRandomNumber();
        const float U3 = generateRandomNumber();
        const float h_wi = m->sampleHeight(wi, h, U2);
        const float h_wo = m->sampleHeight(wo, h, U3);
        // both leave microsurface
        if( (h_wi == FLT_MAX) && (h_wo == FLT_MAX))
            V += 1.0 / (double)N;
        
    }
        cout << "analytic␣=␣\t" << G2 << endl;
        cout << "stochastic␣=␣\t" << V << endl;
        delete m;
}



#endif /* test_masking_shadowing_h */
