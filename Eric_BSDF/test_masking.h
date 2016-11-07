//
//  test_masking.h
//  eric_original
//
//  Created by Administrator on 07/11/2016.
//  Copyright © 2016 Chi Wang. All rights reserved.
//

#ifndef test_masking_h
#define test_masking_h


void test_masking()
{

    const vec3 wo = generateRandomDirectionUp();
    Microsurface * m = generateRandomMicrosurface();
    // analytic
    float G1 = 1.0f / (1.0f + m->m_microsurfaceslope->Lambda(wo));
    // sampling
    const int N = 100000; double V = 0;
    for(int n=0 ; n<N ; ++n)
    {
        // generate random height
        const float U1 = generateRandomNumber();
        const float h = m->m_microsurfaceheight->invC1(U1);
        // next height
        const float U2 = generateRandomNumber();
        const float h_wo = m->sampleHeight(wo, h, U2);
        // leave microsurface
        if(h_wo == FLT_MAX)
        V += 1.0 / (double)N;

    }
    cout << "analytic␣=␣\t" << G1 << endl;
    cout << "stochastic␣=␣\t" << V << endl;
    
    delete m;
}


#endif /* test_masking_h */
