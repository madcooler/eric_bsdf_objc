//
//  test_normalization_bsdf.h
//  eric_original
//
//  Created by Administrator on 07/11/2016.
//  Copyright © 2016 Chi Wang. All rights reserved.
//

#ifndef test_normalization_bsdf_h
#define test_normalization_bsdf_h

// test that \int f(wi, wo) dwo = 1
void test_normalization_bsdf()
{
    const vec3 wi = generateRandomDirectionUp();
    Microsurface * m = generateRandomMicrosurface();

    // quadrature \int f_p(wi, wo) dwo
    double value_quadrature = 0;
    for(double theta_o=0 ; theta_o < M_PI ; theta_o += 0.005)
        for(double phi_o=0 ; phi_o < 2.0*M_PI ; phi_o += 0.005)
        {
            const vec3 wo(cos(phi_o)*sin(theta_o), sin(phi_o)*sin(theta_o), cos(theta_o));
            // stochastic evaluation
            const int N = 10;
            double value_current = 0;
            for(int n=0 ; n<N ; ++n)
            {
                value_current += (double)m->eval(wi, wo) / (double) N;
            }
            value_quadrature += 0.005*0.005*abs(sin(theta_o)) * value_current;
        }
    // display
    cout << "\\int␣f_p(wi,␣wo)␣dwo␣=␣\t\t" << value_quadrature << endl;
    delete m;
}
#endif /* test_normalization_bsdf_h */
