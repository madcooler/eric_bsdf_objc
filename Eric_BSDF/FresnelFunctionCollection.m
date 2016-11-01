/*!
 
 @file      FresnelFunctionCollection.m
 
 @brief     This file implements the Fresnel function collection interface used in imp.
 
 @author    Chi Wang (chi@cgg.mff.cuni.cz), Lukas Novosad (novosad@cgg.mff.cuni.cz)
 @date      13/01/16
 
 */
#import "FresnelFunctionCollection.h"
// imp project imports
#import "BasicFunctionCollection.h"

@implementation FresnelFunctionCollection


+(double) Fresnel_Refraction
            :(double) theta
            :(double) n1
            :(double) k1
            :(double) n2
            :(double) k2
            :(double *) phase
            :(double *) Ts
            :(double *) Tp
            :(double *) TsTp
            :(double *) factor
            ;
{
    double nab = n1 / n2;
    double theta_t = asin( n1*sin(theta) / n2 );
    
    *Tp = ( 2 * nab * cos(theta) / ( cos(theta) + nab * cos(theta_t))) *
          ( 2 * nab * cos(theta) / ( cos(theta) + nab * cos(theta_t)));
    *Ts = ( 2 * nab * cos(theta) / ( nab * cos(theta) + cos(theta_t))) *
          ( 2 * nab * cos(theta) / ( nab * cos(theta) + cos(theta_t)));
    phase = 0;
    *TsTp = ( 4 * nab * nab * cos(theta) * cos(theta) ) /
            ( ( cos(theta) + nab * cos(theta_t) ) * ( nab * cos(theta) + cos(theta_t) ) );
    
    double nba = n2 / n1;
    *factor = /* nba * nba * */ nba * ( cos(theta_t)/cos(theta) );
    
    return *factor * ( *Ts + *Tp ) / 2;
}


+(double) Fresnel_Reflection
            :(double) theta
            :(double) n1
            :(double) k1
            :(double) n2
            :(double) k2
            :(double *) rs
            :(double *) rp
            :(double *) phaseS
            :(double *) phaseP
            :(double *) Fs
            :(double *) Fp
            ;
{
    double temp = ( n2*n2 - k2*k2 - n1*n1*sin(theta)*sin(theta) );
    double p2   = ( sqrt(temp*temp+4*n2*n2*k2*k2) + temp ) / 2;
    double q2   = ( sqrt(temp*temp+4*n2*n2*k2*k2) - temp ) / 2;
    double p    =   sqrt(p2);
    double q    =   sqrt(q2);
    
    double brewsterAngle = atanf( n2 / n1 );
    
    
    *Fs = ( ( n1 * cos(theta) - p ) * ( n1 * cos(theta) - p ) + q2 ) /
          ( ( n1 * cos(theta) + p ) * ( n1 * cos(theta) + p ) + q2);
    *Fp = ( *Fs ) * ( ( p - n1 * sin(theta) * tan(theta) ) *
                    (   p - n1 * sin(theta) * tan(theta) ) + q2 )
                  / ( ( p + n1 * sin(theta) * tan(theta) )
                    * ( p + n1 * sin(theta) * tan(theta) ) + q2 );
    
    // phase shift for dielectric materials
    // please refer to section 8.2 of < Polarised Light, second edition >,  Dennis Goldstein
    if ( n1 > n2 && k2 == 0 )
    {
        double criticalAngle = asin ( n2 / n1 );
        
        if ( theta < criticalAngle )
            * phaseS  = 0 ;
        else * phaseS = 2 * atan( sqrt( sin( theta ) * sin( theta )
                                    - sin( criticalAngle ) * sin( criticalAngle )
                                      )
                                  / cos(criticalAngle)
                                );
        
        
        if ( theta < brewsterAngle )
            * phaseP = 180 * MATH_DEG_TO_RAD;
        else
            if ( theta < criticalAngle )
                * phaseP  = 0 * MATH_DEG_TO_RAD;
            else * phaseP = 2 * atan( sqrt( sin(theta) * sin(theta)
                                            - sin( criticalAngle ) * sin( criticalAngle )
                                          )
                                      / ( cos(criticalAngle) * sin( criticalAngle ) * sin( criticalAngle )
                                        )
                                    );
    }
    
    if ( n1 < n2 && k2 == 0 )
    {
        * phaseS = 180 * MATH_DEG_TO_RAD;
        * phaseP = ( theta < brewsterAngle ) ? 0
                                             : 180 * MATH_DEG_TO_RAD;
    }
    
    // metal material phase shift
    if ( k2 != 0 )
    {
        * phaseS = atanf( 2 * n1 * q * cos(theta) / ( n1 * n1 * cos(theta) * cos(theta)- p2 - q2 ) );
        * phaseP = atanf(
                        (
                           -2 * n1 * q * cos(theta)
                        * ( p2 + q2 - n1 * n1 * sin(theta) * sin(theta)) )
                        / (
                            ( n2 * n2 + k2 * k2 ) * ( n2 * n2 + k2 * k2) * cos(theta) * cos(theta) - n1 * n1 * ( p2 + q2 )
                          )
                        );
    }
    
    
    *rs = sqrt( ABS (*Fs) );
    *rp = sqrt( ABS (*Fp) );
    
    return ( *Fs + *Fp ) / 2;
}


+(void) imp_fresnel_plain_reflective_attenuation_realvalued_IOR
            :(const ART_GV *)       art_gv
            :(const Vec3D *)        incidentDirection
            :(const Vec3D *)        localNormalDirection
            :(const Vec3D *)        sampledReflectionDirection
            :(const double)         n
            :(const double)         k
            :(ArAttenuation *)      attenuation
            ;
{
    double  attenuation_perpendicular;
    double  attenuation_parallel;
    
    double  cosTheta = ABS( vec3d_vv_dot( sampledReflectionDirection, localNormalDirection));
    
    fresnel_plain_ddd_attenuation_dd(
                                     cosTheta,
                                     n,
                                     k,
                                     & attenuation_perpendicular,
                                     & attenuation_parallel
                                     );
    
    double fresnelterm = 0.5 * attenuation_perpendicular + 0.5 * attenuation_parallel;
    //fresnelterm = 1;
    //    printf("TS: art fresnel: %f\n", fresnelterm);
    
    arattenuation_d_init_a(
                           art_gv,
                           fresnelterm,
                           attenuation
                           );
    
    
}


+(void) imp_fresnel_polarising_reflective_attenuation_realvalued_IOR
            :(const ART_GV *)       art_gv
            :(const Vec3D *)        incidentDirection
            :(const Vec3D *)        localNormalDirection
            :(const Vec3D *)        surfaceNormalDirection
            :(const Vec3D *)        sampledReflectionDirection
            :(const double)         n
            :(const double)         k
            :(ArAttenuation *)      attenuation
            ;
{
    //ART__CODE_IS_WORK_IN_PROGRESS__EXIT_WITH_ERROR
    
    ArSpectrum  *attenuationColourA = spc_alloc(art_gv);
    ArSpectrum  *attenuationColourB = spc_alloc(art_gv);
    ArSpectrum  *attenuationColourC = spc_alloc(art_gv);
    ArSpectrum  *attenuationColourS = spc_alloc(art_gv);
    ArSpectrum  *attenuationColourT = spc_alloc(art_gv);
    
#ifdef FOUNDATION_ASSERTIONS
    // We need to initialize attenuationColourS to meaningful values
    // since it is checked in its full extent when accessed
    // by col_ci before it is completely initialized.                      (ip)
    spc_d_init_s (art_gv, 0.0, attenuationColourS);
#endif
    double cosTheta = ABS(vec3d_vv_dot(incidentDirection, localNormalDirection));
    
    double  attenuation_perpendicular, attenuation_parallel;
    double  retardance_perpendicular, retardance_parallel;
    double  retardance_total;
    
    
    fresnel_ddd_attenuation_dddd(
                                 cosTheta,
                                 n,
                                 k,
                                 & attenuation_perpendicular,
                                 & attenuation_parallel,
                                 & retardance_perpendicular,
                                 & retardance_parallel
                                 );
    
    retardance_total = retardance_perpendicular - retardance_parallel;
    
    double  sqrt_reflectance_term =
    sqrt( attenuation_perpendicular * attenuation_parallel );
    
    spc_d_init_s(art_gv,
                 0.5 * attenuation_perpendicular
                 + 0.5 * attenuation_parallel,
                 attenuationColourA
                 );
    
    spc_d_init_s(art_gv,
                 0.5 * attenuation_perpendicular
                 - 0.5 * attenuation_parallel,
                 attenuationColourB
                 );
    
    spc_d_init_s(art_gv,
                 cos(retardance_total)
                 * sqrt_reflectance_term,
                 attenuationColourC
                 );
    
    spc_d_init_s(art_gv,
                 sin(retardance_total)
                 * sqrt_reflectance_term,
                 attenuationColourS
                 );
    
    spc_ds_mul_s(art_gv,
                 -1,
                 attenuationColourS,
                 attenuationColourT
                 );
    
    ArReferenceFrame refframe_exit, refframe_entry;
    
    ArPathDirection pathDirection = arpathdirection_from_light;
    
    arreframe_vv_pd_init_rf_rf(
                                  art_gv,
                                  incidentDirection,
                                  sampledReflectionDirection,
                                  pathDirection,
                                & refframe_entry,
                                & refframe_exit
                                );
    
    ArMuellerMatrix * fresnelMM = armuellermatrix_alloc( art_gv );
    
    spc_s_init_s( art_gv, attenuationColourA, ARMUELLER_M_I( *fresnelMM, 0) );
    spc_s_init_s( art_gv, attenuationColourB, ARMUELLER_M_I( *fresnelMM, 1) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 2) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 3) );
    spc_s_init_s( art_gv, attenuationColourB, ARMUELLER_M_I( *fresnelMM, 4) );
    spc_s_init_s( art_gv, attenuationColourA, ARMUELLER_M_I( *fresnelMM, 5) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 6) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 7) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 8) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 9) );
    spc_s_init_s( art_gv, attenuationColourC, ARMUELLER_M_I( *fresnelMM, 10) );
    spc_s_init_s( art_gv, attenuationColourS, ARMUELLER_M_I( *fresnelMM, 11) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 12) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 13) );
    spc_s_init_s( art_gv, attenuationColourT, ARMUELLER_M_I( *fresnelMM, 14) );
    spc_s_init_s( art_gv, attenuationColourC, ARMUELLER_M_I( *fresnelMM, 15) );
    
    
    arattenuation_mm_rr_init_polarising_a(
                                            art_gv,
                                            fresnelMM,
                                          & refframe_entry,
                                          & refframe_exit,
                                            attenuation
                                          );
    
    armuellermatrix_free( art_gv, fresnelMM );
    
    
    spc_free( art_gv, attenuationColourA );
    spc_free( art_gv, attenuationColourB );
    spc_free( art_gv, attenuationColourS );
    spc_free( art_gv, attenuationColourT );
    spc_free( art_gv, attenuationColourC );
    
}


+(void) imp_fresnel_plain_refractive_attenuation
            :(const ART_GV *)       art_gv
            :(const Vec3D *)        incidentDirection
            :(const Vec3D *)        localNormalDirection
            :(const Vec3D *)        sampledReflectionDirection
            :(const double)         n
            :(const double)         k
            :(ArAttenuation *)      attenuation
            ;
{
    double attenuation_perpendicular, attenuation_parallel;
    
    double cosTheta = ABS(vec3d_vv_dot(incidentDirection, localNormalDirection));

    fresnel_plain_ddd_attenuation_dd(
                                    cosTheta,
                                    n,
                                    k,
                                    & attenuation_perpendicular,
                                    & attenuation_parallel
                                    );
        
    attenuation_perpendicular =     1.0 - attenuation_perpendicular;
    attenuation_parallel      =     1.0 - attenuation_parallel;

    double fresnelTerm = 0.5 * attenuation_perpendicular + 0.5 * attenuation_parallel;
    //fresnelterm = 1;
    //    printf("TS: art fresnel: %f\n", fresnelterm);
    
    arattenuation_d_init_a(
                           art_gv,
                           fresnelTerm,
                           attenuation
                           );

}


+(void) imp_fresnel_polarising_refractive_attenuation
            :(const ART_GV *)   art_gv
            :(const Vec3D *)    incidentDirection
            :(const Vec3D *)    localNormalDirection
            :(const Vec3D *)    surfaceNormalDirection
            :(const Vec3D *)    sampledReflectionDirection
            :(const Vec3D *)    sampledRefractionDirection
            :(const double)     n
            :(const double)     k
            :(ArAttenuation *)  attenuation
            ;
{
    ArSpectrum  *attenuationColourA = spc_alloc(art_gv);
    ArSpectrum  *attenuationColourB = spc_alloc(art_gv);
    ArSpectrum  *attenuationColourC = spc_alloc(art_gv);
    ArSpectrum  *attenuationColourS = spc_alloc(art_gv);
    ArSpectrum  *attenuationColourT = spc_alloc(art_gv);
    
#ifdef FOUNDATION_ASSERTIONS
    // We need to initialize attenuationColourS to meaningful values
    // since it is checked in its full extent when accessed
    // by spc_si before it is completely initialized.                      (ip)
    spc_d_init_s (art_gv, 0.0, attenuationColourS);
#endif
    double cosTheta = ABS(vec3d_vv_dot(sampledReflectionDirection, localNormalDirection));
    
    double  attenuation_perpendicular, attenuation_parallel;
    double  retardance_perpendicular, retardance_parallel;
    double  retardance_total;
    
    
    fresnel_ddd_attenuation_dddd(
          cosTheta,
          n,
          k,
        & attenuation_perpendicular,
        & attenuation_parallel,
        & retardance_perpendicular,
        & retardance_parallel
        );
    
    attenuation_perpendicular = 1.0 - attenuation_perpendicular;
    attenuation_parallel      = 1.0 - attenuation_parallel;
        
    // TODO: Refraction retardance
    // For refraction, phase retardance should be 0.
        
    retardance_total = 0;//retardance_perpendicular - retardance_parallel;

    spc_d_init_s(art_gv,
                 0.5 * attenuation_perpendicular + 0.5 * attenuation_parallel,
                 attenuationColourA
                 );
    
    spc_d_init_s(art_gv,
                 0.5 * attenuation_perpendicular - 0.5 * attenuation_parallel,
                 attenuationColourB
                 );
        
    double sqrt_reflectance_term = sqrt( attenuation_perpendicular * attenuation_parallel );
    
    spc_d_init_s(art_gv,
                 cos(retardance_total)
                 * sqrt_reflectance_term,
                 attenuationColourC
                 );
    
    spc_d_init_s(art_gv,
                 sin(retardance_total)
                 * sqrt_reflectance_term,
                 attenuationColourS
                 );
    
    spc_ds_mul_s(art_gv,
                 -1,
                 attenuationColourS,
                 attenuationColourT
                 );
    
    ArReferenceFrame  refframe_entry, refframe_exit;
    
    // the reflection function also works for refractions
    
    ArPathDirection pathDirection=arpathdirection_from_light;
    
    arreframe_vvv_pd_init_rf_rf(
                                art_gv,
                                incidentDirection,
                                surfaceNormalDirection,
                                sampledRefractionDirection,
                                pathDirection,
                                & refframe_entry,
                                & refframe_exit
                                );
    
    ArMuellerMatrix * fresnelMM = armuellermatrix_alloc( art_gv );
    
    spc_s_init_s( art_gv, attenuationColourA, ARMUELLER_M_I( *fresnelMM, 0) );
    spc_s_init_s( art_gv, attenuationColourB, ARMUELLER_M_I( *fresnelMM, 1) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 2) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 3) );
    spc_s_init_s( art_gv, attenuationColourB, ARMUELLER_M_I( *fresnelMM, 4) );
    spc_s_init_s( art_gv, attenuationColourA, ARMUELLER_M_I( *fresnelMM, 5) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 6) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 7) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 8) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 9) );
    spc_s_init_s( art_gv, attenuationColourC, ARMUELLER_M_I( *fresnelMM, 10) );
    spc_s_init_s( art_gv, attenuationColourS, ARMUELLER_M_I( *fresnelMM, 11) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 12) );
    spc_d_init_s( art_gv, 0.0               , ARMUELLER_M_I( *fresnelMM, 13) );
    spc_s_init_s( art_gv, attenuationColourT, ARMUELLER_M_I( *fresnelMM, 14) );
    spc_s_init_s( art_gv, attenuationColourC, ARMUELLER_M_I( *fresnelMM, 15) );
    
    
    
    arattenuation_mm_rr_init_polarising_a(
                                          art_gv,
                                          fresnelMM,
                                          & refframe_entry,
                                          & refframe_exit,
                                          attenuation
                                          );
    
    armuellermatrix_free( art_gv, fresnelMM );
    
    
    spc_free( art_gv, attenuationColourA );
    spc_free( art_gv, attenuationColourB );
    spc_free( art_gv, attenuationColourS );
    spc_free( art_gv, attenuationColourT );
    spc_free( art_gv, attenuationColourC );
    
}


@end