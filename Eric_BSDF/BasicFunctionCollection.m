/*!
 
 @file      BasicFunctionCollection.m
 
 @brief     This file implements the basic function collection interface used in imp.
 
 @author    Chi Wang (chi@cgg.mff.cuni.cz), Lukas Novosad (novosad@cgg.mff.cuni.cz)
 @date      08/01/16
 
 */
#import "BasicFunctionCollection.h"


@implementation BasicFunctionCollection


+(void) GetRotationAngles
            :(const Vec3D *) incidentDirection
            :(const Vec3D *) outgoingDirection
            :(const Vec3D *) surfaceNormalDirection
            :(const Vec3D *) localNormalDirection
            :(double *) globalToLocalAngle
            :(double *) localToGlobalAngle
            ;
{
    Vec3D Si, Sr;
    Vec3D Ei, Er;
    
    vec3d_vv_cross_v( surfaceNormalDirection,
                      incidentDirection,
                      &Si
                    );
    vec3d_norm_v( &Si );
    
    vec3d_vv_cross_v( surfaceNormalDirection,
                      outgoingDirection,
                      &Sr
                    );
    vec3d_norm_v( &Sr );
    
    vec3d_vv_cross_v( localNormalDirection,
                      incidentDirection,
                      &Ei
                    );
    vec3d_norm_v( &Ei );
    
    
    vec3d_vv_cross_v( localNormalDirection,
                      outgoingDirection,
                      &Er
                    );
    vec3d_norm_v( &Er );

    
    double cos_rotation_angle_in;
    double cos_rotation_angle_out;
    
    m_ddd_clamp_d(-1.0, 1.0, vec3d_vv_dot(&Si, &Ei), &cos_rotation_angle_in);
    m_ddd_clamp_d(-1.0, 1.0, vec3d_vv_dot(&Sr, &Er), &cos_rotation_angle_out);
    
    double rotation_angle_in  = acos(cos_rotation_angle_in);
    double rotation_angle_out = acos(cos_rotation_angle_out);
    
    Vec3D temp_in, temp_out;
    vec3d_vv_cross_v(&Si, &Ei, &temp_in );
    vec3d_vv_cross_v(&Er, &Sr, &temp_out);
    
    if ( vec3d_vv_dot( &temp_in, incidentDirection)  < 0 )
        rotation_angle_in  = 2 * MATH_PI - rotation_angle_in;
    
    if ( vec3d_vv_dot( &temp_out, outgoingDirection) < 0 )
        rotation_angle_out = 2 * MATH_PI - rotation_angle_out;
    
    * globalToLocalAngle = rotation_angle_in;
    * localToGlobalAngle = rotation_angle_out;
    
}

// calculate the rotation angle from rf1 to rf2
+(void) GetRotationAngleFromTwoReferneceFrames
            :(const ART_GV           *) art_gv
            :(const ArReferenceFrame *) rf1
            :(const ArReferenceFrame *) rf2
            :(      double           *) angle
            ;
{
//    if ( arrefframe_rf_rf_d_coaxial(
//            art_gv,
//            rf1,
//            rf2,
//            3.0 DEGREES)
//        )
    {
        double cos_angle = fabs(vec3d_vv_dot(
                & rf1->c[0],
                & rf2->c[0]
                ));
        
        m_dd_clamp_d(-1,1,&cos_angle);
        
        double temp = vec3d_vv_dot(
                & rf1->c[0],
                & rf2->c[1]
                );
        if (temp > 0)
        {
            // rotate a clockwise angle
            * angle = - acos(cos_angle);
        }
        else
            // rotate a counter-clockwise angle
            * angle = acos(cos_angle);
    }
}



+(void) RotateMuellerMatrixBetweenTwoReferneceFrames
            :(const ART_GV                  *) art_gv
            :(const ArReferenceFrame        *) rf1
            :(const ArReferenceFrame        *) rf2
            :(      ArMuellerMatrixSample   *) mm
            ;
{
    // calculate rotation angle
    double rotate_angle = 0;
    [ BasicFunctionCollection  GetRotationAngleFromTwoReferneceFrames
                    :   art_gv
                    :   rf1
                    :   rf2
                    : & rotate_angle
        ];
        
    if ( isnan(rotate_angle) )
        return;
        
    // calculate rotation muelle matrix
    ArMuellerMatrixSample * rotationMatrix;
    rotationMatrix = armuellermatrixsample_alloc(art_gv);
        
    armuellermatrixsample_d_init_rotator_m(
                    art_gv,
                    rotate_angle,
                    rotationMatrix
                    );
    
    //armuellermatrixsample_mi_debugprintf(art_gv, rotationMatrix, 3);
    //armuellermatrixsample_mi_debugprintf(art_gv, mm, 3);
    // do rotation
    armuellermatrixsample_mm_mul_m(
                    art_gv,
                    rotationMatrix,
                    mm,
                    mm
                    );
    
    armuellermatrixsample_free(art_gv, rotationMatrix);
}


+(void) GetThetaAndPhiFromVec3D
        :(const Vec3D *) vec
        :(double *)      theta
        :(double *)      phi
        ;
{
    double x,y,z;
    x = VEC3D_I( *vec, 0);
    y = VEC3D_I( *vec, 1);
    z = VEC3D_I( *vec, 2);
    
    m_ddd_clamp_d(0, MATH_PI, acos(z), theta);
    
    double ABSPhi ;
    m_ddd_clamp_d(0, MATH_PI_DIV_2, ABS(atan( y / x)), & ABSPhi);
    
    if (x >= 0 && y >= 0)
        * phi = ABSPhi;
    if (x < 0 && y > 0)
        * phi = MATH_PI - ABSPhi;
    if (x < 0 && y < 0)
        * phi = MATH_PI + ABSPhi;
    if (x > 0 && y < 0)
        * phi = MATH_2_MUL_PI - ABSPhi;
    
}


+(void) GetVec3DFromThetaAndPhi
        :(const double *)   theta
        :(const double *)   phi
        :(Vec3D *)          vec
{
    
    double x = sin(*theta) * cos(*phi);
    double y = sin(*theta) * sin(*phi);
    double z = cos(*theta);
    * vec = VEC3D(x,y,z);
    
}


+(void) GetXAndYCoordinatesFromThetaAndPhi
        :(double) theta
        :(double) phi
        :(IVec2D) imageSize
        :(int *)  x
        :(int *)  y
{
    // Determine the center of the image.
    int centerX = XC(imageSize) / 2;
    int centerY = YC(imageSize) / 2;
    
    // Pick the smaller radius.
    int minRadius = 0;
    if (centerX < centerY)
        minRadius = centerX;
    else
        minRadius = centerY;
    
    // Calculate the shift via spherical coordinates (e.g. projection into x-y-plane).
    double shiftX = minRadius * sin(theta) * cos(phi);
    double shiftY = minRadius * sin(theta) * sin(phi);
    
    // Add the shift to the center values.
    double xx = centerX;
    double yy = centerY;
    xx = xx + shiftX;
    yy = yy - shiftY;
    
    // Return the pixel position via the two pointers to int.
    *x = (int)xx;
    *y = (int)yy;
    
}

+(void) GetXAndYCoordinatesFromThetaAndPhi_AveragedMapping
        :(double) theta
        :(double) phi
        :(IVec2D) imageSize
        :(int *)  x
        :(int *)  y
{
    // Determine the center of the image.
    int centerX = XC(imageSize) / 2;
    int centerY = YC(imageSize) / 2;
    
    // Pick the smaller radius.
    double minRadius = 0;
    if (centerX < centerY)
        minRadius = centerX;
    else
        minRadius = centerY;
    
    // Calculate the shift via spherical coordinates (e.g. projection into x-y-plane).
    double thetaDeg = theta * MATH_RAD_TO_DEG;
    
    if ( thetaDeg > 90)
        thetaDeg = thetaDeg - (thetaDeg - 90);
    
    minRadius       = minRadius * thetaDeg / 90;
    
    double shiftX   = minRadius * cos(phi);
    double shiftY   = minRadius * sin(phi);
    
    // Add the shift to the center values.
    double xx       = centerX;
    double yy       = centerY;
    xx              = xx + shiftX;
    yy              = yy - shiftY;
    
    // Return the pixel position via the two pointers to int.
    *x = (int)xx;
    *y = (int)yy;
    
}



+(void) GetSpecularRefractionDirection
        :(const Vec3D *)    incidentDirection
        :(const Vec3D *)    normalDirection
        :(const double)     n_i
        :(const double)     n_t
        :(Vec3D *)          refractionDirection
{
    
    //debugprintf(VEC3D_FORMAT( "%10.9f" ), VEC3D_V_PRINTF( *incidentDirection ))
    
    double eta = n_i/n_t;
    
    // Make sure the local normal is in the same side of the incident direction.
    Vec3D localNormal = *normalDirection;
    
    double cos_angle_in;
    m_ddd_clamp_d(-1.0, 1.0, vec3d_vv_dot( incidentDirection, normalDirection ), & cos_angle_in);
    
    if( cos_angle_in > 0 )
        vec3d_v_negate_v( normalDirection, &localNormal );
    
    double idotN = ABS(cos_angle_in);
    
    double k = 1 - eta*eta*(1 - idotN*idotN); //cos_t * cos_t

	if (k < 0)
	{
        // total internal reflection
		//NSLog(@"TIR");
	} 
	else
	{
		
        Vec3D temp1;
        vec3d_dv_mul_v(eta, incidentDirection, &temp1);
        
        Vec3D temp2;
        vec3d_dv_mul_v(eta*idotN-sqrt(k), &localNormal, &temp2);
        
        vec3d_vv_add_v(&temp1, &temp2, refractionDirection);
        
        vec3d_norm_v( refractionDirection );
	}
    
    
    //^^^^^^^ test ^^^^^^^
    /*
    double cost = vec3d_vv_dot( refractionDirection ,
                                normalDirection
                              );
    double cosi = idotN;
    double sint = sqrt(1-cost*cost);
    double sini = sqrt(1-cosi*cosi);
//    printf(" theta_i %f   theta_t  %f",asin(sini),asin(sint));

    if( ABS( n_i*sini - n_t*sint) > 0.0001 )
        printf(" refraction direction error ");
    
    */
    //debugprintf(VEC3D_FORMAT( "%10.9f" ), VEC3D_V_PRINTF( *incidentDirection ))
    //debugprintf(VEC3D_FORMAT( "%10.9f" ), VEC3D_V_PRINTF( *refractionDirection ))
    
}


+(void) GetSpecularReflectionDirection
        :(const Vec3D *) incidentDirection
        :(const Vec3D *) normalDirection
        :(Vec3D *)       reflectionDirection
{
    
    Vec3D localNormal = *normalDirection;

    double cos_angle_in;
    
    m_ddd_clamp_d(-1.0, 1.0, vec3d_vv_dot( incidentDirection, normalDirection ), & cos_angle_in);
    
    // Flip the angle in case the normal direction and the incoming direction do not match.
    // This is the case for incoming directions from "below" the surface; used for layered materials.
    if( cos_angle_in > 0 )
        vec3d_v_negate_v( normalDirection, &localNormal );

    vec3d_vv_reflect_v(incidentDirection, &localNormal, reflectionDirection);
    
}

+(BOOL) IsNormalReflection
        :(const Vec3D *) incomingDirection
        :(const Vec3D *) outgoingDirection
{
    double zofincoming = ZC( *incomingDirection );
    double zofoutgoing = ZC( *outgoingDirection );
    
    if( zofincoming*zofoutgoing < 0 )
        return YES;
    else return NO;
}

+(BOOL) IsNormalRefraction
        :(const Vec3D *) incomingDirection
        :(const Vec3D *) outgoingDirection
{
    double zofincoming = ZC( *incomingDirection );
    double zofoutgoing = ZC( *outgoingDirection );
    
    if( zofincoming*zofoutgoing > 0 )
        return YES;
    else return NO;
}

+(void) GetStokesValueFromArLightAlphaSample
        :(const ART_GV             *) art_gv
        :(const ArLightAlphaSample *) lightSample
        :(      double             *) valueI
        :(      double             *) valueQ
        :(      double             *) valueU
        :(      double             *) valueV
{
    ArStokesVectorSample stokes;
    arlightsample_l_to_sv(art_gv, lightSample->light, & stokes);
        
    ArSpectralSample spectralSampleI = ARSTOKESVECTORSAMPLE_SV_I(stokes, 0);
    ArSpectralSample spectralSampleQ = ARSTOKESVECTORSAMPLE_SV_I(stokes, 1);
    ArSpectralSample spectralSampleU = ARSTOKESVECTORSAMPLE_SV_I(stokes, 2);
    ArSpectralSample spectralSampleV = ARSTOKESVECTORSAMPLE_SV_I(stokes, 3);
        
    *valueI = C1_C_PRINTF(SS_C(spectralSampleI));
    *valueQ = C1_C_PRINTF(SS_C(spectralSampleQ));
    *valueU = C1_C_PRINTF(SS_C(spectralSampleU));
    *valueV = C1_C_PRINTF(SS_C(spectralSampleV));
    
    // There will be a error when free the 'stokes' here
//    arstokesvectorsample_free(art_gv, & stokes);

}

+(void) GetStokesValueFromArStokesVectorSample
        :(const ART_GV               *) art_gv
        :(const ArStokesVectorSample *) stokes
        :(      double               *) valueI
        :(      double               *) valueQ
        :(      double               *) valueU
        :(      double               *) valueV
{
    //ArStokesVectorSample stokes;
    //arlightsample_l_to_sv(art_gv, lightSample->light, & stokes);
        
    ArSpectralSample spectralSampleI = ARSTOKESVECTORSAMPLE_SV_I(*stokes, 0);
    ArSpectralSample spectralSampleQ = ARSTOKESVECTORSAMPLE_SV_I(*stokes, 1);
    ArSpectralSample spectralSampleU = ARSTOKESVECTORSAMPLE_SV_I(*stokes, 2);
    ArSpectralSample spectralSampleV = ARSTOKESVECTORSAMPLE_SV_I(*stokes, 3);
        
    *valueI = C1_C_PRINTF(SS_C(spectralSampleI));
    *valueQ = C1_C_PRINTF(SS_C(spectralSampleQ));
    *valueU = C1_C_PRINTF(SS_C(spectralSampleU));
    *valueV = C1_C_PRINTF(SS_C(spectralSampleV));
    
    // There will be a error when free the 'stokes' here
//    arstokesvectorsample_free(art_gv, & stokes);

}

+(void) GetIntensityValueFromArLightAlphaSample
        :(const ART_GV             *) art_gv
        :(const ArLightAlphaSample *) lightSample
        :(      double             *) valueI
{
    ArStokesVectorSample stokes;
    arlightsample_l_to_sv(art_gv, lightSample->light, & stokes);
        
    ArSpectralSample spectralSampleI = ARSTOKESVECTORSAMPLE_SV_I(stokes, 0);
    
    *valueI = C1_C_PRINTF(SS_C(spectralSampleI));
}


+(double) DistanceBetweenPnt3D
        :(const Pnt3D *) p1
        :(const Pnt3D *) p2
{
    double x = PNT3D_I(*p1, 0) - PNT3D_I(*p2, 0);
    double y = PNT3D_I(*p1, 1) - PNT3D_I(*p2, 1);
    double z = PNT3D_I(*p1, 2) - PNT3D_I(*p2, 2);
    
    return sqrt( x*x + y*y + z*z);
}


+(void) ApplyPolarisationFilter
        :(const ART_GV         *)     art_gv
        :(const double          )     angle
        :(const double          )     strength
        :(ArStokesVectorSample *)     sv
{
    ArMuellerMatrixSample * mm;
    mm = armuellermatrixsample_alloc(art_gv);
        
    armuellermatrixsample_dd_init_linear_polariser_m(
                art_gv,
                angle,
                strength,
                mm);
    ArStokesVectorSample * sv_out = arstokesvectorsample_alloc(art_gv);
    
    //armuellermatrixsample_mi_debugprintf(art_gv, mm, 3);
    arstokesvectorsample_sv_mm_mul_sv(art_gv, sv, mm, sv_out);
    //arstokesvectorsample_sv_debugprintf(art_gv, sv_out);
    
    
    *sv = *sv_out;
    arstokesvectorsample_free(art_gv, sv_out);
    armuellermatrixsample_free(art_gv, mm);
    //arstokesvectorsample_sv_debugprintf(art_gv, sv);
}

+(void) NormaliseStokesVector
        :(const ART_GV         *)     art_gv
        :(ArStokesVectorSample *)     sv
{
    double I = ARSTOKESVECTORSAMPLE_SV_I(*sv, 0).c.x[0];
    double Q = ARSTOKESVECTORSAMPLE_SV_I(*sv, 1).c.x[0];
    double U = ARSTOKESVECTORSAMPLE_SV_I(*sv, 2).c.x[0];
    double V = ARSTOKESVECTORSAMPLE_SV_I(*sv, 3).c.x[0];
    
    arstokesvectorsample_dddd_init_sv(art_gv, 1, Q/I, U/I, V/I, sv);
}


#define INIT_MM_ID(__mm,__i,__d) \
    ss_d_init_s( art_gv, (__d), ARMUELLER_S_M_I(__mm,__i) );

// retarder
+(void) ApplyCompensatorFilter
        :(const ART_GV         *) art_gv
        :(const double)           phase_shift
        :(const double)           angle
        :(ArStokesVectorSample *) sv
{
    ArMuellerMatrixSample * mm;
    mm = armuellermatrixsample_alloc(art_gv);
    

    double cos_phase = cos(phase_shift);
    double sin_phase = sin(phase_shift);
    
    if(phase_shift == M_PI / 2)
    {
        cos_phase = 0;
        sin_phase = 1;
    }
    
    INIT_MM_ID( *mm,  0, 1.0 );
    INIT_MM_ID( *mm,  1, 0.0 );
    INIT_MM_ID( *mm,  2, 0.0 );
    INIT_MM_ID( *mm,  3, 0.0 );

    INIT_MM_ID( *mm,  4, 0.0 );
    INIT_MM_ID( *mm,  5, cos( 2 * angle) * cos( 2 * angle ) + cos_phase * sin( 2 * angle) * sin( 2 * angle ) );
    INIT_MM_ID( *mm,  6, cos( 2 * angle) * sin( 2 * angle ) - cos_phase * cos( 2 * angle) * sin( 2 * angle )  );
    INIT_MM_ID( *mm,  7, sin( 2 * angle ) * sin_phase );

    INIT_MM_ID( *mm,  8, 0 );
    INIT_MM_ID( *mm,  9, cos( 2 * angle) * sin( 2 * angle ) - cos_phase * cos( 2 * angle) * sin( 2 * angle ) );
    INIT_MM_ID( *mm, 10, cos_phase * cos( 2 * angle) * cos( 2*angle ) + sin( 2 * angle) * sin( 2 * angle ) );
    INIT_MM_ID( *mm, 11, -cos( 2 * angle ) * sin_phase );

    INIT_MM_ID( *mm, 12, 0.0 );
    INIT_MM_ID( *mm, 13, -sin( 2 * angle ) * sin_phase );
    INIT_MM_ID( *mm, 14, cos( 2 * angle ) * sin_phase );
    INIT_MM_ID( *mm, 15, cos_phase );
    
    ArStokesVectorSample * sv_out = arstokesvectorsample_alloc(art_gv);
    
    //armuellermatrixsample_mi_debugprintf(art_gv, mm, 3);
    arstokesvectorsample_sv_mm_mul_sv(art_gv, sv, mm, sv_out);
    //arstokesvectorsample_sv_debugprintf(art_gv, sv_out);
    
    *sv = *sv_out;
    arstokesvectorsample_free(art_gv, sv_out);
    armuellermatrixsample_free(art_gv, mm);
    //arstokesvectorsample_sv_debugprintf(art_gv, sv);

}

+(void) GetNormalisedStokesVectorAndDOLP
        :(const ART_GV              *) art_gv
        :(const ArLightAlphaSample  *) light
        :(      double              *) dolp
        :(ArStokesVectorSample      *) sv
{

    double I,Q,U,V;
    
    [self GetStokesValueFromArLightAlphaSample
                :   art_gv
                :   light
                : & I
                : & Q
                : & U
                : & V
                ];
    *dolp = sqrt( Q*Q + U*U )/I;
    
    Q = Q/I;
    U = U/I;
    V = V/I;
    I = 1;
    
    arstokesvectorsample_dddd_init_sv( art_gv, I, Q, U, V, sv );
}

+(void) ApplyPolarisationAnalyser
        :(const ART_GV              *) art_gv
        :(const ArLightAlphaSample  *) light
        :(      double               ) angle
        :(      double              *) intensity
{

    ArStokesVectorSample * sv = arstokesvectorsample_alloc(art_gv);
    
    arlightalphasample_l_to_sv(art_gv, light, sv);
    
    [self ApplyPolarisationFilter
        : art_gv
        : angle
        : 1
        : sv
        ];
    
    *intensity = ARSTOKESVECTORSAMPLE_SV_I(*sv, 0).c.x[0];
}


+(void) ApplyPolarisationAnalyserForStokes
        :(const ART_GV              *) art_gv
        :(const ArStokesVectorSample*) sv
        :(      double               ) angle
        :(      double              *) intensity
{
    ArStokesVectorSample * sv_temp = arstokesvectorsample_alloc(art_gv);

    arstokesvectorsample_sv_init_sv(art_gv, sv, sv_temp);
    
    [self ApplyPolarisationFilter
        : art_gv
        : angle
        : 1
        : sv_temp
        ];
    
    *intensity = ARSTOKESVECTORSAMPLE_SV_I(*sv_temp, 0).c.x[0];
}


+(void) GetNormalisedStokesVectorSampleAndDOLP
        :(const ART_GV              *) art_gv
        :(const ArStokesVectorSample*) light
        :(      double              *) dolp
        :(ArStokesVectorSample      *) sv
{
    ArSpectralSample spectralSampleI = ARSTOKESVECTORSAMPLE_SV_I(*light, 0);
    ArSpectralSample spectralSampleQ = ARSTOKESVECTORSAMPLE_SV_I(*light, 1);
    ArSpectralSample spectralSampleU = ARSTOKESVECTORSAMPLE_SV_I(*light, 2);
    ArSpectralSample spectralSampleV = ARSTOKESVECTORSAMPLE_SV_I(*light, 3);
        
    double I = C1_C_PRINTF(SS_C(spectralSampleI));
    double Q = C1_C_PRINTF(SS_C(spectralSampleQ));
    double U = C1_C_PRINTF(SS_C(spectralSampleU));
    double V = C1_C_PRINTF(SS_C(spectralSampleV));
    
    *dolp = sqrt( Q*Q + U*U )/I;
    
    Q = Q/I;
    U = U/I;
    V = V/I;
    I = 1;
    
    arstokesvectorsample_dddd_init_sv( art_gv, I, Q, U, V, sv );
    
}


// uniform hemisphere sampling
+(void) SampleDiffuseScatteringDirection
    : ( ArcObject <ArpRandomGenerator>  *)  randomGenerator
    : ( Vec3D                           *)  outgoing_dir_local
    : ( double                          *)  PDF
{
    double u1 , u2 ;
    
	[ randomGenerator getValuesFromNewSequences :& u1 :& u2 ];
    
    // from PBRT second edition 13.6.1
    double z    = u1;
    double r    = sqrt(1 - z*z);
    double phi  = 2 * M_PI * u2;
    double x    = r * cosf(phi);
    double y    = r * sinf(phi);
    
	* outgoing_dir_local = VEC3D( x, y, z);
    * PDF =  MATH_1_DIV_2_PI;
}


@end
