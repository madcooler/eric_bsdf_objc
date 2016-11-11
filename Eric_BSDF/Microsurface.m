//
//  MicrosurfaceHeight.h
//  imp
//
//  Created by Chi Wang on 30/10/16.
//
//

#import "Microsurface.h"


@implementation Microsurface

-(id) init
    :(BOOL)     height_uniform
    :(BOOL)     slope_beckmann
    :(float)    alpha_x
    :(float)    alpha_y
{
    self = [super init];
    
    if(height_uniform)
        m_microsurfaceheight = [ALLOC_INIT_OBJECT (MicrosurfaceHeightUniform)];
    else
        m_microsurfaceheight = [ALLOC_INIT_OBJECT (MicrosurfaceHeightGaussian)];
    
    if(slope_beckmann)
        m_microsurfaceslope = [ ALLOC_INIT_OBJECT(MicrosurfaceSlopeBeckmann)
                                : alpha_x
                                : alpha_y
                                ];
    else
        m_microsurfaceslope = [ALLOC_INIT_OBJECT(MicrosurfaceSlopeGGX)
                                : alpha_x
                                : alpha_y
                                ];

//    randomGenerator = ARCRANDOMGENERATOR_NEW(
//                        RANDOM_SEQUENCE,
//                        100000000,
//                        ART_GLOBAL_REPORTER
//                        );

    return self;
    
}

-(void) init_material
    :(float)    n
    :(float)    k
{
    m_eta = n;
    m_k   = k;
}

-(void) dealloc
{
    m_microsurfaceslope     = nil;
    m_microsurfaceheight    = nil;
    [super dealloc];
}

-(void) initRandomNumberGenerator
    :(ArcObject <ArpRandomGenerator>  *) randomGen
{
    randomGenerator = randomGen;
}

-(double) generateRandomNumber
{
    double u1;
    u1 = [ randomGenerator valueFromNewSequence];
//    u1 = 0.4;
    return u1;
}

-(float) evalBSDF
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    : (const int )    scatteringOrder
{
    struct vec3 wo = Vec3D_to_vec3(wo_art);
    wo = normalize(wo);
    double cosine_weight = wo.z;
    return [self eval:wi_art :wo_art :scatteringOrder]/cosine_weight;
}

-(float) eval
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    : (const int )    scatteringOrder
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    struct vec3 wo = Vec3D_to_vec3(wo_art);
    
    if(wo.z < 0)
		return 0;
	// init
	struct vec3 wr = negtive(wi);
	float hr = 1.0f + [ m_microsurfaceheight invC1:0.999f];

	float sum = 0;
	
    Vec3D wr_art = vec3_to_Vec3D( wr );
    
	// random walk
	int current_scatteringOrder = 0;	
	while(scatteringOrder==0 || current_scatteringOrder <= scatteringOrder)
	{
		// next height
		float U = [ self generateRandomNumber ];
        
		hr = [self sampleHeight: & wr_art: hr : U];

		// leave the microsurface?
		if( hr == FLT_MAX )
			break;
		else
			current_scatteringOrder++;

		// next event estimation
        
        Vec3D neg_wr_art = wr_art;
        vec3d_negate_v(&neg_wr_art);
        
		float phasefunction = [self evalPhaseFunction
                                : & neg_wr_art
                                :   wo_art
                               ];
		float shadowing = [self G_1 :wo_art : hr];
        
		float I = phasefunction * shadowing;

		if ( IsFiniteNumber(I) && (scatteringOrder==0 || current_scatteringOrder==scatteringOrder) )
			sum += I;
		
		// next direction
        
        wr_art = [self samplePhaseFunction: & neg_wr_art];
        
		wr = Vec3D_to_vec3(& wr_art);
			
		// if NaN (should not happen, just in case)
		if( (hr != hr) || (wr.z != wr.z) ) 
			return 0.0f;
	}

	return sum;
}

-(Vec3D) sample
    : (const Vec3D *) wi_art
    : (      int    ) scatteringOrder
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    
    // init
	struct vec3 wr = negtive(wi) ;
	float hr = 1.0f + [ m_microsurfaceheight invC1: 0.999f];
	
    Vec3D wr_art = vec3_to_Vec3D( wr );
    
    
	// random walk
	scatteringOrder = 0;	
	while( true )
	{
		// next height
		float U = [ self generateRandomNumber ];
        
        
		hr = [self sampleHeight: & wr_art: hr : U];
		

		// leave the microsurface?
		if( hr == FLT_MAX )
			break;
		else
			scatteringOrder++;

		// next direction
		//wr = samplePhaseFunction(-wr);
        
        Vec3D neg_wr_art = wr_art;
        vec3d_negate_v(&neg_wr_art);
        wr_art = [self samplePhaseFunction: & neg_wr_art];
        wr = Vec3D_to_vec3(&wr_art);
        
        
		// if NaN (should not happen, just in case)
		if( (hr != hr) || (wr.z != wr.z) ) 
			return VEC3D(0, 0, 1) ;

	}

	return wr_art;
}

-(Vec3D) sample
    : (const Vec3D *) wi_art
{
    int scatteringOrder = 0;
    return [self sample : wi_art : scatteringOrder];
}

-(float) G_1
    : (const Vec3D *) wi_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    
    if(wi.z > 0.9999f)
		return 1.0f;
	if(wi.z <= 0.0f)
		return 0.0f;

	// Lambda
	const float Lambda = [m_microsurfaceslope Lambda: wi_art];
	// value
	const float value = 1.0f / (1.0f + Lambda);
	return value;
}

// masking function at height h0
-(float) G_1
    : (const Vec3D *) wi_art
    : (const float  ) h0
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    
    if(wi.z > 0.9999f)
		return 1.0f;
	if(wi.z <= 0.0f)
		return 0.0f;

	// height CDF
	const float C1_h0  = [m_microsurfaceheight C1: h0];
	// Lambda
	const float Lambda = [m_microsurfaceslope Lambda: wi_art];
	// value
	const float value = powf(C1_h0, Lambda);
	return value;

}

// sample height in outgoing direction
-(float) sampleHeight
    : (const Vec3D *) wr_art
    : (const float  ) hr
    : (const float  ) U
{
    struct vec3 wr = Vec3D_to_vec3(wr_art);
    
    if(wr.z > 0.9999f)
		return FLT_MAX;
	if(wr.z < -0.9999f)
	{
		const float value = [ m_microsurfaceheight invC1
                : ( U * [m_microsurfaceheight C1: (hr)] )];
		return value;
	}
	if(fabs(wr.z) < 0.0001f)
		return hr;

	// probability of intersection
	const float G_1_ = [self G_1: wr_art: hr];
		
	if (U > 1.0f - G_1_) // leave the microsurface
		return FLT_MAX;

	const float h = [m_microsurfaceheight invC1
                :   [ m_microsurfaceheight C1:(hr)]
                    / powf(
                            (1.0f - U) ,
                            1.0f / [m_microsurfaceslope Lambda:(wr_art)]
                        )
                ];
	return h;
    
}

-(double) Fresnel_Reflection_plain
            :(double) theta
            :(double) n1
            :(double) k1
            :(double) n2
            :(double) k2
            ;
{
    double temp = ( n2*n2 - k2*k2 - n1*n1*sin(theta)*sin(theta) );
    double p2   = ( sqrt(temp*temp+4*n2*n2*k2*k2) + temp ) / 2;
    double q2   = ( sqrt(temp*temp+4*n2*n2*k2*k2) - temp ) / 2;
    double p    =   sqrt(p2);
    double q    =   sqrt(q2);
    
    double Fs,Fp;
    
    Fs = ( ( n1 * cos(theta) - p ) * ( n1 * cos(theta) - p ) + q2 ) /
          ( ( n1 * cos(theta) + p ) * ( n1 * cos(theta) + p ) + q2);
    Fp = ( Fs ) * ( ( p - n1 * sin(theta) * tan(theta) ) *
                    (   p - n1 * sin(theta) * tan(theta) ) + q2 )
                  / ( ( p + n1 * sin(theta) * tan(theta) )
                    * ( p + n1 * sin(theta) * tan(theta) ) + q2 );
    
    return ( Fs + Fp ) / 2;
}

-(double) Fresnel_Refraction_plain
            :(double) theta
            :(double) n1
            :(double) k1
            :(double) n2
            :(double) k2
            ;
{
    double nab = n1 / n2;
    double theta_t = asin( n1*sin(theta) / n2 );
    
    double Tp,Ts,factor;
    
    Tp = ( 2 * nab * cos(theta) / ( cos(theta) + nab * cos(theta_t))) *
          ( 2 * nab * cos(theta) / ( cos(theta) + nab * cos(theta_t)));
    Ts = ( 2 * nab * cos(theta) / ( nab * cos(theta) + cos(theta_t))) *
          ( 2 * nab * cos(theta) / ( nab * cos(theta) + cos(theta_t)));
    
    double nba = n2 / n1;
    factor = /* nba * nba * */ nba * ( cos(theta_t)/cos(theta) );
    
    return factor * ( Ts + Tp ) / 2;
}

-(void) reflect
    :(const Vec3D *) wi_art
    :(const Vec3D *) wm_art
    :(      Vec3D *) wo_art
{
    // in our implementation, we assume wi points to the surface

    Vec3D neg_wi_art = *wi_art;
    vec3d_negate_v(&neg_wi_art);
    
    // this function failed here sometimes
    // vec3d_vv_reflect_v(&neg_wi_art,&wm_art,&wo_art);
    
    double factor = 2 * fabs(vec3d_vv_dot(wi_art,wm_art));
    Vec3D wm_art_scaled = *wm_art;
    vec3d_d_mul_v(factor,&wm_art_scaled);
    vec3d_vv_add_v(&neg_wi_art,&wm_art_scaled,wo_art);
    vec3d_norm_v(wo_art);

}

-(void) refract
    : (const Vec3D *) wi_art
    : (const Vec3D *) wm_art
    : (const float)   eta
    : (      Vec3D *) wo_art
{
    
	const float cos_theta_i = vec3d_vv_dot(wi_art, wm_art);
	const float cos_theta_t2 = 1.0f - (1.0f-cos_theta_i*cos_theta_i) / (eta*eta);
	const float cos_theta_t = -sqrtf(fmaxf(0.0f,cos_theta_t2));
//
//    
//    
//	return wm * (dot(wi, wm) / eta + cos_theta_t) - wi / eta;
    
    Vec3D temp1 = *wm_art;
    double scaled_fator = (cos_theta_i/eta) + cos_theta_t;
    vec3d_d_mul_v(scaled_fator, &temp1);
    
    Vec3D temp2 = *wi_art;
    vec3d_negate_v(&temp2);
    vec3d_d_mul_v(1/eta, &temp2);
    
    vec3d_vv_add_v(&temp1, &temp2, wo_art);
    vec3d_norm_v(wo_art);
}

@end



@implementation MicrosurfaceConductor


-(id) init
    :(BOOL)     height_uniform
    :(BOOL)     slope_beckmann
    :(float)    alpha_x
    :(float)    alpha_y
{
    self = [super init
                : height_uniform
                : slope_beckmann
                : alpha_x
                : alpha_y
                ];
    
    return self;
}


// evaluate local phase function
-(float) evalPhaseFunction
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    struct vec3 wo = Vec3D_to_vec3(wo_art);
    struct vec3 wh = vec3_add(wi,wo);
    wh = normalize(wh);
    
    Vec3D wh_art = VEC3D(0,0,0);
    wh_art = vec3_to_Vec3D(wh);
    
    // half vector 
	//const vec3 wh = normalize(wi+wo);
  
	if(wh.z < 0.0f)
		return 0.0f;
    
    double F = 0;
    
    double cosine = vec3d_vv_dot( wi_art, &wh_art);
    
	F = [self Fresnel_Reflection_plain
                : cosine
                : 1
                : 0
                : m_eta
                : m_k];
	// value
	const float value = 0.25f
                        * F
                        * [ m_microsurfaceslope D_wi
                                    :   wi_art
                                    : & wh_art
                                    ]
                        / vec3d_vv_dot( wi_art, & wh_art );
	return value;
    
    
}

// sample local phase function
-(Vec3D) samplePhaseFunction
    : (const Vec3D *) wi_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    
    const float U1 = [ self generateRandomNumber ];
	const float U2 = [ self generateRandomNumber ];

	Vec3D wm_art = [m_microsurfaceslope sampleD_wi
                    : wi_art
                    : U1
                    : U2
                    ];
    
    
    Vec3D wo_art = VEC3D(0,0,0);
    
    [self reflect : wi_art : & wm_art : & wo_art];
    
    // in this case, wi DO NOT points to the surface,
	// reflect
	//const Vec3D wo = -wi + 2.0f * wm * dot(wi, wm);
    
	return wo_art;
    
}

// evaluate BSDF limited to single scattering 
// this is in average equivalent to eval(wi, wo, 1);
-(float) evalSingleScattering
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
{

    struct vec3 wi = Vec3D_to_vec3(wi_art);
    struct vec3 wo = Vec3D_to_vec3(wo_art);
    
    Vec3D neg_wi_art = *wi_art;
    Vec3D neg_wo_art = *wo_art;
    
    vec3d_negate_v(&neg_wi_art);
    vec3d_negate_v(&neg_wo_art);
    
    Vec3D wh_art = VEC3D(0,0,0);
    
    vec3d_vv_add_v( wi_art, wo_art, & wh_art );
    vec3d_norm_v( & wh_art );
    
    // half vector 
	//const vec3 wh = normalize(wi+wo);
    
    struct vec3 wh = Vec3D_to_vec3( & wh_art );
//    // half-vector
//	const vec3 wh = normalize(wi+wo);

	const float D = [m_microsurfaceslope D: & wh_art];

	// masking-shadowing 
	const float G2 = 1.0f / (1.0f + [m_microsurfaceslope Lambda:(wi_art)] + [m_microsurfaceslope Lambda:(wo_art)]);

	// BRDF * cos
	const float value = D * G2 / (4.0f * wi.z);

	return value;
}


@end


@implementation MicrosurfaceDiffuse


-(id) init
    :(BOOL)     height_uniform
    :(BOOL)     slope_beckmann
    :(float)    alpha_x
    :(float)    alpha_y
{
    self = [super init
                : height_uniform
                : slope_beckmann
                : alpha_x
                : alpha_y
                ];
    
    return self;
}


// evaluate local phase function
-(float) evalPhaseFunction
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
{
    const float U1 = [ self generateRandomNumber ];
	const float U2 = [ self generateRandomNumber ];
    
	Vec3D wm_art = [m_microsurfaceslope sampleD_wi
                        : wi_art
                        : U1
                        : U2
                        ];
	
	// value
	const float value = 1.0f/M_PI * fmaxf(0.0f, vec3d_vv_dot(wo_art, & wm_art));
    
	return value;
}

// sample local phase function
-(Vec3D) samplePhaseFunction
    : (const Vec3D *) wi_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    
    
    const float U1 = [ self generateRandomNumber ];
	const float U2 = [ self generateRandomNumber ];
	const float U3 = [ self generateRandomNumber ];
	const float U4 = [ self generateRandomNumber ];

	//vec3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);
    Vec3D wm_art = [m_microsurfaceslope sampleD_wi
                        : wi_art
                        : U1
                        : U2
                        ];
    struct vec3 wm = Vec3D_to_vec3( &wm_art );
    
	// sample diffuse reflection
	struct vec3 w1, w2;
	buildOrthonormalBasis( &w1, &w2, &wm);

	float r1 = 2.0f*U3 - 1.0f;
	float r2 = 2.0f*U4 - 1.0f;

	// concentric map code from
	// http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html
	float phi, r;
	if (r1 == 0 && r2 == 0) {
		r = phi = 0;
	} else if (r1*r1 > r2*r2) {
		r = r1;
		phi = (M_PI/4.0f) * (r2/r1);
	} else {
		r = r2;
		phi = (M_PI/2.0f) - (r1/r2) * (M_PI/4.0f);
	}
	float x = r*cosf(phi);
	float y = r*sinf(phi);
	float z = sqrtf(fmaxf(0.0f, 1.0f - x*x - y*y));
	
    //vec3 wo = x*w1 + y*w2 + z*wm;
    
    Vec3D w1_art = VEC3D(0,0,0);
    Vec3D w2_art = VEC3D(0,0,0);
    
    w1_art = vec3_to_Vec3D(w1);
    w2_art = vec3_to_Vec3D(w2);
    
    Vec3D wo_art = VEC3D(0,0,0);
    
    vec3d_dv_mul_dv_mul_dv_mul_add3_v(
                x, & w1_art,
                y, & w2_art,
                z, & wm_art,
                   & wo_art
                );

	return wo_art;
}

// evaluate BSDF limited to single scattering 
// this is in average equivalent to eval(wi, wo, 1);
-(float) evalSingleScattering
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    struct vec3 wo = Vec3D_to_vec3(wo_art);
    
    // sample visible microfacet
	const float U1 = [ self generateRandomNumber ];
	const float U2 = [ self generateRandomNumber ];
	//const vec3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);
    
    Vec3D wm_art = [m_microsurfaceslope sampleD_wi
                        : wi_art
                        : U1
                        : U2
                        ];
    struct vec3 wm = Vec3D_to_vec3( & wm_art );


	// shadowing given masking
	const float Lambda_i = [m_microsurfaceslope Lambda:(wi_art)];
	const float Lambda_o = [m_microsurfaceslope Lambda:(wo_art)];
	float G2_given_G1 = (1.0f + Lambda_i) / (1.0f + Lambda_i + Lambda_o);

	// evaluate diffuse and shadowing given masking
	const float value = 1.0f / (float)M_PI * fmaxf(0.0f, vec3d_vv_dot( & wm_art, wo_art)) * G2_given_G1;

	return value;
}


@end


@implementation MicrosurfaceDielectric

-(id) init
    :(BOOL)     height_uniform
    :(BOOL)     slope_beckmann
    :(float)    alpha_x
    :(float)    alpha_y
    :(float)    refratcion_index
{
    self = [super init
                : height_uniform
                : slope_beckmann
                : alpha_x
                : alpha_y
                ];
    
    m_eta = refratcion_index;
    
    return self;
}

-(float) evalBSDF
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    : (const int )    scatteringOrder
{
    struct vec3 wo = Vec3D_to_vec3(wo_art);
    wo = normalize(wo);
    double cosine_weight = wo.z;
    return [self eval
                : wi_art
                : wo_art
                : scatteringOrder
                ] / cosine_weight;
}

-(float) eval
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    : (const int )    scatteringOrder
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    struct vec3 wo = Vec3D_to_vec3(wo_art);
    
    Vec3D neg_wi_art = *wi_art;
    Vec3D neg_wo_art = *wo_art;
    
    vec3d_negate_v(&neg_wi_art);
    vec3d_negate_v(&neg_wo_art);
    
    // init
	//vec3 wr = -wi;
    struct vec3 wr = negtive(wi);
    Vec3D wr_art = vec3_to_Vec3D( wr );
    
	float hr = 1.0f + [ m_microsurfaceheight invC1:0.999f];
	BOOL outside = YES;

	float sum = 0.0f;
	
	// random walk
	int current_scatteringOrder = 0;	
	while(scatteringOrder==0 || current_scatteringOrder <= scatteringOrder)
	{
		// next height
		float U = [ self generateRandomNumber ];
        
		//hr = (outside) ? sampleHeight(wr, hr, U) : -sampleHeight(-wr, -hr, U);
        
        Vec3D neg_wr_art = wr_art;
        vec3d_negate_v( &neg_wr_art );
        
        if(outside)
        {
            hr = [self sampleHeight: & wr_art: hr : U];
        }
        else
        {
            hr = -[self sampleHeight: & neg_wr_art: -hr : U];
        }
        
		// leave the microsurface?
		if( hr == FLT_MAX || hr == -FLT_MAX)
			break;
		else
			current_scatteringOrder++;

		// next event estimation
		//float phasefunction = evalPhaseFunction(-wr, wo, outside, (wo.z>0) );
        
        float phasefunction = [self evalPhaseFunction
                                : & neg_wr_art
                                :   wo_art
                                :   outside
                                :   wo.z>0
                               ];
        
		//float shadowing = (wo.z>0) ? G_1(wo, hr) : G_1(-wo, -hr);
        
        float shadowing = (wo.z>0) ? [self G_1 :wo_art : hr] : [self G_1 : &neg_wo_art : -hr];
        
		float I = phasefunction * shadowing;

		if ( IsFiniteNumber(I) && (scatteringOrder==0 || current_scatteringOrder==scatteringOrder) )
			sum += I;
		
		// next direction
		//wr = samplePhaseFunction(-wr, outside, outside);
        
        wr_art = [self samplePhaseFunction: & neg_wr_art :outside : &outside];
        wr = Vec3D_to_vec3(&wr_art);
        
		// if NaN (should not happen, just in case)
		if( (hr != hr) || (wr.z != wr.z) ) 
			return 0.0f;
	}

	return sum;
}

// sample BSDF with a random walk
// scatteringOrder is set to the number of bounces computed for this sample
-(Vec3D) sample
    : (const Vec3D *) wi_art
    : (      int    ) scatteringOrder
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    
    Vec3D neg_wi_art = *wi_art;
    vec3d_negate_v(&neg_wi_art);

    // init
	//vec3 wr = -wi;
    struct vec3 wr = negtive(wi);
    Vec3D wr_art = vec3_to_Vec3D( wr );
    
	float hr = 1.0f + [ m_microsurfaceheight invC1:0.999f];
	BOOL outside = YES;
	
	// random walk
	scatteringOrder = 0;	
	while(true)
	{
		// next height
		float U = [ self generateRandomNumber ];
        Vec3D neg_wr_art = wr_art;
        vec3d_negate_v(&neg_wr_art);
        
		//hr = (outside) ? sampleHeight(wr, hr, U) : -sampleHeight(-wr, -hr, U);
        
        if(outside)
        {
            hr = [self sampleHeight: & wr_art: hr : U];
        }
        else
        {
            hr = -[self sampleHeight: & neg_wr_art: -hr : U];
        }
        
		// leave the microsurface?
		if( hr == FLT_MAX || hr == -FLT_MAX)
			break;
		else
			scatteringOrder++;

		// next direction
		//wr = samplePhaseFunction(-wr, outside, outside);
        wr_art = [self samplePhaseFunction: & neg_wr_art :outside : &outside];
        wr = Vec3D_to_vec3(&wr_art);
        
		// if NaN (should not happen, just in case)
		if( (hr != hr) || (wr.z != wr.z) ) 
			return VEC3D(0,0,1);
	}

	return wr_art;
}


// evaluate local phase function
-(float) evalPhaseFunction
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
{

    return [self evalPhaseFunction :wi_art : wo_art : YES : YES ] +
           [self evalPhaseFunction :wi_art : wo_art : YES : NO];
}

// evaluate local phase function
-(float) evalPhaseFunction
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    : (const BOOL   ) wi_outside
    : (const BOOL   ) wo_outside
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    struct vec3 wo = Vec3D_to_vec3(wo_art);
    
    Vec3D neg_wi_art = *wi_art;
    Vec3D neg_wo_art = *wo_art;
    
    vec3d_negate_v(&neg_wi_art);
    vec3d_negate_v(&neg_wo_art);
    
    const float eta = wi_outside ? m_eta : 1.0f / m_eta;

	if( wi_outside == wo_outside ) // reflection
	{
        Vec3D wh_art = VEC3D(0,0,0);
    
        vec3d_vv_add_v(wi_art, wo_art, & wh_art );
        vec3d_norm_v( &wh_art );
        
        Vec3D neg_wh_art = wh_art;
        vec3d_negate_v( & neg_wh_art );
        
        // half vector
        //const vec3 wh = normalize(wi+wo);
        
        //struct vec3 wh = Vec3D_to_vec3( & wh_art );
		
		// value
		const float value = (wi_outside) ?
                       (0.25f *
                       [m_microsurfaceslope D_wi : wi_art: &wh_art]
                       / vec3d_vv_dot(wi_art, &wh_art)
                       * [self Fresnel : wi_art :& wh_art :eta] ):
                       (0.25f
                       * [m_microsurfaceslope D_wi : & neg_wi_art : &neg_wh_art]
                       / vec3d_vv_dot( & neg_wi_art,  &neg_wh_art)
                       * [self Fresnel : &neg_wi_art :& neg_wh_art: eta] );
		return value;
	}
	else // transmission
	{
        
        Vec3D wh_art = VEC3D(0, 0, 0);
        
        vec3d_dv_mul_dv_mul_add_v(1, wi_art, eta, wo_art, &wh_art);
        vec3d_norm_v(&wh_art);
        vec3d_negate_v(&wh_art);
        
        if (wi_outside)
        {
            vec3d_d_mul_v(sign(VEC3D_I(wh_art,2)), &wh_art);
        }
        else
            vec3d_d_mul_v(sign(-VEC3D_I(wh_art,2)), &wh_art);
        
		//vec3 wh = -normalize(wi+wo*eta);
		//wh *= (wi_outside) ? (sign(wh.z)) : (-sign(wh.z));

		if(vec3d_vv_dot(&wh_art, wi_art) < 0)
			return 0;
        
        Vec3D neg_wh_art = wh_art;
        vec3d_negate_v( & neg_wh_art );
        
		float value;
		if(wi_outside)
        {
			value = eta*eta *
                    (1.0f - [self Fresnel : wi_art :& wh_art :eta])
                    * [m_microsurfaceslope D_wi : wi_art: &wh_art]
                    * fmaxf(0.0f, -vec3d_vv_dot(wo_art,&wh_art))
                    * 1.0f / powf(vec3d_vv_dot(wi_art, &wh_art)+eta*vec3d_vv_dot(wo_art,&wh_art), 2.0f);
		}
		else
		{
			value = eta*eta * (1.0f - [self Fresnel : &neg_wi_art :& neg_wh_art: eta])
                * [m_microsurfaceslope D_wi : & neg_wi_art : &neg_wh_art]
                * fmaxf(0.0f, -vec3d_vv_dot(& neg_wo_art,  &neg_wh_art)) *
				1.0f / powf(vec3d_vv_dot( & neg_wi_art,  &neg_wh_art)
                + eta * vec3d_vv_dot(& neg_wo_art,  &neg_wh_art), 2.0f);
		}

		return value;	
	}
}


// sample local phase function
-(Vec3D) samplePhaseFunction
    : (const Vec3D *) wi_art
{
    BOOL wo_outside;
    Vec3D vec_art = [self samplePhaseFunction: wi_art : YES : &wo_outside];
	return vec_art;
}

-(Vec3D) samplePhaseFunction
    : (const Vec3D *) wi_art
    : (const BOOL   ) wi_outside
    : (      BOOL  *) wo_outside
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    
    const float U1 = [ self generateRandomNumber ];
	const float U2 = [ self generateRandomNumber ];

	const float eta = wi_outside ? m_eta : 1.0f / m_eta;

    
    Vec3D neg_wi_art = * wi_art;
    vec3d_negate_v(&neg_wi_art);
    
    Vec3D wm_art = VEC3D(0,0,0);
    
    if(wi_outside)
        wm_art = [m_microsurfaceslope sampleD_wi
                        : wi_art
                        : U1
                        : U2
                        ];
    else
    {
        wm_art = [m_microsurfaceslope sampleD_wi
                        : & neg_wi_art
                        : U1
                        : U2
                        ];
        vec3d_negate_v(&wm_art);
    }
    
    
	const float F = [self Fresnel : wi_art : &wm_art: eta];

	if( [self generateRandomNumber] < F )
//    if( 0 < F )
	{
		Vec3D wo_art = VEC3D(0,0,0);
        
        [self reflect : wi_art : & wm_art : & wo_art];
        
        // in this case, wi DO NOT points to the surface,
        // reflect
        //const Vec3D wo = -wi + 2.0f * wm * dot(wi, wm);
        
        return wo_art;
    }
	else
	{
		*wo_outside = !wi_outside;
        
        Vec3D wo_art = VEC3D(0,0,0);
        
//        [BasicFunctionCollection GetSpecularRefractionDirection
//                : & neg_wi_art
//                : & wm_art
//                :   1
//                :   eta
//                : & wo_art
//                ];
//        vec3d_norm_v(&wo_art);
        
        [self refract:wi_art :&wm_art :eta :&wo_art];
        
        
		//const vec3 wo = refract(wi, wm, eta);
        
		return wo_art;
	}
}

// evaluate BSDF limited to single scattering 
// this is in average equivalent to eval(wi, wo, 1);
-(float) evalSingleScattering
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    struct vec3 wo = Vec3D_to_vec3(wo_art);
    
    Vec3D neg_wi_art = *wi_art;
    Vec3D neg_wo_art = *wo_art;
    
    vec3d_negate_v(&neg_wi_art);
    vec3d_negate_v(&neg_wo_art);

    bool wi_outside = true;
	bool wo_outside = wo.z > 0;

	const float eta = m_eta;

	if(wo_outside) // reflection
	{
		// D
        Vec3D wh_art = VEC3D(0,0,0);
        vec3d_vv_add_v( wi_art, wo_art, & wh_art );
        vec3d_norm_v( & wh_art );
        
		//const vec3 wh = normalize(vec3(wi+wo));
		const float D = [m_microsurfaceslope D : & wh_art];
		
		// masking shadowing
		const float Lambda_i = [m_microsurfaceslope Lambda:(wi_art)];
		const float Lambda_o = [m_microsurfaceslope Lambda:(wo_art)];

		const float G2 = 1.0f / (1.0f + Lambda_i + Lambda_o);

		// BRDF
		const float value = [self Fresnel: wi_art : &wh_art : eta] * D * G2 / (4.0f * wi.z);
		return value;
	}
	else // refraction
	{
		// D
        Vec3D wh_art = VEC3D(0,0,0);
        vec3d_dv_mul_dv_mul_add_v( 1, wi_art, eta, wo_art, &wh_art);
        vec3d_norm_v( & wh_art );
        vec3d_negate_v( & wh_art );
        
		//vec3 wh = -normalize(wi+wo*eta);
		if( eta < 1.0f )
			vec3d_negate_v( & wh_art );
        
		const float D = [m_microsurfaceslope D :( & wh_art ) ];

		// G2
		const float Lambda_i = [m_microsurfaceslope Lambda:(wi_art)];
        const float Lambda_o = [m_microsurfaceslope Lambda:(&neg_wo_art)];
        
		
		const float G2 = (float) eric_beta(1.0f+Lambda_i, 1.0f+Lambda_o);

		// BSDF
		const float value = fmaxf(0.0f, vec3d_vv_dot(wi_art, &wh_art))
                            *
                            fmaxf(0.0f, -vec3d_vv_dot(wo_art,&wh_art))
                            *
							1.0f / wi.z * eta * eta
                            *
                            (1.0f - [self Fresnel : wi_art :& wh_art :eta])
                            *
							G2
                            *
                            D / powf(
                                    vec3d_vv_dot(wi_art, &wh_art)
                                    +
                                    eta * vec3d_vv_dot(wo_art,&wh_art)
                                    , 2.0f);
		return value;
	}
}

-(float) Fresnel
    : (const Vec3D *) wi
    : (const Vec3D *) wm
    : (const float)   eta
{
    
	const float cos_theta_i = vec3d_vv_dot(wi, wm);
    
	const float cos_theta_t2 = 1.0f - (1.0f-cos_theta_i*cos_theta_i) / (eta*eta);

	// total internal reflection 
	if (cos_theta_t2 <= 0.0f) return 1.0f;

	const float cos_theta_t = sqrtf(cos_theta_t2);

	const float Rs = (cos_theta_i - eta * cos_theta_t) / (cos_theta_i + eta * cos_theta_t);
	const float Rp = (eta * cos_theta_i - cos_theta_t) / (eta * cos_theta_i + cos_theta_t);

	const float F = 0.5f * (Rs * Rs + Rp * Rp);
	return F;
}


@end

