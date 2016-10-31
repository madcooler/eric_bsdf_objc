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
        m_microsurfaceheight = [MicrosurfaceHeightUniform init];
    else
        m_microsurfaceheight = [MicrosurfaceHeightGaussian init];
    
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

    return self;
    
}

-(void) dealloc
{
    m_microsurfaceslope     = nil;
    m_microsurfaceheight    = nil;
    [super dealloc];
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
	
	// random walk
	int current_scatteringOrder = 0;	
	while(scatteringOrder==0 || current_scatteringOrder <= scatteringOrder)
	{
		// next height
		float U = generateRandomNumber();
        
        Vec3D wr_art = vec3_to_Vec3D( wr );
        
		hr = [self sampleHeight: & wr_art: hr : U];

		// leave the microsurface?
		if( hr == FLT_MAX )
			break;
		else
			current_scatteringOrder++;

		// next event estimation
        
        Vec3D neg_wr_art = vec3_to_Vec3D(negtive(wr));
		float phasefunction = [self evalPhaseFunction
                                : & neg_wr_art
                                :   wo_art
                               ];
		float shadowing = [self G_1 :wo_art : hr];
        
		float I = phasefunction * shadowing;

		if ( IsFiniteNumber(I) && (scatteringOrder==0 || current_scatteringOrder==scatteringOrder) )
			sum += I;
		
		// next direction
        
        Vec3D next_dir = [self samplePhaseFunction: & neg_wr_art];
        
		wr = Vec3D_to_vec3(& next_dir);
			
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
	float hr = 1.0f + [ m_microsurfaceheight invC1:0.999f];
	
	// random walk
	scatteringOrder = 0;	
	while(true)
	{
		// next height
		float U = generateRandomNumber();
        
        Vec3D wr_art = vec3_to_Vec3D( wr );
		hr = [self sampleHeight: & wr_art: hr : U];
		

		// leave the microsurface?
		if( hr == FLT_MAX )
			break;
		else
			scatteringOrder++;

		// next direction
		//wr = samplePhaseFunction(-wr);
        
        Vec3D neg_wr_art = vec3_to_Vec3D(negtive(wr));
        Vec3D next_dir = [self samplePhaseFunction: & neg_wr_art];
		wr = Vec3D_to_vec3(& next_dir);


		// if NaN (should not happen, just in case)
		if( (hr != hr) || (wr.z != wr.z) ) 
			return VEC3D(0, 0, 1) ;
	}

	return vec3_to_Vec3D(wr);
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
    : (const Vec3D *) wo_art
    : (const float  ) h0
    : (const float  ) U
{
    struct vec3 wo = Vec3D_to_vec3(wo_art);
    
    if(wo.z > 0.9999f)
		return FLT_MAX;
	if(wo.z < -0.9999f)
	{
		const float value = [ m_microsurfaceheight invC1
                : ( U * [m_microsurfaceheight C1: (h0)] )];
		return value;
	}
	if(fabs(wo.z) < 0.0001f)
		return h0;

	// probability of intersection
	const float G_1_ = [self G_1: wo_art: h0];
		
	if (U > 1.0f - G_1_) // leave the microsurface
		return FLT_MAX;

	const float h = [m_microsurfaceheight invC1
                :   [ m_microsurfaceheight C1:(h0)]
                    / powf(
                            (1.0f - U) ,
                            1.0f / [m_microsurfaceslope Lambda:(wo_art)]
                        )
                ];
	return h;
    
}

// evaluate local phase function
-(float) evalPhaseFunction
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    struct vec3 wo = Vec3D_to_vec3(wo_art);
}

// sample local phase function
-(Vec3D) samplePhaseFunction
    : (const Vec3D *) wi_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
}

// evaluate BSDF limited to single scattering 
// this is in average equivalent to eval(wi, wo, 1);
-(float) evalSingleScattering
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    struct vec3 wo = Vec3D_to_vec3(wo_art);
}




@end
