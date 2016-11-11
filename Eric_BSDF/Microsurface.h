//
//  MicrosurfaceHeight.h
//  imp
//
//  Created by Chi Wang on 30/10/16.
//
//

#import <AdvancedRenderingToolkit.h>
#import "MicrosurfaceSlope.h"
#import "MicrosurfaceHeight.h"


// in this code, incident light should always points away from surface
@interface Microsurface : ArcObject
{
    float m_eta;
    float m_k;
@public

    // height distribution
	MicrosurfaceHeight  * m_microsurfaceheight;
    
	// slope distribution
	MicrosurfaceSlope   * m_microsurfaceslope;
    
    ArcObject<ArpRandomGenerator> * randomGenerator;
}

-(id) init
    :(BOOL)     height_uniform
    :(BOOL)     slope_beckmann
    :(float)    alpha_x
    :(float)    alpha_y
    ;

-(void) init_material
    :(float)    n
    :(float)    k
    ;

-(void) dealloc
    ;

-(void) initRandomNumberGenerator
    :(ArcObject <ArpRandomGenerator>  *) randomGen
    ;

// evaluate BSDF with a random walk (stochastic but unbiased)
// scatteringOrder=0 --> contribution from all scattering events
// scatteringOrder=1 --> contribution from 1st bounce only
// scatteringOrder=2 --> contribution from 2nd bounce only, etc..
// this method evaluate cosine weighted BSDF
-(float) eval
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    : (const int )    scatteringOrder
    ;

// evaluate unweighted BSDF
-(float) evalBSDF
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    : (const int )    scatteringOrder
    ;

// sample BSDF with a random walk
// scatteringOrder is set to the number of bounces computed for this sample
-(Vec3D) sample
    : (const Vec3D *) wi_art
    : (      int    ) scatteringOrder
    ;

-(Vec3D) sample
    : (const Vec3D *) wi_art
    ;

-(float) G_1
    : (const Vec3D *) wi_art
    ;

// masking function at height h0
-(float) G_1
    : (const Vec3D *) wi_art
    : (const float  ) h0
    ;

// sample height in outgoing direction
-(float) sampleHeight
    : (const Vec3D *) wo_art
    : (const float  ) h0
    : (const float  ) U
    ;

// evaluate local phase function
-(float) evalPhaseFunction
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    ;

// sample local phase function
-(Vec3D) samplePhaseFunction
    : (const Vec3D *) wi_art
    ;

// evaluate BSDF limited to single scattering 
// this is in average equivalent to eval(wi, wo, 1);
-(float) evalSingleScattering
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    ;

-(double) Fresnel_Reflection_plain
            :(double)   theta
            :(double)   n1
            :(double)   k1
            :(double)   n2
            :(double)   k2
            ;

-(double) Fresnel_Refraction_plain
            :(double) theta
            :(double) n1
            :(double) k1
            :(double) n2
            :(double) k2
            ;

-(void) reflect
    :(const Vec3D *) wi_art
    :(const Vec3D *) wm_art
    :(      Vec3D *) wo_art
    ;

-(void) refract
    : (const Vec3D *) wi_art
    : (const Vec3D *) wm_art
    : (const float)   eta
    : (      Vec3D *) wo_art
    ;

@end




@interface MicrosurfaceConductor : Microsurface
{
}

-(id) init
    :(BOOL)     height_uniform
    :(BOOL)     slope_beckmann
    :(float)    alpha_x
    :(float)    alpha_y
    ;

// evaluate local phase function
-(float) evalPhaseFunction
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    ;

// sample local phase function
-(Vec3D) samplePhaseFunction
    : (const Vec3D *) wi_art
    ;

// evaluate BSDF limited to single scattering 
// this is in average equivalent to eval(wi, wo, 1);
-(float) evalSingleScattering
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    ;

@end


@interface MicrosurfaceDiffuse : Microsurface
{
}

-(id) init
    :(BOOL)     height_uniform
    :(BOOL)     slope_beckmann
    :(float)    alpha_x
    :(float)    alpha_y
    ;

// evaluate local phase function
-(float) evalPhaseFunction
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    ;

// sample local phase function
-(Vec3D) samplePhaseFunction
    : (const Vec3D *) wi_art
    ;

// evaluate BSDF limited to single scattering 
// this is in average equivalent to eval(wi, wo, 1);
-(float) evalSingleScattering
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    ;

@end



@interface MicrosurfaceDielectric : Microsurface
{
    
}

-(id) init
    :(BOOL)     height_uniform
    :(BOOL)     slope_beckmann
    :(float)    alpha_x
    :(float)    alpha_y
    :(float)    refratcion_index
    ;

-(float) Fresnel
    : (const Vec3D *) wi
    : (const Vec3D *) wm
    : (const float)   eta
    ;



    
// evaluate BSDF with a random walk (stochastic but unbiased)
// scatteringOrder=0 --> contribution from all scattering events
// scatteringOrder=1 --> contribution from 1st bounce only
// scatteringOrder=2 --> contribution from 2nd bounce only, etc..
// evaluate cosine-weighted BSDF
-(float) eval
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    : (const int )    scatteringOrder
    ;

// evaluate unweighted BSDF
-(float) evalBSDF
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    : (const int )    scatteringOrder
    ;

// sample BSDF with a random walk
// scatteringOrder is set to the number of bounces computed for this sample
-(Vec3D) sample
    : (const Vec3D *) wi_art
    : (      int    ) scatteringOrder
    ;


// evaluate local phase function
-(float) evalPhaseFunction
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    ;

// evaluate local phase function
-(float) evalPhaseFunction
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    : (const BOOL   ) wi_outside
    : (const BOOL   ) wo_outside
    ;


// sample local phase function
-(Vec3D) samplePhaseFunction
    : (const Vec3D *) wi_art
    ;

-(Vec3D) samplePhaseFunction
    : (const Vec3D *) wi_art
    : (const BOOL   ) wi_outside
    : (      BOOL  *) wo_outside
    ;

// evaluate BSDF limited to single scattering 
// this is in average equivalent to eval(wi, wo, 1);
-(float) evalSingleScattering
    : (const Vec3D *) wi_art
    : (const Vec3D *) wo_art
    ;

@end





