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

@interface Microsurface : ArcObject
{
@public

    // height distribution
	MicrosurfaceHeight  * m_microsurfaceheight;
    
	// slope distribution
	MicrosurfaceSlope   * m_microsurfaceslope;
}

-(id) init
    :(BOOL)     height_uniform
    :(BOOL)     slope_beckmann
    :(float)    alpha_x
    :(float)    alpha_y
    ;

-(void) dealloc
    ;

    // evaluate BSDF with a random walk (stochastic but unbiased)
	// scatteringOrder=0 --> contribution from all scattering events
	// scatteringOrder=1 --> contribution from 1st bounce only
	// scatteringOrder=2 --> contribution from 2nd bounce only, etc..
-(float) eval
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






@end
