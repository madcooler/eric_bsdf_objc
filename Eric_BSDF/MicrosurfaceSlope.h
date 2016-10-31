//
//  MicrosurfaceSlope.h
//  art
//
//  Created by Administrator on 30/10/16.
//
//

#import <AdvancedRenderingToolkit.h>
#include "EricBSDF_BasicFunctions.h"

@interface MicrosurfaceSlope : ArcObject
{
@public
    float m_alpha_x;
    float m_alpha_y;
}

-(id)init
    ;
    
-(id)init
    : (const float) alpha_x
    : (const float) alpha_x
    ;

// projected roughness in wi
-(float) alpha_i
    : (const Vec3D *) wi_art
    ;

// distribution of normals (NDF)
-(float) D
    : (const Vec3D *) wm_art
    ;

// distribution of visible normals (VNDF)
-(float) D_wi
    : (const Vec3D *) wi_art
    : (const Vec3D *) wm_art
    ;

// sample the VNDF
-(Vec3D) sampleD_wi
    : (const Vec3D *)         wi_art
    : (const float )    U1
    : (const float )    U2
    ;

// distribution of slopes
-(float) P22
    : (const float )    slope_x
    : (const float )    slope_y
    ;

// Smith's Lambda function
-(float) Lambda
    : (const Vec3D *) wi_art
    ;

// projected area towards incident direction
-(float) projectedArea
    : (const Vec3D *) wi_art
    ;

// sample the distribution of visible slopes with alpha=1.0
-(Vec2D) sampleP22_11
    : (const float )    theta_i
    : (const float )    U1
    : (const float )    U2
    ;


@end


@interface MicrosurfaceSlopeBeckmann : MicrosurfaceSlope
{

}

-(id)init
    ;
    
-(id)init
    : (const float) alpha_x
    : (const float) alpha_x
    ;


// distribution of slopes
-(float) P22
    : (const float )    slope_x
    : (const float )    slope_y
    ;

// Smith's Lambda function
-(float) Lambda
    : (const Vec3D *) wi_art
    ;

// projected area towards incident direction
-(float) projectedArea
    : (const Vec3D *) wi_art
    ;

// sample the distribution of visible slopes with alpha=1.0
-(Vec2D) sampleP22_11
    : (const float )    theta_i
    : (const float )    U_1
    : (const float )    U_2
    ;


@end

