//
//  MicrosurfaceHeight.h
//  imp
//
//  Created by Chi Wang on 30/10/16.
//
//

#import "MicrosurfaceSlope.h"

@implementation MicrosurfaceSlope

-(id)init
    : (const float) alpha_x
    : (const float) alpha_y
{
    self = [super init];
    
    m_alpha_x = alpha_x;
    m_alpha_y = alpha_y;
    
    return self;
}

-(id)init
{
    self = [super init];
    
    m_alpha_x = 1.0f;
    m_alpha_y = 1.0f;
    
    return self;
}

// projected roughness in wi
-(float) alpha_i
    : (const Vec3D *) wi_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    
    const float invSinTheta2 = 1.0f / (1.0f - wi.z*wi.z);
	const float cosPhi2 = wi.x*wi.x*invSinTheta2;
	const float sinPhi2 = wi.y*wi.y*invSinTheta2;
	const float alpha_i = sqrtf( cosPhi2*m_alpha_x*m_alpha_x + sinPhi2*m_alpha_y*m_alpha_y ); 
	return alpha_i;
}

// distribution of normals (NDF)
-(float) D
    : (const Vec3D *) wm_art
{
    struct vec3 wm = Vec3D_to_vec3(wm_art);
    
    if( wm.z <= 0.0f)
		return 0.0f;

	// slope of wm
	const float slope_x = -wm.x/wm.z;
	const float slope_y = -wm.y/wm.z;

	// value
	const float value = [self P22 :slope_x :slope_y] / (wm.z*wm.z*wm.z*wm.z);
	return value;
}

// distribution of visible normals (VNDF)
-(float) D_wi
    : (const Vec3D *) wi_art
    : (const Vec3D *) wm_art
{
    struct vec3 wm = Vec3D_to_vec3(wm_art);
    
    if( wm.z <= 0.0f)
		return 0.0f;

	// normalization coefficient
	const float projectedarea = [self projectedArea:wi_art];
    
	if(projectedarea == 0)
		return 0;
	const float c = 1.0f / projectedarea;

	// value
	const float value = c * fmaxf(0.0f, vec3d_vv_dot(wi_art, wm_art)) * [self D:wm_art];
	return value;
}

// sample the VNDF
-(Vec3D) sampleD_wi
    : (const Vec3D *)   wi_art
    : (const float )    U1
    : (const float )    U2
{
    
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    
	// stretch to match configuration with alpha=1.0	
	const struct vec3 wi_11 = normalize(vec3(m_alpha_x * wi.x, m_alpha_y * wi.y, wi.z));

	// sample visible slope with alpha=1.0
	Vec2D slope_11_art = [self sampleP22_11
                            : acosf(wi_11.z)
                            : U1
                            : U2
                            ];
    
    struct vec2 slope_11 = Vec2D_to_vec2( &slope_11_art );

	// align with view direction
	const float phi = atan2(wi_11.y, wi_11.x);
	struct vec2 slope = vec2(
                cosf(phi)*slope_11.x - sinf(phi)*slope_11.y,
                sinf(phi)*slope_11.x + cos(phi)*slope_11.y
                );

	// stretch back
	slope.x *= m_alpha_x;
	slope.y *= m_alpha_y;

	// if numerical instability
	if( (slope.x != slope.x) || !IsFiniteNumber(slope.x) ) 
	{
		if(wi.z > 0)
            return vec3_to_Vec3D(vec3(0.0f,0.0f,1.0f));
		else
        {
            return vec3_to_Vec3D(normalize(vec3(wi.x, wi.y, 0.0f)));
        }
	}

	// compute normal
	const struct vec3 wm = normalize(vec3(-slope.x, -slope.y, 1.0f));
    Vec3D wm_art = vec3_to_Vec3D(wm);
    
	return wm_art;
}

//// distribution of slopes
//-(float) P22
//    : (const float )    slope_x
//    : (const float )    slope_y
//{
//
//}
//
//// Smith's Lambda function
//-(float) Lambda
//    : (const Vec3D *) wi_art
//{
//}
//
//// projected area towards incident direction
//-(float) projectedArea
//    : (const Vec3D *) wi_art
//{
//}
//
//// sample the distribution of visible slopes with alpha=1.0
//-(Vec2D) sampleP22_11
//    : (const float )    theta_i
//    : (const float )    U1
//    : (const float )    U2
//{
//}


@end





@implementation MicrosurfaceSlopeBeckmann

-(id)init
    : (const float) alpha_x
    : (const float) alpha_y
{
    self =[super init
        : alpha_x
        : alpha_y
        ];
    return self;
}

-(id)init
{
    self =[super init];
    return self;
}

// distribution of slopes
-(float) P22
    : (const float )    slope_x
    : (const float )    slope_y
{
    const float value = 1.0f / (M_PI * m_alpha_x * m_alpha_y) * expf(-slope_x*slope_x/(m_alpha_x*m_alpha_x) - slope_y*slope_y/(m_alpha_y*m_alpha_y) );
	return value;
}

// Smith's Lambda function
-(float) Lambda
    : (const Vec3D *) wi_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    
    if(wi.z > 0.9999f)
		return 0.0f;
	if(wi.z < -0.9999f)
		return -1.0f;

	// a
	const float theta_i = acosf(wi.z);
	const float a = 1.0f/tanf(theta_i)/ [self alpha_i:wi_art];

	// value
	const float value = 0.5f*((float)erf(a) - 1.0f) + INV_2_SQRT_M_PI / a * expf(-a*a);

	return value;
}

// projected area towards incident direction
-(float) projectedArea
    : (const Vec3D *) wi_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    
    if(wi.z > 0.9999f)
		return 1.0f;
	if(wi.z < -0.9999f)
		return 0.0f;

	// a
	const float alphai = [self alpha_i : wi_art];
	const float theta_i = acosf(wi.z);
	const float a = 1.0f/tanf(theta_i)/alphai;

	// value
	const float value = 0.5f*((float)erf(a) + 1.0f)*wi.z + INV_2_SQRT_M_PI * alphai * sinf(theta_i) * expf(-a*a);

	return value;
}

// sample the distribution of visible slopes with alpha=1.0
-(Vec2D) sampleP22_11
    : (const float )    theta_i
    : (const float )    U_1
    : (const float )    U_2
{
    struct vec2 slope;
    Vec2D slope_art;
    
	if(theta_i < 0.0001f)
	{
		const float r = sqrtf(-logf(U_1));
		const float phi = 6.28318530718f * U_2;
		slope.x = r * cosf(phi);
		slope.y = r * sinf(phi);
        slope_art = vec2_to_Vec2D(slope);
		return slope_art;
	}

	// constant
	const float sin_theta_i = sinf(theta_i);
	const float cos_theta_i = cosf(theta_i);

	// slope associated to theta_i
	const float slope_i = cos_theta_i/sin_theta_i;

	// projected area
	const float a = cos_theta_i/sin_theta_i;	
	const float projectedarea = 0.5f*((float)erf(a) + 1.0f)*cos_theta_i + INV_2_SQRT_M_PI * sin_theta_i * expf(-a*a);
	if(projectedarea < 0.0001f || projectedarea!=projectedarea)
		return vec2_to_Vec2D(vec2(0,0));
	// VNDF normalization factor
	const float c = 1.0f / projectedarea;

	// search 
	float erf_min = -0.9999f;
	float erf_max = fmaxf(erf_min, (float)erf(slope_i));
	float erf_current = 0.5f * (erf_min+erf_max);

	while(erf_max-erf_min > 0.00001f)
	{
		if (!(erf_current >= erf_min && erf_current <= erf_max))
			erf_current = 0.5f * (erf_min + erf_max);

		// evaluate slope
		const float slope = eric_erfinv(erf_current);

		// CDF
		const float CDF = (slope>=slope_i) ? 1.0f : c * (INV_2_SQRT_M_PI*sin_theta_i*expf(-slope*slope) + cos_theta_i*(0.5f+0.5f*(float)erf(slope)));
		const float diff = CDF - U_1;

		// test estimate
		if( fabs(diff) < 0.00001f )
			break;

		// update bounds
		if(diff > 0.0f)
		{
			if(erf_max == erf_current)
				break;
			erf_max = erf_current;
		}
		else
		{
			if(erf_min == erf_current)
				break;
			erf_min = erf_current;
		}

		// update estimate
		const float derivative = 0.5f*c*cos_theta_i - 0.5f*c*sin_theta_i * slope;
		erf_current -= diff/derivative;
	}

	slope.x = eric_erfinv(fminf(erf_max, fmaxf(erf_min, erf_current)));
	slope.y = eric_erfinv(2.0f*U_2-1.0f);
    slope_art = vec2_to_Vec2D(slope);
	return slope_art;

}


@end



@implementation MicrosurfaceSlopeGGX

-(id)init
    : (const float) alpha_x
    : (const float) alpha_y
{
    self =[super init
        : alpha_x
        : alpha_y
        ];
    return self;
}

-(id)init
{
    self =[super init];
    return self;
}

// distribution of slopes
-(float) P22
    : (const float )    slope_x
    : (const float )    slope_y
{
    const float tmp = 1.0f + slope_x*slope_x/(m_alpha_x*m_alpha_x) + slope_y*slope_y/(m_alpha_y*m_alpha_y);
	const float value = 1.0f / (M_PI * m_alpha_x * m_alpha_y) / (tmp * tmp);
	return value;
}

// Smith's Lambda function
-(float) Lambda
    : (const Vec3D *) wi_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    
    if(wi.z > 0.9999f)
		return 0.0f;
	if(wi.z < -0.9999f)
		return -1.0f;

	// a
	const float theta_i = acosf(wi.z);
	const float a = 1.0f/tanf(theta_i)/ [self alpha_i:wi_art];

	// value
	const float value = 0.5f * (-1.0f + sign(a) * sqrtf(1 + 1/(a*a)));

	return value;
}

// projected area towards incident direction
-(float) projectedArea
    : (const Vec3D *) wi_art
{
    struct vec3 wi = Vec3D_to_vec3(wi_art);
    
    if(wi.z > 0.9999f)
		return 1.0f;
	if(wi.z < -0.9999f)
		return 0.0f;

	// a
	const float alphai = [self alpha_i : wi_art];
	const float theta_i = acosf(wi.z);
	const float sin_theta_i = sinf(theta_i);

	// value
	const float value = 0.5f * ( wi.z + sqrtf( wi.z * wi.z + sin_theta_i * sin_theta_i * alphai * alphai));

	return value;
}

// sample the distribution of visible slopes with alpha=1.0
-(Vec2D) sampleP22_11
    : (const float )    theta_i
    : (const float )    U_1
    : (const float )    U_2
{
    struct vec2 slope;
    Vec2D slope_art;
    
    
    if(theta_i < 0.0001f)
	{
		const float r = sqrtf( U_1 / ( 1.0f - U_1 ) );
		const float phi = 6.28318530718f * U_2;
		slope.x = r * cosf(phi);
		slope.y = r * sinf(phi);
		slope_art = vec2_to_Vec2D(slope);
		return slope_art;

	}
    
    // constant
	const float sin_theta_i = sinf(theta_i);
	const float cos_theta_i = cosf(theta_i);
	const float tan_theta_i = sin_theta_i/cos_theta_i;

	// slope associated to theta_i
	const float slope_i = cos_theta_i/sin_theta_i;

	// projected area
	const float projectedarea = 0.5f * (cos_theta_i + 1.0f);
	if(projectedarea < 0.0001f || projectedarea!=projectedarea)
		return vec2_to_Vec2D(vec2(0,0));
    
	// normalization coefficient
	const float c = 1.0f / projectedarea;

	const float A = 2.0f * U_1/ cos_theta_i /c - 1.0f;
	const float B = tan_theta_i;
	const float tmp = 1.0f / (A*A-1.0f);

	const float D = sqrtf(fmaxf(0.0f, B*B*tmp*tmp - (A*A-B*B)*tmp));
	const float slope_x_1 = B*tmp - D;
	const float slope_x_2 = B*tmp + D;
	slope.x = (A < 0.0f || slope_x_2 > 1.0f/tan_theta_i) ? slope_x_1 : slope_x_2;

	float U2;
	float S;
	if(U_2 > 0.5f)
	{
        S = 1.0f;
        U2 = 2.0f*(U_2-0.5f);
	}
	else
	{
        S = -1.0f;
        U2 = 2.0f*(0.5f-U_2);
	}
	const float z = ( U2 * ( U2 * ( U2 * 0.27385f-0.73369f)+ 0.46341f)) / ( U2 * ( U2 * ( U2 * 0.093073f + 0.309420f ) - 1.000000f ) + 0.597999f);
	slope.y = S * z * sqrtf(1.0f+slope.x*slope.x);
    
    slope_art = vec2_to_Vec2D(slope);
	return slope_art;
}


@end























