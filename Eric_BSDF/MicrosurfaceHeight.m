//
//  MicrosurfaceHeight.h
//  imp
//
//  Created by Chi Wang on 30/10/16.
//
//

#import "MicrosurfaceHeight.h"

@implementation MicrosurfaceHeight

-(id)init
{
    self = [super init];
    return self;
}

-(float) P1
    :(const float) h
{
    return 0;
}

-(float) C1
    :(const float) h
{
    return 0;
}

-(float) invC1
    :(const float) U
{
    return 0;
}

@end


@implementation MicrosurfaceHeightUniform
-(id)init
{
    self = [super init];
    return self;
}

-(float) P1
    :(const float) h
{
    const float value = (h >= -1.0f && h <= 1.0f) ? 0.5f : 0.0f;
	return value;
}

-(float) C1
    :(const float) h
{
    const float value = fminf(1.0f, fmaxf(0.0f, 0.5f*(h+1.0f)));
	return value;
}

-(float) invC1
    :(const float) U
{
    const float h = fmaxf(-1.0f, fminf(1.0f, 2.0f*U-1.0f));
	return h;
}

@end

@implementation MicrosurfaceHeightGaussian

-(id)init
{
    self = [super init];
    return self;
}

-(float) P1
    :(const float) h
{
    const float value = INV_SQRT_2_M_PI * expf(-0.5f * h*h);
	return value;
}

-(float) C1
    :(const float) h
{
    const float value = 0.5f + 0.5f * (float)eric_erf(INV_SQRT_2*h);
	return value;
}

-(float) invC1
    :(const float) U
{
    const float h = SQRT_2 * eric_erfinv(2.0f*U - 1.0f);
	return h;
}

@end