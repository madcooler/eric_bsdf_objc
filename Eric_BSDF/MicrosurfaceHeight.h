//
//  MicrosurfaceHeight.h
//  imp
//
//  Created by Chi Wang on 30/10/16.
//
//

#import <AdvancedRenderingToolkit.h>
#include "EricBSDF_BasicFunctions.h"


@interface MicrosurfaceHeight : ArcObject
{

}
-(id)init
    ;

-(float) P1
    :(const float) h;

-(float) C1
    :(const float) h;

-(float) invC1
    :(const float) U;

@end



@interface MicrosurfaceHeightUniform : MicrosurfaceHeight
{

}

-(id)init
    ;
    
-(float) P1
    :(const float) h;

-(float) C1
    :(const float) h;

-(float) invC1
    :(const float) U;

@end



@interface MicrosurfaceHeightGaussian : MicrosurfaceHeight
{

}

-(id)init
    ;

-(float) P1
    :(const float) h;

-(float) C1
    :(const float) h;

-(float) invC1
    :(const float) U;
    
@end
