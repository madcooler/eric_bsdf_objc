/*!
 
 @file      FresnelFunctionCollection.h
 
 @brief     This file declares the basic function collection interface used in imp.
            `FresnelFunctionCollection` is a static class which provides static Fresnel functions.
            Please refer to the individual function descriptions for more details.
 
 @author    Chi Wang (chi@cgg.mff.cuni.cz), Lukas Novosad (novosad@cgg.mff.cuni.cz)
 @date      13/01/16
 
 */
// ART library imports
#import <AdvancedRenderingToolkit.h>


@interface FresnelFunctionCollection : NSObject


/*!
 `Fresnel_Refraction` calculates the Fresnel refraction based on the parameters passed to the method.
 
 @param theta The angle of the incoming direction.
 @param n1      The index of refraction for the first medium, e.g. the medium which will be left.
 @param k1      The absorption coefficient for the first medium, e.g. the medium which will be left.
 @param n2      The index of refraction for the second medium, e.g. the medium which will be entered.
 @param k2      The absorption coefficient for the second medium, e.g. the medium which will be entered.
 
 @param phase   The phase shift between the incoming light and the refracted light.
 @param Ts      The orthogonal Fresnel transmittance (e.g. refraction).
 @param Tp      The parallel Fresnel transmittance (e.g. refraction).
 @param TsTp    The ...
 @param factor  The ...
 
 @return Returns the Fresnel refraction.
 
 */
+(double) Fresnel_Refraction
            :(double)   theta
            :(double)   n1
            :(double)   k1
            :(double)   n2
            :(double)   k2
            :(double *) phase
            :(double *) Ts
            :(double *) Tp
            :(double *) TsTp
            :(double *) factor
            ;


/*!
 `Fresnel_Reflection` calculates the Fresnel reflection based on the parameters passed to the method.
 
 @param theta The angle of the incoming direction.
 @param n1      The index of refraction for the first medium, e.g. the medium which will be left.
 @param k1      The absorption coefficient for the first medium, e.g. the medium which will be left.
 @param n2      The index of refraction for the second medium, e.g. the medium which will be entered.
 @param k2      The absorption coefficient for the second medium, e.g. the medium which will be entered.
 
 @param rs      The ...
 @param rp      The ...
 @param phaseS  The phase shift of the orthogonal Fresnel reflection.
 @param phaseP  The phase shift of the parallel Fresnel reflection.
 @param Fs      The orthogonal Fresnel reflectance.
 @param Fp      The parallel Fresnel reflectance.
 
 @return Returns the Fresnel transmittance.
 
 */
+(double) Fresnel_Reflection
            :(double)   theta
            :(double)   n1
            :(double)   k1
            :(double)   n2
            :(double)   k2
            :(double *) rs
            :(double *) rp
            :(double *) phaseS
            :(double *) phaseP
            :(double *) Fs
            :(double *) Fp
            ;


/*!
 `imp_fresnel_plain_reflective_attenuation_realvalued_IOR` calculates the Fresnel reflection based on the parameters passed to the method.
 This function has been inspired by the ART function `fresnel_plain_reflective_attenuation_realvalued_IOR`.
 
 @param art_gv                      The ART global values pointer.
 @param incidentDirection           The incident light direction.
 @param localNormalDirection        The local normal direction, e.g. the (sampled) half-vector.
 @param normalDirection             The normal direction.
 @param sampledReflectionDirection  The sampled reflection direction.
 @param n                           The index of refraction of ...
 @param k                           The absorption coefficient of ...
 @param attenuation                 The attenuation of ...
 
 @return Returns the non-polarising Fresnel reflectance.
 
 */
+(void) imp_fresnel_plain_reflective_attenuation_realvalued_IOR
            :(const ART_GV *)       art_gv
            :(const Vec3D *)        incidentDirection
            :(const Vec3D *)        localNormalDirection
            :(const Vec3D *)        sampledReflectionDirection
            :(const double)         n
            :(const double)         k
            :(ArAttenuation *)      attenuation
            ;


/*!
 `imp_fresnel_polarising_reflective_attenuation_realvalued_IOR` calculates the Fresnel reflection based on the parameters passed to the method.
 This function has been inspired by the ART function `fresnel_polarising_reflective_attenuation_realvalued_IOR`.
 
 @param art_gv                      The ART global values pointer.
 @param incidentDirection           The incident light direction.
 @param localNormalDirection        The local normal direction, e.g. the (sampled) half-vector.
 @param normalDirection             The normal direction.
 @param sampledReflectionDirection  The sampled reflection direction.
 @param n                           The index of refraction of ...
 @param k                           The absorption coefficient of ...
 @param attenuation                 The attenuation of ...
 
 @return Returns the polarising Fresnel reflectance.
 
 */

+(void) imp_fresnel_polarising_reflective_attenuation_realvalued_IOR
            :(const ART_GV *)       art_gv
            :(const Vec3D *)        incidentDirection
            :(const Vec3D *)        localNormalDirection
            :(const Vec3D *)        normalDirection
            :(const Vec3D *)        sampledReflectionDirection
            :(const double)         n
            :(const double)         k
            :(ArAttenuation *)      attenuation
            ;



/*!
 `imp_fresnel_plain_refractive_attenuation` calculates the Fresnel refraction based on the parameters passed to the method.
 This function has been inspired by the ART function `fresnel_plain_refractive_attenuation`.
 
 @param art_gv                      The ART global values pointer.
 @param incidentDirection           The incident light direction.
 @param localNormalDirection        The local normal direction, e.g. the (sampled) half-vector.
 @param normalDirection             The normal direction.
 @param sampledReflectionDirection  The sampled reflection direction.
 @param n                           The index of refraction of ...
 @param k                           The absorption coefficient of ...
 @param attenuation                 The attenuation of ...
 
 @return Returns the non-polarising Fresnel refraction.
 
 */
+(void) imp_fresnel_plain_refractive_attenuation
            :(const ART_GV *)       art_gv
            :(const Vec3D *)        incidentDirection
            :(const Vec3D *)        localNormalDirection
            :(const Vec3D *)        sampledReflectionDirection
            :(const double)         n
            :(const double)         k
            :(ArAttenuation *)      attenuation
            ;


/*!
 `imp_fresnel_polarising_refractive_attenuation` calculates the Fresnel refraction based on the parameters passed to the method.
 This function has been inspired by the ART function `fresnel_polarising_reflractive_attenuation`.
 
 @param art_gv                      The ART global values pointer.
 @param incidentDirection           The incident light direction.
 @param localNormalDirection        The local normal direction, e.g. the (sampled) half-vector.
 @param normalDirection             The normal direction.
 @param sampledReflectionDirection  The sampled reflection direction.
 @param sampledRefractionDireciton  The sampled refraction direction.
 @param n                           The index of refraction of ...
 @param k                           The absorption coefficient of ...
 @param attenuation                 The attenuation of ...
 
 @return Returns the polarising Fresnel refraction.
 
 */
+(void) imp_fresnel_polarising_refractive_attenuation
            :(const ART_GV *)   art_gv
            :(const Vec3D *)    incidentDirection
            :(const Vec3D *)    localNormalDirection
            :(const Vec3D *)    normalDirection
            :(const Vec3D *)    sampledReflectionDirection
            :(const Vec3D *)    sampledRefractionDirection
            :(const double)     n
            :(const double)     k
            :(ArAttenuation *)  attenuation
            ;


@end
