/*!
 
 @file      BasicFunctionCollection.h
 
 @brief     This file declares the basic function collection interface used in imp.
            `BasicFunctionCollection` is a static class which provides static helper functions.
            Please refer to the individual function descriptions for more details.
 
 @author    Chi Wang (chi@cgg.mff.cuni.cz), Lukas Novosad (novosad@cgg.mff.cuni.cz)
 @date      08/01/16
 
 */
// ART library imports
#import <AdvancedRenderingToolkit.h>


@interface BasicFunctionCollection : ArcObject
{

}

/*!
 `GetRotationAngles` calculates the rotation angles for both the global to local transformation and for the local to global transformation.
 
 @param incidentDirection       The incident (light) direction.
 @param outgoingDirection       The outgoing (light) direction.
 @param surfaceNormalDirection  The surface normal direction.
 @param localNormalDirection    The local (e.g. microsurface) normal direction.
 @param globalToLocalAngle      The rotation angle for the global to local transformation.
 @param localToGlobalAngle      The rotation angle for the local to global transformation.
 
 @return The rotation angles are returned via the pointers to globalToLocalAngle and localToGlobalAngle.
 */
+(void) GetRotationAngles
            :(const Vec3D *)    incidentDirection
            :(const Vec3D *)    outgoingDirection
            :(const Vec3D *)    surfaceNormalDirection
            :(const Vec3D *)    localNormalDirection
            :(     double *)    globalToLocalAngle
            :(     double *)    localToGlobalAngle
            ;


// calculate the rotation angle from rf1 to rf2
+(void) GetRotationAngleFromTwoReferneceFrames
            :(const ART_GV           *) art_gv
            :(const ArReferenceFrame *) rf1
            :(const ArReferenceFrame *) rf2
            :(      double           *) angle
            ;

// rotate MM from rf1 to rf2
+(void) RotateMuellerMatrixBetweenTwoReferneceFrames
            :(const ART_GV                  *) art_gv
            :(const ArReferenceFrame        *) rf1
            :(const ArReferenceFrame        *) rf2
            :(      ArMuellerMatrixSample   *) mm
            ;

/*!
 `GetThetaAndPhiFromVec3D` extracts the spherical angles from the vector passed to the method.
 
 @param vec     The vector from which the spherical angles will be extracted.
 @param theta   The elevation angle of the vector.
 @param phi     The azimuthal angle of the vector.
 
 @return        The spherical angles are returned via the pointers to theta and phi.
 */
+(void) GetThetaAndPhiFromVec3D
            :(const Vec3D *)    vec
            :(double *)         theta
            :(double *)         phi
            ;



/*!
 `GetVec3DFromThetaAndPhi` constructs a vector from the spherical angles passed to the method.
 
 @param theta   The elevation angle of the vector.
 @param phi     The azimuthal angle of the vector.
 @param vec     The vectot based on the elvation and azimuthal angles.

 @return        The vector is returned via the pointer to vec.
 */
+(void) GetVec3DFromThetaAndPhi
            :(const double *)   theta
            :(const double *)   phi
            :(Vec3D *)          vec
            ;


/*!
 `GetXAndYCoordinatesFromThetaAndPhi` obtains the corresponding x and y coordinates from the spherical angles theta and phi.
 
 @param theta       The spherical elevation angle.
 @param phi         The spherical azimuthal angle.
 @param imageSize   The size of the image, used to determine the center of the image.
 @param x           The pointer to the x-coordinate (e.g. abscissa) to be returned.
 @param y           The pointer to the y-coordinate (e.g. ordinate) to be returned.
 
 @return            The cartesian coordinates are returned via the pointers to x and y.
 */
+(void) GetXAndYCoordinatesFromThetaAndPhi
            :(double)           theta
            :(double)           phi
            :(IVec2D)           imageSize
            :(int *)            x
            :(int *)            y
            ;

+(void) GetXAndYCoordinatesFromThetaAndPhi_AveragedMapping
            :(double)           theta
            :(double)           phi
            :(IVec2D)           imageSize
            :(int *)            x
            :(int *)            y
            ;


/*!
 `GetSpecularRefractionDirection` calculates the specular refraction based on the data passed to the method.
 
 @param incidentDirection   The incident (light) direction.
 @param normalDirection     The (surface) normal direction.
 @param eta                 The index of refraction between the two media, where the refraction occurs.
 @param refractionDirection The refracted (light) direction.
 
 @return The refracted (light) direction is returned via the pointer to refractionDirection.
 */
+(void) GetSpecularRefractionDirection
            :(const Vec3D *)    incidentDirection
            :(const Vec3D *)    normalDirection
            :(const double)     n_i
            :(const double)     n_t
            :(Vec3D *)          refractionDirection
            ;



/*!
 `GetSpecularReflectionDirection`calculate the specular reflection based on the data passed to the method.
 This method is intended for layered models: In case the normal points away from the surface,
 the normal is flipped before the calculations are performed.
 
 @param incidentDirection       The incident (light) direction.
 @param normalDirection         The (surface) normal direction.
 @param reflectionDirection     The reflected (light) direction.
 
 @return The reflected (light) direction is returned via the pointer to reflectionDirection.
 
 */
+(void) GetSpecularReflectionDirection
            :(const Vec3D *)    incidentDirection
            :(const Vec3D *)    normalDirection
            :(Vec3D *)          reflectionDirection
            ;

/*!
 `IsNormalReflection` determine the status of reflection. 
 
 For a normal reflection, the sign of the z component of the direction of the incoming and outgoing should be reversed.
 Otherwise, the multireflection should happen at the interface.
 
 @param incomingDirection     The incident (light) direction.
 @param outgoingDirection     The reflected (light) direction.
 
 @return YES or NO.
 
 */
+(BOOL) IsNormalReflection
        :(const Vec3D *) incomingDirection
        :(const Vec3D *) outgoingDirection
        ;

/*!
 `IsNormalRefraction` determine the status of reflection.
 
 For a normal reflection, the sign of the z component of the direction of the incoming and outgoing should be the same.
 
 
 @param incomingDirection     The incident (light) direction.
 @param outgoingDirection     The refracted (light) direction.
 
 @return YES or NO.
 
 */

+(BOOL) IsNormalRefraction
        :(const Vec3D *) incomingDirection
        :(const Vec3D *) outgoingDirection
        ;

/*!
 `GetStokesValueFromArLightAlphaSample` extrat the Stokes Vector component from ArLightAlphaSample structure
 
 @param ArLightAlphaSample    The pointer to a ArLightAlphaSample instance.
 @param valueI                The first  component(intensity) of the Stokes Vector.
 @param valueQ                The second component of the Stokes Vector.
 @param valueU                The third  component of the Stokes Vector.
 @param valueV                The fourth component of the Stokes Vector.
 
 @return The the Stokes Vector component value are returned by the pointer to valueI, valueQ, valueU, and valueV.
 */

+(void) GetStokesValueFromArLightAlphaSample
        :(const ART_GV             *) art_gv
        :(const ArLightAlphaSample *) lightSample
        :(      double             *) valueI
        :(      double             *) valueQ
        :(      double             *) valueU
        :(      double             *) valueV
        ;

+(void) GetStokesValueFromArStokesVectorSample
        :(const ART_GV               *) art_gv
        :(const ArStokesVectorSample *) stokes
        :(      double               *) valueI
        :(      double               *) valueQ
        :(      double               *) valueU
        :(      double               *) valueV
        ;

+(void) GetIntensityValueFromArLightAlphaSample
        :(const ART_GV             *) art_gv
        :(const ArLightAlphaSample *) lightSample
        :(      double             *) valueI
        ;

+(double) DistanceBetweenPnt3D
        :(const Pnt3D *) p1
        :(const Pnt3D *) p2
        ;


+(void) ApplyPolarisationFilter
        :(const ART_GV         *)     art_gv
        :(const double          )     angle
        :(const double          )     strength
        :(ArStokesVectorSample *)     sv
        ;

+(void) ApplyCompensatorFilter
        :(const ART_GV         *) art_gv
        :(const double)           phase_shift
        :(const double)           angle
        :(ArStokesVectorSample *) sv
        ;

+(void) GetNormalisedStokesVectorAndDOLP
        :(const ART_GV              *) art_gv
        :(const ArLightAlphaSample  *) light
        :(      double              *) dolp
        :(ArStokesVectorSample      *) sv
        ;
        
+(void) NormaliseStokesVector
        :(const ART_GV         *)     art_gv
        :(ArStokesVectorSample *)     sv
        ;

+(void) ApplyPolarisationAnalyser
        :(const ART_GV              *) art_gv
        :(const ArLightAlphaSample  *) light
        :(      double               ) angle
        :(      double              *) intensity
        ;

+(void) ApplyPolarisationAnalyserForStokes
        :(const ART_GV              *) art_gv
        :(const ArStokesVectorSample*) sv
        :(      double               ) angle
        :(      double              *) intensity
        ;
+(void) GetNormalisedStokesVectorSampleAndDOLP
        :(const ART_GV              *) art_gv
        :(const ArStokesVectorSample*) light
        :(      double              *) dolp
        :(ArStokesVectorSample      *) sv
        ;
// sample Lambertian reflection direction
+(void) SampleDiffuseScatteringDirection
        : ( ArcObject <ArpRandomGenerator>  *)  randomGenerator
        : ( Vec3D                           *)  outgoing_dir_local
        : ( double                          *)  PDF
        ;
@end
