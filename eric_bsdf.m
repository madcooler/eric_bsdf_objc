/*!
 
 @file      imp.m
 
 @brief     This file is the main ART file for the imp project.
 
 @author    Chi Wang (chi@cgg.mff.cuni.cz), Lukas Novosad (novosad@cgg.mff.cuni.cz)
 @date
 
 */
// ART library imports
#import "AdvancedRenderingToolkit.h"

#import "Eric_BSDF/Microsurface.h"

Vec3D generateIncidentLight_Rad(double theta_rad)
{
    float theta = theta_rad ;
    float phi   = 180 * MATH_DEG_TO_RAD ;
    // calculate cartesian coordinates from spherical angles
    float zz  = cosf(theta);
    float tmp = sinf(theta);
    float xx  = tmp * cosf(phi);
    float yy  = tmp * sinf(phi);
    
    // create a downward pointing, incident light direction vector based on the previously generated cartesian coordinates
    Vec3D incidentLightDirection = VEC3D(xx, yy, zz);
    return incidentLightDirection;
}

void test_normalisation(
    ART_GV                          * art_gv,
    ArcObject <ArpRandomGenerator>  *  randomGenerator
    )
{
    struct vec3 wi = vec3(-0.64278757572174072,0, 0.76604443788528442);
    //generateRandomDirectionUp();
    //struct vec3 wo = vec3(0.4903473153283226, 0.20888626573301827, 0.84612412702771999);
    
    
    Vec3D wi_art = generateIncidentLight_Rad(1.01);
    vec3d_norm_v(&wi_art);
    
    
    Microsurface * m = [ALLOC_INIT_OBJECT(MicrosurfaceDielectric)
                        : NO
                        : YES
                        : 0.25
                        : 0.25
                        : 1.5
                        ];
    
    
    [m initRandomNumberGenerator:randomGenerator];
    
    
    // quadrature \int p(wi, wo) dwo
    double value_quadrature = 0;
    double normalised_bsdf  = 0;
    double bsdf_1 = 0;
    double bsdf_2 = 0;
    double bsdf_3 = 0;
    
    int count = 0;
    
    double step = 0.01;
    
//    for(double theta_o=M_PI/2 ; theta_o < M_PI ; theta_o += step)
    for(double theta_o=0 ; theta_o < M_PI/2 ; theta_o += step)
        for(double phi_o=0 ; phi_o < 2.0*M_PI ; phi_o += step) {
            double x = cos(phi_o)*sin(theta_o);
            double y = sin(phi_o)*sin(theta_o);
            double z = cos(theta_o);
            struct vec3 wo = vec3(x,y,z);
            Vec3D wo_art = vec3_to_Vec3D(wo);
            vec3d_norm_v(&wo_art);
        
            const int N = 10;
            double value_current = 0;
            double bounce_1 = 0;
            double bounce_2 = 0;
            double bounce_3 = 0;
            for(int n=0 ; n<N ; ++n)
            {
                count++;
                if(theta_o > 2.45)
                {
                    printf("");
                }

//                value_current += [m eval:&wi_art:&wo_art:0];
                bounce_1 += [m eval:&wi_art:&wo_art:1];
                bounce_2 += [m eval:&wi_art:&wo_art:2];
                bounce_3 += [m eval:&wi_art:&wo_art:3];
            }
                normalised_bsdf += step*step*fabs(sin(theta_o)) * value_current / N;
                bsdf_1 += step*step*fabs(sin(theta_o)) * bounce_1 / N;
                bsdf_2 += step*step*fabs(sin(theta_o)) * bounce_2 / N;
                bsdf_3 += step*step*fabs(sin(theta_o)) * bounce_3 / N;
        
//            value_quadrature += 0.005*0.005*fabs(sin(theta_o)) * (double)[m evalPhaseFunction:&wi_art:&wo_art];
        }
    // display
    printf( "\\int f_p(wi, wo) dwo = \t\t %f \n",  normalised_bsdf );
    printf( "\\int p(wi, wo) dwo = \t\t %f \n",  value_quadrature );
    printf( "1st bounce = \t\t %f \n",  bsdf_1 );
    printf( "2nd bounce = \t\t %f \n",  bsdf_2 );
    printf( "3rd bounce = \t\t %f \n",  bsdf_3 );
    RELEASE_OBJECT(m);
}

void test_reflection_sample(
    ART_GV                          * art_gv,
    ArcObject <ArpRandomGenerator>  *  randomGenerator
    )
{
    struct vec3 wi = vec3(-0.64278757572174072,0, 0.76604443788528442);
    //generateRandomDirectionUp();
//    struct vec3 wo = vec3(0.4903473153283226, 0.20888626573301827, 0.84612412702771999);
    
    
    Vec3D wi_art = vec3_to_Vec3D(wi);
    vec3d_norm_v(&wi_art);
//    Vec3D wo_art = vec3_to_Vec3D(wo);
//    vec3d_norm_v(&wo_art);
    
    Microsurface * m = [ALLOC_INIT_OBJECT(MicrosurfaceConductor)
                        : NO
                        : NO
                        : 0.1
                        : 0.1
                        //: 1.5
                        ];
    
    [m initRandomNumberGenerator:randomGenerator];
    
    Vec3D sampledReflectionDirection = [m sample
                        : & wi_art
                        ];


}



// this test has been passed 
void test_single_scattering(
    ART_GV                          * art_gv,
    ArcObject <ArpRandomGenerator>  *  randomGenerator
    )
{
    struct vec3 wi = vec3(-0.64278757572174072,0, 0.76604443788528442);
    //generateRandomDirectionUp();
    struct vec3 wo = vec3(0.4903473153283226, 0.20888626573301827, 0.84612412702771999);
    
    
    Vec3D wi_art = vec3_to_Vec3D(wi);
    vec3d_norm_v(&wi_art);
    Vec3D wo_art = vec3_to_Vec3D(wo);
    vec3d_norm_v(&wo_art);
    
    Microsurface * m = [ALLOC_INIT_OBJECT(MicrosurfaceDielectric)
                        : NO
                        : NO
                        : 0.1
                        : 0.1
                        : 1.5
                        ];
    
    [m initRandomNumberGenerator:randomGenerator];
    
    const int N = 100000;
    // eval truncated random walk (loop because eval is stochastic)
    double V = 0;
    for(int n=0 ; n<N ; ++n)
    {
        V += [m eval: &wi_art: &wo_art :2] / (double)N;
    }
    
    // eval single (loop because single of diffuse is also stochastic)
    double V_single = 0;
    for(int n=0 ; n<N ; ++n)
    {
//        V_single += [m evalSingleScattering: &wi_art: &wo_art] / (double)N;

    }

    printf ( "random␣walk␣cut␣after␣1st␣bounce␣=␣%f \n",V );
    printf ( "single␣scattering␣=␣%f \n" , V_single );

}

// this test has been passed
void test_sample_bsdf(
    ART_GV                          * art_gv,
    ArcObject <ArpRandomGenerator>  * randomGenerator
    )
{
    struct vec3 wi = vec3(-0.64278757572174072,0, 0.76604443788528442);
//    struct vec3 wi = vec3(2, -3, 6);
    //generateRandomDirectionUp();
    //struct vec3 wo = vec3(1, 3, -7);
    
    
    Vec3D wi_art = vec3_to_Vec3D(wi);
    //vec3d_negate_v(&wi_art);
    vec3d_norm_v(&wi_art);
//    Vec3D wo_art = vec3_to_Vec3D(wo);
//    vec3d_norm_v(&wo_art);
    
    Microsurface * m = [ALLOC_INIT_OBJECT(MicrosurfaceConductor)
                        : NO
                        : NO
                        : 0.1
                        : 0.1
//                        : 1.5
                        ];
    
    [m initRandomNumberGenerator:randomGenerator];
    
    // quadrature with eval p(wi, wo)
    double quadrature_int_x = 0;
    double quadrature_int_y = 0;
    double quadrature_int_z = 0;
    double quadrature_int_x2 = 0;
    double quadrature_int_y2 = 0;
    double quadrature_int_z2 = 0;
    
    for(double theta_o=0 ; theta_o < M_PI ; theta_o += 0.005)
        for(double phi_o=0 ; phi_o < 2.0*M_PI ; phi_o += 0.005)
        {
            Vec3D wo_art = VEC3D(cos(phi_o)*sin(theta_o), sin(phi_o)*sin(theta_o), cos(theta_o));
            
            // stochastic evaluation
            const int N = 10;
            double value_current = 0;
            for(int n=0 ; n<N ; ++n)
            {
                value_current += (double)[m eval:&wi_art: &wo_art:0] / (double) N;
            }
            
            struct vec3 wo = Vec3D_to_vec3(&wo_art);
            const double d = 0.005*0.005*fabs(sin(theta_o)) * value_current;
            quadrature_int_x += d * wo.x;
            quadrature_int_y += d * wo.y;
            quadrature_int_z += d * wo.z;
            quadrature_int_x2 += d * wo.x * wo.x;
            quadrature_int_y2 += d * wo.y * wo.y;
            quadrature_int_z2 += d * wo.z * wo.z;
        }


}


// define some default values
#define     IMP_DEFAULT_INCIDENT_ANGLE          40 DEGREES
#define     IMP_DEFAULT_IMAGE_RESOLUTION        512
#define     IMP_DEFAULT_NUMBER_OF_SAMPLES       10000000


int eric_bsdf(
        int        argc,
        char    ** argv,
        ART_GV   * art_gv
        )
{
    // write your code here

    struct vec3 v1,vv;
    struct vec2 v2;
    
    v1 = vec3(1, 2, -3);
    v2 = vec2(2, 3);

    v1 = normalize(v1);
    
    struct vec3 wi_11 = vec3(4,2,3);
    
    wi_11 = normalize(vec3(5,2,3));
    
    Vec3D v3d_i,v3d_o;
    v3d_i = vec3_to_Vec3D(v1);
    v3d_o = vec3_to_Vec3D(wi_11);
    
    ArcObject <ArpRandomGenerator> * randomGenerator;

    randomGenerator = ARCRANDOMGENERATOR_NEW(
                        RANDOM_SEQUENCE,
                        100000000,
                        ART_GLOBAL_REPORTER
                        );

    
//    MicrosurfaceConductor * mc;
//    mc = [ALLOC_INIT_OBJECT(MicrosurfaceConductor) : YES: YES:0.1:0.1];
//    float a = [mc eval:&v3d_i :&v3d_o :0];
    
//    test_single_scattering(art_gv,randomGenerator);
    test_normalisation(art_gv, randomGenerator);
    
//    test_sample_bsdf(art_gv,randomGenerator);
    
//      test_reflection_sample(art_gv, randomGenerator);
/* ---------------------------------------------------------------------------
    The standard options common to all ART command line applications are
    defined by this macro.
------------------------------------------------------------------------aw- */

    ART_APPLICATION_DEFINE_STANDARD_OPTIONS_WITH_FEATURES(
        "Simulation",
          art_appfeatures_pick_random_seed
        | art_appfeatures_alter_output_filename
        | art_appfeatures_mandatory_output_name
        | art_appfeatures_open_result
        | art_appfeatures_change_isr
        );        
    
/* ---------------------------------------------------------------------------
    The application-specific command line option definitions follow. These 
    appear on the help/usage screen in the order in which they are defined 
    here, after the generic options defined above.
------------------------------------------------------------------------aw- */

    ART_APPLICATION_MAIN_OPTIONS_FOLLOW

    id resOpt =
        [ INTEGER_OPTION
            :   "resolution"
            :   "r"
            :   "<#pixels>"
            :   "override default resolution of image" 
            ];

    id samplesOpt =
        [ INTEGER_OPTION
             :   "samples"
             :   "s"
             :   "<#samples/pixel>"
             :   "override default number of samples"
             ];
    
    id incidentOpt =
        [ FLOAT_OPTION
             :   "incidentAngle"
             :   "ia"
             :   "<degrees>"
             :   "incident light angle"
             ];
    
    // LuNo(24/2/2016): Added for testing/debugging purposes.
    // Chi            : Used for geometry rendering purposes.
    id debugOpt =
        [ INTEGER_OPTION
            :   "debugMode"
            :   "dm"
            :   "<debug mode number>"
            :   "select the debug mode"
            ];
    
    id usePLYOpt =
        [ STRING_OPTION
             :   "ply file"
             :   "ply"
             :   "<filename>"
             :   "set ply geometry file"
             ];
    id gpvOpt =
        [ [ INTEGER_OPTION
            :   "geometryPreview"
            :   "gpv"
            :   "<#samples>"
            :   "normal shaded geometry preview"
            ] withDefaultIntegerValue: 1 ];
    
    

    //   We do the following circus with asprintf only so that we get
    //   the values from the default macros into the usage string.

    char  * description_string;

    asprintf(
	& description_string,
	  "imp -o <outputfile> [options]\n\n"
          "Unless something else is specified via the application options, the\n"
          "result of the simulation will be a raw image %dx%d pixels in size,\n"
          "with %d rays shot at the BSDF sample area at an angle of %3.1f degrees.",
             IMP_DEFAULT_IMAGE_RESOLUTION,
             IMP_DEFAULT_IMAGE_RESOLUTION,
             IMP_DEFAULT_NUMBER_OF_SAMPLES,
             IMP_DEFAULT_INCIDENT_ANGLE * MATH_RAD_TO_DEG
             );
    
    /* ---------------------------------------------------------------------------
     
     Release history:
     
     1.0 alpha 1:
     
     Initial version, empty shell executable imported into git.
     
     ------------------------------------------------------------------------aw- */

    ART_NO_INPUT_FILES_APPLICATION_STARTUP(
        "imp",
        "brute force BRDF simulator, v. 1.0 alpha 1",
        description_string
        );
    
    
    IVec2D  imageSize;
    
    if ( [ resOpt hasBeenSpecified ] )
        imageSize = IVEC2D( [ resOpt integerValue ], [ resOpt integerValue ] );
    else
        imageSize =
        IVEC2D(
//               IMP_DEFAULT_IMAGE_RESOLUTION,
//               IMP_DEFAULT_IMAGE_RESOLUTION
                1024/1,768/1
//                512*4,512*4
               );
    
    //   As always with ART images, the two colour types one specifies when
    //   creating an image have the following meaning:
    //
    //   - imageColourType: the data type the image expects as input
    //   - fileColourType : what we want it to write to disk for us
    
    ArColourType  imageColourType;
    ArColourType  fileColourType;
    
    //   The image colour type is always "native" (i.e. whatever model we
    //   are using), except for RGB images: to avoid issues with different
    //   RGB colour spaces, any RGB results get written to disk as CIE XYZ
    
    if ( art_isr( art_gv ) == arcolour_ut_rgb || art_isr( art_gv ) ==   arcolour_ut_rgb_polarisable )
        imageColourType = arcolour_ciexyz;
    else
        imageColourType = art_isr( art_gv );
    
    //   In the case of a renderer directly writing its output to file, the
    //   image and file colour data types are the same.
    
    fileColourType = imageColourType;
    
    //   Image resolution in DPI. In this day and age of ubiquitious high
    //   DPI displays, you can also put a value of 144.0 there, to get
    //   Retina screen compatible BRDF images. :)
    
    if ( [ gpvOpt hasBeenSpecified ] )
    {
          // this is for gpv mode
        imageColourType = arspectrum_rgba;
        fileColourType  = arspectrum_rgba32;
    }
    
    FVec2D  resolution = FVEC2D( 72.0, 72.0 );
    
    ArnImageInfo * imageInfo =
    [ ALLOC_INIT_OBJECT(ArnImageInfo)
     :   imageSize
     :   imageColourType
     :   fileColourType
     :   resolution
     ];
    
    char * imageFileName = 0;
    
    // this is for gpv mode
    if ( [ gpvOpt hasBeenSpecified ] )
    {
        arstring_pe_copy_add_extension_p(
              ART_APPLICATION_MAIN_FILENAME,
              ARFTIFF_EXTENSION,
            & imageFileName
            );
    }
    else
    {
        arstring_pe_copy_add_extension_p(
              ART_APPLICATION_MAIN_FILENAME,
              ARFARTRAW_EXTENSION,
            & imageFileName
            );
    }
    
    arnoderefdynarray_debugprintf((ArNodeRefDynArray*)ART_APPLICATION_NODESTACK);
    // Push the output (e.g. result) image onto the node stack.
    ART_APPLICATION_NODESTACK_PUSH( [ ALLOC_INIT_OBJECT(ArnFileImage) :imageFileName :imageInfo ] )
    // arnoderefdynarray_debugprintf((ArNodeRefDynArray*)ART_APPLICATION_NODESTACK);
    
    // read in the number of samples to shoot
    unsigned int samplesToShoot = IMP_DEFAULT_NUMBER_OF_SAMPLES;
    if ( [ samplesOpt hasBeenSpecified ] )
        samplesToShoot = [ samplesOpt integerValue ];
    
    // read in the incident angle for theta
    double incidentAngle = IMP_DEFAULT_INCIDENT_ANGLE;
    if ( [ incidentOpt hasBeenSpecified ] )
        incidentAngle = MATH_DEG_TO_RAD * [ incidentOpt doubleValue ];
    
    // read ply filename
    char * ply_filename;
    BOOL   applyPLYfile = NO;
    if( [usePLYOpt hasBeenSpecified] )
    {
        ply_filename = [ usePLYOpt cStringValue ];
        applyPLYfile = YES;
    }
    
    BOOL gpv_MODE = NO;
    if( [gpvOpt hasBeenSpecified] )
    {
        gpv_MODE = YES;
    }
    
    
    
    //LuNo(24/2/2016): Added for testing/debugging purposes.
    // Here a developer may add a new debug mode in order to jumb into his test code.
    if ( [ debugOpt hasBeenSpecified ] )
    {
        switch ( [debugOpt integerValue] )
        {
        
        //------ debug mode 1: Testing heightfield geometry generation and layer intersections. ------//
        case 1:
        {
            debugprintf("imp::main: debug mode %li detected\n", (long)[ debugOpt integerValue ])
            
            

            break;
        }
        
        // The default case: Report the wrong debug mode number and do nothing.
        default:
            debugprintf("imp::main: You have provided an invalid debug mode (%li)\n", (long)[ debugOpt integerValue ])
            break;
                
        }
        
    }
    // This corresponds to the normal operating mode of imp when no debug mode has been specified.
    else
    {

        id imp_actionSequence = 0;

//        imp_actionSequence =
//            ACTION_SEQUENCE(
//            
//            [ ALLOC_INIT_OBJECT(ArnBRDFSampler)
//                :   samplesToShoot
//                :   RANDOM_SEQUENCE
//                :   incidentAngle
//                :   applyPLYfile
//                :   ply_filename
//                :   applyPolarisationFilter
//                :   polarosationFilterAngle
//                :   applyCompensatorFilter
//                :   compensatorFilterAngle
//                :   imageFileName
//                ],
//
//                        
//            [ IMAGECONVERSION_ARTRAW_TO_ARTCSP
//                removeSource: NO
//            ],
//
//            STANDARD_GLOBAL_TONEMAPPING_OPERATOR,
//                        
//            STANDARD_LUMINANCE_CLIPPING,
//
//            [ IMAGECONVERSION_ARTCSP_TO_TIFF
//                removeSource:    YES
//                bitsPerChannel:  8
//            ],
//            
//            ART_VIEWING_ACTION,
//            
//            ACTION_SEQUENCE_END
//
//            );
        
        // Stack machine startup
        //[ imp_actionSequence performOn :ART_APPLICATION_NODESTACK ];
    
    
        RELEASE_OBJECT(imp_actionSequence);
    }

    // Release all instantiated objects.
    RELEASE_OBJECT(imageInfo);
    FREE_ARRAY(description_string);
    FREE_ARRAY(imageFileName);
    
    // If all went well, return zero.
    return 0;   
             
}

ADVANCED_RENDERING_TOOLKIT_MAIN(eric_bsdf)
