/*!
 
 @file      imp.m
 
 @brief     This file is the main ART file for the imp project.
 
 @author    Chi Wang (chi@cgg.mff.cuni.cz), Lukas Novosad (novosad@cgg.mff.cuni.cz)
 @date
 
 */
// ART library imports
#import "AdvancedRenderingToolkit.h"

#import "Eric_BSDF/MicrosurfaceSlope.h"

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
    
    v1 = vec3(1, 2, 3);
    v2 = vec2(2, 3);

    v1 = normalize(v1);
    
    struct vec3 wi_11 = vec3(4,2,3);
    
    wi_11 = normalize(vec3(5,2,3));
    
    Vec3D v3d;
    v3d = vec3_to_Vec3D(v1);
    
    
    MicrosurfaceSlopeBeckmann * msb;
    msb = [MicrosurfaceSlopeBeckmann init];


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
