/** \file ssp_colours.h
 * \ingroup plot
 * \brief Colour definitions for ssp.
 * 
 * Colour definitions for the Hue-Saturation-Lightness/Brightness
 *                        and Hue-Saturation-Value models
 * All colour values are in the range 0.0 to 1.0.
 *
 */

#define  HSB_RED      0.0
#define  HSB_YELLOW   0.167
#define  HSB_GREEN    0.333
#define  HSB_CYAN     0.5
#define  HSB_BLUE     0.667
#define  HSB_MAGENTA  0.833

struct hsv_colour {
   float h;
   float s;
   float v;
}; /* end struct */

struct hls_colour {
   float h;
   float s;
   float l;
}; /* end struct */

struct rgb_colour {
   float r;
   float g;
   float b;
}; /* end struct */

struct ssp_colour {
   struct hsv_colour hsv;
   struct rgb_colour rgb;
}; /* end struct */

int hsv2rgb( struct hsv_colour *hsv, struct rgb_colour *rgb );
int hls2rgb( struct hls_colour *hls, struct rgb_colour *rgb );

int ssp_setColourFromHSBvalues (
    struct ssp_colour *colour,
    float h, float s, float b);
