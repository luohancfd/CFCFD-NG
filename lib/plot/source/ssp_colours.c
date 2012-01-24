/** \file  ssp_colours.c
 * \ingroup plot
 * \brief Utility routines related to colours.
 */

#include <stdio.h>
#include <math.h>
#include "ssp_colours.h"

/*-----------------------------------------------------------------*/

float ssp_RGB( float H, float M1, float M2 ) {
    float Value;

    if (H < 0.0) {
        H += 1.0;
    } else if (H > 1.0) {
        H -= 1.0;
    } /* end if */

    if (H < 1.0/6.0) {
        Value = M1 + (M2 - M1) * H * 6.0;
    } else if (H < 0.5) {
        Value = M2;
    } else if (H < 2.0/3.0) {
        Value = M1 + (M2 - M1) * (2.0/3.0 - H) * 6.0;
    } else {
        Value = M1;
    } /* end if */

    return Value;
} /* end function ssp_RGB */
 
/*-----------------------------------------------------------------*/

int hls2rgb( struct hls_colour *hls, struct rgb_colour *rgb ) {
    /*
     * Purpose:
     * Convert Hue Saturation Lightness/Brightness to Red Green Blue
     * All values are in the range [0.0 .. 1.0]
     *
     * Reference:
     * J Foley et al Principles of Computer Graphics...
     *
     * The version in 
     * D. F. Rogers
     * Procedural Elements for Computer Graphics
     * McGraw Hill 1985
     * appears to have a few errors.
     *
     * Am still not sure of this coding...
     */
    float S, H, L, M1, M2;
   
    S = hls->s;  /* Saturation */
    H = hls->h;  /* Hue */
    L = hls->l;  /* Lightness */
   
    if ( S == 0.0 ) {
        /* 
         * Achromatic case, set level of grey 
         */
        rgb->r = L;
        rgb->g = L;
        rgb->b = L;
    } else {
        if (L <= 0.5) {
            M2 = L * (1.0 + S);
        } else {
            M2 = L + S + L * S;
        } /* end if */
        M1 = 2.0 * L - M2;
        /* 
         * Determine levels of primary colours. 
         */
        if (H > 1.0 || H < 0.0) {
            H = 0.0;
        } /* end if */

        rgb->r = ssp_RGB( H+1.0/3.0, M1, M2 );
        rgb->g = ssp_RGB( H,         M1, M2 );
        rgb->b = ssp_RGB( H-1.0/3.0, M1, M2 );
    } /* end if */

    #if 0
    printf( "hls = %g %g %g, rgb = %g %g %g\n", 
        hls->h, hls->l, hls->s, rgb->r, rgb->g, rgb->b );
    #endif
    return 0;
} /* end function hls2rgb */

/*-----------------------------------------------------------------*/
 
int hsv2rgb( struct hsv_colour *hsv, struct rgb_colour *rgb ) {
    /*
     * Purpose:
     * Convert Hue Saturation Value to Red Green Blue
     * All values are in the range [0.0 .. 1.0]
     *
     * Reference:
     * D. F. Rogers
     * Procedural Elements for Computer Graphics
     * McGraw Hill 1985
     */
    float S, H, V, F, M, N, K;
    int   I;
   
    S = hsv->s;  /* Saturation */
    H = hsv->h;  /* Hue */
    V = hsv->v;  /* value */
   
    if ( S == 0.0 ) {
        /* 
         * Achromatic case, set level of grey 
         */
        rgb->r = V;
        rgb->g = V;
        rgb->b = V;
    } else {
        /* 
         * Determine levels of primary colours. 
         */
        if (H >= 1.0) {
            H = 0.0;
        } else {
            H = H * 6;
        } /* end if */
        I = (int) H;   /* should be in the range 0..5 */
        F = H - I;     /* fractional part */

        M = V * (1 - S);
        N = V * (1 - S * F);
        K = V * (1 - S * (1 - F));

        if (I == 0) { rgb->r = V; rgb->g = K; rgb->b = M; }
        if (I == 1) { rgb->r = N; rgb->g = V; rgb->b = M; }
        if (I == 2) { rgb->r = M; rgb->g = V; rgb->b = K; }
        if (I == 3) { rgb->r = M; rgb->g = N; rgb->b = V; }
        if (I == 4) { rgb->r = K; rgb->g = M; rgb->b = V; }
        if (I == 5) { rgb->r = V; rgb->g = M; rgb->b = N; }
    } /* end if */

    #if 0
    printf( "hsv = %g %g %g, rgb = %g %g %g\n", 
        hsv->h, hsv->s, hsv->v, rgb->r, rgb->g, rgb->b );
    #endif
    return 0;
} /* end function hsv2rgb */

/*------------------------------------------------------------*/

int ssp_setColourFromHSBvalues (
    struct ssp_colour *colour,
    float h, float s, float b) {
    /*
     * Purpose...
     * -------
     * Set the parameters for the HSB colour model as used
     * in postscript.  
     *
     * Input...
     * -----
     * h   : 0.0 <= h <= 1.0, hue
     *       0.0 = pure red, 1/3 = pure green, 2/3 = pure blue
     * s   : 0.0 <= s <= 1.0, saturation
     *       0.0 = grey (or no colour), 1.0 = maximum colour concentration
     * b   : 0.0 <= b <= 1.0, brightness
     *       0.0 = black, 0.5 = colours, 1.0 = white
     *
     */

    /*
     * Bring values to within range; 0.0 <= h < 1.0
     */
    if (h < 0.0) {
       h = h + 1.0 + ((int) h);
    } else if (h >= 1.0) {
       h = h - ((int) h);
    } /* end if */
    
    if (s <= 0.0) {
       s = 0.0;
    } else if (s >= 1.0) {
       s = 1.0;
    } /* end if */

    if (b <= 0.0) {
       b = 0.0;
    } else if (b >= 1.0) {
       b = 1.0;
    } /* end if */

    /*
     * Determine RGB equivalent.
     */
    colour->hsv.h = h;
    colour->hsv.s = s;
    colour->hsv.v = b;
    hsv2rgb( &(colour->hsv), &(colour->rgb) );
    
    return 0;
} /* end ssp_setColourFromHSBvalues */


