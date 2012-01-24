/* playbez.h
 * Header file, function declarations for playbez.c
 * -----------------------------------------------------------
 */

#ifndef COMPILER_H
/* #  error We need compiler.h to be included. */
#  include "compiler.h"
#endif

#ifndef GEOM_H
/* #  error We need geom.h to be included. */
#  include "geom.h"
#endif

#ifndef SSP_H
/* #  error We need ssp.h to be included. */
#  include "ssp.h"
#endif


#if (PROTO)
   /* Full function prototypes... */
   int root_menu (void);
   int file_menu (void);
   int edit_menu (void);
   int plot_menu (void);

   int edit_bezier (void);
   int edit_data_set (void);

   int save_data_set (void);
   int read_data_set (void);
   int save_bezier (void);
   int read_bezier (void);

   initialize_picture (void);
   change_view (void);
   redraw_all (void);

   int fit_bezier (void);
   double penalty_function ( int npar, double par[] );
   int fit_spline (void);

#else
   /* K&R style prototypes... */
   int root_menu ();
   int file_menu ();
   int edit_menu ();
   int plot_menu ();

   int edit_bezier ();
   int edit_data_set ();

   int save_data_set ();
   int read_data_set ();
   int save_bezier ();
   int read_bezier ();

   initialize_picture ();
   change_view ();
   redraw_all ();

   int fit_bezier ();
   double penalty_function ();
   int fit_spline ();

#endif