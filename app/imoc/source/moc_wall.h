/** \file moc_wall.h
 * \ingroup imoc
 * \brief Header file for the wall-definition functions.
 *
 * \author PA Jacobs
 *
 * \version 04-Jan-2000 : Initial hack
 *
 */

/*
 * Global definitions.
 */

/*
 * Function prototypes.
 */
int InitWall( void );
int WallIsPresent( int iw );
int WallGetNumberOfPoints( int iw );
double WallGetPointX( int iw, int ip );
double WallGetPointY( int iw, int ip );
int WallDeletePoints( int iw );
int WallAddPoint( int iw, double x, double y );

double WallPosX( int iw, double t );
double WallPosY( int iw, double t );
double WallSlope( int iw, double t );
double WallFindT( int iw, 
   double x0, double y0, double cos_th, double sin_th );

int SaveWall( int iw, char *FileName );
int LoadWall( int iw, char *FileName );

