/* moc_kernel.i 
 * SWIG interface definition file.
 */

%module imoc
%{
#include "moc_kernel.h"
#include "moc_gasdynamic.h"
#include "moc_unitproc.h"
#include "moc_wall.h"
%}

/*
 * Make the following functions visible to the Tcl interpreter
 */

extern int SetDebugLevel_C( int value );
extern int GetDebugLevel( void );

extern int SetGamma_C( double value );
extern double GetGamma( void );
extern int SetAxiFlag_C( int value );
extern int GetAxiFlag( void );
extern int GetNumberOfNodes( void );
extern int ValidNode( int id );

extern int InitMOC( void );

extern int CreateNode( int id );
extern int DeleteNode( int id );
extern int SetNodeData_C( int id, char *parameter, char *valueString );
extern char *GetNodeData_C( int id, char *parameter );

extern int GetNextNodeId( int idStart );
extern char *ListNodesNear( double x, double y, double tol );

extern int SaveNodes( char *FileName );
extern int LoadNodes( char *FileName );

extern double T0_T(double M, double g);
extern double P0_P(double M, double g);

extern double NuFromM(double M, double g);
extern double MFromNu( double Nu, double g );
extern double MachAngle( double M );

extern int InteriorNode( int node1, int node2, int node4 );
extern int InsertNode( int node1, int node2, int node4, double alpha );
extern int CMinusWallNode( int iw, int node1, int node4 );
extern int CPlusWallNode( int iw, int node2, int node4 );
extern int CPlusFreeBndyNode( int node0, int node2, int node4 );
extern int CMinusFreeBndyNode( int node0, int node1, int node4 );
extern int AddStreamNode( int node0, int node1, int node2, int node4,
                          int test_only );
extern int StepStreamNode( int node0, int node4, double dL );
extern int InterpolateNode( double x, double y, double R, int node4 );

extern int InitWall( void );
extern int WallIsPresent( int iw );
extern int WallGetNumberOfPoints( int iw );
extern double WallGetPointX( int iw, int ip );
extern double WallGetPointY( int iw, int ip );
extern int WallDeletePoints( int iw );
extern int WallAddPoint( int iw, double x, double y );

extern double WallPosX( int iw, double t );
extern double WallPosY( int iw, double t );
extern double WallSlope( int iw, double t );
extern double WallFindT( int iw, 
   double x0, double y0, double cos_th, double sin_th );

extern int SaveWall( int iw, char *FileName );
extern int LoadWall( int iw, char *FileName );

