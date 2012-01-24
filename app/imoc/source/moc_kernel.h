/** \file moc_kernel.h
 * \ingroup imoc
 * \brief Header file for the Method-of-Characteristics kernel code, 1998
 *
 * \author PA Jacobs
 *
 * \version 26-Sep-1998 : Initial hack
 * \version 09-Jan-2000 : Added CZeroUp, CZeroDown fields
 *
 */

/*
 * Global Constants
 */
#define YES 1
#define NO  0

#define MOC_OK     0
#define MOC_ERROR -1

#define NO_NODE       -1
#define INTERNAL_NODE  1
#define WALL_NODE      2
#define AXIS_NODE      3

/*
 * Global data definitions.
 */
struct NodeData {
   int ValidData;   /* 1 = valid data, 0 otherwise */

   double X;        /* axial coordinate */
   double Y;        /* radial coordinate */
   double Nu;       /* Prandtl-Meyer function, radians */
   double Theta;    /* flow angle, radians */
   double Mach;     /* Mach number */

   double P0;       /* Total pressure */
   double T0;       /* Total temperature */

   int CPlusUp;     /* upstream node along C+ characteristic */
   int CMinusUp;    /* upstream node along C- characteristic */
   int CPlusDown;   /* downstream node along C+ characteristic */
   int CMinusDown;  /* downstream node along C- characteristic */
   int CZeroUp;     /* upstream along streamline */
   int CZeroDown;   /* downstream along streamline */
}; /* end struct NodeData */

/*
 * Function prototypes.
 */
int SetDebugLevel_C( int value );
int GetDebugLevel( void );

int SetGamma_C( double value );
double GetGamma( void );
int SetAxiFlag_C( int value );
int GetAxiFlag( void );
int GetNumberOfNodes( void );
int ValidNode( int id );

int InitMOC( void );

int CreateNode( int id );
int DeleteNode( int id );
struct NodeData * GetNodePtr( int id );
int SetNodeData_C( int id, char *parameter, char *valueString );
char *GetNodeData_C( int id, char *parameter );

int GetNextNodeId( int idStart );
char *ListNodesNear( double x, double y, double tol );
int FindNodesNear( double x, double y, double tol, 
   int id_near_array[], int maxCount );

int SaveNodes( char *FileName );
int LoadNodes( char *FileName );

