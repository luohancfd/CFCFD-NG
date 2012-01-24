/** \file moc_unitproc.h
 * \ingroup imoc
 * \brief Header file for basic unit processes.
 */

int InteriorNode( int node1, int node2, int node4 );

int InsertNode( int node1, int node2, int node4, double alpha );

int CMinusWallNode( int iw, int node1, int node4 );
int CPlusWallNode( int iw, int node2, int node4 );

int CPlusFreeBndyNode( int node0, int node2, int node4 );
int CMinusFreeBndyNode( int node0, int node1, int node4 );

int AddStreamNode( int node0, int node1, int node2, int node4,
                   int test_only );
int StepStreamNode( int node0, int node4, double dL );

int InterpolateNode( double x, double y, double R, int node4 );

