/** \file extract.c
 *  \ingroup plot
 *  \brief Extract columns from a GENERIC plotting file -- no longer used.
 *
 * Read a GENERIC data file containing data in the form
 * (x, f1, f2, ... fn) and produce a file containing
 * specified columns.
 * The data format may be deduced from the section of 
 * code shown below which reads the file.
 * Note that several blocks of data may be included in the
 * one GENERIC file and these are combined for the output.
 *
 * \author PA Jacobs
 */

/*-------------------------------------------------------*/

#include "../../util/source/compiler.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*-------------------------------------------------------*/

int main ()
{
#define  BLOCKS  10
int    ivar, nvar, ix, iy, nnx[BLOCKS], nblock, jb;
#define  DIMVAR  15
char   VarName[DIMVAR][32];
double **Var[BLOCKS];

#define  NCHAR  256
char   title[NCHAR];
char   txt[NCHAR], DataFileName[32], PlotFileName[32];
FILE   *dfp, *pfp;
int    screen, pfile;
int    icol1, icol2;


printf ("--------------------------------------\n");
printf ("Generic Line Plotting: Extract Columns\n");
printf ("--------------------------------------\n\n");
printf ("Enter GENERIC Data File Name : ");
scanf ("%s", DataFileName);
if ((dfp = fopen(DataFileName, "r")) == NULL)
   {
   printf ("Could not open data file: %s\n", DataFileName);
   exit (-1);
   }

/*
 * ***********************
 * * Read the data file. *
 * ***********************
 */

fgets (title, NCHAR, dfp);
printf ("\nTitle: %s\n", title);

fgets (txt, NCHAR, dfp);
sscanf (txt, "%d", &nvar);
printf ("Number of variables: %d\n", nvar);
if (nvar > DIMVAR)
   {
   printf ("DIMVAR is too small; rebuild code.\n");
   exit (-1);
   }

for (ivar = 0; ivar < nvar; ++ivar)
   {
   fgets (txt, NCHAR, dfp);
   sscanf (txt, "%s", VarName[ivar]);
   printf ("Variable[%d]: %s\n", ivar, VarName[ivar]);
   }

/* NOTE : we start reading file without discarding line ends. */
fscanf (dfp, "%d", &nblock);
printf ("Number of blocks: %d\n", nblock);
if (nblock > BLOCKS)
   {
   printf ("BLOCKS is too small; rebuild code.\n");
   exit (-1);
   }

for (jb = 0; jb < nblock; ++jb)
   {
   fscanf (dfp, "%d", &(nnx[jb]) );
   printf ("nnx = %d\n", nnx[jb] );
   
   printf ("Allocating memory for data block[%d]...\n", jb);
   if ((Var[jb] = (double **) malloc (nvar * sizeof(double *)))
      == NULL)
      {
      printf ("Allocation failure\n");
      exit (-1);
      }
   for (ivar = 0; ivar < nvar; ++ivar)
      {
      if ( (Var[jb][ivar] = 
	   (double *) malloc(nnx[jb] * sizeof(double))) 
           == NULL )
         {
         printf ("Allocation failure: %d\n", ivar);
         exit (-1);
         }
      }

   for (ix = 0; ix < nnx[jb]; ++ix)
      for (ivar = 0; ivar < nvar; ++ivar)
         fscanf (dfp, "%lf", &Var[jb][ivar][ix]);

   }   /* for (jb... */

/*
 ********************************************
 * Decide which variables are to be dumped. *
 ********************************************
 */

Make_a_Dump:
printf ("\n-------------------------------\n");

PickVariable:
for (ivar = 0; ivar < nvar; ++ivar)
   printf ("%d=%s ", ivar, VarName[ivar]);
printf ("\nColumn 1 - Which variable (-1(exit), 0 ... %d): ", nvar-1);
scanf ("%d", &icol1);
if (icol1 < 0) exit (0);
if (icol1 >= nvar) goto PickVariable;
printf ("Column 2 - Which variable (-1(exit), 0 ... %d): ", nvar-1);
scanf ("%d", &icol2);
if (icol2 < 0) exit (0);
if (icol2 >= nvar) goto PickVariable;

/*
 * Now dump the selected data.
 */

printf ("Enter Dump File Name : ");
scanf ("%s", PlotFileName);
if ((pfp = fopen(PlotFileName, "w")) == NULL)
   {
   printf ("Could not open DUMP file: %s\n", PlotFileName);
   exit (-1);
   }

printf ("Write data...\n");
for (jb = 0; jb < nblock; ++jb)
   {
   for (ix = 0; ix < nnx[jb]; ++ix)
      {
      fprintf (pfp, "%e %e\n", Var[jb][icol1][ix], 
			       Var[jb][icol2][ix]);
      }
   }


goto Make_a_Dump;

}  /* end of extract.c */

