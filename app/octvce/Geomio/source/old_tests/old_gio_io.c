/**\file Source file for reading (conversion to VTK format) and writing bodies*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "gio_kernel.h"
#include "gio_io.h"

/*------------------------------------------------------------------*/

/**\brief Write body info (just arbitrary collections of hexahedrons for the moment) to VTK mesh as
   polygonal faces*/

int write_bodies(void) 
{
  FILE *fout, *fin;
  int num_bodies = 0;
  int i, j, k, m, vtx_num;
  char buffer[512];
  short int finished_reading = FALSE;
  
  double lb[3]; double ub[3]; /*Bounding box of body - for bounding box == body for hexahedrons*/

  fin = fopen("ov_bodies.par", "r");
  
  if(fin == NULL)
    {
      printf("write_bodies() could not open file ov_bodies.par\n");
      return(ERROR);
    }

  /*Ignore first 3 lines*/
  fgets(buffer, sizeof(buffer), fin);
  fgets(buffer, sizeof(buffer), fin);
  fgets(buffer, sizeof(buffer), fin);

  while(finished_reading == FALSE)
    {
      fgets(buffer, sizeof(buffer), fin);
      if(sscanf(buffer, "lb: %lf", &lb[0]) == 0) /*No more bodies to read*/
	{
	  printf("write_bodies() has no more correct bodies to read - the no. bodies successfully read is %d\n", num_bodies);
	  finished_reading = TRUE;
	}
      else if(sscanf(buffer, "lb: %lf %lf %lf", &lb[0], &lb[1], &lb[2]) < 3)
	{
	  printf("write_bodies() read lower bounds incorrectly\n");
	  fclose(fin);
	  return(ERROR);
	}
      else /*First line read correctly*/
	{
	  fgets(buffer, sizeof(buffer), fin);

	  if(sscanf(buffer, "ub: %lf %lf %lf", &ub[0], &ub[1], &ub[2]) < 3)
	    {
	      printf("write_bodies() read upper bounds incorrectly\n");
	      fclose(fin);
	      return(ERROR);
	    }
	  else if((ub[0] < lb[0]) || (ub[1] < lb[1]) || (ub[2] < lb[2]))
	    {		
	      printf("write_bodies() saw a weird body for body_num %d.  Please check it.\n", num_bodies+1);
	      fclose(fin);
	      return(ERROR);
	    }
	  else /*Body has been read in correctly*/
	    {
	      num_bodies++;
	      fgets(buffer, sizeof(buffer), fin); /*Space in between bodies*/
	    }
	}
    }

  fclose(fin); /*Go back to start to begin writing output*/
  fin = fopen("ov_bodies.par", "r");

  if(num_bodies != 0)
    {
      fout = fopen("body_plot.vtk", "w");
      if(fout == NULL)
	{
	  printf("write_bodies() could not open file body_plot\n");
	  return(ERROR);
	}

      fprintf(fout, "# vtk DataFile Version 2.0\n");
      fprintf(fout, "Hexahedral bodies for plotting\n");
      fprintf(fout, "ASCII\n");
      fprintf(fout, "DATASET POLYDATA\n");

      fprintf(fout, "POINTS %d float\n", num_bodies*8);
      finished_reading = FALSE;

      fgets(buffer, sizeof(buffer), fin);
      fgets(buffer, sizeof(buffer), fin);
      fgets(buffer, sizeof(buffer), fin); /*Ignore first 3 lines again*/
      
      while(finished_reading == FALSE)
	{
	  fgets(buffer, sizeof(buffer), fin);
	  if(sscanf(buffer, "lb: %lf", &lb[0]) == 0) /*No more bodies to read*/
	    finished_reading = TRUE;
	  else 
	    {
	      sscanf(buffer, "lb: %lf %lf %lf", &lb[0], &lb[1], &lb[2]);
	      fgets(buffer, sizeof(buffer), fin);
	      sscanf(buffer, "ub: %lf %lf %lf", &ub[0], &ub[1], &ub[2]);
	      fgets(buffer, sizeof(buffer), fin);

	      /*Get verticies (numbered like a VTK hexahedron) from bounding box*/
	      fprintf(fout, "%20.12e %20.12e %20.12e\n", lb[0], lb[1], lb[2]);
	      fprintf(fout, "%20.12e %20.12e %20.12e\n", ub[0], lb[1], lb[2]);
	      fprintf(fout, "%20.12e %20.12e %20.12e\n", ub[0], ub[1], lb[2]);
	      fprintf(fout, "%20.12e %20.12e %20.12e\n", lb[0], ub[1], lb[2]);
	      fprintf(fout, "%20.12e %20.12e %20.12e\n", lb[0], lb[1], ub[2]);
	      fprintf(fout, "%20.12e %20.12e %20.12e\n", ub[0], lb[1], ub[2]);
	      fprintf(fout, "%20.12e %20.12e %20.12e\n", ub[0], ub[1], ub[2]);
	      fprintf(fout, "%20.12e %20.12e %20.12e\n", lb[0], ub[1], ub[2]); 
	    }
	}

      fprintf(fout, "POLYGONS %d %d\n", num_bodies*6, num_bodies*6*5); 
      vtx_num = 0;
      for(i = 0; i < num_bodies; i++)
	{
	  for(j = 0; j < 6; j++) /*Cycle through each face*/
	    {
	      fprintf(fout, "%d ", 4);
	      
	      for(k = 0; k < 4; k++) /*Cycle through each face vertex*/
		{
		  switch(j)
		    {
		    case 0:
		      m = 0; break;
		    case 1:
		      m = 4; break;
		    case 2: 
		      if(k <= 1)
			m = 0;
		      else if(k == 2) 
			m = 3;
		      else m = 1;
		      break;
		    case 3:
		      if(k <= 1)
			m = 2;
		      else if(k == 2)
			m = 5;
		      else m = 3;
		      break;
		    case 4:
		      if(k == 0)
			m = 0;
		      else if(k == 1)
			m = 3;
		      else if(k == 2)
			m = 5;
		      else m = 0;
		      break;
		    case 5:
		      if(k <= 1)
			m = 1;
		      else if(k == 2)
			m = 4;
		      else m = 2;
		      break;
		    }
		  fprintf(fout, "%d ", k + m + vtx_num); /*All points connected anti-clockwise (from face exterior)*/
		}
	      
	      fprintf(fout, "\n");
	    }
	  vtx_num = vtx_num + 8; /*Get next set of points/verticies belonging to next body*/
	}

      fprintf(fout, "\n");
      fclose(fout);
    }
  else if(num_bodies == 0)
    printf("write_bodies() didn't read any bodies\n");

  fclose(fin);
  
  return(num_bodies);
}

/*------------------------------------------------------------------*/

/**\brief A 3D ASCII XML VTK reader to read a parallel polydata file for the list of bodies.  Error checking
   for all possible XML coding errors isn't done; the reader assumes the VTK file has already been properly written.
   Returns number of files successfully processed (added to body ADT)*/

int read_body_files(char pvtp_name[])
{
  int c = 'A'; /*Initialize to some value not EOF*/
  int i;
  int count = 0;
  int strlen1 = 56;
  int strlen2 = 56;
  int strlen3 = 256;
  char temp1[56], temp2[56], temp3[256];
  FILE *fin;

  fin = fopen(pvtp_name, "r");
  if(fin == NULL) {
    printf("read_body_files(): Couldn't open file to read.  Hope this was intended.  Exiting ...\n");
    return(NOERROR);
  }  

  while(c != EOF) {
    
    c = getc(fin);
      
    if(c == '<') { /*After the left bracket we want to search for some keywords*/
	
      while(c != EOF) { 
	c = getc(fin);
	  
	if(isalpha(c)) { 

	  /*To start with, we want to look for word "VTKFile", then "type", then ""PPolyData"".  Then afterward
	    we want to get the piece sources*/

	  i = 0;
	  strncpy(temp1, "\0", strlen1); /*Get next word - it seems safter to copy over the whole array*/
	  while(isalpha(c)) {
	    temp1[i] = (char) c;
	    c = getc(fin);
	    i++;
	    if(i > strlen1) { /*Most likely the result of a coding error*/
	      printf("read_body_files(): String overflow.  Exiting ...\n");
	      return(ERROR);
	    }
	  }
	    
	  if(count >= 1) { /*Do 2nd set of tests*/
	    
	    if(strncmp(temp1, "Piece", strlen(temp1)) == 0) { /*So 1st keyword "Piece" is found*/
		
	      while(isspace(c))
		c = getc(fin);
		
	      i = 0;
	      strncpy(temp2, "\0", strlen2);
	      if(isalpha(c)) {
		while(isalpha(c)) {
		  temp2[i] = (char) c;
		  c = getc(fin);
		  i++;
		}
		if(i > strlen2) { 
		  printf("read_body_files(): String overflow.  Exiting ...\n");
		  return(ERROR);
		}
	      }
		
	      if(strncmp(temp2, "Source", strlen(temp2)) == 0) { /*So 2nd keyword "Source" found*/
	
		while(c != EOF) {
		  if(c == '=')
		    break;
		  else c = getc(fin);
		}
		    
		while(c != EOF) {
		  if(c == '"')
		    break;
		  else c = getc(fin);
		}

		c = getc(fin);
		while(isspace(c))
		  c = getc(fin);

		/*Now read in the files*/
		i = 0;
		strncpy(temp3, "\0", strlen3);
		while(c != EOF) {
		  if(c == '"')
		    break;

		  temp3[i] = (char) c;
		  c = getc(fin);
		  if(c == EOF) {
		    printf("read_body_files(): Can't seem to read a body.  Exiting ...\n");
		    return(ERROR);
		  }

		  i++;
		  if(i >= strlen3) {
		    printf("read_body_files(): File name too long.  Exiting ...\n");
		    return(ERROR);
		  }
		}
		  
		if(strlen(temp3) > 0) { /*Now process this body*/
		  if(process_body(temp3) == ERROR) 
		    return(ERROR);		  
		  else {
		    count++;
		    break;
		  }
		}
	      }
	    } 
	  } else if(count == 0) { /*Do 1st set of tests*/
	    
	    if(strncmp(temp1, "VTKFile", strlen(temp1)) == 0) { /*So the 1st keyword "VTKFile" is found*/
		
	      while(isspace(c))
		c = getc(fin);  /*Want to get to next word, which should be "type"*/		
	      
	      i = 0;
	      strncpy(temp2, "\0", strlen2); /*Wipe the 2nd string array and look for the 2nd keyword*/
	      if(isalpha(c)) {
		while(isalpha(c)) {
		  temp2[i] = (char) c;
		  c = getc(fin);
		  i++;
		}
		if(i > strlen2) { 
		  printf("read_body_files(): String overflow.  Exiting ...\n");
		  return(ERROR);
		}
	      }
	      
	      if(strncmp(temp2, "type", strlen(temp2)) == 0) { /*2nd keyword is found*/
		
		while(c != EOF) {
		  if(c == '=')
		    break;
		  else c = getc(fin);
		}
		    
		while(c != EOF) {
		  if(c == '"')
		    break;
		  else c = getc(fin);
		}

		c = getc(fin);
		while(isspace(c))
		  c = getc(fin);
		
		/*Now finally search for 3rd keyword, which should be "PPolyData"*/
		
		i = 0; 
		strncpy(temp3, "\0", strlen3);
		if(isalpha(c)) {
		  while(isalpha(c)) {
		    temp3[i] = (char) c;  
		    c = getc(fin);
		    i++;
		  }
		  if(i > strlen3) { 
		    printf("read_body_files(): String overflow.  Exiting ...\n");
		    return(ERROR);
		  }
		}
			
		break; /*We're done with the first set of comparisons*/
	      }
	    }	
	  }    
	}
      }          

      if(count == 0) {	
	if((strncmp(temp1, "VTKFile", strlen(temp1)) != 0) || (strncmp(temp2, "type", strlen(temp2)) != 0) || 
	   (strncmp(temp3, "PPolyData", strlen(temp3)) != 0)) {
	  printf("read_body_files(): Error - file doesn't appear to be a vtkPolyData file - check format definition.  Exiting ...\n");
	  return(ERROR);
	}
	else count++; /*Move onto next set of tests*/	
      }
    }
  }

  fclose(fin);
  
  return(count - 1);
}

/*------------------------------------------------------------------*/

/**\brief A 3D ASCII XML VTK reader to read a vtp file and initialize a body structure (faces, points and edges)
   as defined in gio_kernel.h corresponding to the description in the file.  Error checking for all possible XML coding
   errors is not done as it is assumed that the file is already properly written.  As dimensionality=3 assumed, won't look
   at NumberOfComponents data*/

short int process_body(char vtp_name[])
{
  int c = 'A';
  int i, j, k, num_points, num_polys, connect_length, avg_face_size;
  short int exit, points_now, polys_now, connect_now, offset_now;
  int count = 0;
  char temp1[72], temp2[72], temp3[72];
  int string_length = 72;
  
  double **points;
  int *face_array;
  int **faces;
  int *offsets;
  FILE *fin;
  
  exit = FALSE;
  num_points = 0;
  num_polys = 0;
  avg_face_size = 4; /*For the moment assume each face is a quadrilateral*/

  fin = fopen(vtp_name, "r");
  if(fin == NULL) {
    printf("Process_body(): Can't open file %s.  Exiting ...\n", vtp_name);
    return(ERROR);
  }

  while(c != EOF) {
    c = getc(fin);
    
    if(c == '<') { /*After left bracket we want to search for some keywords*/
	
      while(c != EOF) { 
	c = getc(fin);
	  
	if(isalpha(c)) { 

	  /*To start with, we want to look for word "VTKFile", then "type", then ""PolyData"".  Then afterward
	    we want to get geometry data*/

	  i = 0;
	  strncpy(temp1, "\0", string_length); /*Get next word - it seems safter to copy over the whole array*/
	  while(isalpha(c)) {
	    temp1[i] = (char) c;
	    c = getc(fin);
	    i++;
	    if(i > string_length) { /*Most likely the result of a coding error*/
	      printf("Process_body(): String overflow.  Exiting ...\n");
	      return(ERROR);
	    }
	  }

	  if(count >= 2) {
	    if((strncmp(temp1, "Points", strlen(temp1)) == 0) || (strncmp(temp1, "Polys", strlen(temp1)) == 0)) { 
	      /*Find either Points or Polys module*/

	      if(strncmp(temp1, "Points", strlen(temp1)) == 0) {
		points_now = TRUE;
		polys_now = FALSE;
		points = malloc(sizeof(double *)*num_points);
		if(points == NULL) {
		  printf("Process_body(): Can't allocate points array.  Exiting ...\n");
		  return(ERROR);
		}
		for(i = 0; i < num_points; i++) {
		  points[i] = malloc(sizeof(double)*3);
		  if(points[i] == NULL) {
		    printf("Process_body(): Can't allocate points array.  Exiting ...\n");
		    return(ERROR);
		  }
		}		  
	      }
	      else {
		points_now = FALSE;
		polys_now = TRUE;
		face_array = malloc(sizeof(int)*num_polys*avg_face_size); /*For the moment assume each face is a quadrilateral*/
		if(face_array == NULL) {
		  printf("Process_body(): Can't allocate faces array.  Exiting ...\n");
		  return(ERROR);
		}		
		offsets = malloc(sizeof(int)*num_polys);
		if(offsets == NULL) {
		  printf("Process_body(): Can't allocate offsets array.  Exiting ...\n");
		  return(ERROR);
		}
	      }		
	      
	      while(c != EOF) {
		c = getc(fin);
		if(c == EOF) {
		  printf("Process_body(): Can't find end of <%s> module.  Exiting ...\n",temp1);
		  return(ERROR);
		}
		if(c == '>') /*1st right bracket found*/
		  break;
	      }

	      if(points_now == TRUE) {
		while(c != EOF) {
		  c = getc(fin);
		  if(c == EOF) {
		    printf("Process_body(): Can't find end of <%s> module.  Exiting ...\n",temp1);
		    return(ERROR);
		  }
		  if(c == '>') /*2nd and last right bracket found - all values after this before another left bracket should be
				 input into array*/
		    break;
		}

		/*Now start putting numbers into the array*/

		i = 0; /*Start at 0th array position*/
		j = 0;
		while((c != EOF) && (c != '<')) { /*Look for next left bracket as signal termination*/
		  c = getc(fin);
		  if(c == EOF) {
		    printf("Process_body(): Whilst looking through <%s> module, found inappropriate end of file.  Exiting ...\n", temp1);
		    return(ERROR);
		  }
		  else if(!isspace(c)) {
		    strncpy(temp2, "\0", string_length);
		    
		    k = 0;
		    while((!isspace(c)) && (c != EOF)) {
		      temp2[k] = (char) c;
		      c = getc(fin);
		      if(c == EOF) {
			printf("Process_body(): Looked through <%s> module, and found inappropriate end of file.  Exiting ...\n",temp1);
			return(ERROR);
		      }
		      k++;
		      if(k > string_length) { 
			printf("Process_body(): String overflow.  Exiting ...\n");
			return(ERROR);
		      }		      
		    }
		    
		    points[i][j] = strtod(temp1, NULL);
		    j++;
		    if(j == 3) {
		      j = 0; /*Start with next point*/
		      i++; 
		    }		    
		  }
		}	  
	      }
	      else if(polys_now == TRUE) {
		while((c != EOF) && (c != '>')) {
		  c = getc(fin);
		  if(c == EOF) {
		    printf("Process_body(): Can't find end of <%s> module.  Exiting ...\n",temp1);
		    return(ERROR);
		  }

		  if(isalpha(c)) {
		    i = 0;
		    strncpy(temp2, "\0", string_length);
		    while(isalpha(c)) {
		      temp2[i] = (char) c;
		      if(i > string_length) {
			printf("Process_body(): String overflow.  Exiting ...\n");
			return(ERROR);
		      }
		      i++;
		    }
		    
		    if(strncmp(temp2, "Name", strlen(temp2)) == 0) {
		      while(isspace(c))
			c = getc(fin);
		      if(c == '=')
			c = getc(fin);
		      while(isspace(c))
			c = getc(fin);
		      if(c == '"')
			c = getc(fin);
		      while(isspace(c))
			c = getc(fin);
		      if(c == EOF) {
			printf("Process_body(): Can't find end of <%s> module.  Exiting ...\n",temp1);
			return(ERROR);
		      }

		      if(isalpha(c)) {
			i = 0;
			strncpy(temp2, "\0", string_length);
			while(c != EOF) {
			  temp2[i] = (char) c;
			  c = getc(fin);
			  if(c == EOF) {
			    printf("Process_body(): At module <%s> found inappropriate end of file\n",temp1);
			    return(ERROR);
			  }
			  i++;
			  if(i > string_length) {
			    printf("Process_body(): String overflow.  Exiting ...\n");
			    return(ERROR);
			  }
			}

			if(strncmp(temp2, "connectivity", strlen(temp2)) == 0) { 
			  connect_now = TRUE;
			  offset_now = FALSE;
			}
			else if(strncmp(temp2, "offsets", strlen(temp2)) == 0) {
			  offset_now = TRUE;
			  connect_now = FALSE;
			}
			else {
			  printf("Process_body(): At module <%s> found inappropriate name\n",temp1);
			  return(ERROR);
			}

			while((c != EOF) && (c != '>'))
			  c = getc(fin);
			if(c == EOF) {
			  printf("Process_body(): At module <%s> found inappropriate end of file\n",temp1);
			  return(ERROR);
			}

			if(connect_now == TRUE) { /*So now do connectivity array - may have to extend it*/

			  i = 0;			  
			  while((c != EOF) && (c != '<')) {
			    c = getc(fin);
			    if(c == EOF) {
			      printf("Process_body(): At module <%s> found inappropriate end of file\n",temp1);
			      return(ERROR);
			    }

			    if(!isspace(c)) {
			      j = 0;
			      strncpy(temp2, "\0", string_length);
			      while(!isspace(c)) {
				temp2[j] = (char) c;
				j++;
				if(c == EOF) {
				  printf("Process_body(): At module <%s> found inappropriate end of file\n",temp1);
				  return(ERROR);
				}
				if(j > string_length) {
				  printf("Process_body(): String overflow.  Exiting ...\n");
				  return(ERROR);
				}
			      }

			      face_array[i] = (int) strtol(temp2, NULL, 10);
			      i++;

			      if(i >= num_polys*avg_face_size) { /*We need to extend the array*/
				face_array = realloc(face_array, (i+1)*sizeof(int));
				if(face_array == NULL) {
				  printf("Process_body(): Can't allocate face_array.  Exiting ...\n");
				  return(ERROR);
				}
				connect_length = i+1;
			      }
			    }
			  }
			}
			else if(offset_now == TRUE) { /*So now do offset array*/
			  
			  i = 0;			  
			  while((c != EOF) && (c != '<')) {
			    c = getc(fin);
			    if(c == EOF) {
			      printf("Process_body(): At module <%s> found inappropriate end of file\n",temp1);
			      return(ERROR);
			    }

			    if(!isspace(c)) {
			      j = 0;
			      strncpy(temp2, "\0", string_length);
			      while(!isspace(c)) {
				temp2[j] = (char) c;
				j++;
				if(c == EOF) {
				  printf("Process_body(): At module <%s> found inappropriate end of file\n",temp1);
				  return(ERROR);
				}
				if(j > string_length) {
				  printf("Process_body(): String overflow.  Exiting ...\n");
				  return(ERROR);
				}
			      }

			      offsets[i] = (int) strtol(temp2, NULL, 10);
			      i++;
			    }
			  }
			}

			/*Now we've got to look at last data array (offsets or connectivity)*/

			while((c != EOF) && (c != '>'))
			  c = getc(fin);
			if(c == EOF) {
			  printf("Process_body(): Can't find end of <%s> module.  Exiting ...\n",temp1);
			  return(ERROR);
			}

			while((c != EOF) && (c != '>')) {
			  c = getc(fin);
			  if(c == EOF) {
			    printf("Process_body(): Can't find end of <%s> module.  Exiting ...\n",temp1);
			    return(ERROR);
			  }

			  if(isalpha(c)) {
			    i = 0;
			    strncpy(temp2, "\0", string_length);
			    while(isalpha(c)) {
			      temp2[i] = (char) c;
			      if(i > string_length) {
				printf("Process_body(): String overflow.  Exiting ...\n");
				return(ERROR);
			      }
			      i++;
			    }
		    
			    if(strncmp(temp2, "Name", strlen(temp2)) == 0) {
			      while(isspace(c))
				c = getc(fin);
			      if(c == '=')
				c = getc(fin);
			      while(isspace(c))
				c = getc(fin);
			      if(c == '"')
				c = getc(fin);
			      while(isspace(c))
				c = getc(fin);
			      if(c == EOF) {
				printf("Process_body(): Can't find end of <%s> module.  Exiting ...\n",temp1);
				return(ERROR);
			      }

			      if(isalpha(c)) {
				i = 0;
				strncpy(temp2, "\0", string_length);
				while(c != EOF) {
				  temp2[i] = (char) c;
				  c = getc(fin);
				  if(c == EOF) {
				    printf("Process_body(): At module <%s> found inappropriate end of file\n",temp1);
				    return(ERROR);
				  }
				  i++;
				  if(i > string_length) {
				    printf("Process_body(): String overflow.  Exiting ...\n");
				    return(ERROR);
				  }
				}

				if(strncmp(temp2, "connectivity", strlen(temp2)) == 0) { 
				  connect_now = TRUE;
				  offset_now = FALSE;
				}
				else if(strncmp(temp2, "offsets", strlen(temp2)) == 0) {
				  offset_now = TRUE;
				  connect_now = FALSE;
				}
				else {
				  printf("Process_body(): At module <%s> found inappropriate name\n",temp1);
				  return(ERROR);
				}

				while((c != EOF) && (c != '>'))
				  c = getc(fin);
				if(c == EOF) {
				  printf("Process_body(): At module <%s> found inappropriate end of file\n",temp1);
				  return(ERROR);
				}

				if(connect_now == TRUE) { /*So now do connectivity array - may have to extend it*/

				  i = 0;			  
				  while((c != EOF) && (c != '<')) {
				    c = getc(fin);
				    if(c == EOF) {
				      printf("Process_body(): At module <%s> found inappropriate end of file\n",temp1);
				      return(ERROR);
				    }

				    if(!isspace(c)) {
				      j = 0;
				      strncpy(temp2, "\0", string_length);
				      while(!isspace(c)) {
					temp2[j] = (char) c;
					j++;
					if(c == EOF) {
					  printf("Process_body(): At module <%s> found inappropriate end of file\n",temp1);
					  return(ERROR);
					}
					if(j > string_length) {
					  printf("Process_body(): String overflow.  Exiting ...\n");
					  return(ERROR);
					}
				      }

				      face_array[i] = (int) strtol(temp2, NULL, 10);
				      i++;

				      if(i >= num_polys*avg_face_size) { /*We need to extend the array*/
					face_array = realloc(face_array, (i+1)*sizeof(int));
					if(face_array == NULL) {
					  printf("Process_body(): Can't allocate face_array.  Exiting ...\n");
					  return(ERROR);
					}
					connect_length = i+1;
				      }
				    }
				  }
				}
				else if(offset_now == TRUE) { /*So now do offset array*/
			  
				  i = 0;			  
				  while((c != EOF) && (c != '<')) {
				    c = getc(fin);
				    if(c == EOF) {
				      printf("Process_body(): At module <%s> found inappropriate end of file\n",temp1);
				      return(ERROR);
				    }

				    if(!isspace(c)) {
				      j = 0;
				      strncpy(temp2, "\0", string_length);
				      while(!isspace(c)) {
					temp2[j] = (char) c;
					j++;
					if(c == EOF) {
					  printf("Process_body(): At module <%s> found inappropriate end of file\n",temp1);
					  return(ERROR);
					}
					if(j > string_length) {
					  printf("Process_body(): String overflow.  Exiting ...\n");
					  return(ERROR);
					}
				      }

				      offsets[i] = (int) strtol(temp2, NULL, 10);
				      i++;
				    }
				  }
				}			
			      }		      
			    }
			  }
			}
		      }		      
		    }
		  }
		}
	      }  
	    }
	  }
	  else if(count == 1) { /*Look for PolyData module*/
	    if(strncmp(temp1, "PolyData", strlen(temp1)) == 0) { /*Found the module*/ 
	      while(c != EOF) {
		c = getc(fin); 
		if(c == EOF) {
		  printf("Process_body(): Can't find <PolyData> module end.  Exiting ...\n");
		  return(ERROR);
		}
		else if(c == '<')
		  break;
	      }
	      
	      while(c != EOF) { /*Now find NumberOfPoints and NumberOfPolys piece - ignore the rest*/
		c = getc(fin); 
		while(isspace(c))
		  c = getc(fin);

		if(c == EOF) {
		  printf("Process_body(): Can't find <PolyData> piece.  Exiting ...\n");
		  return(ERROR);
		}

		printf("%c\n",c);

		i = 0; /*Look for "Piece" keyword*/
		strncpy(temp1, "\0", string_length); 
		while(isalpha(c)) {
		  temp1[i] = (char) c; 
		  c = getc(fin);
		  i++;
		  if(i > string_length) { 
		    printf("Process_body(): String overflow.  Exiting ...\n");
		    return(ERROR);
		  }
		} 
		
		if(strncmp(temp1, "Piece", strlen(temp1)) == 0) { /*Found the right keyword - next look for NumberOfPoints
								    and NumberOfPolys*/ 
		  while(isspace(c))
		    c = getc(fin);

		  i = 0;
		  strncpy(temp1, "\0", string_length);
		  if(isalpha(c)) {
		    while(isalpha(c)) {
		      temp1[i] = (char) c; 
		      c = getc(fin);
		      i++;
		    } 
		    if(i > string_length) { 
		      printf("Process_body(): String overflow.  Exiting ...\n");
		      return(ERROR);
		    }

		    if(strncmp(temp1, "NumberOfPoints", strlen(temp1)) == 0) { 
		      while(isspace(c))
			c = getc(fin);
		      if(c == '=')
			c = getc(fin);
		      while(isspace(c))
			c = getc(fin);
		      if(c == '"')
			c = getc(fin);
		      while(isspace(c))
			c = getc(fin);

		      if(isdigit(c)) {
			i = 0;
			strncpy(temp1, "\0", string_length);
			while(isdigit(c)) {
			  temp1[i] = (char) c; 
			  c = getc(fin);
			  i++;
			}
			if(i > string_length) { 
			  printf("Process_body(): String overflow.  Exiting ...\n");
			  return(ERROR);
			}
		      }

		      num_points = (int) strtol(temp1, NULL, 10);
		      if(num_points <= 0) {
			printf("Process_body(): NumberOfPoints <= 0.  Exiting ...\n");
			return(ERROR);
		      } 
		    }
		    else if(strncmp(temp1, "NumberOfPolys", strlen(temp1)) == 0) { 
		      while(isspace(c))
			c = getc(fin);
		      if(c == '=')
			c = getc(fin);
		      while(isspace(c))
			c = getc(fin);
		      if(c == '"')
			c = getc(fin);
		      while(isspace(c))
			c = getc(fin);

		      if(isdigit(c)) {
			i = 0;
			strncpy(temp1, "\0", string_length);
			while(isdigit(c)) {
			  temp1[i] = (char) c;
			  c = getc(fin);
			  i++;
			}
			if(i > string_length) { 
			  printf("Process_body(): String overflow.  Exiting ...\n");
			  return(ERROR);
			}
		      }

		      num_polys = (int) strtol(temp1, NULL, 10);
		      if(num_points <= 0) {
			printf("Process_body(): NumberOfPolys <= 0.  Exiting ...\n");
			return(ERROR);
		      }
		    }
		    else {
		      while(!isalpha(c))
			c = getc(fin);
		    }

		    /*Want to get the next one - either NumberOfPolys or NumberOfPoints*/

		    while(c != EOF) {
		      c = getc(fin);

		      if(isalpha(c)) { /*Should start with an 'N'*/
			i = 0; 
			strncpy(temp1, "\0", string_length); 
			while(isalpha(c)) {
			  temp1[i] = (char) c; 
			  c = getc(fin);
			  i++;
			  if(i > string_length) { 
			    printf("Process_body(): String overflow.  Exiting ...\n");
			    return(ERROR);
			  }
			} printf("\n"); 

			if(strncmp(temp1, "NumberOfPoints", strlen(temp1)) == 0) {
			  while(isspace(c))
			    c = getc(fin);
			  if(c == '=')
			    c = getc(fin);
			  while(isspace(c))
			    c = getc(fin);
			  if(c == '"')
			    c = getc(fin);
			  while(isspace(c))
			    c = getc(fin);

			  if(isdigit(c)) {
			    i = 0;
			    strncpy(temp1, "\0", string_length);
			    while(isdigit(c)) {
			      temp1[i] = (char) c;
			      c = getc(fin);
			      i++;
			    }
			    if(i > string_length) { 
			      printf("Process_body(): String overflow.  Exiting ...\n");
			      return(ERROR);
			    }
			  }

			  num_points = (int) strtol(temp1, NULL, 10);
			  if(num_points <= 0) {
			    printf("Process_body(): NumberOfPoints <= 0.  Exiting ...\n");
			    return(ERROR);
			  }
			}
			else if(strncmp(temp1, "NumberOfPolys", strlen(temp1)) == 0) {
			  while(isspace(c))
			    c = getc(fin);
			  if(c == '=')
			    c = getc(fin);
			  while(isspace(c))
			    c = getc(fin);
			  if(c == '"')
			    c = getc(fin);
			  while(isspace(c))
			    c = getc(fin);

			  if(isdigit(c)) {
			    i = 0;
			    strncpy(temp1, "\0", string_length);
			    while(isdigit(c)) {
			      temp1[i] = (char) c;
			      c = getc(fin);
			      i++;
			    }
			    if(i > string_length) { 
			      printf("Process_body(): String overflow.  Exiting ...\n");
			      return(ERROR);
			    }
			  }

			  num_polys = (int) strtol(temp1, NULL, 10);
			  if(num_points <= 0) {
			    printf("Process_body(): NumberOfPolys <= 0.  Exiting ...\n");
			    return(ERROR);
			  }
			}
			else {
			  while(isalpha(c))
			    c = getc(fin);
			}

			exit = TRUE; /*Found PolyData module and got some points and polygons*/
			break;
		      }
		    }		      
		  }		    
		}
		
		if(exit == TRUE) {
		  exit = FALSE;
		  break;
		}
	      }
	    }
	  } /*Have now hopefully gotten no. points and polys*/
	  else if(count == 0) { /*First test if this file is in the correct format*/
	    if(strncmp(temp1, "VTKFile", strlen(temp1)) == 0) { /*So the 1st keyword "VTKFile" is found*/
		
	      while(isspace(c))
		c = getc(fin);  /*Want to get to next word, which should be "type"*/		
	      
	      i = 0;
	      strncpy(temp2, "\0", string_length); /*Wipe the 2nd string array and look for the 2nd keyword*/
	      if(isalpha(c)) {
		while(isalpha(c)) {
		  temp2[i] = (char) c;
		  c = getc(fin);
		  i++;
		}
		if(i > string_length) { 
		  printf("Process_body(): String overflow.  Exiting ...\n");
		  return(ERROR);
		}
	      }
	      
	      if(strncmp(temp2, "type", strlen(temp2)) == 0) { /*2nd keyword is found*/
		
		while(c != EOF) {
		  if(c == '=')
		    break;
		  else c = getc(fin);
		}
		    
		while(c != EOF) {
		  if(c == '"')
		    break;
		  else c = getc(fin);
		}

		c = getc(fin);
		while(isspace(c))
		  c = getc(fin);
		
		/*Now finally search for 3rd keyword, which should be "PolyData"*/
		
		i = 0; 
		strncpy(temp3, "\0", string_length);
		if(isalpha(c)) {
		  while(isalpha(c)) {
		    temp3[i] = (char) c;  
		    c = getc(fin);
		    i++;
		  }
		  if(i > string_length) { 
		    printf("Process_body(): String overflow.  Exiting ...\n");
		    return(ERROR);
		  }
		}
			
		break; /*We're done with the first set of comparisons*/
	      }
	    }
	  }
	}
      }

      if(count == 0) {	
	if((strncmp(temp1, "VTKFile", strlen(temp1)) != 0) || (strncmp(temp2, "type", strlen(temp2)) != 0) || 
	   (strncmp(temp3, "PolyData", strlen(temp3)) != 0)) {
	  printf("Process_body(): Error - file doesn't appear to be a vtkPolyData file - check format definition.  Exiting ...\n");
	  return(ERROR);
	}
	else count++; /*Move onto next set of tests*/
      }
      else if(count == 1) {
	if((num_points <= 0) || (num_polys <= 0)) {
	  printf("Process_body(): Error - NumberOfPoints %d, NumberOfPolys %d\n",num_points, num_polys);
	  return(ERROR);
	}
	else count++;
      }      
    }
  }


  fclose(fin);

  printf("Points are - \n");
  for(i = 0; i < num_points; i++)
    printf("%g %g %g\n", points[i][0],points[i][1],points[i][2]);

  printf("Connectivities are - ");
  for(i = 0; i < connect_length; i++)
    printf("%d ", face_array[i]);
  printf("\n");

  printf("Offsets are - \n");
  for(i = 0; i < num_polys; i++)
    printf("%d ", offsets[i]);
  printf("\n");

  for(i = 0; i < num_points; i++) {
    free(points[i]);
    free(faces[i]);
  }

  free(points);
  free(faces);

  /*We may not need to free points and faces*/

  free(face_array);
  free(offsets);
  
  return(NOERROR);
}

/*------------------------------------------------------------------*/

