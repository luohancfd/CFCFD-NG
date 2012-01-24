/*Copyright (C) Joseph Tang, 2007

    This file is part of OctVCE.

    OctVCE is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OctVCE is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OctVCE.  If not, see <http://www.gnu.org/licenses/>.
*/

/**\file Source file for reading (conversion to VTK version 4.2 ascii format) and writing bodies*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "gio_kernel.h"
#include "gio_io.h"
#include "gio_adts.h"

extern Body_tnode ADTbody;
extern double root_3D_bbox[2][3];
extern int body_count;

/*------------------------------------------------------------------*/

/**\brief A 3D ASCII XML VTK reader to read a parallel polydata file for the list of bodies.  Error checking
   for all possible XML coding errors isn't done; the reader assumes the VTK file has already been properly written.
   Returns number of files successfully processed (added to body ADT).  I don't recommend filenames with spaces.*/

int read_body_files(char pvtp_name[])
{
  int c = 'A'; /*Initialize to some value not EOF*/
  int i;
  int count = 0;
  int strlen1 = 56;
  int strlen2 = 56;
  int strlen3 = 128;
  char temp1[56], temp2[56], temp3[128]; /*Maximum file length*/
  FILE *fin;

  fin = fopen(pvtp_name, "r");
  if(fin == NULL) {
    printf("read_body_files(): Warning. Couldn't open file %s.  Hope this was intended.  Exiting ...\n", pvtp_name);
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
		  if((c == '"') || (isspace(c)))
		    break;

		  temp3[i] = (char) c;
		  c = getc(fin);
		  if(c == EOF) {
		    printf("read_body_files(): Can't seem to read a body.  Exiting ...\n");
		    return(ERROR);
		  }
		  else if((c == '"') || (isspace(c)))
		    break;

		  i++;
		  if(i >= strlen3) {
		    printf("read_body_files(): File name too long.  Exiting ...\n");
		    return(ERROR);
		  }
		}
		  
		if(strlen(temp3) > 0) { /*Now process this body*/
		  if(process_body(temp3, FALSE, NULL, NULL) == ERROR) 
		    return(ERROR);		  
		  else {
		    count++;
		    break;
		  }
		}
	      }
	    } 
	  } 
	  else if(count == 0) { /*Do 1st set of tests*/
	    
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
   errors is not done as it is assumed that the file is already properly written.  

   Importantly, all face connectivities must be defined using right hand system with outward pointing normal.  Sparse geometrical
   checking is performed to determine if polyhedron is proper i.e. if faces declared properly, if some are missing etc. As 
   dimensionality = 3 assumed, won't look at NumberOfComponents data*/

short int process_body(char vtp_name[], short int do_IC, Body *body_ptr, double **IC_bbox)
{ 
  int c = 'A';
  int i, j, k, num_points, num_polys, connect_length, avg_face_size;
  short int vtkfile, vtktype, vtkdata, dealt_points, dealt_connect, dealt_offsets;
  int count;
  char temp1[128], temp2[128];
  int string_length = 128; 
  
  double bbox[6];
  double plane_bbox[3][4];
  Body Solid;
  
  double **points;
  int *face_array;
  int **faces;
  int *face_size;
  int *offsets;
  FILE *fin;
  k = 0;
  vtkfile = FALSE;
  vtktype = FALSE;
  vtkdata = FALSE;
  dealt_points = FALSE;
  dealt_connect = FALSE;
  dealt_offsets = FALSE;
  num_points = 0;
  num_polys = 0;
  avg_face_size = 4; /*For the moment assume each face is a quadrilateral*/

  fin = fopen(vtp_name, "r");
  if(fin == NULL) {
    printf("Process_body(): Can't open file \"%s\".  Exiting ...\n", vtp_name);
    return(ERROR);
  }

  /*First search for the keywords "VTKFile", "type" and "PolyData"*/

  while(c != EOF) {
    c = getc(fin);
    
    if(c == '<') {
      while(!isalpha(c)) 
	c = getc(fin);
      if(c == EOF) {
	printf("Process_body(): Inappropriate end of file.  Exiting ...\n");
	return(ERROR);
      }	
      
      i = 0;
      strncpy(temp1, "\0", string_length);
      while(isalpha(c)) { /*Hopefully this word is 'VTKFile'*/
	temp1[i] = (char) c;
	i++;
	c = getc(fin);
	if(c == EOF) {
	  printf("Process_body(): Inappropriate end of file.  Exiting ...\n");
	  return(ERROR);
	}
	if(i >= string_length) {
	  printf("Process_body(): String overflow.  Exiting ...\n");
	  return(ERROR);
	} 
      }

      if(strncmp(temp1, "VTKFile", strlen(temp1)) == 0) { 
	vtkfile = TRUE;
	break;
      }
    }
  }

  /*Now deduce the vtk file type*/

  while(c != EOF) {
    c = getc(fin);
    
    while((!isalpha(c)) && (c != EOF))
      c = getc(fin);
    if(c == EOF) {
      printf("Process_body(): Inappropriate end of file.  Exiting ...\n");
      return(ERROR);
    }

    i = 0;
    strncpy(temp1, "\0", string_length);
    while(isalpha(c)) { /*Hopefully this word is 'type'*/
      temp1[i] = (char) c;
      i++;
      c = getc(fin);
      if(c == EOF) {
	printf("Process_body(): Inappropriate end of file.  Exiting ...\n");
	return(ERROR);
      }
      if(i >= string_length) {
	printf("Process_body(): String overflow.  Exiting ...\n");
	return(ERROR);
      } 
    }

    if(strncmp(temp1, "type", strlen(temp1)) == 0) { 
      vtktype = TRUE;
      break;
    }    
  }

  /*The type must be "PolyData"*/

  while(c != EOF) {
    c = getc(fin);
    
    if((c == '=') && (c != EOF))
      c = getc(fin);
    if((c == '"') && (c != EOF))
      c = getc(fin);

    while((!isalpha(c)) && (c != EOF))
      c = getc(fin);

    if(c == EOF) {
      printf("Process_body(): Inappropriate end of file.  Exiting ...\n");
      return(ERROR);
    }

    i = 0;
    strncpy(temp1, "\0", string_length);
    while(isalpha(c)) {
      temp1[i] = (char) c;
      i++;
      c = getc(fin);
      if(c == EOF) {
	printf("Process_body(): Inappropriate end of file.  Exiting ...\n");
	return(ERROR);
      }
      if(i >= string_length) {
	printf("Process_body(): String overflow.  Exiting ...\n");
	return(ERROR);
      }
    }

    if(strncmp(temp1, "PolyData", strlen(temp1)) == 0) { 
      vtkdata = TRUE;
      break;
    }    
  }

  if(((vtkfile == FALSE) || (vtktype == FALSE)) && (vtkdata == FALSE)) {
    printf("Process_body(): File %s doesn't seem to be vtkPolyData file.  Check format definitions.  Exiting ...\n",vtp_name);
    return(ERROR);
  }
  
  /*Now we must get piece data, specifically NumberOfPoints and NumberOfPolys*/

  while(c != EOF) {
    c = getc(fin);
    if(c == '<') {
      c = getc(fin);
      break;
    }
  }

  while(!isalpha(c) && (c != EOF))
    c = getc(fin);

  if(isalpha(c)) {
    i = 0;
    strncpy(temp1, "\0", string_length);
    while(isalpha(c)) {
      temp1[i] = (char) c;
      i++;
      c = getc(fin);
      if(c == EOF) {
	printf("Process_body(): Inappropriate end of file.  Exiting ...\n");
	return(ERROR);
      }
      if(i >= string_length) {
	printf("Process_body(): String overflow.  Exiting ...\n");
	return(ERROR);
      }      
    }    
  }

  /*The string should be 'PolyData'*/

  if(strncmp(temp1, "PolyData", strlen(temp1)) != 0) {
    printf("Process_body(): Can't find <PolyData> module.  Exiting ...\n");
    return(ERROR);
  }
  
  while(c != EOF) {
    c = getc(fin);
    if(c == '<') {
      c = getc(fin);
      break;
    }
  }
  
  while((!isalpha(c)) && (c != EOF))
    c = getc(fin);

  i = 0;
  strncpy(temp1, "\0", string_length);
  while(isalpha(c)) {
    temp1[i] = (char) c;
    i++;
    c = getc(fin);
    if(c == EOF) {
      printf("Process_body(): Inappropriate end of file.  Exiting ...\n");
      return(ERROR);
    }
    if(i >= string_length) {
      printf("Process_body(): String overflow.  Exiting ...\n");
      return(ERROR);
    }
  }

  if(strncmp(temp1, "Piece", strlen(temp1)) != 0) {
    printf("Process_body(): Can't find <PolyData> piece.  Exiting ...\n");
    return(ERROR);
  }
  
  while(((num_points == 0) || (num_polys == 0)) && (c != EOF)) { /*Want to find next 2 kewords being "NumberOfPoints" and "NumberOfPolys"*/
  
    while((!isalpha(c)) && (c != EOF)) 
      c = getc(fin);
    
    i = 0;
    strncpy(temp1, "\0", string_length);
    while(isalpha(c)) {
      temp1[i] = (char) c;
      i++;
      c = getc(fin);
      if(c == EOF) {
	printf("Process_body(): Inappropriate end of file.  Exiting ...\n");
	return(ERROR);
      }
      if(i >= string_length) {
	printf("Process_body(): String overflow.  Exiting ...\n");
	return(ERROR);
      }
    }
    
    if((strncmp(temp1, "NumberOfPoints", strlen(temp1)) == 0) || (strncmp(temp1, "NumberOfPolys", strlen(temp1)) == 0)) {
      while((!isdigit(c)) && (c != EOF)) /*Find the next number after this*/
	c = getc(fin);
      
      i = 0;
      strncpy(temp2, "\0", string_length);
      while(isdigit(c)) {
	temp2[i] = (char) c;
	i++;
	c = getc(fin);
	if(c == EOF) {
	  printf("Process_body(): Inappropriate end of file.  Exiting ...\n");
	  return(ERROR);
	}
	if(i >= string_length) {
	  printf("Process_body(): String overflow.  Exiting ...\n");
	  return(ERROR);
	}
      }
      
      if(strncmp(temp1, "NumberOfPoints", strlen(temp1)) == 0) 
	num_points = (int) strtol(temp2, NULL, 10);
      else num_polys = (int) strtol(temp2, NULL, 10);
    }
  } /*Finished looping over finding no. points and no. of polys*/
  
  if((num_points <= 0) || (num_polys <= 0)) {
    printf("Process_body(): Can't seem to find appropriate NumberOfPoints or NumberOfPolys.  Exiting ...\n");
    return(ERROR);
  }
  
  points = malloc(sizeof(double *)*num_points);
  if(points == NULL) {
    printf("Process_body(): Can't allocate points array\n");
    return(ERROR);
  }

  for(i = 0; i < num_points; i++) {
    points[i] = malloc(sizeof(double)*3); /*xyz components*/
    if(points[i] == NULL) {
      printf("Process_body(): Can't allocate points array\n");
      return(ERROR);
    }
  }
  
  offsets = malloc(sizeof(int)*num_polys);
  if(offsets == NULL) {
    printf("Process_body(): Can't allocate offsets array\n");
    return(ERROR);
  }

  face_array = malloc(sizeof(int)*num_polys*avg_face_size); /*Initialize array - may need to lengthen it later*/
  if(face_array == NULL) {
    printf("Process_body(): Can't allocate connectivity array\n");
    return(ERROR);
  }
  
  connect_length = num_polys*avg_face_size;

  /*Now finally we want to obtain connectivity, offset and point info*/
  
  while(((dealt_points == FALSE) || (dealt_offsets == FALSE) || (dealt_connect == FALSE)) && (c != EOF)) {
    
    while((c != '<') && (c != EOF)) 
      c = getc(fin); 
    while((!isalpha(c)) && (c != EOF)) 
      c = getc(fin); 

    i = 0;
    strncpy(temp1, "\0", string_length);
    while(isalpha(c)) {
      temp1[i] = (char) c; 
      i++;
      c = getc(fin);
      if(c == EOF) {
	printf("Process_body(): Inappropriate end of file.  Exiting ...\n");
	return(ERROR);
      }
      if(i >= string_length) {
	printf("Process_body(): String overflow.  Exiting ...\n");
	return(ERROR);
      }      
    } 
    
    if((strncmp(temp1, "Points", strlen(temp1)) == 0) && (dealt_points == FALSE)) {
      
      count = 0;
      while(count <= 2) {
	if(c == '>')
	  count++;
	c = getc(fin);
	if(count == 2)
	  break;
      }
      
      if(c == EOF) {
	printf("Process_body(): Inappropriate end of file in <%s> module.  Exiting ...\n", temp1);
	return(ERROR);	
      }
      
      /*After points piece should expect 2 closing brackets at most.  Next strings after this should be all points*/

      i = 0; j = 0; 

      while(c != EOF) { /*Keep inputting points until left bracket found*/
	if(c == '<')
	  break;
	
	while(isspace(c)) { /*Ignore any whitespace*/
	  c = getc(fin);
	  if(c == '<')
	    break;
	}
	
	if(c == EOF) {
	  printf("Process_body(): Inappropriate end of file in <points> module.  Exiting ...\n");
	  return(ERROR);
	}

	k = 0;
	strncpy(temp1, "\0", string_length);

	if(c != '<') {
	  while(!isspace(c)) {
	    temp1[k] = (char) c;
	    k++;
	    c = getc(fin);
	    if(c == EOF) {
	      printf("Process_body(): Inappropriate end of file in <points> module.  Exiting ...\n");
	      return(ERROR);
	    }
	    else if(c == '<')
	      break;
	    if(k >= string_length) {
	      printf("Process_body(): String overflow.  Exiting ...\n");
	      return(ERROR);
	    } 
	  }
	}
	else break;

	points[i][j] = strtod(temp1, NULL); 
	j++;
	if(j == 3) {
	  j = 0; /*Go to next point*/
	  i++;
	  if(i > num_points) {
	    printf("Process_body(): More written than declared.  Exiting ...\n");
	    return(ERROR);
	  }
	}

	if(c == '<')
	  break;
	c = getc(fin);
      }  
      dealt_points = TRUE;
      
    }
    else if(strncmp(temp1, "Polys", strlen(temp1)) == 0) {
      
      while(c != EOF) {
	if(c == '>')
	  break;
	c = getc(fin);	
      }

      if(c == EOF) {
	printf("Process_body(): Inappropriate end of file in <%s> module.  Exiting ...\n", temp1);
	return(ERROR);	
      }
      
      while(((dealt_offsets == FALSE) || (dealt_connect == FALSE)) && (c != EOF)) {
	
	/*Try to find either the "connectivity" or "offsets" pieces*/
	
	while(c != EOF) {
	  c = getc(fin); 
	  
	  while(!isalpha(c) && (c != EOF)) 
	    c = getc(fin); 
	  
	  if(c == EOF) {
	    printf("Process_body(): Inappropriate end of file in <%s> module.  Exiting ...\n", temp1);
	    return(ERROR);
	  }
	  
	  i = 0;
	  strncpy(temp2, "\0", string_length);
	  while(isalpha(c)) {
	    temp2[i] = (char) c; 
	    c = getc(fin);
	    i++;
	    if(c == EOF) {
	      printf("Process_body(): Inappropriate end of file in <Polys> module.  Exiting ...\n");
	      return(ERROR);
	    }
	    if(k >= string_length) {
	      printf("Process_body(): String overflow.  Exiting ...\n");
	      return(ERROR);
	    }
	  } 
	  
	  if((strncmp(temp2, "connectivity", strlen(temp2)) == 0) || (strncmp(temp2, "offsets", strlen(temp2)) == 0))
	    break;
	}

	if(strncmp(temp2, "connectivity", strlen(temp2)) == 0) { 
	  while(c != EOF) {
	    if(c == '>') { /*After closing right bracket we're ready to input numbers*/
	      c = getc(fin);
	      break;
	    }
	    c = getc(fin); 
	  }
	  
	  if(c == EOF) {
	    printf("Process_body(): Inappropriate end of file in <Polys> module.  Exiting ...\n");
	    return(ERROR);
	  }

	  i = 0;

	  while(c != EOF) { /*Now get the integers*/
	    if(c == '<')
	      break;

	    while(isspace(c)) {
	      c = getc(fin);
	      if(c == '<')
		break;
	    }

	    if(c == EOF) {
	      printf("Process_body(): Inappropriate end of file in <Polys> module.  Exiting ...\n");
	      return(ERROR);
	    }

	    k = 0;
	    strncpy(temp1, "\0", string_length);

	    if(c != '<') {
	      while(!isspace(c)) {
		temp1[k] = (char) c;
		k++;
		c = getc(fin);
		if(c == EOF) {
		  printf("Process_body(): Inappropriate end of file in <Polys> module.  Exiting ...\n");
		  return(ERROR);
		}
		else if(c == '<')
		  break;
		if(k >= string_length) {
		  printf("Process_body(): String overflow.  Exiting ...\n");
		  return(ERROR);
		}
	      }
	    }
	    else break;
	    
	    face_array[i] = (int) strtol(temp1, NULL, 10); 
	    i++; 

	    if(c == '<')
	      break;
	    c = getc(fin);

	    if(i >= num_polys*avg_face_size) { /*Need to extend array size*/
	      face_array = realloc(face_array, sizeof(int)*(i+1)); 
	      connect_length = i+1;
	      if(face_array == NULL) {
		printf("Process_body(): Can't allocate connectivity array.  Exiting ...\n");
		return(ERROR);
	      }		
	    }
	  }

	  dealt_connect = TRUE; 
	}
	else if(strncmp(temp2, "offsets", strlen(temp2)) == 0) {

	  while(c != EOF) {
	    if(c == '>') {
	      c = getc(fin);
	      break;
	    }
	    c = getc(fin);	    
	  }

	  if(c == EOF) {
	    printf("Process_body(): Inappropriate end of file in <Polys> module.  Exiting ...\n");
	    return(ERROR);
	  }

	  i = 0;

	  while(c != EOF) { /*Now get the integers*/
	    if(c == '<')
	      break;
	  
	    while(isspace(c)) {
	      c = getc(fin);
	      if(c == '<')
		break;
	    }

	    if(c == EOF) {
	      printf("Process_body(): Inappropriate end of file in <Polys> module.  Exiting ...\n");
	      return(ERROR);
	    }
	    else if(c == '<')
	      break;

	    k = 0;
	    strncpy(temp1, "\0", string_length);

	    if(c != '<') {
	      while(!isspace(c)) {
		temp1[k] = (char) c; 
		k++;
		c = getc(fin); 

		if(c == EOF) {
		  printf("Process_body(): Inappropriate end of file in <Polys> module.  Exiting ...\n");
		  return(ERROR);
		}
		else if(c == '<')  
		  break;
	      
		if(k >= string_length) {
		  printf("Process_body(): String overflow.  Exiting ...\n");
		  return(ERROR);
		}
	      }
	    }
	    else break;
	    
	    offsets[i] = (int) strtol(temp1, NULL, 10); 
	    if(offsets[i] == 0) {
	      printf("Process_body(): Abnormal offset found.  Exiting ...\n");
	      return(ERROR);
	    }

	    i++;

	    if(c == '<')
	      break;
	    c = getc(fin);	  
	  }

	  dealt_offsets = TRUE;
	}        
      }
    }
  }

  fclose(fin);
  
  /*Now deal with actual faces from connectivity and offset info*/

  faces = malloc(sizeof(int *)*num_polys);
  if(faces == NULL) {
    printf("Process_body(): Can't allocate faces array.  Exiting ...\n");
    return(ERROR);
  }

  face_size = malloc(sizeof(int)*num_polys);
  if(face_size == NULL) {
    printf("Process_body(): Can't allocate face_size array.  Exiting ...\n");
    return(ERROR);
  }

  j = 0;
  for(i = 0; i < num_polys; i++) {
    k = offsets[i]; /*A measure of face size*/
    
    face_size[i] = k - j;
    if(face_size[i] == 0) {
      printf("Process_body(): Found a face of zero size.  Exiting ...\n");
      return(ERROR);
    }

    faces[i] = malloc(sizeof(int)*(k-j));
    if(faces[i] == NULL) {
      printf("Process_body(): Can't allocate faces array.  Exiting ...\n");
      return(ERROR);
    }
    
    for(count = 0; count < (k-j); count++) 
      faces[i][count] = face_array[j + count];
        
    j = offsets[i];
  }

  /*Finally finish processing the body data structure*/

  /*Get the bounding box (as point in 6D)*/

  bbox[0] = points[0][0];
  bbox[1] = points[0][1];
  bbox[2] = points[0][2];
  bbox[3] = bbox[0];
  bbox[4] = bbox[1];
  bbox[5] = bbox[2];

  for(i = 1; i < num_points; i++) {
    if(points[i][0] < bbox[0])
      bbox[0] = points[i][0];
    else if(points[i][0] > bbox[3])
      bbox[3] = points[i][0];

    if(points[i][1] < bbox[1])
      bbox[1] = points[i][1];
    else if(points[i][1] > bbox[4])
      bbox[4] = points[i][1];

    if(points[i][2] < bbox[2])
      bbox[2] = points[i][2];
    else if(points[i][2] > bbox[5])
      bbox[5] = points[i][2];
  }

  if(((bbox[0] <= root_3D_bbox[1][0]) && (bbox[1] <= root_3D_bbox[1][1]) && (bbox[2] <= root_3D_bbox[1][2])) &&
     ((bbox[3] >= root_3D_bbox[0][0]) && (bbox[4] >= root_3D_bbox[0][1]) && (bbox[5] >= root_3D_bbox[0][2]))) {
    /*Bounding box intersects octree root's bounding box*/

    for(i = 0; i < 6; i++) { /*Trim bounding boxes that exceed octree root's domain to domain limits*/      			
      if(i < 3) {
	if(bbox[i] < root_3D_bbox[0][i])
	  bbox[i] = root_3D_bbox[0][i];
      }
      else if(bbox[i] > root_3D_bbox[1][i%3])
	bbox[i] = root_3D_bbox[1][i%3];      
    }

    /*Bounding box on x plane - is YZ plane, extract Y and Z co-ords of 3D bounding box on plane*/
    plane_bbox[0][0] = bbox[1];
    plane_bbox[0][1] = bbox[2];
    plane_bbox[0][2] = bbox[4];
    plane_bbox[0][3] = bbox[5];

    /*Bounding box on y plane*/
    plane_bbox[1][0] = bbox[0];
    plane_bbox[1][1] = bbox[2];
    plane_bbox[1][2] = bbox[3];
    plane_bbox[1][3] = bbox[5];

    /*Bounding box on z plane*/
    plane_bbox[2][0] = bbox[0];
    plane_bbox[2][1] = bbox[1];
    plane_bbox[2][2] = bbox[3];
    plane_bbox[2][3] = bbox[4];

    Solid = malloc(sizeof(struct body));
    if(Solid == NULL) {
      printf("Process_body(): Can't allocate Body data structure.  Exiting ...\n");
      return(ERROR);
    }

    Solid -> num_points = num_points;
    Solid -> num_faces = num_polys;
    Solid -> points = points;
    Solid -> faces = faces;
    Solid -> face_size = face_size;

    Solid -> xy_bound[0] = bbox[0]; Solid -> xy_bound[2] = bbox[3]; /*We can do the xy_bound here rather than recompute it*/
    Solid -> xy_bound[1] = bbox[1]; Solid -> xy_bound[3] = bbox[4];

    /*First build ADT of faces*/

    if(build_face_ADT(Solid, FALSE) == ERROR) {
      printf("Process_body(): Couldn't build face ADT.  Exiting ...\n");
      return(ERROR);
    }

    /*Now add body to body ADTs - thus process_body() ultimately builds the ADT of bodies*/
		  
    if(do_IC != TRUE) { 
      if(add_to_body_ADT(Solid, bbox, plane_bbox, &ADTbody, body_count) == NOERROR)
	body_count++;
    }
    else { 
      *body_ptr = Solid; /*Used when reading in initial conditions for OctVCE*/ 
      IC_bbox[0][0] = bbox[0]; IC_bbox[0][1] = bbox[1]; IC_bbox[0][2] = bbox[2]; 
      IC_bbox[1][0] = bbox[3]; IC_bbox[1][1] = bbox[4]; IC_bbox[1][2] = bbox[5]; 
    }

    Solid = NULL; 
  }
  
  free(face_array);
  free(offsets);
  
  return(NOERROR);
}

/*------------------------------------------------------------------*/

