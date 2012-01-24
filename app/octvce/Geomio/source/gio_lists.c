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

/**\file Source file for list operations*/

#include <stdio.h>
#include <stdlib.h>
#include "gio_kernel.h"
#include "gio_lists.h"

/*------------------------------------------------------------------*/

/**\brief Add node to list of bounding boxes*/

short int add_to_list_bbox(Body_tnode tnode, List_bbox *L) 
{ 
  List_bbox node = malloc(sizeof(struct list_bbox)); 
  
  if(node != NULL)
    {
      node -> tnode = tnode; /*Address of bounding box's body now stored on list*/
      node -> next = (*L); 

      (*L) = node; /*Front of list is updated*/

      return(NOERROR);
    }
  else return(ERROR);
}

/*------------------------------------------------------------------*/

/**\brief Destroy a list - deallocate all nodes and makes pointer to list front NULL*/

void destroy_list_bbox(List_bbox *L) 
{
  List_bbox temp;
  temp = (*L);  

  while(temp != NULL)
    {
      temp = temp -> next;
      free(*L);
      *L = temp;
    }
}

/*------------------------------------------------------------------*/
