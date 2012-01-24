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

/**\file Source file on list functions*/

#include <stdio.h>
#include <stdlib.h>
#include "ov_kernel.h"
#include "ov_lists.h"
#include "ov_setgeom.h"

extern int number_of_threads;
extern int stepnow;

/*------------------------------------------------------------------*/

/**\brief Add node to front of list of leaf nodes*/

short int add_to_leaf_list_front(Cart_cell C, List_leaf *L, List_leaf *Head, List_leaf *Tail)
{                                                
  List_leaf node = malloc(sizeof(struct list_leaf)); 
  
  if(node != NULL)
    {
      node -> cell_loc = C; /*Cell's address now stored on list*/

      node -> thread_num = 0; /*Default thread no. (in serial)*/

      C -> Leaf_list_loc = node; /*Let cell store node's address*/
      
      node -> next = *L; 

      if((*L) != NULL)
	{ 
	  node -> prev = (*L) -> prev;

	  if((*L) -> prev != NULL)
	    {
	      (*L) -> prev -> next = node;
	      (*L) -> prev = node;
	    }
	  else 
	    {
	      (*L) -> prev = node; 
	      *Head = node;
	    }	      
	}
      else 
	{
	  *Head = node;
	  *Tail = node;
	  node -> prev = NULL;
	}
      
      return(NOERROR);
    }
  else return(ERROR); /*Can't allocate list node*/
}

/*------------------------------------------------------------------*/

/**\brief Add node to front of list of linked cells - same as add_to_list_leaf_front()*/

short int add_to_merge_list_front(Cart_cell C, List_merge *L)
{                                                
  List_merge node = malloc(sizeof(struct list_merge)); 
  
  if(node != NULL)
    {
      node -> cell_loc = C; 
      
      node -> next = *L; 
      
      if((*L) != NULL)
	(*L) -> prev = node;
      
      node -> prev = NULL;
      C -> Merge_list_loc = node; 
      
      *L = node; 
      
      return(NOERROR);
    }
  else return(ERROR); 
}

/*------------------------------------------------------------------*/

/**\brief Delete node from list of leaf nodes*/

void delete_from_leaf_list(Cart_cell C, List_leaf *Head, List_leaf *Tail)
{                                      
  List_leaf node = C -> Leaf_list_loc;

  if(node -> prev == NULL) /*At beginning of list*/
    {
      if(node -> next == NULL) /*Only item on list*/
	{
	  (*Head) = NULL; /*So now Head becomes empty list*/
	  (*Tail) = NULL;

	  free(node);
	}
      else 
	{
	  node -> next -> prev = NULL;
	  (*Head) = node -> next; /*Update pointer to beginning of list*/
	  
	  free(node);
	}
    }  
  else if (node -> next == NULL) /*At end of list*/
    {
      node -> prev -> next = NULL;
      (*Tail) = node -> prev;

      free(node);
    }
  else /*In the middle of the list*/
    {
      node -> prev -> next = node -> next;
      node -> next -> prev = node -> prev;

      free(node);
    }

  C -> Leaf_list_loc = NULL; /*Let deleted cell's pointer to list of leaf nodes be NULL*/
}

/*------------------------------------------------------------------*/

/**\brief Delete node from list of linked cells*/

void delete_from_merge_list(Cart_cell C, List_merge *L)
{                                      
  List_merge node = C -> Merge_list_loc;

  if(node -> prev == NULL) 
    {
      if(node -> next == NULL) 
	{
	  (*L) = NULL;
	  free(node);
	}
      else 
	{
	  node -> next -> prev = NULL;
	  (*L) = node -> next; 	  
	  free(node);
	}
    }  
  else if (node -> next == NULL)
    {
      node -> prev -> next = NULL;
      free(node);
    }
  else 
    {
      node -> prev -> next = node -> next;
      node -> next -> prev = node -> prev;
      free(node);
    }

  C -> Merge_list_loc = NULL; 
}

/*------------------------------------------------------------------*/

/**\brief Add a vertex to global list of verticies*/

short int add_to_Vtxlist(Vertex V, List_vtx_glob *L, List_vtx_glob *Tail)
{                                                
  List_vtx_glob node = malloc(sizeof(struct list_vtx_glob)); 
  
  if(node != NULL)
    {
      node -> vtx_loc = V; 
      
      node -> delete = UNKNOWN;
      
      node -> next = *L; 
      
      if((*L) != NULL)
	(*L) -> prev = node;
      else *Tail = node;
      
      node -> prev = NULL;
      V -> Vtxlist_loc = node; 
      
      *L = node; 
      
      return(NOERROR);
    }
  else return(ERROR); 
}

/*------------------------------------------------------------------*/

/**\brief Delete vertex from global list of verticies*/

void delete_from_Vtxlist(Vertex V, List_vtx_glob *L, List_vtx_glob *Tail)
{  
  List_vtx_glob node = V -> Vtxlist_loc;
  
  if(node -> prev == NULL) 
    {
      if(node -> next == NULL) 
	{
	  (*L) = NULL;
	  (*Tail) = NULL;

	  free(node);
	}
      else 
	{
	  node -> next -> prev = NULL;
	  (*L) = node -> next; 	  
	  free(node);
	}
    }  
  else if (node -> next == NULL)
    {
      node -> prev -> next = NULL;
      (*Tail) = node -> prev;

      free(node);
    }
  else 
    {
      node -> prev -> next = node -> next;
      node -> next -> prev = node -> prev;

      free(node);
    }

  V -> Vtxlist_loc = NULL;
}

/*------------------------------------------------------------------*/

/**\brief Go through sublists and ensure no cells to be coarsened are split across boundaries*/

void modify_lists(List_leaf Heads[], List_leaf Tails[])
{
  int i;
  List_leaf node;
  Cart_cell Parent;

  for(i = 0; i < (number_of_threads-1); i++)
    {
      if(Tails[i] -> cell_loc -> adapt_flag == REALLY_COARSEN)
	{
	  Parent = Tails[i] -> cell_loc -> parent;
	  node = Tails[i];
	  while(node -> next -> cell_loc -> parent == Parent) /*Keep going until next cell isn't sibling*/
	    node = node -> next;

	  Tails[i] = node;
	  Heads[i+1] = Tails[i] -> next;
	}
    }
}

/*------------------------------------------------------------------*/

/**\brief 'Break' up global list into sublists for parallel adaptation*/

void break_list(List_leaf Heads[], List_leaf Tails[])
{
  int i;
  
  for(i = 0; i < (number_of_threads - 1); i++)
    {
      Tails[i] -> next = NULL;
      Heads[i+1] -> prev = NULL;
    }
}

/*------------------------------------------------------------------*/

/**\brief Break up global vertex list into sublists*/

void break_vtx_list(List_vtx_glob Heads[], List_vtx_glob Tails[])
{
  int i;
  
  for(i = 0; i < (number_of_threads - 1); i++)
    { 
      Tails[i] -> next = NULL; 
      Heads[i+1] -> prev = NULL; 
    }
}

/*------------------------------------------------------------------*/

/**\brief Merge leaf cell sublists into global list again*/

void merge_lists(List_leaf Heads[], List_leaf Tails[])
{
  int i;

  for(i = 0; i < (number_of_threads - 1); i++)
    {
      Tails[i] -> next = Heads[i+1];
      Heads[i+1] -> prev = Tails[i];
    }  
}

/*------------------------------------------------------------------*/

/**\brief Merge vertex sublists into global list again*/

void merge_vtx_lists(List_vtx_glob Heads[], List_vtx_glob Tails[])
{
  int i;

  for(i = 0; i < (number_of_threads - 1); i++)
    {
      Tails[i] -> next = Heads[i+1];
      Heads[i+1] -> prev = Tails[i];
    } 
}

/*------------------------------------------------------------------*/

/**\brief Go through vertex sublists and free any verticies that are marked for deletion*/

void delete_verticies(List_vtx_glob *Vtxhead, List_vtx_glob *Vtxtail)
{
  List_vtx_glob node = *Vtxhead;
  Vertex Vtx;

  while(node != NULL)
    {
      Vtx = node -> vtx_loc;

      if(node -> delete == TRUE)
	{
	  node = node -> next;
	  delete_from_Vtxlist(Vtx, Vtxhead, Vtxtail);
	  free(Vtx -> vtx_nums);
	  free(Vtx);
	}
      else node = node -> next;
    }  
}

/*------------------------------------------------------------------*/

/**\brief Below are some useful functions for debugging when things really get hoary*/

void count_deletable_verticies(List_vtx_glob *Vtxhead, List_vtx_glob *Vtxtail, int thread_num, int step)
{
  int count, count2;
  double cen[3];

#if 0
  FILE *fout;
  sprintf(fname,"%dvtxdel%d",step,thread_num);
  fout = fopen(fname, "w");
#endif

  List_vtx_glob node = *Vtxhead;
  List_vtx_glob node2;
  Vertex Vtx, Vtx2;  

  count = 0;
  count2 = 0;
#if 0
  while(node != NULL)
    {
      Vtx = node -> vtx_loc;

      if(node -> delete == TRUE)
	{
	  fprintf(fout,"vtx at %e %e %e\n",Vtx->loc[0],Vtx->loc[1],Vtx->loc[2]);
	  count++;	 
	}
      count2++;
      
      node = node -> next;
    }
#endif  

#if 0
  fprintf(fout,"No. vtxs deleted thread %d is %d\n",thread_num,count);
  fprintf(fout,"Tot. no. vtxs in this thread %d is %d\n",thread_num,count2);
#endif

  node = *Vtxhead; /*Check for repetitions*/  
  while(node != NULL) {
    Vtx = node -> vtx_loc;
    
    count = 0;
    cen[0] = Vtx->loc[0]; cen[1]=Vtx->loc[1]; cen[2] = Vtx->loc[2];

    node2 = *Vtxhead;      
    while(node2 != NULL) {
      Vtx2 = node2 -> vtx_loc;

      if((Vtx2->loc[0] == cen[0]) && (Vtx2->loc[1] == cen[1]) && (Vtx2->loc[2] == cen[2])) {
	count++;
      }
      
      node2 = node2 -> next;
    }

    if(count > 1) 
      printf("ATTENTION - Vtx %e %e %e is repeated\n",cen[0],cen[1],cen[2]);

    node = node -> next;
  }

  
}

void check_duplicate_cells(List_leaf Head, List_leaf Tail, int thread_num, int step)
{
  int count;
  short int level;
  double cen[3];

#if 0
  FILE *fout;
  sprintf(fname,"%ddupcell%d",step,thread_num);
  fout = fopen(fname, "w");
#endif

  List_leaf node, node2;
  Cart_cell C, C2;

  node = Head; /*Check for repetitions*/  
  while(node != NULL) {
    C = node -> cell_loc;
    
    count = 0;
    cen[0] = C->centroid[0]; cen[1]=C->centroid[1]; cen[2] = C->centroid[2];
    level = C -> cell_level;

    node2 = Head;      
    while(node2 != NULL) {
      C2 = node2 -> cell_loc;

      if((C2->centroid[0] == cen[0]) && (C2->centroid[1] == cen[1]) && (C2->centroid[2] == cen[2]) && (C2 -> cell_level == level)) {
	count++;
      }
      
      node2 = node2 -> next;
    }

    if(count > 1) 
      printf("ATTENTION - Cell level %hd @ (%e %e %e) is repeated\n",C->cell_level,cen[0],cen[1],cen[2]);

    node = node -> next;
  }
}

void check_flux_connect(List_leaf Head, List_leaf Tail, int thread_num)
{
List_leaf node = Head;
Cart_cell C;
short int i;

while(node != NULL) {
C = node->cell_loc;
for(i = 0; i < 6; i++) {
if((C->face_neighbours[i][0]!= NULL) && (C->face_neighbours[i][1]==NULL)) {
          if(((C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]) > 0) &&
             (C->Flow_data.Face_fluxes[i][0] == NULL) &&
             ((C->Merge_data == NULL) || (C->face_neighbours[i][0]->Merge_data != C->Merge_data))) {
            printf("T %d,  Strangely, dir %hd, area %e here but no flux vector\n",thread_num,i,
                    (C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]));
          }
}
}
node = node -> next;
}
}

void printout_adapt_cells(short int flag, short int just_created, List_leaf Head, List_leaf Tail, int thread_num, int step)
{
  List_leaf node = Head;
  Cart_cell C;
  short int i,j,opp;
  FILE *fout;
  char fname[64];
  double loc[3];
  int count = 0;
  opp = 0;
  if(just_created == 1)
    sprintf(fname, "%dprintout_just_created%d",step,thread_num);
  else if(flag == REFINE)
    sprintf(fname, "%dprintout_refine_cells%d",step,thread_num);
  else if(flag == REALLY_COARSEN)
    sprintf(fname, "%dprintout_coarsen_cells%d",step,thread_num);
  
  
  fout = fopen(fname, "w");

  while(node != NULL) {
    C = node -> cell_loc;

    if(((just_created == FALSE) && (C -> adapt_flag == flag)) || ((just_created==TRUE) && (C->just_created==TRUE))) {
      fprintf(fout,"No. %d, this cell child %hd, level %hd @ (%9.8e %9.8e %9.8e) , type %hd, vol %e (just_created %hd)\n",count+1,C->child_num,
	      C->cell_level,C->centroid[0],C->centroid[1],C->centroid[2],C->cell_type,C->cell_volume,C->just_created);
      fprintf(fout,"This cell's parent flag %hd, L%hd @ (%9.8e %9.8e %9.8e), type %hd, vol %e\n",
	      C->parent->adapt_flag,C->parent->cell_level,C->parent->centroid[0],C->parent->centroid[1],C->parent->centroid[2],
	      C->parent->cell_type,C->parent->cell_volume);

      if(C->Merge_data != NULL)
	fprintf(fout,"Is merged cell\n");
      if(C->is_small_cell == TRUE)
	fprintf(fout,"Is small cell\n");
      if(C -> Leaf_list_loc == NULL)
	fprintf(fout, "Strangely, not on list\n");

      fprintf(fout,"Flux areas - \n");
      for(i = 0; i < 6; i++)
	fprintf(fout, "dir %hd, %e %e %e %e (tot %e)\n",i,C->flux_area[i][0],C->flux_area[i][1],C->flux_area[i][2],C->flux_area[i][3],
		(C->flux_area[i][0])+(C->flux_area[i][1])+(C->flux_area[i][2])+(C->flux_area[i][3]));
      fprintf(fout,"\n");

      fprintf(fout,"Neighbours are - \n");
      for(i = 0; i < 6; i++) {

	switch(i)
	  {
	  case 0:
	    opp = 1; break;
	  case 1:
	    opp = 0; break;
	  case 2:
	    opp = 3; break;
	  case 3:
	    opp = 2; break;
	  case 4:
	    opp = 5; break;
	  case 5:
	    opp = 4; break;
	  }

	if((C->face_neighbours[i][0] == NULL)) {

	  fprintf(fout,"Dir %hd, no neighb\n",i);
	  if(((C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]) > 0) && 
	     (C->Flow_data.Face_fluxes[i][0] == NULL) && (flag != REALLY_COARSEN))	  
	    fprintf(fout, "  Strangely, adapt_flag %hd have not face flux border dir %hd though area %e\n",
		    C->adapt_flag, i, C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]);


	  if(((C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]) == 0) && 
	     (C->Flow_data.Face_fluxes[i][0] != NULL) && (flag != REALLY_COARSEN))	  
	    fprintf(fout, "  Strangely, adapt_flag %hd have face flux border dir %hd though area %e\n",
		    C->adapt_flag,i,C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]);

	}
	else if(C -> face_neighbours[i][1] == NULL) {

	  fprintf(fout,"Dir %hd, 1 neighb (L%hd @ (%9.8e %9.8e %9.8e) type %hd, vol %g (small_cell %hd) flag %hd, just_created %hd\n ",i,
		  C->face_neighbours[i][0]->cell_level,C->face_neighbours[i][0]->centroid[0],C->face_neighbours[i][0]->centroid[1],
		  C->face_neighbours[i][0]->centroid[2],C->face_neighbours[i][0]->cell_type,
		  C->face_neighbours[i][0]->cell_volume,C->face_neighbours[i][0]->is_small_cell,
		  C->face_neighbours[i][0]->adapt_flag,C->face_neighbours[i][0]->just_created);

	  fprintf(fout, "Neighb flux areas (opp_dir %hd) - %e %e %e %e\n", opp,
		  C->face_neighbours[i][0]->flux_area[opp][0],C->face_neighbours[i][0]->flux_area[opp][1],
		  C->face_neighbours[i][0]->flux_area[opp][2],C->face_neighbours[i][0]->flux_area[opp][3]);

	  if(C->face_neighbours[i][0] -> Merge_data != NULL) {
	    fprintf(fout, "  Merged neighbour");
	    if(C->face_neighbours[i][0]->Merge_data == C->Merge_data) {
	      fprintf(fout, " and merged with this cell\n");
	      if(just_created == TRUE)
		if(C->Flow_data.Face_fluxes[i][0] != NULL)
		  fprintf(fout,"  Strangely, flag %hd have face flux here dir %hd\n",C->adapt_flag,i);
	    }
	    else fprintf(fout, "\n");
	  }

	  if(C->Flow_data.Face_fluxes[i][0] == NULL)
	    fprintf(fout, "  No flux vector shared here\n");
	  if((C->face_neighbours[i][0]->cell_type == SOLID) && 
	     ((C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]) > 0)) {
	    fprintf(fout, "  Strangely, solid neighb here but can flux thru\n");
	  }

	  if((C->face_neighbours[i][0] ->cell_type == SOLID) && (C->Flow_data.Face_fluxes[i][0] != NULL))
	     fprintf(fout, "  Strangely, solid neighb here but have flux vector\n");
	  if(C->face_neighbours[i][0]->cell_level > C->cell_level)
	    fprintf(fout,"  Strangely, 1 neighb here finer than cell\n");

	  if((flag!=REALLY_COARSEN) && ((C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]) > 0) &&
	     (C->Flow_data.Face_fluxes[i][0] == NULL) && 
	     ((C->Merge_data == NULL) || (C->face_neighbours[i][0]->Merge_data != C->Merge_data))) {
	    fprintf(fout, "  Strangely, dir %hd, area %e here but no flux vector\n",i,
		    (C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]));
	  }

	  if((flag!=REALLY_COARSEN) && ((C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]) == 0) &&
	     (C->Flow_data.Face_fluxes[i][0] != NULL) && 
	     ((C->Merge_data == NULL) || (C->face_neighbours[i][0]->Merge_data != C->Merge_data))) {
	    fprintf(fout, "  Strangely, dir %hd, area %e here but flux vector\n",i,
		    (C->flux_area[i][0]+C->flux_area[i][1]+C->flux_area[i][2]+C->flux_area[i][3]));
	  }


	}
	else if(C -> face_neighbours[i][1] != NULL) {
	  for(j = 0; j < 4; j++) {
	    fprintf(fout, "Dir %hd quad %hd, neighb (L%hd @ (%9.8e %9.8e %9.8e) type %hd, volume %g (small cell %hd), flag %hd, just_created %hd\n",i,j,
		    C->face_neighbours[i][j]->cell_level,C->face_neighbours[i][j]->centroid[0],
		    C->face_neighbours[i][j]->centroid[1],C->face_neighbours[i][j]->centroid[2],
		    C->face_neighbours[i][j]->cell_type,C->face_neighbours[i][j]->cell_volume,
		    C->face_neighbours[i][j]->is_small_cell,C->face_neighbours[i][j]->adapt_flag,
		    C->face_neighbours[i][j]->just_created);

	    fprintf(fout, "Neighb tot flux areas (opp_dir %hd, cell quad %hd) - %e\n",opp,j,
		    ((C->face_neighbours[i][j]->flux_area[opp][0])+(C->face_neighbours[i][j]->flux_area[opp][1])+
		     (C->face_neighbours[i][j]->flux_area[opp][2])+(C->face_neighbours[i][j]->flux_area[opp][3])));


	    if(C->face_neighbours[i][j] -> Merge_data != NULL) {
	      fprintf(fout, "  Merged neighbour");
	      if(C->face_neighbours[i][j]->Merge_data == C->Merge_data) {
		fprintf(fout, " and merged with this cell\n");
		if(just_created == TRUE)
		  if(C->Flow_data.Face_fluxes[i][0] != NULL)
		    fprintf(fout,"  Strangely, adapt_flag %hd have face flux here dir %hd quad %hd\n",C->adapt_flag,i,j);
	      }
	      else fprintf(fout, "\n");
	    }

	    if(C->Flow_data.Face_fluxes[i][j] == NULL)
	      fprintf(fout, "  No flux vector shared here\n");
	    if((C->face_neighbours[i][j]->cell_type == SOLID) && (C->flux_area[i][j] > 0)) {
	      fprintf(fout, "  Strangely, solid neighb here but can flux thru\n");
	    }

	    if((C->face_neighbours[i][j]->cell_type==SOLID) && (C->Flow_data.Face_fluxes[i][j] != NULL))
	      fprintf(fout, "  Strangely, solid neighb here but have flux vector\n");
	    if(C->face_neighbours[i][j]->cell_level < C->cell_level)
	      fprintf(fout,"  Strangely, finer neighb here actually not as fine as I thought\n");


	    if((flag!=REALLY_COARSEN) && (C->flux_area[i][j] > 0) && (C->Flow_data.Face_fluxes[i][j] == NULL) &&
	       ((C->Merge_data == NULL) || (C->Merge_data != C->face_neighbours[i][j]->Merge_data))) {
	      fprintf(fout,"  Strangely, dir %hd, quad %hd, area %e here but no flux vector\n",i,j,C->flux_area[i][j]);
	    }

	    if((flag!=REALLY_COARSEN) && (C->flux_area[i][j] == 0) && (C->Flow_data.Face_fluxes[i][j] != NULL) &&
	       ((C->Merge_data == NULL) || (C->Merge_data != C->face_neighbours[i][j]->Merge_data))) {
	      fprintf(fout,"  Strangely, dir %hd, quad %hd, area %e here but flux vector\n",i,j,C->flux_area[i][j]);
	    }

	  }
	}
      }
      
      /*Check if verticies should be where they should be*/
      for(i = 0; i < 8; i++) {
#if 1
	get_vtx_loc(C, i, loc);
	if((C->verticies[i]->loc[0] != loc[0]) || (C->verticies[i]->loc[1] != loc[1]) || (C->verticies[i]->loc[2] != loc[2])) {
	  fprintf(fout, "  Strangely, vtx %hd should be %e %e %e\n",i,loc[0],loc[1],loc[2]);
	}
#endif
	if(C -> verticies[i] -> Vtxlist_loc == NULL)
	  fprintf(fout, "  Strangely, vtx doesn't exist on vtxlist\n");
	if(C->verticies[i]->Vtxlist_loc->vtx_loc == NULL)
	  fprintf(fout, "  Strangely, vtx->vtxlist_loc->vtx_loc is NULL\n");
      }

      count++;
      fprintf(fout,"\n\n\n");
    }

    if(node == Tail)
      break;
    else node = node -> next;
  }

  fclose(fout);
}

void reset_adapt_flag(List_leaf L, short int flag)
{
  List_leaf node = L;
  Cart_cell C;
  while(node != NULL) {
    C = node -> cell_loc;

    if(flag == REFINE)
      C -> parent -> adapt_flag = UNKNOWN;
    else C -> adapt_flag = UNKNOWN;

    node = node -> next;
  }
}

void printout_cells_2b_refined(List_leaf L, int step)
{
  List_leaf node = L;
  FILE *fout = fopen("refine_cells","a");
  int count = 0;

  fprintf(fout,"Cells 2b refined for step %d\n",step);
  while(node != NULL) {
    if(node -> cell_loc -> adapt_flag == REFINE) {
      fprintf(fout, "Want to refine cell %hd (%e %e %e)\n",node->cell_loc->cell_level,node->cell_loc->centroid[0],
	      node->cell_loc->centroid[1],node->cell_loc->centroid[2]);
      count++;
    }

    node = node -> next;
  }
  fprintf(fout, "Total - %d cells\n",count);
  fprintf(fout,"\n\n");

  fclose(fout);
}

void printout_cells_2b_coarsened(List_leaf L, int step)
{
  List_leaf node = L;
  Cart_cell Parent;
  FILE *fout = fopen("coarsen_cells","a");
  int i,count = 0;

#if 0
  while(node != NULL) {
    if((node->cell_loc->parent->cell_level==5) && (node->cell_loc->parent->centroid[0]==0.015625) && 
       (node->cell_loc->parent->centroid[1]==-0.046875) && (node->cell_loc->parent->centroid[2]==0.046875)) {
      fprintf(fout,"FOUND THIS CELL\n");
      
      for(i = 0; i < 8; i++) {
	if(node->cell_loc->parent->children[i]->cell_type != SOLID)
	  fprintf(fout,"Child %hd's flags %hd (type %hd, jc %hd)\n",i,node->cell_loc->parent->children[i]->adapt_flag,
		  node->cell_loc->parent->children[i]->cell_type,node->cell_loc->parent->children[i]->just_created);
      }
      break;      
    }
    node = node -> next;
  }
#endif

  node = L;

  fprintf(fout,"Cells to be coarsened for step %d\n",step);
  while(node != NULL) {

    if(node -> cell_loc -> adapt_flag == REALLY_COARSEN) {
      Parent = node -> cell_loc -> parent;
      fprintf(fout, "Want to coarsen cell %hd (%e %e %e)\n",node->cell_loc->parent->cell_level,
	      node->cell_loc->parent->centroid[0],node->cell_loc->parent->centroid[1],node->cell_loc->parent->centroid[2]);

#if 1
      fprintf(fout, "Children params - \n");
      for(i = 0; i < 8; i++) {
	if(Parent->children[i]->cell_type != SOLID)
	  fprintf(fout, "Child %hd, type %hd, flag %hd\n",i,Parent->children[i]->cell_type,Parent->children[i]->adapt_flag);		  
      }
#endif

      count++;

      while((node != NULL) && (node->cell_loc->parent == Parent))
	node = node -> next;      
    }
    else node = node -> next;
  }
  fprintf(fout, "Total - %d cells\n",count);
  fprintf(fout,"\n\n");

  fclose(fout);
}

