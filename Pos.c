#include "Pos.h"

//Produce a warning if debug is enabled since the print statements are now enabled
#ifdef _DEBUG_
#warning debug enabled
#endif

/*
 * srand (time(NULL)); init the randdom number gen
 * int a = rand(); 
 * Max rand number == RAND_MAX
 * TODO Remove the findMinMaxOAtom function since it is not being called by anyone
 */

int FindAtoms(FILE *file_p)
{

   //var
   int output = 0;
   char string[BUFFER_SIZE] = {0,};
   char *string_p = &string[0];
   
   //find the number of atoms in the file
   while(!output)
   {
      fgets(string_p, BUFFER_SIZE, file_p);
      output = atoi(string_p);
   }
   
   fgets(string_p, BUFFER_SIZE, file_p);
   max_bonds = atoi(string_p);
   
   fgets(string_p, BUFFER_SIZE, file_p);
   max_angles = atoi(string_p);
   
   //return the number of atoms in the file
   return output;

}

void ReadTill(FILE *file_p,char *str2, int compare_number)
{

   //var
   char string1[BUFFER_SIZE] = {0,};
   const char *STR2 = str2;
   const char *STR1 = &string1[0];

   //read the file until the word str2 is found
   while(strncmp(STR1,STR2,compare_number) != 0)
   {
      fgets(string1, BUFFER_SIZE, file_p);
   }
 
   //Read the blank line after the word str2
   fgets(string1, BUFFER_SIZE, file_p);

}


/**
* Reads the infomation about the atoms from the dump file.
*/
void ReadAtoms(FILE *file_p)
{

	//variables
	int id,mol,type,ix,iy,iz;
	double q,x,y,z;
	unsigned int index = 0;
	char string[BUFFER_SIZE];
   Atom *p;

   //Loop though the atom list reading in the data
	for (index = 0; index < numAtoms; index++)
	{
      //Reset the pointer for the array
      p = list_atom_g;

		//read the formated lammps input file
		fscanf (file_p, "%i %i %i %lf %lf %lf %lf %i %i %i",\
			&id, &mol, &type, &q, &x, &y, &z, &ix, &iy, &iz);

      //Set the id in terms of the index of the array
      id--;
	   
      //Go to the indes in the array
      p += id;

      //save the important parts of the dump file
      p->id = id + 1;
      p->x = x;
      p->y = y;
      p->z = z;
      p->type = type;
      p->mol = mol;
      p->q = q;
      p->ix = ix;
      p-> iy = iy;
      p-> iz = iz;

		//Get rid of the LF char
      fgets(string, BUFFER_SIZE, file_p);

	}

   //clear out the pointer
   p = NULL;

}

void ReadBonds(FILE *file_p)
{

   //var
   char str[BUFFER_SIZE] = {0,};
   int atom1,atom2,id,type,index;
   Bond *list = list_bond_g;

   //read the bond list
   for (index = 0; index < max_bonds; index++, list++)
   {

      //read the file
      fscanf(file_p, "%i %i %i %i", &id, &type, &atom1, &atom2);
 
      //save the information from the bond list
      list->type = type;
      list->id = id;
      list->atom1 = atom1;
      list-> atom2 = atom2;

      //read the LF from the line
      fgets(str, BUFFER_SIZE, file_p);

   }
   
   list = NULL;

}

void ReadAngles(FILE *file_p)
{

   //var
   int type, atom1, atom2, atom3, id, index;
   Angle *list = list_angle_g;
   char str[BUFFER_SIZE] = {0,};

   //Read the angle list
   for(index = 0; index < max_angles; index++, list++)
   {

       //read the file
      fscanf (file_p, "%i %i %i %i %i", &id, &type, &atom1, &atom2, &atom3);
 
      //save the information from the bond list
      list->type = type;
      list->id = id;
      list->atom1 = atom1;
      list->atom2 = atom2;
      list->atom3 = atom3;

      //read the LF from the line
      fgets(str, BUFFER_SIZE, file_p);

   }

}

void populateSurface( void )
{//start of function

   //var
   Atom *o_atom_surface, *HatomPair;
   Atom *list = list_atom_g;
   Bond *bond_list = list_bond_g;
   Angle *angle_list = list_angle_g;
   int *pair_list = list_pair_g;
   int index_list[MAX_Y_VALUES] = {0,};
   double list_y_values[MAX_Y_VALUES] = {0.0,};
   int count = 0;
   double x,y;
   const double RANDOM_MAX_VALUE_D = (double) RAND_MAX;
   int random_num, index;
   double percentage_populated = 0.5;
   int number_of_atoms_to_populate = top_O_atoms * percentage_populated;

   //Change the pointer pos to the first open slot in the array.
   //numAtoms is the number of atoms and not a indexing array value
   list += numAtoms;
   bond_list += max_bonds;
   angle_list += max_angles;

   //increase the number of angle types and bond types by 1
   max_angle_type++;
   max_atom_types++;
   max_bond_type++;

   //find where the O atoms of the surface are in the atom list
   //Not needed since this is now called in the main function
   //FindStartOAtom();

   //Map out the Y vales for the O atoms of the surface
   mapOutYValueOfAtomList(index_list, &count);

   //Set up the list_y_values array
   setUpYValueList(index_list, count, list_y_values);

   //Find the starting index for the H atoms of the surface
   FindHAtomSurface();

   //The below lines will be in a for loop depending on the percentage that the surface
   //is to be covered

   //Do random number generator to find the x and y coord of the O atom to populate
   srand(time(NULL));
   
      index = 0;
//   for (index = 0; index < number_of_atoms_to_populate; index++)
//   {//Start of for loop
   
      //Find a random Y value
      random_num = rand();
      y = (double) random_num / RANDOM_MAX_VALUE_D;
      y *= MAX_Y;
   
      //Find a random X value
      random_num = rand();
      x = (double) random_num / RANDOM_MAX_VALUE_D;
      x *= MAX_X;

#ifdef _DEBUG_
      printf("Random Values X: %lf and y: %lf\n",x,y);
#endif

      //Find the O atom of the surface to populate
      o_atom_surface = FindOAtom(x, y, index_list, count, list_y_values);


      //check to see if this O atom has already RCVD another H atom.
      //If the O atom already has a H atom attached to it,
      //The program will try again from start.
      if(CheckOAtomBond(o_atom_surface))
      {
         index--;
         //continue;
      }
      else
      {
         *pair_list = o_atom_surface->id;
         pair_list++;
      }

#ifdef _DEBUG_
      //For debugging only
      printf("x: %lf \t y: %lf type: %i\n",o_atom_surface->x,o_atom_surface->y,\
            o_atom_surface->type);
#endif
      
      //Find the H atom that is already bonded with the O atom of the surface
      HatomPair = FindHAtomPair(o_atom_surface);

#ifdef _DEBUG_
      //For debugging only
      puts("The atom id for the H atom that is ponded with The O atom is");
      printf("id = %i that is bonded to id = %i\n", HatomPair->id,\
      o_atom_surface->id);
#endif
      
      //attach new H atom to O atom of the surface
      
      //Set the new H atom to have the same mol id number as the O/H atom
      list->mol = o_atom_surface->mol;

      //increase the number of atoms by one since we are adding one H atom
      numAtoms++;

      //Set the id number to the next id number in the list;
      list->id = numAtoms;
      
      /*change the location of the old H atom to make the required angle for
       *the angle type to work. 109 degs is required or somewhere close to this
       *This value
       */
      HatomPair->x = o_atom_surface->x + 0.25;
      HatomPair->y = o_atom_surface->y + 0.125;
      HatomPair->z = o_atom_surface->z + 0.915481752958518;

      //Give the new H atom x y and z coords based on the O atom
      list->x = o_atom_surface->x + 0.0625;
      list->y = o_atom_surface->y + 0.945788152138317;
      list->z = o_atom_surface->z - 0.133455465511114;

      //give the new H atom a type
      list->type = H_ATOM_SURFACE_POS;

      //Add to the bond list
      //increase the number of bonds by one
      max_bonds++;

      //Set the new bond id as the total number of bonds
      bond_list->id = max_bonds;

      //set the type of the bond 
      bond_list->type = NEW_BOND_TYPE_H_ATOM_SURFACE;

      //set the two atoms that are bonded together
      bond_list->atom1 = o_atom_surface->id;
      bond_list->atom2 = list->id;

      //Add one to the total number of angles
      max_angles++;

      //Add to the angle list
      angle_list->id = max_angles;

      //Add to the angle list
      //Updated the angle type to the new one for the surface
      angle_list->type = NEW_ANGLE_TYPE_SURFACE;

      //Added the information for what atoms that make up the angle
      angle_list->atom1 = HatomPair->id;
      angle_list->atom2 = o_atom_surface->id;
      angle_list->atom3 = list->id;

//   } //End of for loop
} //End of function

void PrintFile(char *Filename)
{
   //vars



}


int CheckOAtomBond(Atom *o_atom)
{

   //var
   int* list = list_pair_g;

   while(*list)
   {
      if(*list == o_atom->id)
         return TRUE;
      //go to the next slot in the array
      list++;
   }

   return FALSE;

}

/**
*  Find the H atom that is bonded with the o_atom_surface varaible in the list.
*  @param Atom* the O atom that the H atom is bonded to.
*  @return Atom* of H atom that is bonded to o_atom_surface
*/
Atom* FindHAtomPair(Atom* o_atom_surface)
{

   //vars
   Atom* output = list_atom_g;
   output += h_surface_index;
   while(output->mol != o_atom_surface->mol)
   {
      output++;
   }
   
   return output;

}

void setUpYValueList(int *index_list, int count, double *list_y_values)
{
   
   //vars
   int index, index2;
   Atom* atom_list = list_atom_g;
  
   for (index = 0; index < count; index++)
   {
      //Find the next index of the atom list 
      index2 = start_id_o_atom + index_list[index];

      //save the y coord value to the list to make searching easier 
      list_y_values[index] = atom_list[index2].y;
   }
   
}

void FindHAtomSurface( void )
{
   //var 
   int h_atom_index = 0;
   Atom *list = list_atom_g;
   
   //start at the end of the O atoms to find th H atoms
   //The next atom should be the H atom
   list += start_id_o_atom + top_O_atoms;
   
   //This is the starting index for H surface atom
   h_atom_index = start_id_o_atom + top_O_atoms;

   //Find the H atom index for the surface
   while (list->type != H_ATOM_SURFACE)
   {
      list++;
      h_atom_index++;
   }

   //set the goboal varaible to the index value for the H atom of the surface
   h_surface_index = h_atom_index;
}

void FindStartOAtom( void )
{ 

   //var
   int index;
   Atom *list = list_atom_g;
   const double Q = -0.95; //charge of O atom we are looking for
   top_O_atoms = 0;

   //look for the O atoms of the surface
   for (index = 0; index < numAtoms; index++, list++)
   {
      //check to see if this is the top O atom
      if(list->q == Q)
      {
         top_O_atoms++;
         start_id_o_atom = index;
         list->type = O_ATOM_SURFACE;
      }
   }
  
   //Find the starting point of the O atoms of the surface based on the 
   //total number of O atoms
   start_id_o_atom -= top_O_atoms - 1;

}

void roateAlongZ(double theta, Atom* atom)
{
}

Atom* FindOAtom(double x, double y, int* index_list, int count, double* y_list)
{
   
   //var
   Atom *output = list_atom_g;
   int yindex;
   
#ifdef _DEBUG_   
   //For debugging only
   puts("From FindOAtom");
   printf("The starting of the O atoms of the surface: %i\n",start_id_o_atom);
   printf("x: %lf \t y: %lf type: %i\n",output->x,output->y,\
         output->type);
   puts("End of FindOAtom");
#endif

   //Find the Y coord of the O atom range
   yindex = findYIndexStart(&y,y_list,&count, index_list);
   
   //Find the X coord of the O atom

   //output += findXIndex(yindex, x);
   yindex= findXIndex(yindex,x);
   //printf("the x value is found to be ");
   output += yindex;
   
   return output;
}

int findYIndexStart(double *y, double *y_list, int *count, int *index_list)
{
   //vars
   double distance_max_sq ,distance_min_sq;
   int index, max = 0;

   //find the location in the array where y is less then the vales stored in the array
   for(index = 0; index < *count; index++)
   {
      if(*y <= y_list[index])
      {
         max = index;
         break;
      }
   }

   //find the distance between the max value and the next lower value
   //If max is 0 then there is nothing that needs to be compared and it is returned
   if(max > 0)
   {
      distance_max_sq = *y - y_list[max];
      distance_max_sq *= distance_max_sq;

      max--;
      distance_min_sq = *y - y_list[max];
      distance_min_sq *= distance_min_sq;
   }
   else
      return index_list[max];
   
   //See which index the y value is closer to and pick that value
   if(distance_min_sq <= distance_max_sq)
      return index_list[max];
   else
      return index_list[max + 1];

}

int findXIndex(int y_index, double x)
{
   //vars
   double distance_min_sq, distance_max_sq;
   int index = 0;
   Atom* atom_list = list_atom_g;
   index = y_index + start_id_o_atom;
   double y = atom_list[index].y;
   char broke = 0;
   int index_minus_one;

   //Loop though the atom list looking for the max x value index
   while(atom_list[index].y == y)
   {

      if(x <= atom_list[index].x)
      {
         broke = 1;
         index_minus_one = index - 1;
         break;  
      }
      //increase the vars for the next itteration of the loop
      index++;
   }

   //check to see if the while loop was broken out of
   //If it was not, then return the largest x value for the system
   if(broke)
   {
      
      //check to see if the x value is not the smallest
      if(atom_list[index_minus_one].y != y)
         return index;
      
      //Find the distnace of the x values
      distance_max_sq = x - atom_list[index].x;
      distance_max_sq *= distance_max_sq;

      distance_min_sq = x - atom_list[index_minus_one].x;
      distance_min_sq *= distance_min_sq;

      //return the smaller of the two distances as the chosen x value
      if(distance_min_sq <= distance_max_sq)
         return index_minus_one;
      else
         return index;
   }
   else 
      return index - 1;

}

void mapOutYValueOfAtomList(int *index_list, int *count)
{
   
   //vars
   Atom* atom_list = list_atom_g;
   *count = 0; //To ensure that the starting pos is always the same
   atom_list += start_id_o_atom;
   double old_y = -1;
   int index = 0;

   //loop though the atom list looking for all posible 
   while(atom_list->type == O_ATOM_SURFACE)
   {
      //If a new y value is found, mark that location and continue on
      if(old_y != atom_list->y)
      {
         //mark the index value with an offset
         index_list[*count] = index;

         //increase the number of values in the list
         *count += 1;
         
         //change the old_y value to the new y value
         old_y = atom_list->y;
      }

      //move to the next atom in the list
      atom_list++;
      index++;
   }


}

//Old function and will be removed at a later time
void findMinMaxOAtom(int *input, double y)
{
   //var
   Atom *list = list_atom_g;
   int index;
   list += start_id_o_atom;

   for (index = 0; index < top_O_atoms; index++, list++)
   {
      if (y <= list->y)
      {
         input[1] = index;
         break;
      }
   }

   input[0] = -1;

   while(list->type == O_ATOM_SURFACE)
   {
      list--;
      input--;
      if(list->y > y)
      {
         input[0] = index;
         break;
      }
   }

   if(input[0] == -1)
   {
      list = list_atom_g;
      list += input[1];
      list--;
      if(list->type == O_ATOM_SURFACE)
         input[0] = input[1] - 1; 
   }
   
}

/**
* Main function for the H bonding file.
*/
int main( void )
{

   //The input file name
   char *filename = "/home/I1nk/NewSurface/List.dat";
   
   //Open the file to be a read only
   FILE *file_p = fopen(filename,"r");

   //Check to see if the file could be opened as a read
   if(file_p == NULL)
   {
      puts("File could not be opened.");
      return ERROR_OPENING_FILE;
   }
   
   //Begin the work
   //Set up the arrays for the atom bond and angle list
   numAtoms = FindAtoms(file_p);
   Atom *list_atom = (Atom*) calloc(numAtoms * 2, sizeof(Atom));
   Bond *list_bond = (Bond*) calloc(numAtoms * 2, sizeof(Bond));
   Angle *list_angle = (Angle*) calloc(numAtoms * 2, sizeof(Angle));
   int *list_pair = (int*) calloc(numAtoms * 2, sizeof(int));

   //Set the global pointers up
   list_atom_g = list_atom;
   list_bond_g = list_bond;
   list_angle_g = list_angle;
   list_pair_g = list_pair;
 
   //Read the input file until the atom list is found
   ReadTill(file_p,"Atoms",5);

   //Read the atoms data and save it to memory
   ReadAtoms(file_p);
   
   //Read the file until Bonds are found
   ReadTill(file_p,"Bonds",5);

   //Read the Bonds data and save it to memory
   ReadBonds(file_p);

   //Read the file until Angles are found
   ReadTill(file_p,"Angles",6);

   //Read the angles data and save it to memory
   ReadAngles(file_p);

   //find where the O atoms of the surface are in the atom list
   FindStartOAtom();

   //Populate the surface with an extra H atom to make the surface positve
   populateSurface();
   
   //Clean the memory
   list_angle_g = NULL;
   list_atom_g = NULL;
   list_bond_g = NULL;
   list_pair_g = NULL;

   //Free the local memory that was set up using calloc
   free(list_bond);
   free(list_angle);
   free(list_atom);
   free(list_pair);


   return 0;
}
