//include protection
#ifndef HBONDING_H
#define HBONDING_H 

//includes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

//struct for holding atom information and other kinds of information
typedef struct ATOMS
{
   double x;
   double y;
   double z;
   unsigned char type;
   unsigned int mol;
   double q;
   int ix;
   int iy;
   int iz;
   int id;

}Atom;

typedef struct BONDS
{

   unsigned int id;
   unsigned int type;
   unsigned int atom1;
   unsigned int atom2;

}Bond;

typedef struct ANGLES
{

   unsigned int id;
   unsigned int type;
   unsigned int atom1;
   unsigned int atom2;
   unsigned int atom3;

}Angle;


//Function prototypes
//Void returns
void ReadAtoms(FILE*);
void ReadTill(FILE*, char*, int);
void ReadBonds(FILE*);
void ReadAngles(FILE*);
void FindStartOAtom(void);
void populateSurface(void);
void mapOutYValueOfAtomList(int*, int*);
void setUpYValueList(int*, int, double*);
void FindHAtomSurface(void);
void PrintFile(char*);

//Non void return functions
int findYIndexStart(double*, double*, int*, int*);
int findXIndex(int, double);
int FindAtoms(FILE*);
int CheckOAtomBond(Atom*);
Atom* FindOAtom(double, double, int*, int, double*);
Atom* FindHAtomPair(Atom*);

//Pound defines
#define TRUE 1
#define FALSE 0
#define ERROR_OPENING_FILE -1
#define O_ATOM_SURFACE 4
#define H_ATOM_SURFACE 3
#define AL_ATOM_SUBSTRIGT 1
#define O_ATOM_SUBSTRIGHT 2
#define H_ATOM_SURFACE_POS 5
#define MAX_Y 49.4362
#define MAX_X 199.794
#define BUFFER_SIZE 100
#define MAX_Y_VALUES 100
#define NEW_BOND_TYPE_H_ATOM_SURFACE 2
#define NEW_ANGLE_TYPE_SURFACE 2
//#define _DEBUG_

//Globals
Atom *list_atom_g;
Bond *list_bond_g;
Angle *list_angle_g;
int *list_pair_g;
int numAtoms = 0;
int max_bonds = 0;
int max_angles = 0;
int max_atom_types = 5;
int max_bond_type = 0;
int max_angle_type = 0;
int top_O_atoms = 0;
int start_id_o_atom = 0; //NOTE this is the index value and not the actual id number
double percentage = 0.25;
double _xmin = 0;
double _xmax = 0;
double _ymin = 0;
double _ymax = 0;
double _zmin = 0;
double _zmax = 0;
int h_surface_index = 0;


#endif
