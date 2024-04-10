#define FTT_CELLS     4
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_amd_printf : enable

#define NULL (global void*)0
#define NULLPrvt (void *)0

#define FTT_2D
//#define DEBUG_MODE

typedef enum
{
  FTT_X = 0,
  FTT_Y,
#ifndef FTT_2D
  FTT_Z,
#endif /* FTT_3D */
  FTT_DIMENSION,
  FTT_XY,
#ifdef FTT_2D
  FTT_XYZ = FTT_XY
#else  /* FTT_3D */
  FTT_XYZ
#endif /* FTT_3D */
} FttComponent;

typedef enum
{
  FTT_RIGHT = 0,
  FTT_LEFT,
  FTT_TOP,
  FTT_BOTTOM,
#if !defined(FTT_2D)
  FTT_FRONT,
  FTT_BACK,
#endif 
  FTT_NEIGHBORS
} FttDirection;

#define FTT_NEIGHBORS_2D (FTT_BOTTOM + 1)

typedef enum {
  FTT_FLAG_ID        = 7,
  FTT_FLAG_DESTROYED = 1 << 3,
  FTT_FLAG_LEAF      = 1 << 4,        /* used only for I/O operations */
  FTT_FLAG_TRAVERSED = FTT_FLAG_LEAF, /* used for face traversal */
  FTT_FLAG_USER      =      5         /* user flags start here */
} FttCellFlags;

struct _FttVector {
  double x, y, z;
};

struct _FttCell {
  unsigned int  flags;
  __global void * data;
  __global struct _FttOct * parent, * children;
};

struct _FttCellNeighbors {
  __global  struct _FttCell * c[FTT_NEIGHBORS];
};

struct _FttCellChildren {
  __global struct _FttCell * c[FTT_CELLS];
};

struct _FttRootCell {
  struct _FttCell cell;

  struct _FttCellNeighbors neighbors;
  struct _FttVector pos;
  unsigned int level;
  __global void * parent;
};

struct _FttOct {
  unsigned int level;
  __global struct _FttCell * parent;
  struct _FttCellNeighbors neighbors;
  struct _FttVector pos;

  struct _FttCell cell[FTT_CELLS];
};

struct _FttCellFace {
  __global struct _FttCell * cell, * neighbor;
  FttDirection d;
};

struct _GfsFaceStateVector {
  double un;
  double v;
};

struct _GfsStateVector {
  /* temporary face variables */
  struct _GfsFaceStateVector f[FTT_NEIGHBORS];

  /* solid boundaries */
  __global struct _GfsSolidVector * solid;

  double place_holder;
};

struct _GfsSolidVector {
  double s[FTT_NEIGHBORS];
  double a, fv;
  __global struct _FttCell * merged;
  struct _FttVector cm, ca, v;
};

struct _GfsGradient {
  double a, b;
};

struct Data_Of_Cell
{
  struct _GfsStateVector a;
  double numbers[20];
};

typedef enum {
  FTT_BOUNDARY,
  FTT_FINE_FINE,
  FTT_FINE_COARSE
} FttFaceType;

#define              ftt_cell_level(c)  ((c)->parent ?\
                                         (c)->parent->level + 1 :\
                                         ((global struct _FttRootCell *) c)->level)

#define  FTT_CELL_IS_LEAF(c)      ((c)->children == NULL)
#define  FTT_CELL_IS_DESTROYED(c) (((c)->flags & FTT_FLAG_DESTROYED) != 0)
#define  FTT_CELL_ID(c)           ((c)->flags & FTT_FLAG_ID)
#define  FTT_CELL_IS_ROOT(c)      ((c)->parent == NULL)

#define GFS_STATE(cell)               ((global struct _GfsStateVector *) (cell)->data)

#define GFS_VALUE(cell,v)            ((&GFS_STATE (cell)->place_holder)[v])

#define GFS_VALUEI(cell, index)     ((&GFS_STATE (cell)->place_holder)[index])

#define GFS_IS_MIXED(cell)      ((cell) != NULL &&\
                                 GFS_STATE (cell)->solid != NULL)

#define GFS_FACE_FRACTION(fa) (GFS_IS_MIXED ((fa)->cell) ?\
                               GFS_STATE ((fa)->cell)->solid->s[(fa)->d] : 1.)

#define GFS_NODATA  1.79769e+308

#define return_if_fail(expr)\
        if(expr != 1) {\
           printf("return_if_fail: %s false\n", expr);\
           return;\
       }

#define return_val_if_fail(expr, val)\
         if(!expr) {\
           printf("return_val_if_fail: %s false\n", expr);\
           return val;\
         }

#define assert(expression)\   
      if(!expression){\
            printf("assert: %s false\n", expression);\ 
            return;\
        }


         

__constant int ftt_opposite_direction[FTT_NEIGHBORS] = {1, 0, 3, 2};//mallon dn xreiazetai pali dhlwsh

#define FTT_OPPOSITE_DIRECTION(d)     (ftt_opposite_direction[d])

__constant double coords[FTT_CELLS * 3] = {-1.0, 1.0, 0.0, 1.0, 1.0, 0.0, -1.0, -1.0, 0.0, 1.0, -1.0, 0.0};

__constant unsigned int index[FTT_NEIGHBORS_2D * FTT_CELLS/2] = {1, 3, 0, 2, 0, 1, 2, 3};

#ifdef FTT_2D
__constant int perpendicular[FTT_NEIGHBORS][FTT_CELLS] = 
  {{-1,  2, -1,  3},
   { 2, -1,  3, -1},
   { 1,  0, -1, -1},
   {-1, -1,  1,  0}};

__constant int neighbor_index[FTT_NEIGHBORS][FTT_CELLS]
    = {{1,-1,3,-3},
       {-2,0,-4,2},
       {-3,-4,0,1},
       {2,3,-1,-2}};
#endif
