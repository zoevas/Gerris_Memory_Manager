#include "/srv/homes/zovasili/Gerris/gerrisopencl/bin/kernel.h"


/*static inline void ftt_cell_neighbors (const __global struct _FttCell  * cell,
			  struct _FttCellNeighbors * neighbors)
{
  //g_return_if_fail (cell != NULL);
  //g_return_if_fail (neighbors != NULL);

  if (!FTT_CELL_IS_LEAF (cell) && neighbors != &cell->children->neighbors) {
    memcpy (neighbors, &cell->children->neighbors, sizeof (struct _FttCellNeighbors));
    return;
  }
};*/


static inline void ftt_cell_neighbors (const global struct _FttCell  * cell,
			  struct _FttCellNeighbors * neighbors)
{ 
  int i;

 return_if_fail(cell != NULL);

 return_if_fail(neighbors != NULLPrvt);

  if (!FTT_CELL_IS_LEAF (cell)) 
    for(i = 0; i < FTT_NEIGHBORS; i++) {
      neighbors->c[i] =  cell->children->neighbors.c[i];
    }
  return;
};

static inline void ftt_cell_children (const global struct _FttCell * cell,
			 struct _FttCellChildren * children)
{
  global struct _FttOct * oct;
  unsigned int i;

  return_if_fail(cell != NULL); 

  return_if_fail(!FTT_CELL_IS_LEAF (cell)); 
 
  return_if_fail(children != NULLPrvt); 

  oct = cell->children;
  for (i = 0; i < FTT_CELLS; i++) {
    children->c[i] = FTT_CELL_IS_DESTROYED (&(oct->cell[i])) ? 
      NULL : &(oct->cell[i]);
  }
};

/**
 * ftt_face_type:
 * @face: a #FttCellFace.
 *
 * Returns: the type of @face.
 */
static inline
FttFaceType ftt_face_type (const struct _FttCellFace * face)
{
  return_val_if_fail (face != NULLPrvt, 0);

  if (face->neighbor == NULL)
    return FTT_BOUNDARY;
  if (ftt_cell_level (face->cell) > ftt_cell_level (face->neighbor))
    return FTT_FINE_COARSE;
  assert (ftt_cell_level (face->cell) == ftt_cell_level (face->neighbor));
  return FTT_FINE_FINE;
}


struct _Gradient {
  double a, b, c;
};


/**
 * ftt_cell_children_direction:
 * @cell: a #FttCell.
 * @d: a direction.
 * @children: a #FttCellChildren.
 *
 * Fills @children with the children (2 in 2D, 4 in 3D)
 * of @cell in direction @d.
 *  200

 * This function fails if @cell is a leaf.
 *
 * Returns: the number of children in direction @d.
 */
static inline unsigned int ftt_cell_children_direction (const __global struct _FttCell * cell,
				   FttDirection d,
				    struct _FttCellChildren * children)
{
  __global struct _FttOct * oct;
  unsigned int i;
 
  return_val_if_fail (cell != NULL, 0);
  return_val_if_fail (!FTT_CELL_IS_LEAF (cell), 0);
  return_val_if_fail (d < FTT_NEIGHBORS, 0);
  return_val_if_fail (children != NULLPrvt, 0);

  oct = cell->children;

  for (i = 0; i < FTT_CELLS/2; i++)
    children->c[i] = FTT_CELL_IS_DESTROYED (&(oct->cell[index[d * FTT_CELLS/2 + i]])) ? 
      NULL : &(oct->cell[index[d * FTT_CELLS/2 + i]]);
  return FTT_CELLS/2;
}


/* face_weighted_gradient_stencil() needs to be updated whenever
   this function is modified                                   */
static void face_weighted_gradient (const struct _FttCellFace * face,
				    struct _Gradient * g,
				    unsigned int v,
				    int max_level,
				    unsigned int dimension)
{
  unsigned int level;

  return_if_fail(face != NULLPrvt); 

  g->a = g->b = 0.;
  if (face->neighbor == NULL)
    return;

  level = ftt_cell_level (face->cell);
  if (ftt_cell_level (face->neighbor) < level) {
    /* neighbor is at a shallower level */
    struct _Gradient  gcf;
    double w = GFS_STATE (face->cell)->f[face->d].v;

   // gcf = gradient_fine_coarse (face, v);
    g->a = w*gcf.a;
    g->b = w*(gcf.b*GFS_VALUEI (face->neighbor, v) + gcf.c);
  }
  else {
    if (level == max_level || FTT_CELL_IS_LEAF (face->neighbor)) {
      /* neighbor is at the same level */
      double w = GFS_STATE (face->cell)->f[face->d].v;

      g->a = w;
      g->b = w*GFS_VALUEI (face->neighbor, v);
    }
    else {
      /* neighbor is at a deeper level */
      struct _FttCellChildren  children;
      struct _FttCellFace f;
      unsigned int i, n;
      
      f.d = ftt_opposite_direction[face->d];
      n = ftt_cell_children_direction (face->neighbor, f.d, &children);
      f.neighbor = face->cell;
      for (i = 0; i < n; i++) 
	if ((f.cell = children.c[i])) {
	  struct _Gradient  gcf;
	  double w = GFS_STATE (f.cell)->f[f.d].v;
	
	  //gcf = gradient_fine_coarse (&f, v);
	  g->a += w*gcf.b;
	  g->b += w*(gcf.a*GFS_VALUEI (f.cell, v) - gcf.c);
	}
      if (dimension > 2) {
	g->a /= n/2.;
	g->b /= n/2.;
      }
    }
  }
}


/**
 * ftt_cell_neighbor_not_cached:
 * @cell: a #FttCell.
 * @d: a direction.
 *
 * Returns: the neighbor of @cell in direction @d or %NULL if @cell
 * has no neighbor in this direction (does not use saved values even
 * if available).  
 */
static inline
__global struct _FttCell * ftt_cell_neighbor_not_cached (const __global struct _FttCell * cell,
					FttDirection d)
{

//#else  /* FTT_3D */
//    = {{1,-1,3,-3,5,-5,7,-7},
//       {-2,0,-4,2,-6,4,-8,6},
//       {-3,-4,0,1,-7,-8,4,5},
//       {2,3,-1,-2,6,7,-5,-6},
//       {-5,-6,-7,-8,0,1,2,3},
//       {4,5,6,7,-1,-2,-3,-4}};
//#endif /* FTT_3D */
  int n;
  __global struct _FttCell * c;

  return_val_if_fail (cell != NULL, NULL);
  return_val_if_fail (d < FTT_NEIGHBORS, NULL);

  if (FTT_CELL_IS_ROOT (cell))
    return ((__global struct _FttRootCell *) cell)->neighbors.c[d];

  n = neighbor_index[d][FTT_CELL_ID (cell)];
  if (n >= 0) /* neighbor belongs to same Oct */
    c = &(cell->parent->cell[n]);
  else {      /* neighbor belongs to neighboring Cell or Oct */
    c = cell->parent->neighbors.c[d];
    if (c != NULL && c->children != NULL)
      c = &(c->children->cell[- n - 1]);
  }
  if (c == NULL || FTT_CELL_IS_DESTROYED (c))
    return NULL;
  else
    return c;
}


/**
 * ftt_cell_neighbor:
 * @cell: a #FttCell.
 * @d: a direction.
 *
 * Returns: the neighbor of @cell in direction @d or %NULL if @cell
 * has no neighbor in this direction.  
 */
static inline
__global struct _FttCell * ftt_cell_neighbor (const __global struct _FttCell * cell,
			     FttDirection d)
{
  return_val_if_fail (cell != NULL, NULL);
  return_val_if_fail (d < FTT_NEIGHBORS, NULL);

  if (!FTT_CELL_IS_LEAF (cell))
    return cell->children->neighbors.c[d];

  return ftt_cell_neighbor_not_cached (cell, d);
}


/**
 * gfs_cell_face:
 * @cell: a #FttCell.
 * @d: a direction.
 *
 * This function is different from ftt_cell_face() because it takes
 * into account the solid fractions.
 *
 * Returns: the face of @cell in direction @d.
 */
struct _FttCellFace gfs_cell_face (__global struct _FttCell * cell,
			   FttDirection d)
{
  struct _FttCellFace f = {cell, NULL, d};

  return_val_if_fail (cell != NULL, f);

  if (!GFS_IS_MIXED (cell) || GFS_STATE (cell)->solid->s[d] > 0.)
    f.neighbor = ftt_cell_neighbor (cell, d);
  return f;
}

/* get_average_neighbor_value_stencil() needs to be updated whenever 
 * this function is modified
 */
static double average_neighbor_value (const struct _FttCellFace * face,
				       unsigned int v,
				       double * x)
{
  /* check for corner refinement violation (topology.fig) */
  assert (ftt_cell_level (face->neighbor) == ftt_cell_level (face->cell));
  
  if (FTT_CELL_IS_LEAF (face->neighbor))
    return GFS_VALUEI (face->neighbor, v);
  else {
    struct _FttCellChildren children;
    double av = 0., a = 0.;
    FttDirection od = FTT_OPPOSITE_DIRECTION (face->d);
    unsigned int i, n;
    
    n = ftt_cell_children_direction (face->neighbor, od, &children);
    for (i = 0; i < n; i++)
      if (children.c[i] && GFS_VALUEI (children.c[i], v) != GFS_NODATA) {
	double w = GFS_IS_MIXED (children.c[i]) ? GFS_STATE (children.c[i])->solid->s[od] : 1.;
	a += w;
	av += w*GFS_VALUEI (children.c[i], v);
      }
    if (a > 0.) {
      *x = 3./4.;
      return av/a;
    }
    else
      return GFS_VALUEI (face->cell, v);
  }
}

/* v = a*v(cell) + b 
 * 
 * First order 1D interpolation.
 *
 * get_interpolate_1D1_stencil() needs to be updated whenever 
 * this function is modified
 */
static struct _GfsGradient interpolate_1D1 (__global struct _FttCell * cell,
				    FttDirection d,
				    double x,
				    unsigned int v)
{
  struct _GfsGradient p = { 1., 0. };
  struct _FttCellFace f;

  f = gfs_cell_face (cell, d);
  if (f.neighbor) {
     double x2 = 1.;
     double p2 = average_neighbor_value (&f, v, &x2);
    if (p2 != GFS_NODATA) {
      double a2 = x/x2;
      p.b += a2*p2;
      p.a -= a2;
    }
  }

  return p;
}

static struct _Gradient gradient_fine_coarse (const struct _FttCellFace * face, unsigned int v)
{
  struct _Gradient g;
  struct _GfsGradient p;

#ifdef FTT_2D
  int dp;
#endif
//#else  
//  gint * dp;
//#endif 
 
  assert(face != NULLPrvt);

  assert (ftt_face_type (face) == FTT_FINE_COARSE);

  dp = perpendicular[face->d][FTT_CELL_ID (face->cell)];
#ifdef FTT_2D
  assert (dp >= 0);
  p = interpolate_1D1 (face->neighbor, dp, 1./4., v);
#endif
//#else  
//  g_assert (dp[0] >= 0 && dp[1] >= 0);
//  p = interpolate_2D1 (face->neighbor, dp[0], dp[1], 1./4., 1./4., v);
//#endif 

  g.a = 2./3.;
  g.b = 2.*p.a/3.;
  g.c = 2.*p.b/3.;

  return g;
}


void ftt_cell_relative_pos (const global struct _FttCell * cell,
			     struct _FttVector * pos)
{
  unsigned int n;

  return_if_fail(cell != NULL); 

  return_if_fail(pos != NULLPrvt); 

  return_if_fail(!FTT_CELL_IS_ROOT (cell)); 

  n = FTT_CELL_ID (cell);
  pos->x = coords[n * 3 + 0]/4.;
  pos->y = coords[n * 3 + 1]/4.;
  pos->z = coords[n * 3 + 2]/4.;
}



void gfs_face_gradient (const struct _FttCellFace * face,
			struct _Gradient * g,
			unsigned int v,
			int max_level)
{
  unsigned int level;

  return_if_fail(face != NULLPrvt); 

  g->a = g->b = 0.;
  if (face->neighbor == NULL || GFS_FACE_FRACTION (face) == 0.)
    return;

  level = ftt_cell_level (face->cell);
  if (ftt_cell_level (face->neighbor) < level) {
    //neighbor is at a shallower level
    struct _Gradient gcf;

    gcf = gradient_fine_coarse (face, v);
    g->a = gcf.a;
    g->b = gcf.b*GFS_VALUEI (face->neighbor, v) + gcf.c;
  }
  else {
    if (level == max_level || FTT_CELL_IS_LEAF (face->neighbor)) {
      //neighbor is at the same level 
      g->a = 1.;
      g->b = GFS_VALUEI (face->neighbor, v);
    }
    else {
     // neighbor is at a deeper level
      struct _FttCellChildren children;
      struct _FttCellFace f;
      unsigned int i, n;
      double s;
      
      f.d = FTT_OPPOSITE_DIRECTION (face->d);
      n = ftt_cell_children_direction (face->neighbor, f.d, &children);
      f.neighbor = face->cell;
      for (i = 0; i < n; i++)
	if ((f.cell = children.c[i])) {
	  struct _Gradient gcf;
	  gcf = gradient_fine_coarse (&f, v);
	  s = GFS_FACE_FRACTION (&f);
	  g->a += s*gcf.b;
	  g->b += s*(gcf.a*GFS_VALUEI (f.cell, v) - gcf.c);
	}
      s = GFS_FACE_FRACTION (face)*n/2.;
      g->a /= s;
      g->b /= s;
    }
  }
}


// to v to pernaw ws v->i
static void get_from_above (global struct _FttCell * parent, unsigned int v)
{
  unsigned int level = ftt_cell_level (parent);
  struct _FttCellNeighbors n;
  struct _FttCellChildren child;
  FttComponent c;
  struct _FttVector h;
  unsigned int i;

  ftt_cell_neighbors (parent, &n);
  for (c = 0; c < FTT_DIMENSION; c++) {
    struct _FttCellFace f;
    struct _Gradient g;
    double g1, g2;
    
    f.cell = parent;
    f.d = 2*c;
    f.neighbor = n.c[f.d];
    gfs_face_gradient (&f, &g, v, level);
    g1 = g.b - g.a*GFS_VALUE (parent, v);
    f.d = 2*c + 1;
    f.neighbor = n.c[f.d];
    gfs_face_gradient (&f, &g, v, level);
    g2 = g.b - g.a*GFS_VALUE (parent, v);
    (&h.x)[c] = (g1 - g2)/2.;
  }

  ftt_cell_children (parent, &child);
  for (i = 0; i < FTT_CELLS; i++) 
    if (child.c[i]) {
      struct _FttVector p;
      
      GFS_VALUE (child.c[i], v) = GFS_VALUE (parent, v);
      ftt_cell_relative_pos (child.c[i], &p);
      for (c = 0; c < FTT_DIMENSION; c++)
	GFS_VALUE (child.c[i], v) += (&p.x)[c]*(&h.x)[c];
    // printf("GFS_VALUE (child.c[%d], v) is %lf\n", i, GFS_VALUE (child.c[i],v));
   }
  
}



__kernel void cell_traverse_level_non_leafs_kernel_func(int max_depth, int i, __global void * RootD,  __global void * treeD, unsigned int v, char print, __global unsigned long * sizes) 
{
  int id;
  unsigned int n;
  struct _GfsStateVector a;//for debug
  struct _FttOct b;

  id = get_global_id(0);
  __global struct _FttRootCell * root = (__global struct _FttRootCell *) RootD;
  __global struct _FttOct * tree = (__global struct _FttOct *) treeD;


  if(print) {
   if(max_depth == 0){
       if(id ==0) {
         #ifdef DEBUG_MODE
         int n = get_global_size(0);
         printf("\nglobal size is %d\n", n); 

         printf("kernel : sizeof(FttOct) = %ld\n", sizeof(struct _FttOct));
         printf("kernel : sizeof(FttOct *) = %ld\n", sizeof(struct _FttOct *));
         printf("kernel : sizeof(struct FttCell) = %ld\n", sizeof(struct _FttCell));
         printf("kernel : sizeof(struct FttCell *) = %ld\n", sizeof(struct _FttCell *));
         printf("kernel : sizeof(struct _GfsStateVector) = %ld\n", sizeof(struct _GfsStateVector));
         printf("kernel : sizeof(struct _GfsStateVector *) = %ld\n", sizeof(struct _GfsStateVector *));
         printf("kernel : sizeof(struct _GfsFaceStateVector) = %ld\n", sizeof(struct _GfsFaceStateVector));
         printf("kernel : sizeof(struct _FttCellNeighbors) = %ld\n", sizeof( struct _FttCellNeighbors));
         printf("kernel : sizeof(struct _FttVector) = %ld\n", sizeof(struct _FttVector));
         printf("kernel : sizeof(struct _FttRootCell) = %ld\n", sizeof(struct _FttRootCell));
         printf("kernel : sizeof(double) =  %ld\n", sizeof(double));
         printf("kernel : sizeof(int) =  %ld\n\n\n", sizeof(int));
         #endif 
         
         //for GPU Debugging
         sizes[0] = sizeof(struct _FttOct);
         sizes[1] = sizeof(struct _FttOct *);
         sizes[2] = sizeof(struct _FttCell);
         sizes[3] = sizeof(struct _FttCell *);
         sizes[4] = sizeof(struct _GfsStateVector);
         sizes[5] = sizeof(struct _GfsStateVector *);
         sizes[6] = sizeof(struct _GfsFaceStateVector);   
         sizes[7] = sizeof( struct _FttCellNeighbors);
         sizes[8] = sizeof(struct _FttVector);
         sizes[9] = sizeof(struct _FttRootCell);     
         sizes[10] = sizeof(double);
         sizes[11] = sizeof(int);    
    
         #ifdef DEBUG_MODE
         printf("######kernel : GfsStateVector########\n");
         printf("difference f[1] - f[0] is %ld\n", (unsigned long) &(a.f[1]) - (unsigned long)&a.f[0]);
         printf("difference f[2] - f[1] is %ld\n", (unsigned long) &(a.f[2]) - (unsigned long)&a.f[1]);
         printf("difference solid - f[3] is %ld\n", (unsigned long) &(a.solid) - (unsigned long)&a.f[3]);
         printf("difference place_holder - solid is %ld\n", (unsigned long) &(a.place_holder) - (unsigned long)&a.solid);
         
         printf("###########kernel : FttOct##############\n");
         printf("difference parent -level is %ld\n", (unsigned long)&(b.parent) - (unsigned long)&(b.level));
         printf("difference neighbors  - parent is %ld\n", (unsigned long)&(b.neighbors) - (unsigned long)&(b.parent));
         printf("difference pos - neighbors is %ld\n", (unsigned long)&(b.pos) - (unsigned long)&(b.neighbors));
         printf("difference cell[0] - pos is %ld\n", (unsigned long)&(b.cell[0]) - (unsigned long)&(b.pos));
         printf("difference cell[1] - cell[0] is %ld\n", (unsigned long)&(b.cell[1]) - (unsigned long)&(b.cell[0]));        
         printf("difference cell[2] - cell[1] is %ld\n", (unsigned long)&(b.cell[2]) - (unsigned long)&(b.cell[1]));  
         printf("difference cell[3] - cell[2] is %ld\n", (unsigned long)&(b.cell[3]) - (unsigned long)&(b.cell[2])); 
         #endif

         //for GPU Debugging
         sizes[12] = (unsigned long) &(a.f[1]) - (unsigned long)&a.f[0];
         sizes[13] = (unsigned long) &(a.f[2]) - (unsigned long)&a.f[1];
         sizes[14] = (unsigned long) &(a.solid) - (unsigned long)&a.f[3];
         sizes[15] = (unsigned long) &(a.place_holder) - (unsigned long)&a.solid;
         sizes[16] = (unsigned long)&(b.parent) - (unsigned long)&(b.level);
         sizes[17] = (unsigned long)&(b.neighbors) - (unsigned long)&(b.parent);
         sizes[18] = (unsigned long)&(b.pos) - (unsigned long)&(b.neighbors);   
         sizes[19] = (unsigned long)&(b.cell[0]) - (unsigned long)&(b.pos);
         sizes[20] = (unsigned long)&(b.cell[1]) - (unsigned long)&(b.cell[0]);
         sizes[21] = (unsigned long)&(b.cell[2]) - (unsigned long)&(b.cell[1]);     
         sizes[22] = (unsigned long)&(b.cell[3]) - (unsigned long)&(b.cell[2]);
     }//id ==0 
   }
  }
 
   
  

  if(id == i) {
         __global struct _FttCell * f;
         f = (__global struct _FttCell *)&root[id]; 
         if (!FTT_CELL_IS_DESTROYED (f))
            if (ftt_cell_level (f) == max_depth && !FTT_CELL_IS_LEAF(f)) { 
             get_from_above(f, v);
           } 
  }

  if(id < 12000) {
      for(n = 0; n < FTT_CELLS; n++) {
          if (!FTT_CELL_IS_DESTROYED (&(tree[id].cell[n])) )
            if (ftt_cell_level (&(tree[id].cell[n])) == max_depth && !FTT_CELL_IS_LEAF (&(tree[id].cell[n])) ) {
              get_from_above(&(tree[id].cell[n]), v);
            } 
      } 
      
   }

};


unsigned long check_address_range(unsigned long ptr, unsigned long treeH, unsigned long rootH, unsigned long boundaryH, unsigned long treeHEnd, unsigned long rootHEnd, unsigned long boundaryHEnd)
{
	   if(rootH <= ptr  && ptr <= rootHEnd)
	        return rootH;
	   else if(treeH <= ptr && ptr <= treeHEnd)
	        return treeH;
	   else if (boundaryH <= ptr && ptr <= boundaryHEnd)
	        return boundaryH; 
          
};

/**/
__kernel void convertTreeAddresses(__global void * treeD, __global void * dataD, __global void * rootD, __global void * boundaryD, unsigned long treeH, unsigned long dataH, unsigned long rootH, unsigned long boundaryH, unsigned long treeHEnd, unsigned long rootHEnd, unsigned long boundaryHEnd, int number_of_roots)
{
    int id = get_global_id(0);
    int j, i;
    
    __global struct _FttOct * td = (__global struct _FttOct *)treeD;
    __global struct _FttRootCell * rd = (__global struct _FttRootCell *)rootD;
    __global struct _FttRootCell * bd = (__global struct _FttRootCell *)boundaryD;

    unsigned long DeviceBase;
    unsigned long HostBase;
    
   #ifdef DEBUG_MODE
    if(id ==0) {
       printf("KERNEL BEFORE: dataD is %p\n", dataD);
       printf("KERNEL BEFORE: boundaryD is %p\n", boundaryD);
       printf("KERNEL BEFORE: rootD is %p\n", rootD);
       printf("number_of_roots is %d\n", number_of_roots);
       
       for(i = 0; i < number_of_roots; i++)
         printf("KERNEL_BEFORE : rootD[%d] is %p\n", i, rd + i);
       
       printf("KERNEL BEFORE: rd[%d].cell.children = %p\n", id, rd[id].cell.children);
       printf("KERNEL BEFORE: rd[%d].cell.data = %p\n", id, rd[id].cell.data);
       for(j = 0; j < FTT_CELLS; j++) {
           printf("KERNEL BEFORE: rd[%d].neighbors.c[%d] is %p\n", id, j, rd[id].neighbors.c[j]);
      }
    
      for(i=0; i < 4; i++)
        for(j = 0; j < FTT_CELLS; j++) {
              printf("tree[%d].cell[%d].children = %p, ", id, j, td[id].cell[j].children);
              printf("tree[%d].neighbors.c[%d] = %p\n", id, j, td[id].neighbors.c[j]);
              printf("tree[%d].cell[%d].data = %p, ", id, j, td[id].cell[j].data);
              if(td[id].cell[j].data != NULL)
                printf("GFS_STATE(&(td[%d].cell[%d]))->solid = %p\n", id, j, GFS_STATE(&(td[id].cell[j]))->solid);
        }
    }
    #endif
 
    //Change of addresses 
    if(id < number_of_roots) {
          if(rd[id].cell.parent != NULL)
           rd[id].cell.parent = (__global struct _FttOct *)((unsigned long)td + (unsigned long)rd[id].cell.parent - treeH);
          if(rd[id].cell.children != NULL)
           rd[id].cell.children = (__global struct _FttOct *)((unsigned long)td + (unsigned long)rd[id].cell.children - treeH);
          if(rd[id].cell.data != NULL)
           rd[id].cell.data = (__global void *) ((unsigned long)dataD + ((unsigned long)rd[id].cell.data - dataH));           
          for(j = 0; j < FTT_NEIGHBORS; j++) {
            
            if(rd[id].neighbors.c[j] !=NULL) {
              HostBase = check_address_range((unsigned long)rd[id].neighbors.c[j] , treeH, rootH, boundaryH, treeHEnd, rootHEnd, boundaryHEnd);
              if(HostBase == treeH) {
                DeviceBase = (unsigned long)td;
               //printf("tree base\n");
              }
            else if(HostBase == rootH) {
                DeviceBase = (unsigned long)rd;
                //printf("Root: root base\n");
             }
            else if(HostBase == boundaryH)  {
               DeviceBase = (unsigned long)boundaryD;
               //printf("Root: boundary base\n");
            }
            else {
              
               //printf("Root: neighbor pointer: no range found\n");
             }
              rd[id].neighbors.c[j] = (__global struct _FttCell *)(DeviceBase + (unsigned long)rd[id].neighbors.c[j] - HostBase);
          } 
        }
    }

  if(id < 12000) {
     for(j = 0; j < FTT_NEIGHBORS; j++) {
          if(td[id].neighbors.c[j] != NULL)  {
             HostBase = check_address_range((unsigned long)td[id].neighbors.c[j], treeH, rootH, boundaryH, treeHEnd, rootHEnd, boundaryHEnd);
             if(HostBase == treeH)  {
               DeviceBase = (unsigned long)td;
               //printf("tree base\n");
             }
             else if(HostBase == rootH) {
                DeviceBase = (unsigned long)rd;
               #ifdef DEBUG_MODE
                //printf("Tree: root base\n");
               #endif
             }
             else if(HostBase == boundaryH)  {
                DeviceBase = (unsigned long)boundaryD;
               #ifdef DEBUG_MODE
                //printf("Tree: boundary base\n");
               #endif 
             }
             else {
               #ifdef DEBUG_MODE
                //printf("Tree: neighbor pointer: no range found\n");
               #endif
             }
             td[id].neighbors.c[j] = (__global struct _FttCell *)(DeviceBase + ((unsigned long)(td[id].neighbors.c[j]) - HostBase));
       }
     }
 
     for(j = 0; j < FTT_CELLS; j++) {
         if(td[id].cell[j].parent != NULL)
           td[id].cell[j].parent = (__global struct _FttOct *)((unsigned long)td + ((unsigned long)td[id].cell[j].parent - treeH));
       
         if(td[id].cell[j].children != NULL) 
           td[id].cell[j].children = (__global struct _FttOct *)((unsigned long)td + ((unsigned long)(td[id].cell[j].children) - treeH));
            
         if(td[id].cell[j].data != NULL) {
           td[id].cell[j].data = (__global void *) ((unsigned long)dataD + ((unsigned long)td[id].cell[j].data - dataH));
           if(GFS_STATE(&(td[id].cell[j]))->solid != NULL)
             GFS_STATE(&(td[id].cell[j]))->solid = (__global void *) ((unsigned long)dataD + ((unsigned long)(GFS_STATE(&(td[id].cell[j]))->solid) - dataH));
         }
     }
    }

    //Printing
    #ifdef DEBUG_MODE 
    if(id == 0) {
     for(i = 0; i < 4; i++)
       for(j = 0; j < FTT_CELLS; j++) {    
            printf("KERNEL AFTER: td[%d].cell[%d].children = %p, ", id, j, td[id].cell[j].children);
            printf("td[%d].neighbors.c[%d] = %p,  \n", id, j,  td[id].neighbors.c[j]);
            printf("KERNEL AFTER: td[%d].cell[%d].data = %p, ", id, j, td[id].cell[j].data);
            if(td[id].cell[j].data !=NULL) 
              printf("GFS_STATE(&(td[%d].cell[%d]))->solid = %p\n", id, j, GFS_STATE(&(td[id].cell[j]))->solid);
        }         
    }
   #endif

};

//####################################################################################################################################

int power(int base, int expot)
{
    int result = 1;
    while(expot) { result *= base; expot--; }
    return result;
};

/*__kernel void cell_traverse_level_non_leafs_kernel_func(int max_depth, int i, __global struct _FttRootCell * RootD, __global void * v, __global void *temp)
{

  int id, parent_id;
  unsigned int level;

  id = get_global_id(0);

  parent_id = id/4;
  if(parent_id % 4 == 0)
     parent_id = parent_id - 1;
   
  __global struct _FttCell ** carray;
  carray = (__global struct _FttCell **)&temp;
  carray[parent_id] = NULL;
  carray[0] = &(RootD[i].cell); 

 if(id<5) { 
   if(ftt_cell_level(carray[0]) == max_depth && !FTT_CELL_IS_LEAF (carray[0])) {
        if(id ==0) {
          //call function
        }
   }
   else if(!FTT_CELL_IS_LEAF (carray[0])) {
        if(id != 0) {
             // while(carray[parent_id] == NULL){} 
             // if(!FTT_CELL_IS_LEAF (carray[parent_id])) {
             // __global struct _FttOct * children = carray[parent_id]->children;
             // carray[id] = (__global struct _FttCell *)&(children->cell[id % 4]);
            
         //     if(!FTT_CELL_IS_DESTROYED(carray[id])) 
         //       if(ftt_cell_level(carray[id]) == max_depth && !(FTT_CELL_IS_LEAF(carray[id]))) {
                  //call function
         //       }
         //     } 
        }
  }
 }
};
*/


