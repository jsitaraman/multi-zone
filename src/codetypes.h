typedef struct 
{
  int nnodes;  /* number of nodes */
  int ntypes;  /* number of types of polyhedra */ 
  int *nv;     /* number of vertices per polyhedra */
  int *ncells; /* number of cells per polyhedra type */
  double *x;   /* vertex coordinates */
  long *vconn;  /* connectivity of all polyhedra in cumulative format */
  int *vcft;   /* cumulative frequency table to index into vconn */
  int ntri,nquad; /* number of triangles and quads in the exposed cell faces */
  long *eface;   /* connectivity of triangles followed by quads of exposed cell faces */
  int npatch;
  int *patchid;
} GRID;

void generate_grid_(double *, double *, double *, int *, int *, int *);
void get_elem_count_(int *, int *);
void get_tess_(double *, int *);
void parseInputs(char *inputfile,
		 double *outer_d1,double *outer_d2,
		 double *outer_h1, double *outer_h2,
		 double *outer_dx,
		 double *inner_d1,
		 double *inner_h1,
		 double *inner_dx);
void writegrid_tecplot(GRID *g, char *fname);
void writegridsurface_tecplot(GRID *g, char *fname);
void get_annulus_grid(double outer_d1,
		      double outer_d2,
		      double outer_h1,
		      double outer_dx,
		      double trans[3],
		      GRID *g);
void uniquenodes_octree_(double *x,int *itag,int *nn);
void uniquenodes(double *x,int *itag,int *nn);
void merge_grids(GRID *g1, GRID *g2);
void separate_patches(GRID *g, double outer_d1,double outer_d2, double outer_h1, double outer_h2);
