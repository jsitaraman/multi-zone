#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "codetypes.h"
void main(int argc, char *argv[])
{
  double outer_d1,outer_d2,outer_h1,outer_h2,outer_dx;
  double inner_d1,inner_h1,inner_dx;  
  GRID gannulus,gcyl_top,gcyl_bottom,gcyl_inner;
  double trans[3];
    
  parseInputs("input.mz",&outer_d1,&outer_d2,&outer_h1,&outer_h2,
	      &outer_dx,&inner_d1,&inner_h1,&inner_dx);


  trans[0]=0.0;trans[1]=0.0;trans[2]=0.0;
  get_annulus_grid(outer_d1,outer_d2,outer_h1,outer_dx,trans,&gannulus);
  trans[0]=0.0;trans[1]=0.0;trans[2]=(outer_h1+outer_h2)*0.25;
  get_annulus_grid(outer_d2,outer_d2,(outer_h1-outer_h2)*0.5,outer_dx,trans,&gcyl_top);
  trans[0]=0.0;trans[1]=0.0;trans[2]=-(outer_h1+outer_h2)*0.25;
  get_annulus_grid(outer_d2,outer_d2,(outer_h1-outer_h2)*0.5,outer_dx,trans,&gcyl_bottom);

  trans[0]=0.0;trans[1]=0.0;trans[2]=0.0;
  get_annulus_grid(inner_d1,inner_d1,inner_h1,inner_dx,trans,&gcyl_inner);

  /*
  writegrid_tecplot(&gannulus, "annulus.dat");
  writegrid_tecplot(&gcyl_top, "topcyl.dat");
  writegrid_tecplot(&gcyl_bottom,"bottomcyl.dat");
  */
  merge_grids(&gannulus,&gcyl_top);
  merge_grids(&gannulus,&gcyl_bottom);
  //writegrid_tecplot(&gannulus,"merged.dat");
  extract_faces(&gannulus);
  extract_faces(&gcyl_inner);
  separate_patches(&gannulus,outer_d1,outer_d2,outer_h1,outer_h2);
  //writegridsurface_tecplot(&gcyl_inner,"exposed.dat");
  write_cgns("test.cgns",2,&gannulus,&gcyl_inner);
}
