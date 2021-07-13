#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "cgnslib.h"
#include "codetypes.h"
void write_cgns(char *fname,int gridcount,...)
{
  int z,j,ntri,nquad,p,i;
  long isize[3];
  int index_file,index_base,index_zone,index_cx,index_cy,index_cz,index_ele,index_section;
  long nelem_start,nelem_end,nbdyelem,elems_sofar;
  double *rxyz;
  long *vtmp,*ielem;
  int *pcount;
  char basename[100],zonename[100],solname[100],intstr[10],sectionname[100];
  va_list gridlist;
  GRID *g;


  sprintf(basename, "Base");
  if (cg_open(fname,CG_MODE_WRITE,&index_file)) cg_error_exit();
  if (cg_base_write(index_file, basename, 3,3, &index_base)) cg_error_exit();

  va_start(gridlist,gridcount);
  for(z=0;z<gridcount;z++)
    {
      g=(GRID *)va_arg(gridlist, void *);	    
      isize[0] = g->nnodes;
      isize[1] = g->ncells[0];
      isize[2] = 0;
      rxyz=(double *)malloc(sizeof(double)*g->nnodes*3);
      vtmp=(long *)malloc(sizeof(long) *g->ncells[0]*6); 
      sprintf(zonename, "Zone %d",z);
      if(cg_zone_write(index_file,index_base,zonename,isize, Unstructured, &index_zone))
	cg_error_exit();

      // write x-grid coordinates
      for (j = 0; j < g->nnodes; j++) rxyz[j] = g->x[3*j  ];
      if(cg_coord_write(index_file,index_base,index_zone,RealDouble,"CoordinateX",
			rxyz,&index_cx)) cg_error_exit();
      
      // write y-grid coordinates
      for (j = 0; j < g->nnodes; j++) rxyz[j] = g->x[3*j+1];
      if(cg_coord_write(index_file,index_base,index_zone,RealDouble,"CoordinateY",
			rxyz,&index_cy)) cg_error_exit();
      
      // write z-grid coordinates
      for (j = 0; j < g->nnodes; j++) rxyz[j] = g->x[3*j+2];
      if(cg_coord_write(index_file,index_base,index_zone,RealDouble,"CoordinateZ",
			rxyz,&index_cz)) cg_error_exit();

      nelem_start=0;
      nelem_end=g->ncells[0]-1;
      nbdyelem=0;
      for(j=0;j<6*g->ncells[0];j++) vtmp[j]=g->vconn[j]+1;
      if (cg_section_write(index_file,index_base,index_zone,"Elem",PENTA_6,
			  nelem_start,nelem_end,nbdyelem,vtmp,&index_ele)) cg_error_exit();

      pcount=(int *)calloc(g->npatch,sizeof(int));
      for(i=0;i<g->ntri+g->nquad;i++)  pcount[g->patchid[i]]++;
      ielem=(long *) malloc(sizeof(long)*(3*g->ntri+4*g->nquad));
      for(p=0;p<g->npatch;p++) {
	ntri=nquad=0;
	for(i=0;i<g->ntri;i++)
	  if (g->patchid[i]==p) {
	    ielem[3*ntri]=g->eface[3*i]+1;
	    ielem[3*ntri+1]=g->eface[3*i+1]+1;
	    ielem[3*ntri+2]=g->eface[3*i+2]+1;
	    ntri++;
	  }
	if (ntri > 0) {
	  elems_sofar=nelem_end;
	  nelem_start=nelem_end+1;
	  nelem_end=nelem_start+ntri-1;
	  sprintf(sectionname,"TRI%d",p);
	  if (cg_section_write(index_file,index_base,index_zone,sectionname,TRI_3,
			  nelem_start,nelem_end,nbdyelem,ielem,&index_section)) cg_error_exit();
	}
	for(i=0;i<g->nquad;i++)
	  if (g->patchid[i+g->ntri]==p) {
	    ielem[4*nquad]=g->eface[3*g->ntri+4*i]+1;
	    ielem[4*nquad+1]=g->eface[3*g->ntri+4*i+1]+1;
	    ielem[4*nquad+2]=g->eface[3*g->ntri+4*i+2]+1;
	    ielem[4*nquad+3]=g->eface[3*g->ntri+4*i+3]+1;
	    nquad++;
	  }
	if (nquad > 0) {
	  elems_sofar=nelem_end;
	  nelem_start=nelem_end+1;
	  nelem_end=nelem_start+nquad-1;
	  sprintf(sectionname,"QUAD%d",p);
	  if (cg_section_write(index_file,index_base,index_zone,sectionname,QUAD_4,
			  nelem_start,nelem_end,nbdyelem,ielem,&index_section)) cg_error_exit();
	}
      }
      free(rxyz);
      free(vtmp);
      free(pcount);
      free(ielem);
    } 
    va_end(gridlist);
    if(cg_close(index_file)) cg_error_exit();     
}
