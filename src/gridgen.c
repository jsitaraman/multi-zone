#include <math.h>
#include <stdlib.h>
#include "codetypes.h"

void get_annulus_grid(double outer_d1,
		      double outer_d2,
		      double outer_h1,
		      double outer_dx,
		      double trans[3],
		      GRID *g)
{
  int i,j,m,n,offset,nh,nnodes2d,ntri2d,*tri2d,ndom,bb[2],nbnodes;
  double psi,dh,*xb,*yb,*nsf,*x2d;
  double TOL=1e-8;

  ndom=1;
  if (fabs(outer_d1-outer_d2) > TOL) ndom=2;
  
  bb[0]=round(M_PI*outer_d1/outer_dx);
  nbnodes=bb[0];
  if (ndom > 1) {
    bb[1]=round(M_PI*outer_d2/outer_dx);
    nbnodes+=bb[1];
  }
  
  xb=(double *)malloc(sizeof(double)*nbnodes);
  yb=(double *)malloc(sizeof(double)*nbnodes);
  nsf=(double *)malloc(sizeof(double)*nbnodes);

  for(i=0;i<bb[0];i++)
    {
      psi=i*2*M_PI/bb[0];
      xb[i]=outer_d1*0.5*cos(psi);
      yb[i]=outer_d1*0.5*sin(psi);
      if (i > 0) {
	nsf[i]=sqrt((xb[i]-xb[i-1])*(xb[i]-xb[i-1])+
		    (yb[i]-yb[i-1])*(yb[i]-yb[i-1]));
      }
    }
  nsf[0]=nsf[1];
  if (ndom > 1) {
    for(i=bb[0];i<nbnodes;i++)
      {
	psi=2*M_PI-i*2*M_PI/bb[1];
	xb[i]=outer_d2*0.5*cos(psi);
	yb[i]=outer_d2*0.5*sin(psi);
	if (i > bb[0]) {
	  nsf[i]=sqrt((xb[i]-xb[i-1])*(xb[i]-xb[i-1])+
		      (yb[i]-yb[i-1])*(yb[i]-yb[i-1]));
	}
      }
    nsf[bb[0]]=nsf[bb[0]+1];
  }
  /* perform 2-D triangulation using ugrid */
  
  generate_grid_(xb,yb,nsf,bb,&nbnodes,&ndom);
  get_elem_count_(&nnodes2d,&ntri2d);
  x2d=(double *)malloc(sizeof(double)*nnodes2d*2);
  tri2d=(int *)malloc(sizeof(int)*ntri2d*3);
  get_tess_(x2d,tri2d);
  //for(i=0;i<nbnodes;i++) { x2d[2*i]=xb[i];x2d[2*i+1]=yb[i];}

  /* create prismatic mesh by extruding the annulus */
  
  nh=outer_h1/outer_dx;
  dh=outer_h1/nh;

  g->nnodes=(nnodes2d)*(nh+1);
  g->ntypes=1;
  g->nv=(int *)malloc(sizeof(int)*g->ntypes);
  g->ncells=(int *)malloc(sizeof(int)*g->ntypes);
  g->nv[0]=6;
  g->ncells[0]=ntri2d*nh;
  g->x=(double *)malloc(sizeof(double)*g->nnodes*3);
  g->vcft=(int *)malloc(sizeof(int)*(g->ntypes+1));
  g->vconn=(long *)malloc(sizeof(long)*g->ncells[0]*6);
  g->vcft[0]=0;
  g->vcft[1]=g->ncells[0];

  m=n=0;
  offset=-1;
  for(j=0;j<nh+1;j++)
    {
      for(i=0;i<nnodes2d;i++)
	{	      
	  g->x[m++]=x2d[2*i]   + trans[0];
	  g->x[m++]=x2d[2*i+1] + trans[1];
	  g->x[m++]=j*dh-outer_h1*0.5 + trans[2];
	}
      if (j < nh)
	{
	  for(i=0;i<ntri2d;i++)
	    {
	      g->vconn[n++]=(long)(tri2d[3*i]+offset);
	      g->vconn[n++]=(long)(tri2d[3*i+1]+offset);
	      g->vconn[n++]=(long)(tri2d[3*i+2]+offset);
	      g->vconn[n++]=(long)(tri2d[3*i]+offset+nnodes2d);
	      g->vconn[n++]=(long)(tri2d[3*i+1]+offset+nnodes2d);
	      g->vconn[n++]=(long)(tri2d[3*i+2]+offset+nnodes2d);
	    }
	  offset+=nnodes2d;
	}
    }

    
  free(x2d);
  free(tri2d);
  free(xb);
  free(yb);
  free(nsf);
}

