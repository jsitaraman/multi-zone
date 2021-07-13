#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "codetypes.h"

/* interfaces to boundary_interface.f90 */
void get_exposed_faces_prizms_(int *,int *);
void get_face_count_(int *,int *);
void get_faces_(int *, int *);

# define MIN(x,y)  (x) < (y) ? (x) : (y)
# define MAX(x,y)  (x) > (y) ? (x) : (y)


/* merge g2 into to g1 */
/* assuming only one type (prizms) now */
/* will fix it to be general later if required */

void merge_grids(GRID *g1, GRID *g2)
{
  int i,j,k,*fmap,*itag;
  long *vtmp;
  double *xtmp,a[3],b[3],c[3],n[3];
  printf("Merging grids..\n"); 
  xtmp=(double *)malloc(sizeof(double)*g1->nnodes*3);
  vtmp=(long *)malloc(sizeof(long)*g1->ncells[0]*6);
  
  for(i=0;i<3*g1->nnodes;i++) xtmp[i]=g1->x[i];
  for(i=0;i<6*g1->ncells[0];i++) vtmp[i]=g1->vconn[i];

  free(g1->x);
  free(g1->vconn);  

  g1->x=(double *)malloc(sizeof(double)*3*(g1->nnodes+g2->nnodes));
  g1->vconn=(long *)malloc(sizeof(long)*6*(g1->ncells[0]+g2->ncells[0]));

  for(i=0;i<3*g1->nnodes;i++) g1->x[i]=xtmp[i];
  for(i=0;i<3*g2->nnodes;i++) g1->x[i+3*g1->nnodes]=g2->x[i];
  for(i=0;i<6*g1->ncells[0];i++) g1->vconn[i]=vtmp[i];
  for(i=0;i<6*g2->ncells[0];i++) g1->vconn[i+6*g1->ncells[0]]=g2->vconn[i]+g1->nnodes;

  g1->nnodes+=g2->nnodes;
  g1->ncells[0]+=g2->ncells[0];

  itag=(int *)malloc(sizeof(int)*g1->nnodes);
  fmap=(int *)malloc(sizeof(int)*g1->nnodes);

  /* remove common nodes using a floating point check
     ok, this is a lazy implementation now, will do better
     if required */
  
  uniquenodes(g1->x,itag,&g1->nnodes);
  j=0;
  for(i=0;i<g1->nnodes;i++)
    if (itag[i]==i) {
      for(k=0;k<3;k++) g1->x[3*j+k]=g1->x[3*i+k];
      fmap[i]=j;
      j++;
    }
  printf("Number of nodes before and after merge (%d %d)\n",g1->nnodes,j);
  g1->nnodes=j;

  for(i=0;i<6*g1->ncells[0];i++) g1->vconn[i]=fmap[itag[g1->vconn[i]]];

  for(i=0;i<g1->ncells[0];i++) {
    for(j=0;j<3;j++) {
      a[j]=g1->x[3*g1->vconn[6*i]+j];
      b[j]=g1->x[3*g1->vconn[6*i+1]+j];
      c[j]=g1->x[3*g1->vconn[6*i+2]+j];
      b[j]-=a[j];
      c[j]-=a[j];
    }
    n[0]=b[1]*c[2]-b[2]*c[1];
    n[1]=b[2]*c[0]-b[0]*c[2];
    n[2]=b[0]*c[1]-b[1]*c[0];
    if (n[2] < 0) printf("normal is negative\n");
  }

  free(fmap);
  free(itag);
  free(vtmp);
  free(xtmp);
}


void extract_faces(GRID *g)
{
  int i,*ndc6,*nbf,nfaces;

  ndc6=(int *)malloc(sizeof(int)*g->ncells[0]*6);
  for(i=0;i<g->ncells[0]*6;i++) ndc6[i]=(int)g->vconn[i];

  /* leaving f90 as 32 bit integers for now 
     will change this code to C when I get time */
  get_exposed_faces_prizms_(ndc6,&g->ncells[0]);
  get_face_count_(&g->ntri,&g->nquad);
  nbf=(int *)malloc(sizeof(int)*(3*g->ntri+4*g->nquad));
  g->eface=(long *)malloc(sizeof(long)*(3*g->ntri+4*g->nquad));
  nfaces=3*g->ntri+4*g->nquad;
  get_faces_(nbf,&nfaces);

  for(i=0;i<nfaces;i++) g->eface[i]=(long)nbf[i];
  g->patchid=(int *) calloc((g->ntri+g->nquad),sizeof(int));
  g->npatch=1; 
  free(ndc6);
  free(nbf);  
}

void separate_patches(GRID *g, double outer_d1,double outer_d2, double outer_h1, double outer_h2)
{
  int i,j,k;
  double xmin[3],xmax[3];
  double TOL=1e-8;
  
  for(i=0;i<g->ntri;i++)
    {
      xmin[0]=xmin[1]=xmin[2]=1E15;
      xmax[0]=xmax[1]=xmax[2]=-1E15;      
      for(j=0;j<3;j++)
	for(k=0;k<3;k++)
	  {
	    xmin[k]=MIN(xmin[k],g->x[3*g->eface[3*i+j]+k]);
	    xmax[k]=MAX(xmax[k],g->x[3*g->eface[3*i+j]+k]);
	  }
      if (fabs(xmax[2]+outer_h1*0.5) < TOL) g->patchid[i]=0;
      if (fabs(xmin[2]-outer_h1*0.5) < TOL) g->patchid[i]=1;
      if (fabs(xmax[2]+outer_h2*0.5) < TOL) g->patchid[i]=2;
      if (fabs(xmax[2]-outer_h2*0.5) < TOL) g->patchid[i]=3;
    }
  
  for(i=0;i<g->nquad;i++)
    {
      xmin[0]=xmin[1]=xmin[2]=1E15;
      xmax[0]=xmax[1]=xmax[2]=-1E15;      
      for(j=0;j<4;j++)
	for(k=0;k<3;k++)
	  {
	    xmin[k]=MIN(xmin[k],g->x[3*g->eface[3*g->ntri+4*i+j]+k]);
	    xmax[k]=MAX(xmax[k],g->x[3*g->eface[3*g->ntri+4*i+j]+k]);
	  }
      if (sqrt(xmin[0]*xmin[0]+
	       xmin[1]*xmin[1]) > outer_d1*0.3) {
	g->patchid[i+g->ntri]=4;
      } else {
	g->patchid[i+g->ntri]=5;
      }
    }
  g->npatch=6;
}
