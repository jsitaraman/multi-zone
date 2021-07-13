#include<stdio.h>
#include "codetypes.h"
void parseInputs(char *inputfile,
		 double *outer_d1,double *outer_d2,
		 double *outer_h1, double *outer_h2,
		 double *outer_dx,
		 double *inner_d1,
		 double *inner_h1,
		 double *inner_dx)
{
  FILE *fp;
  char line[256];
  char comments[100];
  fp=fopen(inputfile,"r");
  fgets(line,256,fp);  sscanf(line,"outer_d1=%lf",outer_d1);
  fgets(line,256,fp);  sscanf(line,"outer_d2=%lf",outer_d2);
  fgets(line,256,fp);  sscanf(line,"outer_h1=%lf",outer_h1);
  fgets(line,256,fp);  sscanf(line,"outer_h2=%lf",outer_h2);
  fgets(line,256,fp);  sscanf(line,"outer_dx=%lf",outer_dx);
  fgets(line,256,fp);  sscanf(line,"inner_d1=%lf",inner_d1);
  fgets(line,256,fp);  sscanf(line,"inner_h1=%lf",inner_h1);
  fgets(line,256,fp);  sscanf(line,"inner_dx=%lf",inner_dx);
  fclose(fp);
}

void writegrid_tecplot(GRID *g, char *fname)
{
  int i;
  FILE *fp;
  
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"triangle file\"\n");
  fprintf(fp,"VARIABLES =\"X\", \"Y\", \"Z\", \"\n");
  fprintf(fp,"ZONE T=\"PRISMS\",N= %d ,E= %d ,ET=BRICK F=FEPOINT\n",g->nnodes,g->ncells[0]);
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f %f %f\n",g->x[3*i],g->x[3*i+1],g->x[3*i+2]);
  for(i=0;i<g->ncells[0];i++)
    fprintf(fp,"%d %d %d %d %d %d %d %d\n",
	    g->vconn[6*i]+1,g->vconn[6*i+1]+1,g->vconn[6*i+2]+1,g->vconn[6*i+2]+1,
	    g->vconn[6*i+3]+1,g->vconn[6*i+4]+1,g->vconn[6*i+5]+1,g->vconn[6*i+5]+1);	    
  fclose(fp);
}


void writegridsurface_tecplot(GRID *g, char *fname)
{
  int i,m;
  FILE *fp;
  int *btag;
  
  fp=fopen(fname,"w");
  fprintf(fp,"TITLE =\"triangle file\"\n");
  fprintf(fp,"VARIABLES =\"X\", \"Y\", \"Z\",\"BTAG\",\"\n");
  fprintf(fp,"ZONE T=\"EXP\",N= %d ,E= %d ,ET=QUADRILATERAL F=FEPOINT\n",
	  g->nnodes,g->ntri+g->nquad);
  btag=(int *)calloc(g->nnodes,sizeof(int));

  for(i=0;i<g->ntri;i++)
    for(m=0;m<3;m++)
      btag[g->eface[3*i+m]]=g->patchid[i];
  for(i=0;i<g->nquad;i++)
    for(m=0;m<4;m++)
      btag[g->eface[3*g->ntri+4*i+m]]=g->patchid[i+g->ntri];

  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f %f %f %d\n",g->x[3*i],g->x[3*i+1],g->x[3*i+2],btag[i]);
  m=0;
  for(i=0;i<g->ntri;i++)
    {
      fprintf(fp,"%ld %ld %ld %ld\n",g->eface[m]+1,g->eface[m+1]+1,g->eface[m+2]+1,g->eface[m+2]+1);
      m+=3;
    }

  for(i=0;i<g->nquad;i++)
    {
      fprintf(fp,"%ld %ld %ld %ld\n",g->eface[m]+1,g->eface[m+1]+1,g->eface[m+2]+1,g->eface[m+3]+1);
      m+=4;
    }
  fclose(fp);

}
