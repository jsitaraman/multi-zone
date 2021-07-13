/* This file is part of the Tioga software library */

/* Tioga  is a tool for overset grid assembly on parallel distributed systems */
/* Copyright (C) 2015 Jay Sitaraman */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA */
# include <math.h>
# include <stdlib.h>
# include <stdio.h>
# define TIOGA_MIN(x,y)  (x) < (y) ? (x) : (y)
# define TIOGA_MAX(x,y)  (x) > (y) ? (x) : (y)
# define TIOGA_FREE(a1)  {free(a1);a1=NULL;}
# define TOL 1e-8
# define BIGVALUE 1E15
# define BASE 0
//
// Use an octree to find and remove common nodes 
//
void uniqNodesTree(double *coord,
		   int *itag,
		   int *elementsAvailable,
		   int ndim, int nav)
{  
  int nd=ndim;
  int i,j,ibox,k;
  int p1,p2;
  int *tmpint;
  int npts[8],iv[3],cft[9];
  double xmax[3],xmin[3],xmid[3],dx[3],xp[3],dist;
  int icheck=1;
  //
  // if there are more than 20 elements divide the tree
  //
  if (nav > 20) {
    //
    // find the bound of the boxes
    //
    icheck=0;
    xmin[0]=xmin[1]=xmin[2]=BIGVALUE;
    xmax[0]=xmax[1]=xmax[2]=-BIGVALUE;
    for(i=0;i<nav;i++)
      for(j=0;j<nd;j++)
	{
	  xmin[j]=TIOGA_MIN(xmin[j],coord[ndim*elementsAvailable[i]+j]);
	  xmax[j]=TIOGA_MAX(xmax[j],coord[ndim*elementsAvailable[i]+j]);
	}
    for(j=0;j<nd;j++) { 
      xmid[j]=(xmax[j]+xmin[j])*0.5;
      dx[j]=(xmax[j]-xmin[j])*0.5+TOL;
    }
    for(j=0;j<8;j++) npts[j]=0;
    for(i=0;i<nav;i++)
      {
	for(j=0;j<3;j++) {
	  xp[j]=coord[ndim*elementsAvailable[i]+j]-xmid[j];
	  iv[j]=floor(xp[j]/dx[j])+1;
          iv[j]=TIOGA_MIN(1,TIOGA_MAX(0,iv[j]));
	}
	ibox= 4*iv[0]+2*iv[1]+iv[2];
	npts[ibox]++;
      }
    for(j=0;j<8;j++) if(npts[j]==nav) icheck=1;
    if (!icheck) {
    cft[0]=0;
    for(j=0;j<8;j++)
      cft[j+1]=cft[j]+npts[j];
    tmpint=(int *)malloc(sizeof(int)*nav);
    for(i=0;i<nav;i++)
      {
	for(j=0;j<3;j++){
	  xp[j]=coord[ndim*elementsAvailable[i]+j]-xmid[j];
	  iv[j]=floor(xp[j]/dx[j])+1;
          iv[j]=TIOGA_MIN(1,TIOGA_MAX(0,iv[j]));
	  }
	  ibox= 4*iv[0]+2*iv[1]+iv[2];
	  npts[ibox]=npts[ibox]-1;
	  tmpint[npts[ibox]+cft[ibox]]=elementsAvailable[i];
	}
      for(i=0;i<nav;i++)
	elementsAvailable[i]=tmpint[i];
      TIOGA_FREE(tmpint);
      for(j=0;j<8;j++)
	if (cft[j+1] > cft[j])
	  uniqNodesTree(coord,itag,&(elementsAvailable[cft[j]]),
			ndim,cft[j+1]-cft[j]);
    }
  }
  if (icheck) {
    for(i=0;i<nav;i++)
      {
	p1=elementsAvailable[i];
	for(j=i+1;j<nav;j++)
	  {
	    p2=elementsAvailable[j];	    
	    dist=0.0;
            for(k=0;k<nd;k++) dist+=fabs(coord[ndim*p1+k]-coord[ndim*p2+k]);
            if (p1==25456) printf("%d %d %f\n",p1,p2,dist);
	    if (dist < TOL)
	      {
		if (p1 > p2) {
		  itag[p1]=itag[p2];
		}
		else {
		  itag[p2]=itag[p1];
		}
	      }
	  }
      }
  }
}
/*
 * Create a unique hash for list of coordinates with duplicates in 
 * them. Find the rtag as max of all duplicate samples. itag contains
 * the hash to the real point
 */
void uniquenodes_octree_(double *x,int *itag,int *nn)
{
  int nelem=*nn;
  int *elementsAvailable;
  int i;
  printf("nn=%d\n",*nn);
  i=23175;
  printf("%f %f %f\n",x[3*i],x[3*i+1],x[3*i+2]);
  i=25456;
  printf("%f %f %f\n",x[3*i],x[3*i+1],x[3*i+2]);

  elementsAvailable=(int *)malloc(sizeof(int)*nelem);
  for(i=0;i<nelem;i++)
    {
     elementsAvailable[i]=i;
     itag[i]=i+BASE;
    }

  uniqNodesTree(x,itag,elementsAvailable,3,nelem);
  TIOGA_FREE(elementsAvailable);
}

/*
 * Create a unique hash for list of coordinates with duplicates in 
 * them. Find the rtag as max of all duplicate samples. itag contains
 * the hash to the real point
 */
void uniquenodes(double *x,int *itag,int *nn)
{
  int NSUB=101;
  int i,j,k,m,ij,i3,jj,kk,ll,p1,p2,indx,jmax,kmax,lmax,nsblks,jkmax;
  double xmax[3],xmin[3],ds,dsi,dsx,dsxi,dsy,dsyi,dsz,dszi;
  int *cft,*numpts,*ilist;
  int nnodes=*nn;

  for(j=0;j<3;j++) xmax[j]=-1E15;
  for(j=0;j<3;j++) xmin[j]=1E15;
  
  for(i=0;i<nnodes;i++)
    for(j=0;j<3;j++) {
      xmax[j]=TIOGA_MAX(xmax[j],x[3*i+j]);
      xmin[j]=TIOGA_MIN(xmin[j],x[3*i+j]);
    }

  ds=(xmax[0]-xmin[0]+xmax[1]-xmin[1]+xmax[2]-xmin[2])/3.0/NSUB;
  dsi=1.0/ds;
  for(j=0;j<3;j++) xmax[j]+=ds;
  for(j=0;j<3;j++) xmin[j]-=ds;
  
  jmax=TIOGA_MIN(round((xmax[0]-xmin[0])*dsi),NSUB);
  jmax=TIOGA_MAX(jmax,1);
  dsx=(xmax[0]-xmin[0]+TOL)/jmax;
  dsxi=1./dsx;    
  kmax=TIOGA_MIN(round((xmax[1]-xmin[1])*dsi),NSUB);
  kmax=TIOGA_MAX(kmax,1);
  dsy=(xmax[1]-xmin[1]+TOL)/kmax;
  dsyi=1./dsy;
  lmax=TIOGA_MIN(round((xmax[2]-xmin[2])*dsi),NSUB);
  lmax=TIOGA_MAX(lmax,1);
  dsz=(xmax[2]-xmin[2]+TOL)/lmax;
  dszi=1./dsz;
  nsblks=jmax*kmax*lmax;
  jkmax=jmax*kmax;
  cft=(int *)malloc(sizeof(int)*(nsblks+1));
  numpts=(int *)malloc(sizeof(int)*nsblks);
  ilist=(int *)malloc(sizeof(int)*nnodes);

  for(i=0;i<nsblks;i++) numpts[i]=0;
  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      jj=(int)((x[i3]-xmin[0])*dsxi);
      kk=(int)((x[i3+1]-xmin[1])*dsyi);
      ll=(int)((x[i3+2]-xmin[2])*dszi);
      indx=ll*jkmax+kk*jmax+jj;
      numpts[indx]=numpts[indx]+1;
    }

  cft[0]=0;
  for(i=0;i<nsblks;i++) cft[i+1]=cft[i]+numpts[i];
  
  for(i=0;i<nnodes;i++)
    {
      i3=3*i;
      jj=(int)((x[i3]-xmin[0])*dsxi);
      kk=(int)((x[i3+1]-xmin[1])*dsyi);
      ll=(int)((x[i3+2]-xmin[2])*dszi);
      indx=ll*jkmax+kk*jmax+jj;
      ilist[cft[indx]+numpts[indx]-1]=i;
      numpts[indx]--;
      itag[i]=i;
    }
  
  for(i=0;i<nsblks;i++)
    for(j=cft[i];j<cft[i+1];j++)
      {
	p1=ilist[j];
	for(k=j+1;k<cft[i+1];k++)
	  {
	    p2=ilist[k];
	    if ( fabs(x[3*p1  ]-x[3*p2  ])+
		 fabs(x[3*p1+1]-x[3*p2+1])+
		 fabs(x[3*p1+2]-x[3*p2+2]) < TOL)
	      {
		if (p1 > p2) {
		  itag[p1]=itag[p2];
		}
		else {
		  itag[p2]=itag[p1];
		}
	      }
	  }
      }
  /*
  m=0;
  for(i=0;i<nnodes;i++)
    if (itag[i]==i+1) {
     m++;
   }
  */
  TIOGA_FREE(ilist);
  TIOGA_FREE(cft);
  TIOGA_FREE(numpts);
}
