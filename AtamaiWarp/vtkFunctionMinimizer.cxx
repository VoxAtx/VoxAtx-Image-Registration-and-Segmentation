
/*=========================================================================

  Program:   AtamaiRegistration for VTK
  Module:    $RCSfile: vtkFunctionMinimizer.cxx,v $
  Language:  C++
  Date:      $Date: 2007/08/24 20:02:25 $
  Version:   $Revision: 1.2 $

Copyright (c) 2006 Atamai, Inc.
All rights reserved.

THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE SOFTWARE "AS IS"
WITHOUT EXPRESSED OR IMPLIED WARRANTY INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE.  IN NO EVENT SHALL ANY COPYRIGHT HOLDER OR OTHER PARTY WHO MAY
MODIFY AND/OR REDISTRIBUTE THE SOFTWARE UNDER THE TERMS OF THIS LICENSE
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, LOSS OF DATA OR DATA BECOMING INACCURATE
OR LOSS OF PROFIT OR BUSINESS INTERRUPTION) ARISING IN ANY WAY OUT OF
THE USE OR INABILITY TO USE THE SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.
=========================================================================*/

#include "vtkFunctionMinimizer.h"
#include "vtkObjectFactory.h"

//----------------------------------------------------------------------------
static double amotry(double **p, double *y, double *ptry, double *psum, 
		     int ndim, void (*funk)(void *data), void *data,
		     double *result, int ihi, double fac)
{
  int j;
  double fac1,fac2,ytry;

  fac1 = (1.0-fac)/(double)ndim;
  fac2 = fac1-fac;
  for (j = 0; j < ndim; j++)
    {
    ptry[j] = psum[j]*fac1 - p[ihi][j]*fac2;
    }
  (*funk)(data);
  ytry = *result; 
  if (ytry < y[ihi]) 
    {
    y[ihi] = ytry;
    for (j = 0; j < ndim; j++)
      {
      psum[j] += ptry[j]-p[ihi][j];
      p[ihi][j] = ptry[j];
      }
    }
  return ytry;
}

//----------------------------------------------------------------------------
static int amoeba(double **p, double *y, double *ptry, int ndim, double ftol,
	   void (*funk)(void *data), void *data, double *result,
	   int *nfunk, int maxnfunk)
{
  int i,ihi,ilo,inhi,j,mpts;
  double rtol,sum,swap,ysave,ytry;
  double *psum = new double[ndim];
 
  mpts = ndim+1;
  *nfunk = 0;

  for (j = 0; j < ndim; j++)
    {
    sum = 0.0;
    for (i = 0; i < mpts; i++)
     {
     sum += p[i][j];
     }
    psum[j] = sum;
    }

  for (;;) 
    {
    ilo = 0;
    if (y[0] > y[1]) 
      {
      ihi = 0;
      inhi = 1;
      }
    else 
      {
      ihi = 1;
      inhi = 0;
      }
    for (i = 0; i < mpts; i++)
      {
      if (y[i] <= y[ilo]) 
	ilo = i;
      if (y[i] > y[ihi]) 
	{
	inhi = ihi;
	ihi = i;
	}
      else if (y[i] > y[inhi] && i != ihi)
	inhi = i;
      }
    
    if (fabs(y[ihi])+fabs(y[ilo]) < ftol)
      rtol = double(2.0*fabs(y[ihi]-y[ilo]));
    else
      rtol = double(2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])));

    if (rtol < ftol) {
	swap = y[1];
	y[1] = y[ilo];
	y[ilo] = swap;
	for (i = 0; i < ndim; i++) {
	  swap = p[0][i];
	  p[0][i] = p[ilo][i];
	  p[ilo][i] = swap;
	}
	break;
    }
    if (*nfunk >= maxnfunk) {
      delete [] psum;
      return -1;      /* break if greater than max number of func evals */
    }

    *nfunk += 2;
    ytry = amotry(p,y,ptry,psum,ndim,funk,data,result,ihi,double(-1.0));
    if (ytry <= y[ilo])
      ytry = amotry(p,y,ptry,psum,ndim,funk,data,result,ihi,double(2.0));
    else if (ytry >= y[inhi]) {
      ysave = y[ihi];
      ytry = amotry(p,y,ptry,psum,ndim,funk,data,result,ihi,double(0.5));
      if (ytry >= ysave) {
	for (i = 0; i < mpts; i++) {
	  if (i != ilo) {
	    for (j = 0; j < ndim; j++)
	      p[i][j] = ptry[j] = psum[j] = (p[i][j] + p[ilo][j])/double(2.0);
	    (*funk)(data);
	    y[i] = *result;
	  }
	}
	*nfunk += ndim;
	
	for (j = 0; j < ndim; j++) {
	  sum = 0;
	  for (i = 0; i < mpts; i++)
	    sum += p[i][j];
	  psum[j] = sum;
	}
      }
    }
    else
      --(*nfunk);
  }
  
  delete [] psum;
  
  return 0;
}

//----------------------------------------------------------------------------
static double minimize(double *parameters, double **vertices, int ndim, 
		       void (*funk)(void *data), void *data, double *result,
		       double tolerance, int maxiterations, int *iterations)
{
  double *y = new double[ndim+1];

  for (int k = 0; k < ndim+1; k++) 
    {
    for (int l = 0; l < ndim; l++)
      {
      parameters[l] = vertices[k][l];
      }

    (*funk)(data);
    y[k] = *result;
    }

  amoeba(vertices,y,parameters,ndim,tolerance,funk,data,result,
	 iterations,maxiterations);