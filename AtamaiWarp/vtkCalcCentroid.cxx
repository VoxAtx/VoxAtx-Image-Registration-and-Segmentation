/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkCalcCentroid.cxx,v $
  Language:  C++
  Date:      $Date: 2007/08/24 20:02:25 $
  Version:   $Revision: 1.8 $
  Thanks:    Thanks to Yves who developed this class.

Copyright (c) 1993-2000 Ken Martin, Will Schroeder, Bill Lorensen 
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

 * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
   of any contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

 * Modified source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#include "vtkCalcCentroid.h"
#include "vtkObjectFactory.h"

//--------------------------------------------------------------------------
vtkCalcCentroid* vtkCalcCentroid::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkCalcCentroid");
  if(ret)
    {
    return (vtkCalcCentroid*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkCalcCentroid;
}

//--------------------------------------------------------------------------
// Constructs with initial 0 values.
vtkCalcCentroid::vtkCalcCentroid()
{
  this->Input = NULL;
  this->Centroid[0] = 0.0;
  this->Centroid[1] = 0.0;
  this->Centroid[2] = 0.0;
  for (int i=0; i<=8; i++) 
    this->CovarianceMatrix[i]=0.0;
}

//--------------------------------------------------------------------------
vtkCalcCentroid::~vtkCalcCentroid()
{
}

//--------------------------------------------------------------------------
// Function to set up the covariance matrix
template <class T>
static int vtkCalculateCovarianceMatrix(vtkImageData * input, 
                                        T *inPtr,
                                        double *centroid,
                                        int *inputExtent, 
                                        double *covar)
{
  vtkIdType inInc0, inInc1, inInc2;
  int idx0, idx1, idx2;
  T *inPtr0, *inPtr1, *inPtr2;

  double sxx = 0.0, sxy = 0.0, sxz = 0.0;
  double syx = 0.0, syy = 0.0, syz = 0.0;
  double szx = 0.0, szy = 0.0, szz = 0.0, si = 0.0;
  float dataCentroid[3];  //centroid in data coordinates
  vtkFloatingPointType *spacing = input->GetSpacing();
  vtkFloatingPointType *origin = input->GetOrigin();

  dataCentroid[0] = (centroid[0] - origin[0])/spacing[0];
  dataCentroid[1] = (centroid[1] - origin[1])/spacing[1];
  dataCentroid[2] = (centroid[2] - origin[2])/spacing[2];
  
  input->GetIncrements(inInc0, inInc1, inInc2);

  inPtr2 = inPtr;
  for (idx2 = inputExtent[4]; idx2 <= inputExtent[5]; ++idx2)
    {
      inPtr1 = inPtr2;
      for (idx1 = inputExtent[2]; idx1 <= inputExtent[3]; ++idx1)
	{
	  inPtr0 = inPtr1;
	  for (idx0 = inputExtent[0]; idx0 <= inputExtent[1]; ++idx0)
	    {
	      sxx += (idx0-centroid[0]) * (idx0-centroid[0]) * *inPtr0;
	      sxy += (idx0-centroid[0]) * (idx1-centroid[1]) * *inPtr0;
	      sxz += (idx0-centroid[0]) * (idx2-centroid[2]) * *inPtr0;

	      syx += (idx1-centroid[1]) * (idx0-centroid[0]) * *inPtr0;
	      syy += (idx1-centroid[1]) * (idx1-centroid[1]) * *inPtr0;
	     