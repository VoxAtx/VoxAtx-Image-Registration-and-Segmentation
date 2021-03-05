
/*=========================================================================

  Module: vtkImageNeighborhoodCorrelation.cxx

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageNeighborhoodCorrelation.h"

#include "vtkImageSimilarityMetricInternals.h"

#include <vtkObjectFactory.h>
#include <vtkImageData.h>
#include <vtkImageStencilData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkTemplateAliasMacro.h>
#include <vtkVersion.h>

// turn off 64-bit ints when templating over all types
# undef VTK_USE_INT64
# define VTK_USE_INT64 0
# undef VTK_USE_UINT64
# define VTK_USE_UINT64 0

#include <math.h>

vtkStandardNewMacro(vtkImageNeighborhoodCorrelation);

//----------------------------------------------------------------------------
// Data needed for each thread.
class vtkImageNeighborhoodCorrelationThreadData
{
public:
  vtkImageNeighborhoodCorrelationThreadData() : Result(0.0), Count(0) {}

  double Result;
  vtkIdType Count;
};

class vtkImageNeighborhoodCorrelationTLS
  : public vtkImageSimilarityMetricTLS<vtkImageNeighborhoodCorrelationThreadData>
{
};

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageNeighborhoodCorrelation::vtkImageNeighborhoodCorrelation()
{
  this->NeighborhoodRadius[0] = 7;
  this->NeighborhoodRadius[1] = 7;
  this->NeighborhoodRadius[2] = 7;
}

//----------------------------------------------------------------------------
vtkImageNeighborhoodCorrelation::~vtkImageNeighborhoodCorrelation()
{
}

//----------------------------------------------------------------------------
void vtkImageNeighborhoodCorrelation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "NeighborhoodRadius: " << this->NeighborhoodRadius[0] << " "
     << this->NeighborhoodRadius[1] << " " << this->NeighborhoodRadius[2] << "\n";
}

// begin anonymous namespace
namespace {

//----------------------------------------------------------------------------