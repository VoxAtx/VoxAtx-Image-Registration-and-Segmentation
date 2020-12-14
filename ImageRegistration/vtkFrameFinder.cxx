
/*=========================================================================

  Program:   Atamai Classes for VTK
  Module:    vtkFrameFinder.cxx

Copyright (c) 2005 Atamai, Inc.
All rights reserved.

Use, modification and redistribution of the software, in source or
binary forms, are permitted provided that the following terms and
conditions are met:

1) Redistribution of the source code, in verbatim or modified
   form, must retain the above copyright notice, this license,
   the following disclaimer, and any notices that refer to this
   license and/or the following disclaimer.

2) Redistribution in binary form must include the above copyright
   notice, a copy of this license and the following disclaimer
   in the documentation or with other materials provided with the
   distribution.

3) Modified copies of the source code must be clearly marked as such,
   and must not be misrepresented as verbatim copies of the source code.

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

#include "vtkFrameFinder.h"

// VTK header files
#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkMatrix4x4.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkVersion.h"
#include "vtkMath.h"
#include "vtkSmartPointer.h"

#include "vtkTemplateAliasMacro.h"
// turn off 64-bit ints when templating over all types
# undef VTK_USE_INT64
# define VTK_USE_INT64 0
# undef VTK_USE_UINT64
# define VTK_USE_UINT64 0

#include <math.h>
#include <string.h>

#include <stack>
#include <map>
#include <vector>
#include <utility>
#include <limits>

vtkStandardNewMacro(vtkFrameFinder);
vtkCxxSetObjectMacro(vtkFrameFinder,DICOMPatientMatrix,vtkMatrix4x4);

//----------------------------------------------------------------------------
vtkFrameFinder::vtkFrameFinder()
{
  this->ImageToFrameMatrix = vtkMatrix4x4::New();
  this->DICOMPatientMatrix = 0;
  this->Success = 0;
  this->UseAnteriorFiducial = 1;
  this->UsePosteriorFiducial = 1;

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(2);
}

//----------------------------------------------------------------------------
vtkFrameFinder::~vtkFrameFinder()
{
  if (this->DICOMPatientMatrix)
    {
    this->DICOMPatientMatrix->Delete();
    }
  if (this->ImageToFrameMatrix)
    {
    this->ImageToFrameMatrix->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkFrameFinder::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "UseAnteriorFiducial: "
     << (this->UseAnteriorFiducial ? "On\n" : "Off\n");

  os << indent << "UsePosteriorFiducial: "
     << (this->UsePosteriorFiducial ? "On\n" : "Off\n");

  os << indent << "ImageToFrameMatrix: "
     << this->ImageToFrameMatrix << "\n";

  os << indent << "DICOMPatientMatrix: "
     << this->DICOMPatientMatrix << "\n";

  os << indent << "Success: " << this->Success << "\n";
}

//----------------------------------------------------------------------------
int vtkFrameFinder::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkFrameFinder::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
  return 1;
}

//----------------------------------------------------------------------------
void vtkFrameFinder::SetInputData(vtkDataObject* input)
{
#if VTK_MAJOR_VERSION <= 5
  if (input)
    {
    this->SetInputConnection(0, input->GetProducerPort());
    }
  else