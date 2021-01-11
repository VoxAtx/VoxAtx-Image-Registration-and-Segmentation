/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkITKXFMWriter.cxx

Copyright (c) 2006 Atamai, Inc.

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

#include "vtkITKXFMWriter.h"

#include "vtkObjectFactory.h"

#include "vtkMath.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkCollection.h"
#include "vtkTransform.h"
#include "vtkHomogeneousTransform.h"
#include "vtkGeneralTransform.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPoints.h"

#include <ctype.h>
#include <stdio.h>

#if !defined(_WIN32) || defined(__CYGWIN__)
# include <unistd.h> /* unlink */
#else
# include <io.h> /* unlink */
#endif

#include <stack>

//--------------------------------------------------------------------------
vtkStandardNewMacro(vtkITKXFMWriter);

//-------------------------------------------------------------------------
vtkITKXFMWriter::vtkITKXFMWriter()
{
  this->FileName = 0;
  this->Transform = 0;
  this->TransformCenter[0] = 0.0;
  this->TransformCenter[1] = 0.0;
  this->TransformCenter[2] = 0.0;
  this->Transforms = vtkCollection::New();
}

//-------------------------------------------------------------------------
vtkITKXFMWriter::~vtkITKXFMWriter()
{
  if (this->Transforms)
    {
    this->Transforms->Delete();
    }
  if (this->Transform)
    {
    this->Transform->Delete();
    }
  if (this->FileName)
    {
    delete [] this->FileName;
    }
}

//-------------------------------------------------------------------------
void vtkITKXFMWriter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "FileName: "
     << (this->FileName ? this->FileName : "none") << "\n";
  os << indent << "Transform: " << this->Transform << "\n";
  if (this->Transform)
    {
    this->Transform->PrintSelf(os, indent.GetNextIndent());
    }
  os << indent << "TransformCenter: "
     << this-