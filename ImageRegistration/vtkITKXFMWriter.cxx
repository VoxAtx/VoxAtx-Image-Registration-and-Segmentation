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
     << this->TransformCenter[0] << " "
     << this->TransformCenter[1] << " "
     << this->TransformCenter[2] << "\n";
  os << indent << "NumberOfTransforms: "
     << this->Transforms->GetNumberOfItems() << "\n";
}

//-------------------------------------------------------------------------
int vtkITKXFMWriter::WriteLinearTransform(
  ostream &outfile, vtkHomogeneousTransform *transform)
{
  vtkMatrix4x4 *matrix = transform->GetMatrix();

  if (matrix->GetElement(3,0) != 0.0 ||
      matrix->GetElement(3,1) != 0.0 ||
      matrix->GetElement(3,2) != 0.0 ||
      matrix->GetElement(3,3) != 1.0)
    {
    vtkErrorMacro("WriteLinearTransform: The transform is not linear");
    return 0;
    }

  double c[4] = { 0.0, 0.0, 0.0, 1.0 };
  this->GetTransformCenter(c);
  double t[4];
  matrix->MultiplyPoint(c, t);
  t[0] -= c[0];
  t[1] -= c[1];
  t[2] -= c[2];

  double p[12];
  for (int i = 0; i < 3; i++)
    {
    for (int j = 0; j < 3; j++)
      {
      p[3*i + j] =  matrix->GetElement(i, j);
      }
    p[9 + i] = t[i];
    }

  if (vtkITKXFMWriter::IsMatFile(this->FileName))
    {
    // write the transform as a Matlab Level 4 file, little endian
    char head[20] = {
      0x00,0x00,0x00,0x00, // type is IEEE little endian
      0x0C,0x00,0x00,0x00, // 12 rows
      0x01,0x00,0x00,0x00, // 1 column
      0x00,0x00,0x00,0x00, // double-precision
      0x1b,0x00,0x00,0x00  // array name is 27 bytes (with null)
    };

    for (int j = 0; j < 2; j++)
      {
      // first iteration writes the parameters,
      // second iteration writes the fixed parameters
      double *dp = (j == 0 ? p : c);
      const char *name = (j == 0 ? "AffineTransform_double_3_3" : "fixed");
      int n = (j == 0 ? 12 : 3);
      int m = static_cast<int>(strlen(name) + 1);
      head[4] = n;
      head[16] = m;

      outfile.write(head, 20);
      outfile.write(name, m);

      for (int k = 0; k < n; k++)
        {
        // write little-endian double-precision
        union { double d; unsigned long long l; } u;
        u.d = dp[k];
        unsigned long long l = u.l;
        char op[8];
        op[0] = static_cast<unsigned char>(l);
        op[1] = static_cast<unsigned char>(l >> 8);
        op[2] = static_cast<unsigned char>(l >> 16);
        op[3] = static_cast<unsigned char>(l >> 24);
        l >>= 32;
        op[4] = static_cast<unsigned char>(l);
        op[5] = static_cast<unsigned char>(l >> 8);
        op[6] = static_cast<unsigned char>(l >> 16);
        op[7] = static_cast<unsigned char>(l >> 24);
        outfile.write(op, 8);
        }
      }
    }
  else
    {
    // write the transform as text
    outfile << "Transform: MatrixOffsetTransformBase_double_3_3\n";

    outfile << "Parameters:";

    outfile.precision(15);

    for (int k = 0; k < 12; k++)
      {
      outfile << " " << p[k];
      }
    outfile << "\n";

    outfile << "FixedParameters:";
    outfile << " " << c[0] << " " << c[1] << " " << c[2];

    outfile << "\n";
    }

  return 1;
}

//-------------------------------------------------------------------------
int vtkITKXFMWriter::WriteTransform(
  ostream &outfile, vtkAbstractTransform *transform)
{
  if (transform->IsA("vtkHomogeneousTransform"))
    {
    return this->WriteLinearTransform(
      outfile, (vtkHomogeneousTransform *)transform);
    }

  vtkErrorMacro("Unsupported transform type "
                << transform->GetClassName());

  return 0;
}

//-------------------------------------------------------------------------
int vtkITKXFMWriter::WriteFile()
{
  // Check that a transform has been set.
  if (!this->Transform)
    {
    vtkErrorMacro("WriteFile: No input transform has been set.");
    return 0;
    }
  // Check that the file name has been set.
  if (!this->FileName)
    {
    vtkErrorMacro("WriteFile: No file name has been set.");
    return 0;
    }

  // Is this a matlab file?
  bool isMat = vtkITKXFMWriter::IsMatFile(this->FileName);

  // Open the file.
  ofstr