
/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMorphologicalInterpolator.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkMorphologicalInterpolator.h"
#include "vtkImageInterpolatorInternals.h"
#include "vtkImageData.h"
#include "vtkDataArray.h"
#include "vtkObjectFactory.h"

#include "vtkTemplateAliasMacro.h"
// turn off 64-bit ints when templating over all types, because
// they cannot be faithfully represented by doubles
# undef VTK_USE_INT64
# define VTK_USE_INT64 0
# undef VTK_USE_UINT64
# define VTK_USE_UINT64 0

vtkStandardNewMacro(vtkMorphologicalInterpolator);

//----------------------------------------------------------------------------
vtkMorphologicalInterpolator::vtkMorphologicalInterpolator()
{
  this->Operation = VTK_IMIOPERATION_DILATE;
  this->Radius[0] = 0.5;
  this->Radius[1] = 0.5;
  this->Radius[2] = 0.5;
}

//----------------------------------------------------------------------------
vtkMorphologicalInterpolator::~vtkMorphologicalInterpolator()
{
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "Operation: "
     << this->GetOperationAsString() << "\n";
  os << indent << "Radius: " << this->Radius[0] << " "
     << this->Radius[1] << " " << this->Radius[2] << "\n";
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::ComputeSupportSize(
  const double vtkNotUsed(matrix)[16], int size[3])
{
  this->InternalUpdate();

  size[0] = 1 + 2*static_cast<int>(this->InternalRadius[0]);
  size[1] = 1 + 2*static_cast<int>(this->InternalRadius[1]);
  size[2] = 1 + 2*static_cast<int>(this->InternalRadius[2]);
}

//----------------------------------------------------------------------------
bool vtkMorphologicalInterpolator::IsSeparable()
{
  return true;
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::SetOperation(int mode)
{
  static int minmode = VTK_IMIOPERATION_DILATE;
  static int maxmode = VTK_IMIOPERATION_ERODE;
  mode = ((mode > minmode) ? mode : minmode);
  mode = ((mode < maxmode) ? mode : maxmode);
  if (this->Operation != mode)
    {
    this->Operation = mode;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
const char *vtkMorphologicalInterpolator::GetOperationAsString()
{
  const char *result = "";

  switch (this->Operation)
    {
    case VTK_IMIOPERATION_DILATE:
      result = "Dilate";
      break;
    case VTK_IMIOPERATION_ERODE:
      result = "Erode";
      break;
    }

  return result;
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::SetRadius(
  double x, double y, double z)
{
  if (this->Radius[0] != x ||
      this->Radius[1] != y ||
      this->Radius[2] != z)
    {
    this->Radius[0] = x;
    this->Radius[1] = y;
    this->Radius[2] = z;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::ComputeInternalRadius(
  const double radius[3])
{
  const double rmin = 1e-17;
  const double rmax = (VTK_IMI_KERNEL_SIZE_MAX-1.0)*0.5;

  // clamp the radius to reasonable values
  for (int i = 0; i < 3; i++)
    {
    double r = radius[i];
    // clamp the radius to a reasonable range
    r = (r >= 0 ? r : 0);
    r = (r <= rmax ? r : rmax);
    this->InternalRadius[i] = r;
    // compute the inverse of the radius
    r = (r >= rmin ? r : rmin); 
    this->InternalRadius[i+3] = 1.0/r;
    }
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::InternalDeepCopy(
  vtkAbstractImageInterpolator *a)
{
  vtkMorphologicalInterpolator *obj =
    vtkMorphologicalInterpolator::SafeDownCast(a);
  if (obj)
    {
    this->SetOperation(obj->Operation);
    this->SetRadius(obj->Radius);
    for (int i = 0; i < 6; i++) // InternalRadius has six elements
      {
      this->InternalRadius[i] = obj->InternalRadius[i];
      }
    }
}

//----------------------------------------------------------------------------
void vtkMorphologicalInterpolator::InternalUpdate()
{
  this->ComputeInternalRadius(this->Radius);
  this->InterpolationInfo->InterpolationMode = this->Operation;
  this->InterpolationInfo->ExtraInfo = this->InternalRadius;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//  Interpolation subroutines and associated code
//----------------------------------------------------------------------------

namespace {

//----------------------------------------------------------------------------
template<class F, class T>
struct vtkImageMorphInterpolate
{
  static void General(
    vtkInterpolationInfo *info, const F point[3], F *outPtr);
};

//----------------------------------------------------------------------------
template <class F, class T>
void vtkImageMorphInterpolate<F, T>::General(
  vtkInterpolationInfo *info, const F point[3], F *outPtr)
{
  const T *inPtr = static_cast<const T *>(info->Pointer);
  int *inExt = info->Extent;
  vtkIdType *inInc = info->Increments;
  int numscalars = info->NumberOfComponents;
  int operation = info->InterpolationMode;
  double *radius = static_cast<double *>(info->ExtraInfo);
  double *invRadius = &radius[3];

  // find the closest point and offset from closest point
  double fx, fy, fz;
  int inIdX0 = vtkInterpolationMath::Floor(point[0] + 0.5, fx);
  int inIdY0 = vtkInterpolationMath::Floor(point[1] + 0.5, fy);
  int inIdZ0 = vtkInterpolationMath::Floor(point[2] + 0.5, fz);
  fx -= 0.5;
  fy -= 0.5;
  fz -= 0.5;

  // change arrays into locals
  vtkIdType inIncX = inInc[0];
  vtkIdType inIncY = inInc[1];
  vtkIdType inIncZ = inInc[2];

  int minX = inExt[0];
  int maxX = inExt[1];
  int minY = inExt[2];
  int maxY = inExt[3];
  int minZ = inExt[4];
  int maxZ = inExt[5];

  // arrays for the memory offsets
  vtkIdType factX[VTK_IMI_KERNEL_SIZE_MAX];
  vtkIdType factY[VTK_IMI_KERNEL_SIZE_MAX];
  vtkIdType factZ[VTK_IMI_KERNEL_SIZE_MAX];

  // handle the borders, compute the block size
  int rx = static_cast<int>(radius[0]);
  int ry = static_cast<int>(radius[1]);
  int rz = static_cast<int>(radius[2]);
  int xi = inIdX0 - rx;
  int xm = 2*rx+1;
  int yi = inIdY0 - ry;
  int ym = 2*ry+1;
  int zi = inIdZ0 - rz;
  int zm = 2*rz+1;
  int mm = xm;
  mm = ((mm >= ym) ? mm : ym);
  mm = ((mm >= zm) ? mm : zm);

  switch (info->BorderMode)
    {
    case VTK_IMAGE_BORDER_REPEAT:
      {
      int l = 0;
      int m = mm;
      do
        {
        factX[l] = vtkInterpolationMath::Wrap(xi, minX, maxX)*inIncX;
        factY[l] = vtkInterpolationMath::Wrap(yi, minY, maxY)*inIncY;
        factZ[l] = vtkInterpolationMath::Wrap(zi, minZ, maxZ)*inIncZ;
        l++; xi++; yi++; zi++;
        }
      while (--m);
      }
      break;

    case VTK_IMAGE_BORDER_MIRROR:
      {
      int l = 0;
      int m = mm;
      do
        {
        factX[l] = vtkInterpolationMath::Mirror(xi, minX, maxX)*inIncX;
        factY[l] = vtkInterpolationMath::Mirror(yi, minY, maxY)*inIncY;
        factZ[l] = vtkInterpolationMath::Mirror(zi, minZ, maxZ)*inIncZ;
        l++; xi++; yi++; zi++;
        }
      while (--m);
      }
      break;

    default:
      {
      int l = 0;
      int m = mm;
      do
        {
        factX[l] = vtkInterpolationMath::Clamp(xi, minX, maxX)*inIncX;
        factY[l] = vtkInterpolationMath::Clamp(yi, minY, maxY)*inIncY;
        factZ[l] = vtkInterpolationMath::Clamp(zi, minZ, maxZ)*inIncZ;
        l++; xi++; yi++; zi++;
        }
      while (--m);
      }
      break;
    }

  // the squared distances to each pixel
  F fX[VTK_IMI_KERNEL_SIZE_MAX];
  F fY[VTK_IMI_KERNEL_SIZE_MAX];
  F fZ[VTK_IMI_KERNEL_SIZE_MAX];

  double x = -rx - fx;
  double y = -ry - fy;
  double z = -rz - fz;
  double d;
  int ll = 0;
  do
    {
    d = fabs(x) - 0.5; d = (d < 0 ? 0 : d); d *= invRadius[0]; fX[ll] = d*d;
    d = fabs(y) - 0.5; d = (d < 0 ? 0 : d); d *= invRadius[1]; fY[ll] = d*d;
    d = fabs(z) - 0.5; d = (d < 0 ? 0 : d); d *= invRadius[2]; fZ[ll] = d*d;
    ll++; x++; y++; z++;
    }
  while (--mm);

  // check if only one slice in a particular direction
  int multipleY = (minY != maxY);
  int multipleZ = (minZ != maxZ);

  // the limits to use when doing the interpolation
  int k1 = rz - rz*multipleZ;
  int k2 = rz + rz*multipleZ;
  int j1 = ry - ry*multipleY;
  int j2 = ry + ry*multipleY;

  do // loop over components
    {
    F val = inPtr[factZ[rz] + factY[ry] + factX[rx]];
    int k = k1;
    do // loop over z
      {
      F ifz = fZ[k] - (1.0 + VTK_INTERPOLATE_FLOOR_TOL);
      vtkIdType factz = factZ[k];
      int j = j1;
      do // loop over y
        {
        F ify = fY[j];
        F fzy = ifz + ify;
        vtkIdType factzy = factz + factY[j];
        // loop over x
        const T *tmpPtr = inPtr + factzy;
        const F *tmpfX = fX;
        const vtkIdType *tmpfactX = factX;
        int l = xm;