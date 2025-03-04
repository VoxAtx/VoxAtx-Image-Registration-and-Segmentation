
/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGaussianInterpolator.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkGaussianInterpolator.h"
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

// masks for storing window and size in a single integer
#define VTK_INTERPOLATION_WINDOW_MASK        0x0000007f
#define VTK_INTERPOLATION_WINDOW_XBLUR_MASK  0x00008000
#define VTK_INTERPOLATION_WINDOW_XSIZE_MASK  0x00007f00
#define VTK_INTERPOLATION_WINDOW_XSIZE_SHIFT 8
#define VTK_INTERPOLATION_WINDOW_YBLUR_MASK  0x00800000
#define VTK_INTERPOLATION_WINDOW_YSIZE_MASK  0x007f0000
#define VTK_INTERPOLATION_WINDOW_YSIZE_SHIFT 16
#define VTK_INTERPOLATION_WINDOW_ZBLUR_MASK  0x80000000
#define VTK_INTERPOLATION_WINDOW_ZSIZE_MASK  0x7f000000
#define VTK_INTERPOLATION_WINDOW_ZSIZE_SHIFT 24

// kernel lookup table size must be 256*n where n is kernel half-width
// in order to provide sufficient precision for 16-bit images
#define VTK_GAUSS_KERNEL_TABLE_DIVISIONS 256

vtkStandardNewMacro(vtkGaussianInterpolator);

//----------------------------------------------------------------------------
vtkGaussianInterpolator::vtkGaussianInterpolator()
{
  this->KernelType = VTK_GAUSSIAN_INTERPOLATION;
  this->RadiusFactors[0] = 3;
  this->RadiusFactors[1] = 3;
  this->RadiusFactors[2] = 3;
  this->KernelLookupTable[0] = NULL;
  this->KernelLookupTable[1] = NULL;
  this->KernelLookupTable[2] = NULL;
  this->KernelSize[0] = 6;
  this->KernelSize[1] = 6;
  this->KernelSize[2] = 6;
  this->Antialiasing = 0;
  this->Renormalization = 0;
  this->BlurFactors[0] = 1.0;
  this->BlurFactors[1] = 1.0;
  this->BlurFactors[2] = 1.0;
  this->LastBlurFactors[0] = 1.0;
  this->LastBlurFactors[1] = 1.0;
  this->LastBlurFactors[2] = 1.0;
}

//----------------------------------------------------------------------------
vtkGaussianInterpolator::~vtkGaussianInterpolator()
{
  if (this->KernelLookupTable[0])
    {
    this->FreeKernelLookupTable();
    }
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "KernelType: " << this->GetKernelTypeAsString() << "\n";
  os << indent << "RadiusFactors: " << this->RadiusFactors[0] << " "
     << this->RadiusFactors[1] << " " << this->RadiusFactors[2] << "\n";
  os << indent << "BlurFactors: " << this->BlurFactors[0] << " "
     << this->BlurFactors[1] << " " << this->BlurFactors[2] << "\n";
  os << indent << "Antialiasing: "
     << (this->Antialiasing ? "On\n" : "Off\n");
  os << indent << "Renormalization: "
     << (this->Renormalization ? "On\n" : "Off\n");
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::ComputeSupportSize(
  const double matrix[16], int size[3])
{
  // compute the default support size for when matrix is null
  if (this->Antialiasing)
    {
    size[0] = VTK_GAUSS_KERNEL_SIZE_MAX;
    size[1] = VTK_GAUSS_KERNEL_SIZE_MAX;
    size[2] = VTK_GAUSS_KERNEL_SIZE_MAX;
    }
  else
    {
    for (int i = 0; i < 3; i++)
      {
      // use blur factors to compute support size
      size[i] = 2*static_cast<int>(
        this->RadiusFactors[i] + 1.0 - VTK_INTERPOLATE_FLOOR_TOL);
      double rowscale = this->BlurFactors[i];
      if (rowscale > (1.0 + VTK_INTERPOLATE_FLOOR_TOL))
        {
        size[i] = 2*static_cast<int>(
          rowscale*this->RadiusFactors[i] + 1.0 - VTK_INTERPOLATE_FLOOR_TOL);
        }
      }
    }

  if (matrix == NULL)
    {
    return;
    }

  if (this->Antialiasing)
    {
    // if antialiasing is on, initialize blur factors to 1
    for (int i = 0; i < 3; i++)
      {
      this->BlurFactors[i] = 1.0;
      this->KernelSize[i] = 2*static_cast<int>(
        this->RadiusFactors[i] + 1.0 - VTK_INTERPOLATE_FLOOR_TOL);
      }
    }
  else
    {
    // keep blur factors, use kernel size computed from blur factors
    this->KernelSize[0] = size[0];
    this->KernelSize[1] = size[1];
    this->KernelSize[2] = size[2];
    }

  // if matrix does perspective, use the defaults just computed
  if (matrix[12] != 0 || matrix[13] != 0 || matrix[14] != 0 ||
      matrix[15] != 1.0)
    {
    return;
    }

  // use matrix to compute blur factors and kernel size
  for (int i = 0; i < 3; i++)
    {
    double rowscale = 0.0;
    for (int j = 0; j < 3; j++)
      {
      // compute the scale from a row of the matrix
      double x = matrix[4*i + j];
      rowscale += x*x;

      // verify that the element is an integer:
      // check fraction that remains after floor operation
      double f;
      vtkInterpolationMath::Floor(x, f);
      }

    if (this->Antialiasing)
      {
      // rowscale is the subsampling factor in a particular direction
      rowscale = sqrt(rowscale);
      }
    else
      {
      // ignore computed value, use factor provided by SetBlurFactors()
      rowscale = this->BlurFactors[i];
      }

    // if scale is greater than one, expand kernel size
    if (rowscale > (1.0 + VTK_INTERPOLATE_FLOOR_TOL))
      {
      // need extra suport for antialiasing
      this->BlurFactors[i] = rowscale;
      int s = 2*static_cast<int>(
        rowscale*this->RadiusFactors[i] + 1.0 - VTK_INTERPOLATE_FLOOR_TOL);
      size[i] = s;
      this->KernelSize[i] = s;
      }
    }

  // rebuild the kernel lookup tables
  this->InternalUpdate();
}

//----------------------------------------------------------------------------
bool vtkGaussianInterpolator::IsSeparable()
{
  return true;
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::SetRadiusFactors(
  double x, double y, double z)
{
  if (this->RadiusFactors[0] != x ||
      this->RadiusFactors[1] != y ||
      this->RadiusFactors[2] != z)
    {
    this->RadiusFactors[0] = x;
    this->RadiusFactors[1] = y;
    this->RadiusFactors[2] = z;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::SetKernelType(int ktype)
{
  static int mintype = VTK_GAUSSIAN_INTERPOLATION;
  static int maxtype = VTK_APPLEDORN10_INTERPOLATION;
  ktype = ((ktype > mintype) ? ktype : mintype);
  ktype = ((ktype < maxtype) ? ktype : maxtype);
  if (this->KernelType != ktype)
    {
    this->KernelType = ktype;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
const char *vtkGaussianInterpolator::GetKernelTypeAsString()
{
  const char *result = "";

  switch (this->KernelType)
    {
    case VTK_GAUSSIAN_INTERPOLATION:
      result = "Gaussian";
      break;
    case VTK_APPLEDORN2_INTERPOLATION:
      result = "Appledorn2";
      break;
    case VTK_APPLEDORN6_INTERPOLATION:
      result = "Appledorn6";
      break;
    case VTK_APPLEDORN10_INTERPOLATION:
      result = "Appledorn10";
      break;
    }

  return result;
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::SetBlurFactors(double x, double y, double z)
{
  if (this->BlurFactors[0] != x ||
      this->BlurFactors[1] != y ||
      this->BlurFactors[2] != z)
    {
    this->BlurFactors[0] = x;
    this->BlurFactors[1] = y;
    this->BlurFactors[2] = z;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::SetAntialiasing(int val)
{
  val = (val != 0);
  if (this->Antialiasing != val)
    {
    this->Antialiasing = val;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::SetRenormalization(int val)
{
  val = (val != 0);
  if (this->Renormalization != val)
    {
    this->Renormalization = val;
    this->Modified();
    }
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::InternalDeepCopy(
  vtkAbstractImageInterpolator *a)
{
  vtkGaussianInterpolator *obj = vtkGaussianInterpolator::SafeDownCast(a);
  if (obj)
    {
    this->SetKernelType(obj->KernelType);
    this->SetRadiusFactors(obj->RadiusFactors);
    this->SetAntialiasing(obj->Antialiasing);
    if (this->Antialiasing)
      {
      // if blur factors were computed, then don't call "modified"
      obj->GetBlurFactors(this->BlurFactors);
      }
    else
      {
      this->SetBlurFactors(obj->BlurFactors);
      }
    }

  this->KernelSize[0] = 6;
  this->KernelSize[1] = 6;
  this->KernelSize[2] = 6;

  if (this->KernelLookupTable[0])
    {
    this->FreeKernelLookupTable();
    }
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::InternalUpdate()
{
  bool blurchange = false;
  int mode = this->KernelType;
  int hsize[3];
  for (int i = 0; i < 3; i++)
    {
    static int minsize = 1;
    static int maxsize = VTK_GAUSS_KERNEL_SIZE_MAX/2;
    int size = this->KernelSize[i]/2;
    size = ((size > minsize) ? size : minsize);
    size = ((size < maxsize) ? size : maxsize);
    hsize[i] = size;
    blurchange |= (fabs(this->BlurFactors[i] - this->LastBlurFactors[i]) >=
                   VTK_INTERPOLATE_FLOOR_TOL);
    }

  if (this->BlurFactors[0] > 1.0 + VTK_INTERPOLATE_FLOOR_TOL)
    {
    mode |= VTK_INTERPOLATION_WINDOW_XBLUR_MASK;
    }
  if (this->BlurFactors[1] > 1.0 + VTK_INTERPOLATE_FLOOR_TOL)
    {
    mode |= VTK_INTERPOLATION_WINDOW_YBLUR_MASK;
    }
  if (this->BlurFactors[2] > 1.0 + VTK_INTERPOLATE_FLOOR_TOL)
    {
    mode |= VTK_INTERPOLATION_WINDOW_ZBLUR_MASK;
    }

  mode |= (hsize[0] << VTK_INTERPOLATION_WINDOW_XSIZE_SHIFT);
  mode |= (hsize[1] << VTK_INTERPOLATION_WINDOW_YSIZE_SHIFT);
  mode |= (hsize[2] << VTK_INTERPOLATION_WINDOW_ZSIZE_SHIFT);

  if (this->InterpolationInfo->InterpolationMode != mode ||
      blurchange ||
      this->KernelLookupTable[0] == NULL)
    {
    this->BuildKernelLookupTable();
    }

  this->InterpolationInfo->InterpolationMode = mode;
  this->InterpolationInfo->ExtraInfo = this->KernelLookupTable;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//  Interpolation subroutines and associated code
//----------------------------------------------------------------------------

namespace {

// Equations were originally taken from the following paper:
//  C.R. Appledorn, "A New Approach to the Interpolation of Sampled Data."
//  IEEE Transactions on Medical Imaging, 15(3), 1996.
// I substituted f = 1/sqrt(v) to simplify some of the math (v = variance)
// and did other minor rearrangement and refactoring, followed by a
// numerical comparison of the results of this code against an exact
// transcription of Appledorn's equations into a Python program.

struct vtkAppledorn
{
  static double G0(double x, double f);
  static double G2(double x, double f);
  static double G6(double x, double f);
  static double G10(double x, double f);

  static double D2(double x);
  static double D6(double x);
  static double D10(double x);
};

inline double vtkAppledorn::G0(double x, double f)
{
  // 1/sqrt(2*pi)*f*exp(1/2*x^2*f^2)
  return 0.3989422804014327*f*exp(-0.5*(x*x)*(f*f));
}

inline double vtkAppledorn::G2(double x, double f)
{
  // f^4*(x^2 - 1/f^2)
  double f2 = f*f;
  return ((f2*f2)*(x*x) - f2);
}

inline double vtkAppledorn::G6(double x, double f)
{
  // f^12*(x^6 - 15*x^4/f^2 + 45*x^2/f^4 - 15/f^6)
  double x2 = x*x;
  double x4 = x2*x2;
  double f2 = f*f;
  double f4 = f2*f2;
  double f6 = f4*f2;
  return (f6*(x4*x2) - 15.0*f4*x4 + 45.0*f2*x2 - 15.0)*f6;
}

inline double vtkAppledorn::G10(double x, double f)
{
  // f^20*(x^10 - 45*x^8/f^2 + 630*x^6/f^4 - 3150*x^4/f^6
  //       + 4725*x^2/f^8 - 945/f^10)
  double x2 = x*x;
  double x4 = x2*x2;
  double x6 = x4*x2;
  double f2 = f*f;
  double f4 = f2*f2;
  double f6 = f4*f2;
  double f10 = f6*f4;
  return (f10*(x6*x4) - 45.0*(f4*f4)*(x4*x4) + 630.0*f6*x6 - 3150.0*f4*x4
          + 4725.0*f2*x2 - 945.0)*f10;
}

inline double vtkAppledorn::D2(double x)
{
  // c1 = sqrt(pi)/(1/sqrt(2) + 1)
  // c2 = sqrt(2)*c1
  // k1 = 1/pi*(1/sqrt(2) + 1)
  // k2 = 1/2*k1
  const double c1 = 1.0382794271800315;
  const double c2 = 1.4683488474509689;
  const double k2 = 0.46381149367711955;

  return vtkAppledorn::G0(x,c1) -
         vtkAppledorn::G0(x,c2)*k2*vtkAppledorn::G2(x,c2);
}

inline double vtkAppledorn::D6(double x)
{
  // c1 = sqrt(pi)/(1/sqrt(2) + 1 + 15/24)
  // c2 = sqrt(2)*c1
  // k1 = 1/pi*(1/sqrt(2) + 1)
  // k2 = 1/2*k1
  // k6 = 1/192*k1^3
  const double c1 = 0.76002259639402658;
  const double c2 = 1.0748342635304455;
  const double k2 = 0.86559949658680813;
  const double k6 = 0.027023384702061636;

  return vtkAppledorn::G0(x,c1) -
         vtkAppledorn::G0(x,c2)*(k2*vtkAppledorn::G2(x,c2) +
                                 k6*vtkAppledorn::G6(x,c2));
}

inline double vtkAppledorn::D10(double x)
{
  // c1 = sqrt(pi)/(1/sqrt(2) + 1 + 15/24 + 945/1920)
  // c2 = sqrt(2)*c1
  // k1 = 1/pi*(1/sqrt(2) + 1)
  // k2 = 1/2*k1
  // k6 = 1/192*k1^3
  // k10 = 1/61440*k1^5
  const double c1 = 0.62757406786975101;
  const double c2 = 0.88752375817505513;
  const double k2 = 1.2695213966757273;
  const double k6 = 0.085252835612840125;
  const double k10 = 0.0017175085033765082;

  return vtkAppledorn::G0(x,c1) -
         vtkAppledorn::G0(x,c2)*(k2*vtkAppledorn::G2(x,c2) +
                                 k6*vtkAppledorn::G6(x,c2) +
                                 k10*vtkAppledorn::G10(x,c2));
}

//----------------------------------------------------------------------------
// Gaussian kernel computation: compute half of the interpolation kernel,
// including n sinc lobes, to make a lookup table of size "size".
// In the table, x=0.0 corresponds to index position zero, and
// x = 1.0 corresponds to index position "size", which is just
// beyond the end of the table and holds an implicit value of zero.

struct vtkGaussKernel
{
  template<class F>
  static void D0(F *kernel, int size, double p);
  template<class F>
  static void D2(F *kernel, int size, double p);
  template<class F>
  static void D6(F *kernel, int size, double p);
  template<class F>
  static void D10(F *kernel, int size, double p);
};

template<class F>
void vtkGaussKernel::D10(F *kernel, int size, double p)
{
  double x = 0.0;
  do
    {
    *kernel++ = vtkAppledorn::D10(x);
    x += p;
    }
  while (--size);
}

template<class F>
void vtkGaussKernel::D6(F *kernel, int size, double p)
{
  double x = 0.0;
  do
    {
    *kernel++ = vtkAppledorn::D6(x);
    x += p;
    }
  while (--size);
}

template<class F>
void vtkGaussKernel::D2(F *kernel, int size, double p)
{
  double x = 0.0;
  do
    {
    *kernel++ = vtkAppledorn::D2(x);
    x += p;
    }
  while (--size);
}

template<class F>
void vtkGaussKernel::D0(F *kernel, int size, double p)
{
  const double f = 2.5066282746310002; // sqrt(2*pi)
  double x = 0.0;
  do
    {
    *kernel++ = vtkAppledorn::G0(x, f);
    x += p;
    }
  while (--size);
}

//----------------------------------------------------------------------------
template<class T, class F>
void vtkGaussInterpWeights(T *kernel, F *fX, F fx, int m)
{
  // table bins per unit
  int p = VTK_GAUSS_KERNEL_TABLE_DIVISIONS;

  // compute table interpolation info
  F f = fx*p;
  int offset = static_cast<int>(f);
  f -= offset;
  F r = 1 - f;

  // interpolate the table
  int n = m;
  int i = (1 - (m >> 1))*p - offset;
  do
    {
    int i0 = i;
    int i1 = i + 1;
    int ni = -i0;
    i0 = ((i0 >= 0) ? i0 : ni);
    ni = -i1;
    i1 = ((i1 >= 0) ? i1 : ni);
    *fX++ = r*kernel[i0] + f*kernel[i1];
    i += p;
    }
  while (--n);
}

//----------------------------------------------------------------------------
// Ensure that the set of n coefficients extracted from the kernel table
// will always sum to unity.  This renormalization is needed to ensure that
// the interpolation will not have a DC offset.  For the rationale, see e.g.
//  NA Thacker, A Jackson, D Moriarty, E Vokurka, "Improved quality of
//  re-sliced MR images usng re-normalized sinc interpolation," Journal of
//  Magnetic Resonance Imaging 10:582-588, 1999.
// Parameters:
//  kernel = table containing half of a symmetric kernel
//  m = table offset between lookup positions for adjacent weights
//  n = number of kernel weights (i.e. size of discrete kernel size)
//  The kernel table size must be (n*m+1)/2
template<class T>
void vtkRenormalizeKernel(T *kernel, int m, int n)
{
  // the fact that we only have half the kernel makes the weight
  // lookup more complex: there will be kn direct lookups and km
  // mirror lookups.
  int kn = (n + 1)/2;
  int km = n - kn;

  if (m == 0 || km == 0)
    {
    return;
    }

  // get sum of weights for zero offset
  T w = - (*kernel)*0.5;
  T *kernel2 = kernel;
  int k = kn;
  do
    {
    w += *kernel2;
    kernel2 += m;
    }
  while (--k);

  // divide weights by their sum to renormalize
  w *= 2;
  kernel2 = kernel;
  k = kn;
  do
    {
    *kernel2 /= w;
    kernel2 += m;
    }
  while (--k);

  // need the opposite end of the kernel array
  kernel2 = kernel + km*m;

  int j = (m - 1)/2;
  if (j) do
    {
    // move to next offset
    kernel++;
    kernel2--;

    // get the sum of the weights at this offset
    w = 0;
    T *kernel1 = kernel2;
    k = km;
    do
      {
      w += *kernel1;
      kernel1 -= m;
      }
    while (--k);
    kernel1 = kernel;
    k = kn;
    do
      {
      w += *kernel1;
      kernel1 += m;
      }
    while (--k);

    // divide the weights by their sum to renormalize
    kernel1 = kernel2;
    k = km;
    do
      {
      *kernel1 /= w;
      kernel1 -= m;
      }
    while (--k);
    kernel1 = kernel;
    k = kn;
    do
      {
      *kernel1 /= w;
      kernel1 += m;
      }
    while (--k);
    }
  while (--j);

  // get sum of weights for offset of 0.5 (only applies when m is even)
  if ((m & 1) == 0)
    {
    w = 0;
    kernel++;
    kernel2 = kernel;
    k = km;
    do
      {
      w += *kernel2;
      kernel2 += m;
      }
    while (--k);

    // divide weights by their sum to renormalize
    w *= 2;
    kernel2 = kernel;
    k = km;
    do
      {
      *kernel2 /= w;
      kernel2 += m;
      }
    while (--k);
    }
}

//----------------------------------------------------------------------------
template<class F, class T>
struct vtkImageGaussInterpolate
{
  static void General(
    vtkInterpolationInfo *info, const F point[3], F *outPtr);
};

//----------------------------------------------------------------------------
template <class F, class T>
void vtkImageGaussInterpolate<F, T>::General(
  vtkInterpolationInfo *info, const F point[3], F *outPtr)
{
  const T *inPtr = static_cast<const T *>(info->Pointer);
  int *inExt = info->Extent;
  vtkIdType *inInc = info->Increments;
  int numscalars = info->NumberOfComponents;

  // kernel lookup table
  float **kernel = static_cast<float **>(info->ExtraInfo);

  // size of kernel
  int mode = info->InterpolationMode;
  int xm = 2*((mode & VTK_INTERPOLATION_WINDOW_XSIZE_MASK)
              >> VTK_INTERPOLATION_WINDOW_XSIZE_SHIFT);
  int ym = 2*((mode & VTK_INTERPOLATION_WINDOW_YSIZE_MASK)
              >> VTK_INTERPOLATION_WINDOW_YSIZE_SHIFT);
  int zm = 2*((mode & VTK_INTERPOLATION_WINDOW_ZSIZE_MASK)
              >> VTK_INTERPOLATION_WINDOW_ZSIZE_SHIFT);

  // index to kernel midpoint position
  int xm2 = ((xm - 1) >> 1);
  int ym2 = ((ym - 1) >> 1);
  int zm2 = ((zm - 1) >> 1);

  F fx, fy, fz;
  int inIdX0 = vtkInterpolationMath::Floor(point[0], fx);
  int inIdY0 = vtkInterpolationMath::Floor(point[1], fy);
  int inIdZ0 = vtkInterpolationMath::Floor(point[2], fz);

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

  // the memory offsets
  vtkIdType factX[VTK_GAUSS_KERNEL_SIZE_MAX];
  vtkIdType factY[VTK_GAUSS_KERNEL_SIZE_MAX];
  vtkIdType factZ[VTK_GAUSS_KERNEL_SIZE_MAX];

  // handle the borders
  int xi = inIdX0 - xm2;
  int yi = inIdY0 - ym2;
  int zi = inIdZ0 - zm2;
  int mm = xm;
  mm = ((mm >= ym) ? mm : ym);
  mm = ((mm >= zm) ? mm : zm);

  switch (info->BorderMode)
    {
    case VTK_IMAGE_BORDER_REPEAT:
      {
      int l = 0;
      do
        {
        factX[l] = vtkInterpolationMath::Wrap(xi, minX, maxX)*inIncX;
        factY[l] = vtkInterpolationMath::Wrap(yi, minY, maxY)*inIncY;
        factZ[l] = vtkInterpolationMath::Wrap(zi, minZ, maxZ)*inIncZ;
        l++; xi++; yi++; zi++;
        }
      while (--mm);
      }
      break;

    case VTK_IMAGE_BORDER_MIRROR:
      {
      int l = 0;
      do
        {
        factX[l] = vtkInterpolationMath::Mirror(xi, minX, maxX)*inIncX;
        factY[l] = vtkInterpolationMath::Mirror(yi, minY, maxY)*inIncY;
        factZ[l] = vtkInterpolationMath::Mirror(zi, minZ, maxZ)*inIncZ;
        l++; xi++; yi++; zi++;
        }
      while (--mm);
      }
      break;

    default:
      {
      int l = 0;
      do
        {
        factX[l] = vtkInterpolationMath::Clamp(xi, minX, maxX)*inIncX;
        factY[l] = vtkInterpolationMath::Clamp(yi, minY, maxY)*inIncY;
        factZ[l] = vtkInterpolationMath::Clamp(zi, minZ, maxZ)*inIncZ;
        l++; xi++; yi++; zi++;
        }
      while (--mm);
      }
      break;
    }

  // compute the kernel weights
  F fX[VTK_GAUSS_KERNEL_SIZE_MAX];
  F fY[VTK_GAUSS_KERNEL_SIZE_MAX];
  F fZ[VTK_GAUSS_KERNEL_SIZE_MAX];

  vtkGaussInterpWeights(kernel[0], fX, fx, xm);
  vtkGaussInterpWeights(kernel[1], fY, fy, ym);
  vtkGaussInterpWeights(kernel[2], fZ, fz, zm);

  // check if only one slice in a particular direction
  int multipleY = (minY != maxY);
  int multipleZ = (minZ != maxZ);

  // the limits to use when doing the interpolation
  int k1 = zm2*(1 - multipleZ);
  int k2 = (zm2 + 1)*(multipleZ + 1) - 1;
  int j1 = ym2*(1 - multipleY);
  int j2 = (ym2 + 1)*(multipleY + 1) - 1;

  do // loop over components
    {
    F val = 0;
    int k = k1;
    do // loop over z
      {
      F ifz = fZ[k];
      vtkIdType factz = factZ[k];
      int j = j1;
      do // loop over y
        {
        F ify = fY[j];
        F fzy = ifz*ify;
        vtkIdType factzy = factz + factY[j];
        // loop over x
        const T *tmpPtr = inPtr + factzy;
        const F *tmpfX = fX;
        const vtkIdType *tmpfactX = factX;
        F tmpval = 0;
        int l = xm;
        do
          {
          tmpval += (*tmpfX++)*tmpPtr[(*tmpfactX++)];
          }
        while (--l);
        val += fzy*tmpval;
        }
      while (++j <= j2);
      }
    while (++k <= k2);

    *outPtr++ = val;
    inPtr++;
    }
  while (--numscalars);
}

//----------------------------------------------------------------------------
// Get the interpolation function for the specified data types
template<class F>
void vtkGaussianInterpolatorGetInterpolationFunc(
  void (**interpolate)(vtkInterpolationInfo *, const F [3], F *),
  int dataType)
{
  switch (dataType)
    {
    vtkTemplateAliasMacro(
      *interpolate =
        &(vtkImageGaussInterpolate<F, VTK_TT>::General)
      );
    default:
      *interpolate = 0;
    }
}

//----------------------------------------------------------------------------
// Interpolation for precomputed weights

template <class F, class T>
struct vtkImageGaussRowInterpolate
{
  static void General(
    vtkInterpolationWeights *weights, int idX, int idY, int idZ,
    F *outPtr, int n);
};


//--------------------------------------------------------------------------
// helper function for high-order interpolation
template<class F, class T>
void vtkImageGaussRowInterpolate<F, T>::General(
  vtkInterpolationWeights *weights, int idX, int idY, int idZ,
  F *outPtr, int n)
{
  int stepX = weights->KernelSize[0];
  int stepY = weights->KernelSize[1];
  int stepZ = weights->KernelSize[2];
  idX *= stepX;
  idY *= stepY;
  idZ *= stepZ;
  const F *fX = static_cast<F *>(weights->Weights[0]) + idX;
  const F *fY = static_cast<F *>(weights->Weights[1]) + idY;
  const F *fZ = static_cast<F *>(weights->Weights[2]) + idZ;
  const vtkIdType *factX = weights->Positions[0] + idX;
  const vtkIdType *factY = weights->Positions[1] + idY;
  const vtkIdType *factZ = weights->Positions[2] + idZ;
  const T *inPtr = static_cast<const T *>(weights->Pointer);

  int numscalars = weights->NumberOfComponents;
  for (int i = n; i > 0; --i)
    {
    const T *inPtr0 = inPtr;
    int c = numscalars;
    do // loop over components
      {
      F val = 0;
      int k = 0;
      do // loop over z
        {
        F ifz = fZ[k];
        vtkIdType factz = factZ[k];
        int j = 0;
        do // loop over y
          {
          F ify = fY[j];
          F fzy = ifz*ify;
          vtkIdType factzy = factz + factY[j];
          // loop over x
          const T *tmpPtr = inPtr0 + factzy;
          const F *tmpfX = fX;
          const vtkIdType *tmpfactX = factX;
          F tmpval = 0;
          int l = stepX;
          do
            {
            tmpval += (*tmpfX++)*tmpPtr[(*tmpfactX++)];
            }
          while (--l);
          val += fzy*tmpval;
          }
        while (++j < stepY);
        }
      while (++k < stepZ);

      *outPtr++ = val;
      inPtr0++;
      }
    while (--c);

    factX += stepX;
    fX += stepX;
    }
}

//----------------------------------------------------------------------------
// get row interpolation function for different interpolation modes
// and different scalar types
template<class F>
void vtkGaussianInterpolatorGetRowInterpolationFunc(
  void (**summation)(vtkInterpolationWeights *weights, int idX, int idY,
                     int idZ, F *outPtr, int n),
  int scalarType, int vtkNotUsed(interpolationMode))
{
  switch (scalarType)
    {
    vtkTemplateAliasMacro(
      *summation = &(vtkImageGaussRowInterpolate<F,VTK_TT>::General)
      );
    default:
      *summation = 0;
    }
}

//----------------------------------------------------------------------------
template<class F>
void vtkGaussianInterpolatorPrecomputeWeights(
  const F newmat[16], const int outExt[6], int clipExt[6],
  const F bounds[6], vtkInterpolationWeights *weights)
{
  float **kernel = static_cast<float **>(weights->ExtraInfo);
  weights->WeightType = vtkTypeTraits<F>::VTKTypeID();
  int sizes[3];
  bool blur[3];
  int mode = weights->InterpolationMode;
  sizes[0] = 2*((mode & VTK_INTERPOLATION_WINDOW_XSIZE_MASK)
                >> VTK_INTERPOLATION_WINDOW_XSIZE_SHIFT);
  sizes[1] = 2*((mode & VTK_INTERPOLATION_WINDOW_YSIZE_MASK)
                >> VTK_INTERPOLATION_WINDOW_YSIZE_SHIFT);
  sizes[2] = 2*((mode & VTK_INTERPOLATION_WINDOW_ZSIZE_MASK)
                >> VTK_INTERPOLATION_WINDOW_ZSIZE_SHIFT);
  blur[0] = ((mode & VTK_INTERPOLATION_WINDOW_XBLUR_MASK) != 0);
  blur[1] = ((mode & VTK_INTERPOLATION_WINDOW_YBLUR_MASK) != 0);
  blur[2] = ((mode & VTK_INTERPOLATION_WINDOW_ZBLUR_MASK) != 0);

  // set up input positions table for interpolation
  bool validClip = true;
  for (int j = 0; j < 3; j++)
    {
    // set k to the row for which the element in column j is nonzero,
    // and set matrow to the elements of that row
    int k = 0;
    const F *matrow = newmat;
    while (k < 3 && matrow[j] == 0)
      {
      k++;
      matrow += 4;
      }

    // get the extents
    clipExt[2*j] = outExt[2*j];
    clipExt[2*j + 1] = outExt[2*j + 1];
    int minExt = weights->Extent[2*k];
    int maxExt = weights->Extent[2*k + 1];
    F minBounds = bounds[2*k];
    F maxBounds = bounds[2*k + 1];

    // the kernel size should not exceed the input dimension
    int m = sizes[j];
    int m2 = ((m - 1) >> 1);
    int step = m;
    int inCount = maxExt - minExt + 1;
    step = ((step < inCount) ? step : inCount);

    // allocate space for the weights
    vtkIdType size = step*(outExt[2*j+1] - outExt[2*j] + 1);
    vtkIdType *positions = new vtkIdType[size];
    positions -= step*outExt[2*j];
    F *constants = new F[size];
    constants -= step*outExt[2*j];

    weights->KernelSize[j] = step;
    weights->Positions[j] = positions;
    weights->Weights[j] = constants;
    weights->WeightExtent[2*j] = outExt[2*j];
    weights->WeightExtent[2*j+1] = outExt[2*j+1];

    int region = 0;
    for (int i = outExt[2*j]; i <= outExt[2*j+1]; i++)
      {
      F point = matrow[3] + i*matrow[j];

      F f = 0;
      int idx = vtkInterpolationMath::Floor(point, f);
      int lmax = 1;
      if (step > 1)
        {
        idx -= m2;
        lmax = m;
        }

      int inId[VTK_GAUSS_KERNEL_SIZE_MAX];

      int l = 0;
      switch (weights->BorderMode)
        {
        case VTK_IMAGE_BORDER_REPEAT:
          do
            {
            inId[l] = vtkInterpolationMath::Wrap(idx++, minExt, maxExt);
            }
          while (++l < lmax);
          break;

        case VTK_IMAGE_BORDER_MIRROR:
          do
            {
            inId[l] = vtkInterpolationMath::Mirror(idx++, minExt, maxExt);
            }
          while (++l < lmax);
          break;

        default:
           do
            {
            inId[l] = vtkInterpolationMath::Clamp(idx++, minExt, maxExt);
            }
          while (++l < lmax);
          break;
        }

      // compute the weights and offsets
      vtkIdType inInc = weights->Increments[k];
      if (step == 1)
        {
        positions[step*i] = inId[0]*inInc;
        constants[step*i] = static_cast<F>(1);
        }
      else
        {
        F g[VTK_GAUSS_KERNEL_SIZE_MAX];
        vtkGaussInterpWeights(kernel[j], g, f, m);

        if (step == m)
          {
          int ll = 0;
          do
            {
            positions[step*i + ll] = inId[ll]*inInc;
            constants[step*i + ll] = g[ll];
            }
          while (++ll < step);
          }
        else
          {
          // it gets tricky if the data is thinner than the kernel
          F gg[VTK_GAUSS_KERNEL_SIZE_MAX];
          int ll = 0;
          do { gg[ll] = 0; } while (++ll < m);
          ll = 0;
          do
            {
            int rIdx = inId[ll];
            gg[rIdx] += g[ll];
            }
          while (++ll < m);
          ll = 0;
          do
            {
            positions[step*i + ll] = ll*inInc;
            constants[step*i + ll] = gg[ll];
            }
          while (++ll < step);
          }
        }

      if (point >= minBounds && point <= maxBounds)
        {
        if (region == 0)
          { // entering the input extent
          region = 1;
          clipExt[2*j] = i;
          }
        }
      else
        {
        if (region == 1)
          { // leaving the input extent
          region = 2;
          clipExt[2*j+1] = i - 1;
          }
        }
      }

    if (region == 0 || clipExt[2*j] > clipExt[2*j+1])
      { // never entered input extent!
      validClip = false;
      }
    }

  if (!validClip)
    {
    // output extent doesn't itersect input extent
    for (int j = 0; j < 3; j++)
      {
      clipExt[2*j] = outExt[2*j];
      clipExt[2*j + 1] = outExt[2*j] - 1;
      }
    }
}


//----------------------------------------------------------------------------
} // ends anonymous namespace

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::GetInterpolationFunc(
  void (**func)(vtkInterpolationInfo *, const double [3], double *))
{
  vtkGaussianInterpolatorGetInterpolationFunc(
    func, this->InterpolationInfo->ScalarType);
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::GetInterpolationFunc(
  void (**func)(vtkInterpolationInfo *, const float [3], float *))
{
  vtkGaussianInterpolatorGetInterpolationFunc(
    func, this->InterpolationInfo->ScalarType);
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::GetRowInterpolationFunc(
  void (**func)(vtkInterpolationWeights *, int, int, int, double *, int))
{
  vtkGaussianInterpolatorGetRowInterpolationFunc(
    func, this->InterpolationInfo->ScalarType, this->KernelType);
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::GetRowInterpolationFunc(
  void (**func)(vtkInterpolationWeights *, int, int, int, float *, int))
{
  vtkGaussianInterpolatorGetRowInterpolationFunc(
    func, this->InterpolationInfo->ScalarType, this->KernelType);
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::PrecomputeWeightsForExtent(
  const double matrix[16], const int extent[6], int newExtent[6],
  vtkInterpolationWeights *&weights)
{
  weights = new vtkInterpolationWeights(*this->InterpolationInfo);

  vtkGaussianInterpolatorPrecomputeWeights(
    matrix, extent, newExtent, this->StructuredBoundsDouble, weights);
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::PrecomputeWeightsForExtent(
  const float matrix[16], const int extent[6], int newExtent[6],
  vtkInterpolationWeights *&weights)
{
  weights = new vtkInterpolationWeights(*this->InterpolationInfo);

  vtkGaussianInterpolatorPrecomputeWeights(
    matrix, extent, newExtent, this->StructuredBoundsFloat, weights);
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::FreePrecomputedWeights(
  vtkInterpolationWeights *&weights)
{
  this->Superclass::FreePrecomputedWeights(weights);
}

//----------------------------------------------------------------------------
// build any tables required for the interpolation
void vtkGaussianInterpolator::BuildKernelLookupTable()
{
  if (this->KernelLookupTable[0])
    {
    this->FreeKernelLookupTable();
    }

  float *kernel[3];
  kernel[0] = 0;
  kernel[1] = 0;
  kernel[2] = 0;

  for (int i = 0; i < 3; i++)
    {
    // reuse the X kernel lookup table if possible
    if (i > 0 && this->KernelSize[i] == this->KernelSize[0] &&
        fabs(this->RadiusFactors[i] - this->RadiusFactors[0]) <
          VTK_INTERPOLATE_FLOOR_TOL &&
        fabs(this->BlurFactors[i] - this->BlurFactors[0]) <
          VTK_INTERPOLATE_FLOOR_TOL)
      {
      kernel[i] = kernel[0];
      continue;
      }

    // kernel parameters
    int m = this->KernelSize[i];
    double b = this->BlurFactors[i];

    // blur factor must be restricted to half the max kernel size
    if (b > 0.5*VTK_GAUSS_KERNEL_SIZE_MAX)
      {
      b = 0.5*VTK_GAUSS_KERNEL_SIZE_MAX;
      }

    // compute lookup table size and step size
    int size = m/2*VTK_GAUSS_KERNEL_TABLE_DIVISIONS;
    double p = 1.0/(b*VTK_GAUSS_KERNEL_TABLE_DIVISIONS);

    // allocate and compute the kernel lookup table
    // (add a small safety margin for when the table is interpolated)
    kernel[i] = new float[size + 4];

#if 0 /* print out kernel lookup table, for debugging */
    if (i == 0)
      {
      float *ktemp[4];
      for (int j = 0; j < 4; j++)
        {
        ktemp[j] = new float[size];
        }
      vtkGaussKernel::D0(ktemp[0], size, p);
      vtkGaussKernel::D2(ktemp[1], size, p);
      vtkGaussKernel::D6(ktemp[2], size, p);
      vtkGaussKernel::D10(ktemp[3], size, p);
      double x = 0.0;
      cout << "=============== \n";
      for (int k = 0; k < size; k++)
        {
        cout << x;
        for (int j = 0; j < 4; j++)
          {
          cout << "," << ktemp[j][k];
          }
        cout << "\n";
        x += p;
        }
      cout << "=============== \n";
      for (int j = 0; j < 4; j++)
        {
        delete [] ktemp[j];
        }
      }
#endif

    int cutoff = static_cast<int>(this->RadiusFactors[i]*b/p + 0.5);
    cutoff = (cutoff < size ? cutoff : size);
    cutoff += 1;
    if (this->KernelType == VTK_APPLEDORN10_INTERPOLATION)
      {
      vtkGaussKernel::D10(kernel[i], cutoff, p);
      }
    else if (this->KernelType == VTK_APPLEDORN6_INTERPOLATION)
      {
      vtkGaussKernel::D6(kernel[i], cutoff, p);
      }
    else if (this->KernelType == VTK_APPLEDORN2_INTERPOLATION)
      {
      vtkGaussKernel::D2(kernel[i], cutoff, p);
      }
    else
      {
      vtkGaussKernel::D0(kernel[i], cutoff, p);
      }

    // add a tail of zeros for when table is interpolated
    float *kptr = &kernel[i][cutoff];
    int k = size + 4 - cutoff;
    do
      {
      *kptr++ = 0;
      }
    while (--k);

    // renormalize the table if requested
    if (this->Renormalization)
      {
      vtkRenormalizeKernel(kernel[i], VTK_GAUSS_KERNEL_TABLE_DIVISIONS, m);
      }
    else if (b > 1.0)
      {
      // if kernel stretched to create blur, divide by stretch factor
      float *ktmp = kernel[i];
      float bf = 1.0/b;
      int j = size + 4;
      do
        {
        *ktmp *= bf;
        ktmp++;
        }
      while (--j);
      }
    }

  this->KernelLookupTable[0] = kernel[0];
  this->KernelLookupTable[1] = kernel[1];
  this->KernelLookupTable[2] = kernel[2];

  this->LastBlurFactors[0] = this->BlurFactors[0];
  this->LastBlurFactors[1] = this->BlurFactors[1];
  this->LastBlurFactors[2] = this->BlurFactors[2];
}

//----------------------------------------------------------------------------
void vtkGaussianInterpolator::FreeKernelLookupTable()
{
  float *kernel = this->KernelLookupTable[0];
  if (kernel)
    {
    delete [] kernel;
    for (int i = 1; i < 3; i++)
      {
      if (this->KernelLookupTable[i] != kernel)
        {
        delete [] this->KernelLookupTable[i];
        }
      }
    }
}