
/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkImageAmoebaGrid.cxx,v $
  Language:  C++
  Date:      $Date: 2007/08/24 20:02:25 $
  Version:   $Revision: 1.8 $

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
#include "vtkImageAmoebaGrid.h"

#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkObjectFactory.h"

#if (VTK_MAJOR_VERSION >= 5)
#include "vtkInformation.h"
#include "vtkExecutive.h"
#endif

//----------------------------------------------------------------------------
vtkImageAmoebaGrid* vtkImageAmoebaGrid::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret = vtkObjectFactory::CreateInstance("vtkImageAmoebaGrid");
  if(ret)
    {
    return (vtkImageAmoebaGrid*)ret;
    }
  // If the factory was unable to create the object, then create it here.
  return new vtkImageAmoebaGrid;
}

// on i386 platforms, SetOptimization(2) provides integer-based calculations
// (on other platforms, integer math is slower than float math)
#ifdef i386  // kwang to test the performance gain on altix IA64 system
#define VTK_RESLICE_INTEGER_MATH 1
#endif

//----------------------------------------------------------------------------
vtkImageAmoebaGrid::vtkImageAmoebaGrid()
{
  this->ShrinkFactors[0] = 1.0f;
  this->ShrinkFactors[1] = 1.0f;
  this->ShrinkFactors[2] = 1.0f;
  
  this->KernelRadius[0] = 1;
  this->KernelRadius[1] = 1;
  this->KernelRadius[2] = 1;

  this->ReverseStencil = 0;

  this->LastThreadCount = 0;
  this->VectorLength = NULL;
  this->VectorsMinimized = NULL;
  this->TotalCost = NULL;
  
  this->Tolerance = 0.005;

#if (VTK_MAJOR_VERSION >= 5)
  // we have the image inputs and the optional stencil input
  this->SetNumberOfInputPorts(2);
#endif
}

//----------------------------------------------------------------------------
vtkImageAmoebaGrid::~vtkImageAmoebaGrid()
{
  if (this->VectorLength)
    {
    delete [] this->VectorLength;
    }

  if (this->VectorsMinimized)
    {
    delete [] this->VectorsMinimized;
    }

  if (this->TotalCost)
    {
    delete [] this->TotalCost;
    }
}

//----------------------------------------------------------------------------
#if (VTK_MAJOR_VERSION >= 5)
int vtkImageAmoebaGrid::FillInputPortInformation(int port,
                                                 vtkInformation *info)
{
  if (port == 1)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageStencilData");
    // the stencil input is optional
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }
  else
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    }
  return 1;
}
#endif

//----------------------------------------------------------------------------
void vtkImageAmoebaGrid::SetStencil(vtkImageStencilData *stencil)
{
#if (VTK_MAJOR_VERSION >= 5)
  this->SetNthInputConnection(1, 0, 
    (stencil ? stencil->GetProducerPort() : 0));
#else
  this->vtkProcessObject::SetNthInput(3, stencil);
#endif
}

//----------------------------------------------------------------------------
vtkImageStencilData *vtkImageAmoebaGrid::GetStencil()
{
#if (VTK_MAJOR_VERSION >= 5)
  if (this->GetNumberOfInputConnections(1) < 1) 
    { 
    return NULL;
    }
  return vtkImageStencilData::SafeDownCast(
    this->GetExecutive()->GetInputData(1, 0));
#else
  if (this->NumberOfInputs < 4)
    {
    return NULL;
    }
  else
    {
    return (vtkImageStencilData *)(this->Inputs[3]);
    }    
#endif
}

#ifdef VTK_RESLICE_INTEGER_MATH
//----------------------------------------------------------------------------
// for fixed-point math, define a basic set of macros

// the radix is 14, to support signed values up to 1<<16
#define VTK_FP_RADIX 14
#define VTK_FP_RADIX_MINUS_1 13

// 0.5 in fixed-point
#define VTK_FP_HALF (1<<VTK_FP_RADIX_MINUS_1)

// various integer values in fixed-point
#define VTK_FP_0 0
#define VTK_FP_1 (1<<VTK_FP_RADIX)
#define VTK_FP_2 (2<<VTK_FP_RADIX)
#define VTK_FP_3 (3<<VTK_FP_RADIX)
#define VTK_FP_4 (4<<VTK_FP_RADIX)

// a very nifty hack I discovered at http://www.stereopsis.com/FPU.html,
// it adds (2**(52-radix))*1.5 to a number to get just the right roundoff 
// from double to fixed-point
static inline int vtkCastFloatToFixed(double x)
{
  union { double d; unsigned int i[2]; } dual;
  dual.d = x + 412316860416.0; // (2**(52-radix))*1.5
#ifdef VTK_WORDS_BIGENDIAN
  return dual.i[1];
#else
  return dual.i[0];
#endif
}

//----------------------------------------------------------------------------
// converting the other way is much more straightforward
static inline double vtkCastFixedToFloat(int x)
{
  return x*(1.0/VTK_FP_1);
}

//----------------------------------------------------------------------------
// what follows are a whole bunch of inline functions that provide
// equivalent behaviour for fixed-point numbers vs. floats (note
// that doubles are not supported, but would be easy to add)
static inline int vtkResliceFloor(int x)
{
  return x>>VTK_FP_RADIX;
}

//----------------------------------------------------------------------------
static inline int vtkResliceCeil(int x)
{
  return ((1<<VTK_FP_RADIX) - 1 + x)>>VTK_FP_RADIX;
}

//----------------------------------------------------------------------------
static inline int vtkResliceRound(int x)
{
  return (x + VTK_FP_HALF)>>VTK_FP_RADIX;
}

//----------------------------------------------------------------------------
// convert a fixed-point into an integer plus a fraction  
static inline int vtkResliceFloor(int x, int &f)
{
  int ix = x>>VTK_FP_RADIX;
  f = x - (ix<<VTK_FP_RADIX);

  return ix;
}

//----------------------------------------------------------------------------
// multiplication, the product must be less than 8 for fixed-point
static inline int vtkResliceQuikMul(int xy)
{
  return (xy + VTK_FP_HALF)>>VTK_FP_RADIX;
}

//----------------------------------------------------------------------------
// multiplication of larger fixed-point numbers, does not check overflow
static inline int vtkResliceMultiply(int x, int y)
{
  int hx = x>>VTK_FP_RADIX;
  int hy = y>>VTK_FP_RADIX;
  int lx = x - (hx<<VTK_FP_RADIX);
  int ly = y - (hy<<VTK_FP_RADIX);
  
  return ((lx*ly + VTK_FP_HALF)>>VTK_FP_RADIX) + hx*ly + x*hy;
}

//----------------------------------------------------------------------------
static inline int vtkResliceInverse(int x)
{
  return ((1<<(2*VTK_FP_RADIX + 1))/x + 1)>>1;
}

//----------------------------------------------------------------------------
// cast between float and fixed-point
static inline void vtkResliceCast(float x, int& y)
{
  y = vtkCastFloatToFixed(x);
}

//----------------------------------------------------------------------------
static inline void vtkResliceCast(double x, int& y)
{
  y = vtkCastFloatToFixed(x);
}

//----------------------------------------------------------------------------
static inline void vtkResliceCast(int x, float& y)
{
  y = vtkCastFixedToFloat(x);
}

//----------------------------------------------------------------------------
// 1-x, it gets used a lot
static inline int vtkResliceOneMinusX(int x)
{
  return VTK_FP_1 - x;
}

//----------------------------------------------------------------------------
// check if a number is equal to one
static inline int vtkResliceIsEqualToOne(int x)
{
  return (VTK_FP_1 == x);
}

//----------------------------------------------------------------------------
// check if a number is an integer
static inline int vtkResliceIsInteger(int x)
{
  return (x == ((x>>VTK_FP_RADIX)<<VTK_FP_RADIX));
}

//----------------------------------------------------------------------------
// rounding functions for fixed-point, note that only
// char, unsigned char, short, unsigned short are supported
static inline void vtkResliceRound(int val, char& rnd)
{
  rnd = (char)vtkResliceRound(val);
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(int val, unsigned char& rnd)
{
  rnd = (unsigned char)vtkResliceRound(val);
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(int val, short& rnd)
{
  rnd = (short)vtkResliceRound(val);
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(int val, unsigned short& rnd)
{
  rnd = (unsigned short)vtkResliceRound(val);
}
#endif /* VTK_RESLICE_INTEGER_MATH */

//----------------------------------------------------------------------------
// fast floor() function for converting a float to an int
// (the floor() implementation on some computers is much slower than this,
// because they require some 'exact' behaviour that we don't).

// The 'floor' function on x86 and mips is many times slower than these
// and is used a lot in this code, optimize for different CPU architectures
inline int vtkResliceFloor(double x)
{
#if defined mips || defined sparc || defined __ppc__
  x += 2147483648.0;
  unsigned int i = (unsigned int)(x);
  return (int)(i - 2147483648U);
#elif defined i386 || defined _M_IX86
  union { double d; unsigned short s[4]; unsigned int i[2]; } dual;
  dual.d = x + 103079215104.0;  // (2**(52-16))*1.5
  return (int)((dual.i[1]<<16)|((dual.i[0])>>16));
#elif defined ia64 || defined __ia64__ || defined IA64
  x += 103079215104.0;
  long long i = (long long)(x);
  return (int)(i - 103079215104LL);
#else
  double y = floor(x);
  return (int)(y);
#endif
}

//----------------------------------------------------------------------------
static inline int vtkResliceCeil(double x)
{
  return -vtkResliceFloor(-x - 1.0) - 1;
}

//----------------------------------------------------------------------------
static inline int vtkResliceRound(double x)
{
  return vtkResliceFloor(x + 0.5);
}

//----------------------------------------------------------------------------
static inline int vtkResliceFloor(float x)
{
  return vtkResliceFloor((double)x);
}

//----------------------------------------------------------------------------
static inline int vtkResliceCeil(float x)
{
  return vtkResliceCeil((double)x);
}

//----------------------------------------------------------------------------
static inline int vtkResliceRound(float x)
{
  return vtkResliceRound((double)x);
}

//----------------------------------------------------------------------------
// convert a float into an integer plus a fraction  
static inline int vtkResliceFloor(float x, float &f)
{
  int ix = vtkResliceFloor(x);
  f = x - ix;

  return ix;
}

//----------------------------------------------------------------------------
static inline float vtkResliceQuikMul(float xy)
{
  return xy;
}

//----------------------------------------------------------------------------
static inline float vtkResliceMultiply(float x, float y)
{
  return x*y;
}

//----------------------------------------------------------------------------
// invert a number
static inline float vtkResliceInverse(float x)
{
  return 1.0f/x;
}

//----------------------------------------------------------------------------
static inline void vtkResliceCast(float x, float& y)
{
  y = x;
}

//----------------------------------------------------------------------------
static inline void vtkResliceCast(double x, float& y)
{
  y = (float)x;
}

//----------------------------------------------------------------------------
static inline float vtkResliceOneMinusX(float x)
{
  return 1.0f - x;
}

//----------------------------------------------------------------------------
static inline int vtkResliceIsEqualToOne(float x)
{
  return (1.0f == x);
}

//----------------------------------------------------------------------------
static inline int vtkResliceIsInteger(float x)
{
  return (x == vtkResliceFloor(x));
}


//----------------------------------------------------------------------------
// rounding functions for each type, with some crazy stunts to avoid
// the use of the 'floor' function which is too slow on x86
template<class T>
static inline void vtkResliceRound(float val, T& rnd)
{
  rnd = vtkResliceRound(val);
}

//----------------------------------------------------------------------------
template<class T>
static inline void vtkResliceRound(double val, T& rnd)
{
  rnd = vtkResliceRound(val);
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(float val, float& rnd)
{
  rnd = val;
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(double val, float& rnd)
{
  rnd = val;
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(float val, double& rnd)
{
  rnd = val;
}

//----------------------------------------------------------------------------
static inline void vtkResliceRound(double val, double& rnd)
{
  rnd = val;
}

//----------------------------------------------------------------------------
void vtkImageAmoebaGrid::ComputeInputUpdateExtents(vtkDataObject *output)
{
  this->vtkImageMultipleInputFilter::ComputeInputUpdateExtents(output);

  vtkImageStencilData *stencil = this->GetStencil();
  if (stencil)
    {
    stencil->SetUpdateExtent(output->GetUpdateExtent());
    }
  
  if (this->LastThreadCount != this->GetNumberOfThreads())
    {
    this->LastThreadCount = this->GetNumberOfThreads();
    if (this->VectorLength)
      {
      delete [] this->VectorLength;
      }
    this->VectorLength = new double[this->LastThreadCount];

    if (this->VectorsMinimized)
      {
      delete [] this->VectorsMinimized;
      }
    this->VectorsMinimized = new int[this->LastThreadCount];

    if (this->TotalCost)
      {
      delete [] this->TotalCost;
      }
    this->TotalCost = new double[this->LastThreadCount];

    // these arrays should be initialized otherwise will cause problem 
    // during summing. add by kwang 05/20/2004
    for (int i=0; i<this->LastThreadCount; i++)
      {
        this->VectorLength[i] = 0.0f;
        this->VectorsMinimized[i] = 0;
        this->TotalCost[i] = 0.0f;
      }
    }  
}

//----------------------------------------------------------------------------
// Compute and return the mean vector length
float vtkImageAmoebaGrid::GetMeanVectorLength()
{
  double totalLength = 0;
  for (int x=0; x<this->LastThreadCount; x++)
    {
    totalLength += this->VectorLength[x];
    }
  int totalVectors = this->GetVectorsMinimized();
  if (totalVectors)
    {
    return totalLength / totalVectors;
    }
  else
    {
    return 0.0f;
    }
}

//----------------------------------------------------------------------------
// Compute and return the mean cost function
float vtkImageAmoebaGrid::GetMeanCost()
{
  double totalCost = 0.0f;
  for (int x=0; x<this->LastThreadCount; x++)
    {
    totalCost += this->TotalCost[x];
    }
  int totalVectors = this->GetVectorsMinimized();
  if (totalVectors)
    {
    return totalCost / totalVectors;
    }
  else
    {
    return 0.0f;
    }
}

//----------------------------------------------------------------------------
// Compute and return the total number of vectors minimized in the most recent
// iteration.
int vtkImageAmoebaGrid::GetVectorsMinimized()
{
  int minimized = 0;
  for (int x=0; x<this->LastThreadCount; x++)
    {
    minimized += this->VectorsMinimized[x];
    }
  return minimized;
}

//----------------------------------------------------------------------------