
/*=========================================================================

  Program:   Atamai Image Registration and Segmentation
  Module:    vtkImageConnectivityFilter.cxx

  Copyright (c) 2014 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkImageConnectivityFilter.h"

#include "vtkMath.h"
#include "vtkImageData.h"
#include "vtkDataSet.h"
#include "vtkPointData.h"
#include "vtkImageStencilData.h"
#include "vtkImageStencilIterator.h"
#include "vtkImageIterator.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkIdTypeArray.h"
#include "vtkIntArray.h"
#include "vtkImageStencilData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkTemplateAliasMacro.h"
#include "vtkTypeTraits.h"
#include "vtkSmartPointer.h"
#include "vtkVersion.h"

#include <vector>
#include <stack>
#include <algorithm>

vtkStandardNewMacro(vtkImageConnectivityFilter);

//----------------------------------------------------------------------------
// Constructor sets default values
vtkImageConnectivityFilter::vtkImageConnectivityFilter()
{
  this->LabelMode = SeedScalar;
  this->ExtractionMode = SeededRegions;

  this->ScalarRange[0] = 0.5;
  this->ScalarRange[1] = VTK_DOUBLE_MAX;

  this->SizeRange[0] = 1;
#if VTK_MAJOR_VERSION >= 6
  this->SizeRange[1] = VTK_ID_MAX;
#else
  this->SizeRange[1] = VTK_LARGE_ID;
#endif

  this->LabelConstantValue = 255;

  this->ActiveComponent = 0;

  this->LabelScalarType = VTK_UNSIGNED_CHAR;

  this->GenerateRegionExtents = 0;

  this->ExtractedRegionLabels = vtkIdTypeArray::New();
  this->ExtractedRegionSizes = vtkIdTypeArray::New();
  this->ExtractedRegionSeedIds = vtkIdTypeArray::New();
  this->ExtractedRegionExtents = vtkIntArray::New();
  this->ExtractedRegionExtents->SetNumberOfComponents(6);

  this->SetNumberOfInputPorts(3);
}

//----------------------------------------------------------------------------
vtkImageConnectivityFilter::~vtkImageConnectivityFilter()
{
  if (this->ExtractedRegionSizes)
    {
    this->ExtractedRegionSizes->Delete();
    }
  if (this->ExtractedRegionLabels)
    {
    this->ExtractedRegionLabels->Delete();
    }
  if (this->ExtractedRegionSeedIds)
    {
    this->ExtractedRegionSeedIds->Delete();
    }
  if (this->ExtractedRegionExtents)
    {
    this->ExtractedRegionExtents->Delete();
    }
}

//----------------------------------------------------------------------------
int vtkImageConnectivityFilter::FillInputPortInformation(
  int port, vtkInformation* info)
{
  if (port == 2)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }
  else if (port == 1)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageStencilData");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    }
  else
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    }
  return 1;
}

//----------------------------------------------------------------------------
void vtkImageConnectivityFilter::SetStencilConnection(
  vtkAlgorithmOutput *stencil)
{
  this->SetInputConnection(1, stencil);
}

//----------------------------------------------------------------------------
vtkAlgorithmOutput *vtkImageConnectivityFilter::GetStencilConnection()
{
  return this->GetInputConnection(1, 0);
}

//----------------------------------------------------------------------------
void vtkImageConnectivityFilter::SetStencilData(vtkImageStencilData *stencil)
{
#if VTK_MAJOR_VERSION >= 6
  this->SetInputData(1, stencil);
#else
  this->SetInput(1, stencil);
#endif
}

//----------------------------------------------------------------------------
void vtkImageConnectivityFilter::SetSeedConnection(
  vtkAlgorithmOutput *seeds)
{
  this->SetInputConnection(2, seeds);
}

//----------------------------------------------------------------------------
vtkAlgorithmOutput *vtkImageConnectivityFilter::GetSeedConnection()
{
  return this->GetInputConnection(2, 0);
}

//----------------------------------------------------------------------------
void vtkImageConnectivityFilter::SetSeedData(vtkDataSet *seeds)
{
#if VTK_MAJOR_VERSION >= 6
  this->SetInputData(2, seeds);
#else
  this->SetInput(2, seeds);
#endif
}

//----------------------------------------------------------------------------
const char *vtkImageConnectivityFilter::GetLabelScalarTypeAsString()
{
  const char *result = "Unknown";
  switch (this->LabelScalarType)
    {
    case VTK_UNSIGNED_CHAR:
      result = "UnsignedChar";
      break;
    case VTK_SHORT:
      result = "Short";
      break;
    case VTK_UNSIGNED_SHORT:
      result = "UnsignedShort";
      break;
    case VTK_INT:
      result = "Int";
      break;
    }
  return result;
}

//----------------------------------------------------------------------------
const char *vtkImageConnectivityFilter::GetLabelModeAsString()
{
  const char *result = "Unknown";
  switch (this->LabelMode)
    {
    case SeedScalar:
      result = "SeedScalar";
      break;
    case ConstantValue:
      result = "ConstantValue";
      break;
    case SizeRank:
      result = "SizeRank";
      break;
    }
  return result;
}

//----------------------------------------------------------------------------
const char *vtkImageConnectivityFilter::GetExtractionModeAsString()
{
  const char *result = "Unknown";
  switch (this->ExtractionMode)
    {
    case SeededRegions:
      result = "SeededRegions";
      break;
    case AllRegions:
      result = "AllRegions";
      break;
    case LargestRegion:
      result = "LargestRegion";
      break;
    }
  return result;
}

//----------------------------------------------------------------------------
vtkIdType vtkImageConnectivityFilter::GetNumberOfExtractedRegions()
{
  return this->ExtractedRegionLabels->GetNumberOfTuples();
}

//----------------------------------------------------------------------------
namespace {

// Methods for the connectivity algorithm
class vtkICF
{
public:
  // Simple struct that holds information about a region.
  struct Region;

  // A class that is a vector of regions.
  class RegionVector;

protected:
  // Simple class that holds a seed location and a scalar value.
  class Seed;

  // A functor to assist in comparing region sizes.
  struct CompareSize;

  // Remove all but the largest region from the output image.
  template<class OT>
  static void PruneAllButLargest(
    vtkImageData *outData, OT *outPtr, vtkImageStencilData *stencil,
    int extent[6], const OT& value, vtkICF::RegionVector& regionInfo);

  // Remove the smallest region from the output image.
  // This is called when there are no labels left, i.e. when the label
  // value reaches the maximum allowed by the output data type.
  template<class OT>
  static void PruneSmallestRegion(
    vtkImageData *outData, OT *outPtr, vtkImageStencilData *stencil,
    int extent[6], vtkICF::RegionVector& regionInfo);

  // Remove all islands that aren't in the given range of sizes
  template<class OT>
  static void PruneBySize(
    vtkImageData *outData, OT *outPtr, vtkImageStencilData *stencil,
    int extent[6], vtkIdType sizeRange[2], vtkICF::RegionVector& regionInfo);

  // This is the function that grows a region from a seed.
  template<class OT>
  static vtkIdType Fill(
    OT *outPtr, vtkIdType outInc[3], int outLimits[6],
    unsigned char *maskPtr, int maxIdx[3], int fillExtent[6],
    std::stack<vtkICF::Seed> &seedStack);

  // Add a region to the list of regions.
  template<class OT>
  static void AddRegion(
    vtkImageData *outData, OT *outPtr, vtkImageStencilData *stencil,
    int extent[6], vtkIdType sizeRange[2], vtkICF::RegionVector& regionInfo,
    vtkIdType voxelCount, vtkIdType regionId, int regionExtent[6],
    int extractionMode);

  // Fill the ExtractedRegionSizes and ExtractedRegionLabels arrays.
  static void GenerateRegionArrays(
    vtkImageConnectivityFilter *self, vtkICF::RegionVector& regionInfo,
    vtkDataArray *seedScalars, int extent[6], int minLabel, int maxLabel);

  // Relabel the image, usually the last method called.
  template<class OT>
  static void Relabel(
    vtkImageData *outData, OT *outPtr,
    vtkImageStencilData *stencil, int extent[6],
    vtkIdTypeArray *labelMap);

  // Sort the ExtractedRegionLabels array and the other arrays.
  static void SortRegionArrays(vtkImageConnectivityFilter *self);

  // Finalize the output
  template<class OT>
  static void Finish(
    vtkImageConnectivityFilter *self, vtkImageData *outData,
    OT *outPtr, vtkImageStencilData *stencil, int extent[6],
    vtkDataArray *seedScalars, vtkICF::RegionVector& regionInfo);

  // Subtract the lower extent limit from "limits", and return the
  // extent size subtract 1 in maxIdx.
  static int *ZeroBaseExtent(
    const int wholeExtent[6], int extent[6], int maxIdx[3]);

  // Execute method for when point seeds are provided.
  template <class OT>
  static void SeededExecute(
    vtkImageConnectivityFilter *self,
    vtkImageData *outData, vtkDataSet *seedData, vtkImageStencilData *stencil,
    OT *outPtr, unsigned char *maskPtr, int extent[6],
    vtkICF::RegionVector& regionInfo);

  // Execute method for when no seeds are provided.
  template <class OT>
  static void SeedlessExecute(
    vtkImageConnectivityFilter *self,
    vtkImageData *outData, vtkImageStencilData *stencil,
    OT *outPtr, unsigned char *maskPtr, int extent[6],
    vtkICF::RegionVector& regionInfo);

public:
  // Create a bit mask from the input
  template<class IT>
  static void ExecuteInput(
    vtkImageConnectivityFilter *self, vtkImageData *inData, IT *inPtr,
    unsigned char *maskPtr, vtkImageStencilData *stencil, int extent[6]);

  // Generate the output
  template <class OT>
  static void ExecuteOutput(
    vtkImageConnectivityFilter *self,
    vtkImageData *outData, vtkDataSet *seedData, vtkImageStencilData *stencil,
    OT *outPtr, unsigned char *maskPtr, int extent[6]);

  // Utility method to find the intersection of two extents.
  // Returns false if the extents do not intersect.
  static bool IntersectExtents(
    const int extent1[6], const int extent2[6], int output[6]);
};

//----------------------------------------------------------------------------
// seed struct: structured coordinates and a value,
// the coords can be accessed with [] and the value with *
class vtkICF::Seed
{
public:
  Seed() {
    pos[0] = 0; pos[1] = 0; pos[2] = 0; value = 0; }

  Seed(int i, int j, int k) {
    pos[0] = i; pos[1] = j; pos[2] = k; value = 0; }

  Seed(int i, int j, int k, int v) {
    pos[0] = i; pos[1] = j; pos[2] = k; value = v; }

  Seed(const vtkICF::Seed &seed) {
    pos[0] = seed.pos[0];
    pos[1] = seed.pos[1];
    pos[2] = seed.pos[2];
    value = seed.value; };

  const int &operator[](int i) const { return pos[i]; }
  int &operator[](int i) { return pos[i]; }

  const int &operator*() const { return value; }
  int &operator*() { return value; }

  const vtkICF::Seed &operator = (const vtkICF::Seed seed) {
    pos[0] = seed.pos[0];
    pos[1] = seed.pos[1];
    pos[2] = seed.pos[2];
    value = seed.value;
    return *this; }

private:
  int pos[3];
  int value;
};

//----------------------------------------------------------------------------
// region struct: size and id
struct vtkICF::Region
{
  Region(vtkIdType s, vtkIdType i, const int e[6]) : size(s), id(i) {
    extent[0] = e[0]; extent[1] = e[1]; extent[2] = e[2];
    extent[3] = e[3]; extent[4] = e[4]; extent[5] = e[5]; }
  Region() : size(0), id(0) {
    extent[0] = extent[1] = extent[2] = 0;
    extent[3] = extent[4] = extent[5] = 0; }

  vtkIdType size;
  vtkIdType id;
  int extent[6];
};

//----------------------------------------------------------------------------
class vtkICF::RegionVector : public std::vector<vtkICF::Region>
{
public:
  typedef std::vector<vtkICF::Region>::iterator iterator;

  // get the smallest of the regions in the vector
  iterator smallest()
  {
    // start at 1, because 0 is the background
    iterator small = begin() + 1;
    if (small != end())
      {
      for (iterator i = small + 1; i != end(); ++i)
        {
        if (i->size <= small->size)
          {
          small = i;
          }
        }
      }
    return small;
  }

  // get the largest of the regions in the vector
  iterator largest()
  {
    iterator large = begin() + 1;
    if (large != end())
      {
    // start at 1, because 0 is the background
      for (iterator i = large + 1; i != end(); ++i)
        {
        if (i->size > large->size)
          {
          large = i;
          }
        }
      }
    return large;
  }

};

//----------------------------------------------------------------------------
bool vtkICF::IntersectExtents(
  const int extent1[6], const int extent2[6], int output[6])
{
  bool rval = true;

  int i = 0;
  while (i < 6)
    {
    output[i] = (extent1[i] > extent2[i] ? extent1[i] : extent2[i]);
    i++;
    output[i] = (extent1[i] < extent2[i] ? extent1[i] : extent2[i]);
    rval &= (output[i-1] <= output[i]);
    i++;
    }

  return rval;
}

//----------------------------------------------------------------------------
template<class IT>
void vtkICF::ExecuteInput(
  vtkImageConnectivityFilter *self, vtkImageData *inData, IT *,
  unsigned char *maskPtr, vtkImageStencilData *stencil, int extent[6])
{
  // Get active component (only one component is thresholded)
  int nComponents = inData->GetNumberOfScalarComponents();
  int activeComponent = self->GetActiveComponent();
  if (activeComponent < 0 || activeComponent > nComponents)
    {
    activeComponent = 0;
    }

  // Get the scalar range clamped to the input type range
  double drange[2];
  self->GetScalarRange(drange);
  IT srange[2];
  srange[0] = vtkTypeTraits<IT>::Min();
  srange[1] = vtkTypeTraits<IT>::Max();
  if (drange[0] > static_cast<double>(srange[1]))
    {
    srange[0] = srange[1];
    }
  else if (drange[0] > static_cast<double>(srange[0]))
    {
    srange[0] = static_cast<IT>(drange[0]);
    }
  if (drange[1] < static_cast<double>(srange[0]))
    {