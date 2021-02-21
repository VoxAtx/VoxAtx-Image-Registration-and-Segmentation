/*=========================================================================

  Module: vtkImageMutualInformation.h

  Copyright (c) 2006 Atamai, Inc.
  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageMutualInformation - Mutual information between two images
// .SECTION Description
// vtkImageMutualInformation generates a joint histogram from the two input
// images and uses this joint histogram to compute the mutual information
// between the images.  Each input image must have just a single scalar
// component.  The output of the filter will be the joint histogram, with
// the bins from the first input along the X axis, and the bins from the
// second image along the Y axis.  The number of bins and the bin size
// along each axis must be set before the filter executes.  After the
// filter has executed, the mutual information, normalized mutual information,
// and the value to minimize to register the images can be retrieved.
//
// When using this metric, you must call SetInputRange() to set the range
// for each of the input images.  The values will be clamped to the ranges
// that you specify when the metric is computed.  If the images are likely
// to have some voxel values that are outliers, then a winsorized range
// should be used (that is, a range that excludes the outliers).
//
// References:
//
//  [1] D. Mattes, D.R. Haynor, H. Vesselle, T. Lewellen and W. Eubank,
//      PET-CT Image Registration in the Chest Using Free-form Deformations,
//      IEEE Transactions in Medical Imaging 22:120-128, 2003.
//
//  [2] C. Studholme, D.L.G. Hill and D.J. Hawkes,
//      An Overlap Invariant Measure of 3D Medical Image Alignment,
//      Pattern Recognition 32:71-86, 1999.

#ifndef vtkImageMutualInformation_h
#define vtkImageMutualInformation_h

#include "vtkImageSimilarityMetric.h"

class vtkImageMutualInformationTLS;

class VTK_EXPORT vtkImageMutualInformation : public vtkImageSimilarityMetric
{
public:
  static vtkImageMutualInformation *New();
  vtkTypeMacro(vtkImageMutualInformation, vtkImageSimilarityMetric);

  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the type for the output.  The joint histogram will always be
  // computed using vtkIdType, but since vtkIdType is not directly
  // supported as an image data type, it will be converted to the requested
  // type for use as the output of the filter.  The default type is f