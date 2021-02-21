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
// for each of the input images.  The values will be c