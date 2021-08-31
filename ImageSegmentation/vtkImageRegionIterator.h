/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageRegionIterator.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImageRegionIterator - an image region iterator
// .SECTION Description
// This is an image iterator that can be used to iterate over a
// region of an image.

// .SECTION See also
// vtkImageData vtkImageStencilData vtkImageProgressIterator

#ifndef __WRAP__
#ifndef vtkImageRegionIterator_h
#define vtkImageRegionIterator_h

#include "vtkImageRegionIteratorBase.h"

template<class DType>
class VTK_EXPORT vtkImageRegionIterator :
  public vtkImageRegionIteratorBase
{
public:
  // Description:
  // Default constructor, its use must be followed by Initialize().
  vtkImageRegionIterator()
    {
    this->Increment = 0;
    this->BasePointer = 0;
    this->Pointer = 0;
    this->SpanEndPointer = 0;
    }

  // Description:
  // Create an iterator for the given image, with several options.
  // If a stencil is provided, then the iterator's IsInStencil() method
  // reports whether each span is inside the stencil.  If an extent is
  // provided, it iterates over the extent and ignores the rest of the
  // image (the provided extent must be within the image extent).  If
  // a pointer to the algorithm is provided and threadId 