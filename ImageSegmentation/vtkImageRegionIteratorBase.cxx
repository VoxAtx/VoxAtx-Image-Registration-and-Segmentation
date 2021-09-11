
/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImageRegionIteratorBase.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkImageRegionIteratorBase.h"
#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkAlgorithm.h"

//----------------------------------------------------------------------------
class vtkImageStencilIteratorFriendship
{
public:

  static int *GetExtentListLengths(vtkImageStencilData *stencil)
    {
    return stencil->ExtentListLengths;
    }

  static int **GetExtentLists(vtkImageStencilData *stencil)
    {
    return stencil->ExtentLists;
    }
};

//----------------------------------------------------------------------------
vtkImageRegionIteratorBase::vtkImageRegionIteratorBase()
{
  this->Id = 0;
  this->SpanEnd = 0;
  this->RowEnd = 0;
  this->SliceEnd = 0;
  this->End = 0;

  this->RowEndIncrement = 0;
  this->RowIncrement = 0;
  this->SliceEndIncrement = 0;
  this->SliceIncrement = 0;

  this->Extent[0] = 0;
  this->Extent[1] = 0;
  this->Extent[2] = 0;
  this->Extent[3] = 0;
  this->Extent[4] = 0;
  this->Extent[5] = 0;

  this->Index[0] = 0;
  this->Index[1] = 0;
  this->Index[2] = 0;
  this->StartY = 0;

  this->HasStencil = false;
  this->InStencil = false;
  this->SpanSliceEndIncrement = 0;
  this->SpanSliceIncrement = 0;
  this->SpanIndex = 0;
  this->SpanCountPointer = 0;
  this->SpanListPointer = 0;

  this->Algorithm = 0;
  this->Count = 0;