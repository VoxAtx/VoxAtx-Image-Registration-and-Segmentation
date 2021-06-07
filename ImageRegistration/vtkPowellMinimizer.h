/*=========================================================================

  Module: vtkPowellMinimizer.h

  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkPowellMinimizer - use Powell's method to minimize a function
// .SECTION Description
// vtkPowellMinimizer will modify a set of parameters in order to find
// the minimum of a specified function.  This method conducts a series
// of linear searches and attempts to construct a conjugate set of search
// directions as it goes.

#ifndef vtkPowellMinimizer_h
#define vtkPowellMinimizer_h

#include "vtkFunctionMinimizer.h"

class VTK_EXPORT vtkPowellMinimizer : public vtkFunctionMinimizer
{
public:
  static vtkPowellMinimizer *New();
  vtkTypeMacro(vtkPowellMinimizer,vtkFunctionMinimizer);
  void PrintSelf