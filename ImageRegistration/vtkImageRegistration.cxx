/*=========================================================================

  Module: vtkImageRegistration.cxx

  Copyright (c) 2006 Atamai, Inc.
  Copyright (c) 2016 David Gobbi
  All rights reserved.
  See Copyright.txt or http://dgobbi.github.io/bsd3.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkImageRegistration.h"

// VTK header files
#include <vtkTimerLog.h>
#include <vtkObjectFactory.h>
#include <vtkImageData.h>
#include <vtkImageStencilData.h>
#include <vtkMath.h>
#include <vtkDoubleArray.h>
#include <vtkTransform.h>
#include <vtkMatrixToLinearTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkImageReslice.h>
#include <vtkImageShiftScale.h>
#include <vtkCommand.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkImageHistogramStatistics.h>
#include <vtkImageBSplineCoefficients.h>
#include <vtkImageBSplineInterpolator.h>
#include <vtkImageSincInterpolator.h>
#include <vtkVersion.h>

// Interpolator header files
#include "vtkLabelInterpolator.h"

// Optimizer header files
#include "vtkNelderMeadMinimizer.h"
#include "vtkPowellMinimizer.h"

// Image metric header files
#include "vtkImageSquaredDifference.h"
#include "vtkImageMutualInformation.h"
#include "vtkImageCorrelationRatio.h"
#include "vtkImageCrossCorrelation.h"
#include "vtkImageNeighborhoodCorrelation.h"

// C header files
#include <math.h>

// A macro to assist VTK 5 backwards compatibility
#if VTK_MAJOR_VERSION >= 6
#define SET_INPUT_DATA SetInputData
#define SET_STENCIL_DATA SetStencilData
#else
#define SET_INPUT_DATA SetInput
#define SET_STENCIL_DATA SetStencil
#endif

// A helper class for the optimizer
struct vtkImageRegistrationInfo
{
  vtkTransform *Transform;
  vtkFunctionMinimizer *Optimizer;
  vtkImageSimilarityMetric *Metric;
  vtkMatrix4x4 *InitialMatrix;

  vtkDoubleArray *MetricValues;
  vtkDoubleArray *CostValues;
  vtkDoubleArray *ParameterValues;

  int TransformDimensionality;
  int TransformType;
  int OptimizerType;
  int MetricType;

  double Center[3];

  int NumberOfEvaluations;
};

//----------------------------------------------------------------------------
vtkImageRegistration* vtkImageRegistration::New()
{
  // First try to create the object from the vtkObjectFactory
  vtkObject* ret =
    vtkObjectFactory::CreateInstance("vtkImageRegistration");

  if (ret)
    {
    return (vtkImageRegistration*)ret;
    }

  // If the factory was unable to create the object, then create it here.
  return new vtkImageRegistration;
}

//----------------------------------------------------------------------------
vtkImageRegistration::vtkImageRegistration()
{
  this->OptimizerType = vtkImageRegistration::Powell;
  this->MetricType = vtkImageRegistration::MutualInformation;
  this->InterpolatorType = vtkImageRegistration::Linear;
  this->TransformType = vtkImageRegistration::Rigid;
  this->InitializerType = vtkImageRegistration::None;
  this->TransformDimensionality = 3;

  this->Transform = vtkTransform::New();
  this->Metric = NULL;
  this->Optimizer = NULL;
  this->Interpolator = NULL;

  this->RegistrationInfo = new vtkImageRegistrationInfo;
  this->RegistrationInfo->Transform = NULL;
  this->RegistrationInfo->Optimizer = NULL;
  this->RegistrationInfo->Metric = NULL;
  this->RegistrationInfo->InitialMatrix = NULL;
  this->RegistrationInfo->MetricValues = NULL;
  this->RegistrationInfo->CostValues = NULL;
  this->RegistrationInfo->ParameterValues = NULL;
  this->RegistrationInfo->TransformDimensionality = 0;
  this->RegistrationInfo->TransformType = 0;
  this->RegistrationInfo->OptimizerType = 0;
  this->RegistrationInfo->MetricType = 0;
  this->RegistrationInfo->NumberOfEvaluations = 0;

  this->JointHistogramSize[0] = 64;
  this->JointHistogramSize[1] = 64;
  this->SourceImageRange[0] = 0.0;
  this->SourceImageRange[1] = -1.0;
  this->TargetImageRange[0] = 0.0;
  this->TargetImageRange[1] = -1.0;

  this->InitialTransformMatrix = vtkMatrix4x4::New();
  this->ImageReslice = vtkImageReslice::New();
  this->ImageBSpline = vtkImageBSplineCoefficients::New();
  this->TargetImageTypecast = vtkImageShiftScale::New();
  this->SourceImageTypecast = vtkImageShiftScale::New();

  this->MetricValue = 0.0;
  this->CostValue = 0.0;

  this->CollectValues = false;
  this->MetricValues = vtkDoubleArray::New();
  this->CostValues = vtkDoubleArray::New();
  this->ParameterValues = vtkDoubleArray::New();

  this->CostTolerance = 1e-4;
  this->TransformTolerance = 1e-1;
  this->MaximumNumberOfIterations = 500;
  this->MaximumNumberOfEvaluations = 5000;

  // we have the image inputs and the optional stencil input
  this->SetNumberOfInputPorts(3);
  this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
vtkImageRegistration::~vtkImageRegistration()
{
  // delete vtk objects
  if (this->Optimizer)
    {
    this->Optimizer->Delete();
    }
  if (this->Metric)
    {
    this->Metric->Delete();
    }
  if (this->Interpolator)
    {
    this->Interpolator->Delete();
    }
  if (this->Transform)
    {
    this->Transform->Delete();
    }
  if (this->MetricValues)
    {
    this->MetricValues->Delete();
    }
  if (this->CostValues)
    {
    this->CostValues->Delete();
    }
  if (this->ParameterValues)
    {
    this->ParameterValues->Delete();
    }

  if (this->RegistrationInfo)
    {
    delete this->RegistrationInfo;
    }

  if (this->InitialTransformMatrix)
    {
    this->InitialTransformMatrix->Delete();
    }
  if (this->ImageReslice)
    {
    this->ImageReslice->Delete();
    }
  if (this->SourceImageTypecast)
    {
    this->SourceImageTypecast->Delete();
    }
  if (this->TargetImageTypecast)
    {
    this->TargetImageTypecast->Delete();
    }
  if (this->ImageBSpline)
    {
    this->ImageBSpline->Delete();
    }
}

//----------------------------------------------------------------------------
void vtkImageRegistration::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "OptimizerType: " << this->OptimizerType << "\n";
  os << indent << "MetricType: " << this->MetricType << "\n";
  os << indent << "InterpolatorType: " << this->InterpolatorType << "\n";
  os << indent << "TransformType: " << this->TransformType << "\n";
  os << indent << "TransformDimensionality: "
     << this->TransformDimensionality << "\n";
  os << indent << "InitializerType: " << this->InitializerType << "\n";
  os << indent << "CostTolerance: " << this->CostTolerance << "\n";
  os << indent << "TransformTolerance: " << this->TransformTolerance << "\n";
  os << indent << "MaximumNumberOfIterations: "
     << this->MaximumNumberOfIterations <