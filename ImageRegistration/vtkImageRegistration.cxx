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
     << this->MaximumNumberOfIterations << "\n";
  os << indent << "MaximumNumberOfEvaluations: "
     << this->MaximumNumberOfEvaluations << "\n";
  os << indent << "JointHistogramSize: " << this->JointHistogramSize[0] << " "
     << this->JointHistogramSize[1] << "\n";
  os << indent << "SourceImageRange: " << this->SourceImageRange[0] << " "
     << this->SourceImageRange[1] << "\n";
  os << indent << "TargetImageRange: " << this->TargetImageRange[0] << " "
     << this->TargetImageRange[1] << "\n";
  os << indent << "MetricValue: " << this->MetricValue << "\n";
  os << indent << "CostValue: " << this->CostValue << "\n";
  os << indent << "CollectValues: "
     << (this->CollectValues ? "On\n" : "Off\n");
  os << indent << "MetricValues: " << this->MetricValues << "\n";
  os << indent << "CostValues: " << this->CostValues << "\n";
  os << indent << "ParameterValues: " << this->ParameterValues << "\n";
  os << indent << "NumberOfEvaluations: "
     << this->RegistrationInfo->NumberOfEvaluations << "\n";
}

//----------------------------------------------------------------------------
vtkLinearTransform *vtkImageRegistration::GetTransform()
{
  return this->Transform;
}

//----------------------------------------------------------------------------
int vtkImageRegistration::GetNumberOfEvaluations()
{
  return this->RegistrationInfo->NumberOfEvaluations;
}

//----------------------------------------------------------------------------
void vtkImageRegistration::SetTargetImage(vtkImageData *input)
{
  // Ask the superclass to connect the input.
#if VTK_MAJOR_VERSION >= 6
  this->SetInputDataInternal(1, input);
#else
  this->SetNthInputConnection(1, 0, (input ? input->GetProducerPort() : 0));
#endif
}

//----------------------------------------------------------------------------
vtkImageData* vtkImageRegistration::GetTargetImage()
{
  if (this->GetNumberOfInputConnections(0) < 1)
    {
    return NULL;
    }
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(1, 0));
}

//----------------------------------------------------------------------------
void vtkImageRegistration::SetSourceImage(vtkImageData *input)
{
  // Ask the superclass to connect the input.
#if VTK_MAJOR_VERSION >= 6
  this->SetInputDataInternal(0, input);
#else
  this->SetNthInputConnection(0, 0, (input ? input->GetProducerPort() : 0));
#endif
}

//----------------------------------------------------------------------------
vtkImageData* vtkImageRegistration::GetSourceImage()
{
  if (this->GetNumberOfInputConnections(1) < 1)
    {
    return NULL;
    }
  return vtkImageData::SafeDownCast(this->GetExecutive()->GetInputData(0, 0));
}

//----------------------------------------------------------------------------
void vtkImageRegistration::SetSourceImageStencil(vtkImageStencilData *stencil)
{
  // if stencil is null, then set the input port to null
#if VTK_MAJOR_VERSION >= 6
  this->SetInputDataInternal(2, stencil);
#else
  this->SetNthInputConnection(2, 0,
    (stencil ? stencil->GetProducerPort() : 0));
#endif
}

//----------------------------------------------------------------------------
vtkImageStencilData* vtkImageRegistration::GetSourceImageStencil()
{
  if (this->GetNumberOfInputConnections(2) < 1)
    {
    return NULL;
    }
  return vtkImageStencilData::SafeDownCast(
    this->GetExecutive()->GetInputData(2, 0));
}

//--------------------------------------------------------------------------
namespace {

void vtkTransformRotation(
  vtkTransform *transform, double rx, double ry, double rz)
{
  // angle is the norm of the parameters,
  // axis is the unit vector from the parameters
  double theta2 = rx*rx + ry*ry + rz*rz;
  if (theta2 > 0)
    {
    // compute the quaternion wxyz from angle and vector,
    // then use the quaternion to compute the matrix
    double theta = sqrt(theta2);
    double n = sin(0.5*theta);
    double w = cos(0.5*theta);
    double f = n/theta;
    double x = f*rx;
    double y = f*ry;
    double z = f*rz;

    double ww = w*w;
    double wx = w*x;
    double wy = w*y;
    double wz = w*z;

    double xx = x*x;
    double yy = y*y;
    double zz = z*z;

    double xy = x*y;
    double xz = x*z;
    double yz = y*z;

    double s = ww - xx - yy - zz;

    double matrix[16];
    matrix[0] = xx*2 + s;
    matrix[1] = (xy - wz)*2;
    matrix[2] = (xz + wy)*2;
    matrix[3] = 0.0;
    matrix[4] = (xy + wz)*2;
    matrix[5] = yy*2 + s;
    matrix[6] = (yz - wx)*2;
    matrix[7] = 0.0;
    matrix[8] = (xz - wy)*2;
    matrix[9] = (yz + wx)*2;
    matrix[10] = zz*2 + s;
    matrix[11] = 0.0;
    matrix[12] = 0.0;
    matrix[13] = 0.0;
    matrix[14] = 0.0;
    matrix[15] = 1.0;

    transform->Concatenate(matrix);
    }
}

void vtkSetTransformParameters(vtkImageRegistrationInfo *registrationInfo)
{
  vtkFunctionMinimizer* optimizer = registrationInfo->Optimizer;
  vtkTransform* transform = registrationInfo->Transform;
  vtkMatrix4x4 *initialMatrix = registrationInfo->InitialMatrix;
  int transformType = registrationInfo->TransformType;
  int transformDim = registrationInfo->TransformDimensionality;

  int pcount = 0;

  double tx = optimizer->GetParameterValue(pcount++);
  double ty = optimizer->GetParameterValue(pcount++);
  double tz = 0.0;
  if (transformDim > 2)
    {
    tz = optimizer->GetParameterValue(pcount++);
    }

  double rx = 0.0;
  double ry = 0.0;
  double rz = 0.0;

  if (transformType > vtkImageRegistration::Translation)
    {
    if (transformDim > 2)
      {
      rx = optimizer->GetParameterValue(pcount++);
      ry = optimizer->GetParameterValue(pcount++);
      }
    rz = optimizer->GetParameterValue(pcount++);
    }

  double sx = 1.0;
  double sy = 1.0;
  double sz = 1.0;

  if (transformType > vtkImageRegistration::Rigid)
    {
    sx = exp(optimizer->GetParameterValue(pcount++));
    sy = sx;
    if (transformDim > 2)
      {
      sz = sx;
      }
    }

  if (transformType > vtkImageRegistration::Similarity)
    {
    if (transformDim > 2)
      {
      sx = sz*exp(optimizer->GetParameterValue(pcount++));
      }
    sy = sz*exp(optimizer->GetParameterValue(pcount++));
    }

  bool scaledAtSource =
    (transformType == vtkImageRegistration::ScaleSourceAxes);

  double qx = 0.0;
  double qy = 0.0;
  double qz = 0.0;

  if (transformType >= vtkImageRegistration::Affine)
    {
    if (transformDim > 2)
      {
      qx = optimizer->GetParameterValue(pcount++);
      qy = optimizer->GetParameterValue(pcount++);
      }
    qz = optimizer->GetParameterValue(pcount++);
    }

  double *center = registrationInfo->Center;

  transform->Identity();
  transform->PostMultiply();
  transform->Translate(-center[0], -center[1], -center[2]);
  if (scaledAtSource)
    {
    vtkTransformRotation(transform, -qx, -qy, -qz);
    transform->Scale(sx, sy, sz);
    vtkTransformRotation(transform, qx, qy, qz);
    transform->Concatenate(initialMatrix);
    vtkTransformRotation(transform, rx, ry, rz);
    }
  else
    {
    vtkTransformRotation(transform, rx, ry, rz);
    transform->Concatenate(initialMatrix);
    vtkTransformRotation(transform, -qx, -qy, -qz);
    transform->Scale(sx, sy, sz);
    vtkTransformRotation(transform, qx, qy, qz);
    }
  transform->Translate(center[0], center[1], center[2]);
  transform->Translate(tx,ty,tz);
}

//--------------------------------------------------------------------------
void vtkEvaluateFunction(void * arg)
{
  vtkImageRegistrationInfo *registrationInfo =
    static_cast<vtkImageRegistrationInfo*>(arg);

  vtkFunctionMinimizer *optimizer = registrationInfo->Optimizer;
  vtkImageSimilarityMetric *metric = registrationInfo->Metric;

  vtkSetTransformParameters(registrationInfo);

  registrationInfo->Metric->Update();

  optimizer->SetFunctionValue(metric->GetCost());

  if (registrationInfo->MetricValues)
    {
    registrationInfo->MetricValues->InsertNextValue(metric->GetValue());
    }
  if (registrationInfo->CostValues)
    {
    registrationInfo->CostValues->InsertNextValue(metric->GetCost());
    }
  if (registrationInfo->ParameterValues)
    {
    double parameters[12];
    int n = optimizer->GetNumberOfParameters();
    for (int i = 0; i < n; i++)
      {
      parameters[i] = optimizer->GetParameterValue(i);
      }
    registrationInfo->ParameterValues->InsertNextTuple(parameters);
    }

  registrationInfo->NumberOfEvaluations++;
}

} // end anonymous namespace

//--------------------------------------------------------------------------
void vtkImageRegistration::ComputeImageRange(
  vtkImageData *data, vtkImageStencilData *stencil, double range[2])
{
  vtkImageHistogramStatistics *hist =
    vtkImageHistogramStatistics::New();
  hist->SET_STENCIL_DATA(stencil);
  hist->SET_INPUT_DATA(data);
  hist->SetActiveComponent(0);
  hist->Update();

  range[0] = hist->GetMinimum();
  range[1] = hist->GetMaximum();

  if (range[0] >= range[1])
    {
    range[1] = range[0] + 1.0;
    }

  hist->SET_INPUT_DATA(NULL);
  hist->Delete();
}

//--------------------------------------------------------------------------
void vtkImageRegistration::Initialize(vtkMatrix4x4 *matrix)
{
  // update our inputs
  this->Update();

  int transformDim = this->TransformDimensionality;
  if (transformDim < 2) { transformDim = 2; }
  if (transformDim > 3) { transformDim = 3; }

  vtkImageData *targetImage = this->GetTargetImage();
  vtkImageData *sourceImage = this->GetSourceImage();

  if (targetImage == NULL || sourceImage == NULL)
    {
    vtkErrorMacro("Initialize: Input images are not set");
    return;
    }

  // get the source image center
  double bounds[6];
  double center[3];
  double size[3];
  sourceImage->GetBounds(bounds);
  center[0] = 0.5*(bounds[0] + bounds[1]);
  center[1] = 0.5*(bounds[2] + bounds[3]);
  center[2] = 0.5*(bounds[4] + bounds[5]);
  size[0] = (bounds[1] - bounds[0]);
  size[1] = (bounds[3] - bounds[2]);
  size[2] = (bounds[5] - bounds[4]);

  vtkTransform *transform = this->Transform;
  vtkMatrix4x4 *initialMatrix = this->InitialTransformMatrix;

  // create an initial transform
  initialMatrix->Identity();
  transform->Identity();

  // the initial translation
  double tx = 0.0;
  double ty = 0.0;
  double tz = 0.0;

  // initialize from the supplied matrix
  if (matrix)
    {
    // move the translation into tx, ty, tz variables
    tx = matrix->Element[0][3];
    ty = matrix->Element[1][3];
    tz = matrix->Element[2][3];

    // move rotation/scale/shear into the InitialTransformMatrix
    initialMatrix->DeepCopy(matrix);
    initialMatrix->Element[0][3] = 0.0;
    initialMatrix->Element[1][3] = 0.0;
    initialMatrix->Element[2][3] = 0.0;

    // adjust the translation for the transform centering
    double scenter[4];
    scenter[0] = center[0];
    scenter[1] = center[1];
    scenter[2] = center[2];
    scenter[3] = 1.0;

    initialMatrix->MultiplyPoint(scenter, scenter);

    tx -= center[0] - scenter[0];
    ty -= center[1] - scenter[1];
    tz -= center[2] - scenter[2];
    }

  if (this->InitializerType == vtkImageRegistration::Centered)
    {
    // set an initial translation from one image center to the other image center
    double tbounds[6];
    double tcenter[3];
    targetImage->GetBounds(tbounds);
    tcenter[0] = 0.5*(tbounds[0] + tbounds[1]);
    tcenter[1] = 0.5*(tbounds[2] + tbounds[3]);
    tcenter[2] = 0.5*(tbounds[4] + tbounds[5]);

    tx = tcenter[0] - center[0];
    ty = tcenter[1] - center[1];
    tz = tcenter[2] - center[2];
    }

  if (transformDim <= 2)
    {
    center[2] = 0.0;
    tz = 0.0;
    }

  // do the setup for mutual information
  double sourceImageRange[2];
  double targetImageRange[2];
  sourceImageRange[0] = this->SourceImageRange[0];
  sourceImageRange[1] = this->SourceImageRange[1];
  targetImageRange[0] = this->TargetImageRange[0];
  targetImageRange[1] = this->TargetImageRange[1];

  if (this->MetricType == vtkImageRegistration::MutualInformation ||
      this->MetricType == vtkImageRegistration::NormalizedMutualInformation)
    {
    if (sourceImageRange[0] >= sourceImageRange[1])
      {
      this->ComputeImageRange(sourceImage, this->GetSourceImageStencil(),
        sourceImageRange);
      }
    if (targetImageRange[0] >= targetImageRange[1])
      {
      this->ComputeImageRange(targetImage, NULL,
        targetImageRange);
      }

    if (this->InterpolatorType == vtkImageRegistration::Nearest &&
        this->JointHistogramSize[0] <= 256 &&
        this->JointHistogramSize[1] <= 256)
      {
      // If nearest-neighbor interpolation is used, then the image instensity
      // can be quantized during initialization, instead of being done at each
      // iteration of the registration.

      double sourceScale = ((this->JointHistogramSize[1] - 1)/
        (sourceImageRange[1] - sourceImageRange[0]));
      // The "0.5/sourceScale" causes the vtkImageShiftScale filter to
      // round the value, instead of truncating it.
      double sourceShift = (-sourceImageRange[0] + 0.5/sourceScale);

      vtkImageShiftScale *sourceQuantizer = this->SourceImageTypecast;
      sourceQuantizer->SET_INPUT_DATA(sourceImage);
      sourceQuantizer->SetOutputScalarTypeToUnsignedChar();
      sourceQuantizer->ClampOverflowOn();
      sourceQuantizer->SetShift(sourceShift);
      sourceQuantizer->SetScale(sourceScale);
      sourceQuantizer->Update();
      sourceImage = sourceQuantizer->GetOutput();

      double targetScale = ((this->JointHistogramSize[0] - 1)/
        (targetImageRange[1] - targetImageRange[0]));
      double targetShift = (-targetImageRange[0] + 0.5/targetScale);

      vtkImageShiftScale *targetQuantizer = this->TargetImageTypecast;
      targetQuantizer->SET_INPUT_DATA(targetImage);
      targetQuantizer->SetOutputScalarTypeToUnsignedChar();
      targetQuantizer->ClampOverflowOn();
      targetQuantizer->SetShift(targetShift);
      targetQuantizer->SetScale(targetScale);
      targetQuantizer->Update();
      targetImage = targetQuantizer->GetOutput();

      // the rescaled image range is now the histogram range
      targetImageRange[0] = 0;
      targetImageRange[1] = this->JointHistogramSize[0] - 1;
      sourceImageRange[0] = 0;
      sourceImageRange[1] = this->JointHistogramSize[1] - 1;
      }
    }

  // make sure source range is computed for CorrelationRatio
  if (this->MetricType == vtkImageRegistration::CorrelationRatio)
    {
    if (sourceImageRange[0] >= sourceImageRange[1])
      {
      this->ComputeImageRange(sourceImage, this->GetSourceImageStencil(),
        sourceImageRange);
      }
    }

  // apply b-spline prefilter if b-spline interpolator is used
  if (this->Interpolato