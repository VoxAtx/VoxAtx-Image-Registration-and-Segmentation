/*=========================================================================

Program:   Atamai Image Registration and Segmentation
Module:    DifferenceRegistration.cxx

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

// This example registers two images (assumed to be of the same patient)
// and then subtracts first image from the second.

#include "AIRSConfig.h"

#include <vtkSmartPointer.h>

#include <vtkImageReslice.h>
#include <vtkImageResize.h>
#include <vtkImageSincInterpolator.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkMath.h>
#include <vtkGlobFileNames.h>
#include <vtkIntArray.h>
#include <vtkStringArray.h>

#include <vtkMINCImageReader.h>
#include <vtkDICOMImageReader.h>
#include <vtkMNITransformWriter.h>

// optional readers
#ifdef AIRS_USE_DICOM
#define AIRS_USE_NIFTI
#include <vtkNIFTIReader.h>
#include <vtkNIFTIWriter.h>
#include <vtkDICOMReader.h>
#include <vtkDICOMFileSorter.h>
#include <vtkDICOMMetaData.h>
#endif

#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageSlice.h>
#include <vtkImageStack.h>
#include <vtkImageResliceMapper.h>
#include <vtkImageProperty.h>
#include <vtkImageHistogramStatistics.h>
#include <vtkImageMathematics.h>
#include <vtkImageShiftScale.h>

#include <vtkTimerLog.h>
#include <vtkVersion.h>

#include <vtkImageRegistration.h>

#include <vtksys/SystemTools.hxx>

// A macro to assist VTK 5 backwards compatibility
#if VTK_MAJOR_VERSION >= 6
#define SET_INPUT_DATA SetInputData
#else
#define SET_INPUT_DATA SetInput
#endif

// internal methods for reading images, these methods read the image
// into the specified data object and also provide a matrix for converting
// the data coordinates into patient coordinates.
namespace {

#ifdef AIRS_USE_DICOM

void ReadDICOMImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *directoryName)
{
  vtkDICOMReader *reader = vtkDICOMReader::New();

  bool singleFile = true;
  if (vtksys::SystemTools::FileIsDirectory(directoryName))
    {
    // get all the DICOM files in the directory
    singleFile = false;
    std::string dirString = directoryName;
    vtksys::SystemTools::ConvertToUnixSlashes(dirString);
    vtkSmartPointer<vtkGlobFileNames> glob =
      vtkSmartPointer<vtkGlobFileNames>::New();
    glob->SetDirectory(dirString.c_str());
    glob->AddFileNames("*");

    // sort the files
    vtkSmartPointer<vtkDICOMFileSorter> sorter =
      vtkSmartPointer<vtkDICOMFileSorter>::New();
    sorter->SetInputFileNames(glob->GetFileNames());
    sorter->Update();

    if (sorter->GetNumberOfSeries() == 0)
      {
      fprintf(stderr, "Folder contains no DICOM files: %s\n", directoryName);
      exit(1);
      }
    else if (sorter->GetNumberOfSeries() > 1)
      {
      fprintf(stderr, "Folder contains more than one DICOM series: %s\n",
              directoryName);
      exit(1);
      }
    reader->SetFileNames(sorter->GetFileNamesForSeries(0));
    }
  else
    {
    // was given a single file instead of a directory
    reader->SetFileName(directoryName);
    }

  // For NIfTI coordinate system, use BottomUp
  reader->SetMemoryRowOrderToFileNative();

  reader->UpdateInformation();
  if (reader->GetErrorCode())
    {
    exit(1);
    }

  vtkStringArray *stackArray = reader->GetStackIDs();
  vtkIdType numStacks = stackArray->GetNumberOfValues();
  for (vtkIdType stackId = 0; stackId+1 < numStacks; stackId++)
    {
    // Find the first stack that has more than one image
    if (reader->GetFileIndexArray()->GetNumberOfTuples() > 1)
      {
      break;
      }
    reader->SetDesiredStackID(stackArray->GetValue(stackId+1));
    reader->UpdateInformation();
    }

  if (!singleFile)
    {
    // when reading images, only read 1st component if the
    // image has multiple components or multiple time points
    vtkIntArray *fileArray = reader->GetFileIndexArray();

    // create a filtered list of files
    vtkSmartPointer<vtkStringArray> fileNames =
      vtkSmartPointer<vtkStringArray>::New();
    vtkIdType n = fileArray->GetNumberOfTuples();
    for (vtkIdType i = 0; i < n; i++)
      {
      std::string newFileName =
        reader->GetFileNames()->GetValue(fileArray->GetComponent(i, 0));
      bool alreadyThere = false;
      vtkIdType m = fileNames->GetNumberOfTuples();
      for (vtkIdType j = 0; j < m; j++)
        {
        if (newFileName == fileNames->GetValue(j))
          {
          alreadyThere = true;
          break;
          }
        }
      if (!alreadyThere)
        {
        fileNames->InsertNextValue(newFileName);
        }
      }
    reader->SetFileNames(fileNames);
    }
  reader->SetDesiredTimeIndex(0);

  reader->Update();
  if (reader->GetErrorCode())
    {
    exit(1);
    }

  vtkImageData *image = reader->GetOutput();

  // get the data
  data->CopyStructure(image);
  data->GetPoint