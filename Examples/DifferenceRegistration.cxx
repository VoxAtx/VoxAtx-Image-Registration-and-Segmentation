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
  data->GetPointData()->PassData(image->GetPointData());

  // get the matrix
  matrix->DeepCopy(reader->GetPatientMatrix());

  return;
}

#else /* AIRS_USE_DICOM */

void ReadDICOMImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *directoryName)
{
  // read the image
  vtkSmartPointer<vtkDICOMImageReader> reader =
    vtkSmartPointer<vtkDICOMImageReader>::New();

  reader->SetDirectoryName(directoryName);
  reader->Update();

  // the reader flips the image and reverses the ordering, so undo these
  vtkSmartPointer<vtkImageReslice> flip =
    vtkSmartPointer<vtkImageReslice>::New();

  flip->SetInputConnection(reader->GetOutputPort());
  flip->SetResliceAxesDirectionCosines(
    1,0,0, 0,-1,0, 0,0,-1);
  flip->Update();

  vtkImageData *image = flip->GetOutput();

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());
  data->SetOrigin(0,0,0);

  // generate the matrix
  float *position = reader->GetImagePositionPatient();
  float *orientation = reader->GetImageOrientationPatient();
  float *xdir = &orientation[0];
  float *ydir = &orientation[3];
  float zdir[3];
  vtkMath::Cross(xdir, ydir, zdir);

  for (int i = 0; i < 3; i++)
    {
    matrix->Element[i][0] = xdir[i];
    matrix->Element[i][1] = ydir[i];
    matrix->Element[i][2] = zdir[i];
    matrix->Element[i][3] = position[i];
    }
  matrix->Element[3][0] = 0;
  matrix->Element[3][1] = 0;
  matrix->Element[3][2] = 0;
  matrix->Element[3][3] = 1;
  matrix->Modified();
}

#endif /* AIRS_USE_DICOM */

void ReadMINCImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *fileName)
{
  // read the image
  vtkSmartPointer<vtkMINCImageReader> reader =
    vtkSmartPointer<vtkMINCImageReader>::New();

  reader->SetFileName(fileName);
  reader->Update();

  double spacing[3];
  reader->GetOutput()->GetSpacing(spacing);
  spacing[0] = fabs(spacing[0]);
  spacing[1] = fabs(spacing[1]);
  spacing[2] = fabs(spacing[2]);

  // flip the image rows into a DICOM-style ordering
  vtkSmartPointer<vtkImageReslice> flip =
    vtkSmartPointer<vtkImageReslice>::New();

  flip->SetInputConnection(reader->GetOutputPort());
  flip->SetResliceAxesDirectionCosines(
    -1,0,0, 0,-1,0, 0,0,1);
  flip->SetOutputSpacing(spacing);
  flip->Update();

  vtkImageData *image = flip->GetOutput();

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());

  // generate the matrix, but modify to use DICOM coords
  static double xyFlipMatrix[16] =
    { -1, 0, 0, 0,  0, -1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
  // correct for the flip that was done earlier
  vtkMatrix4x4::Multiply4x4(*reader->GetDirectionCosines()->Element,
                            xyFlipMatrix, *matrix->Element);
  // do the left/right, up/down dicom-to-minc transformation
  vtkMatrix4x4::Multiply4x4(xyFlipMatrix, *matrix->Element, *matrix->Element);
  matrix->Modified();
}

#ifdef AIRS_USE_NIFTI
void ReadNIFTIImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *fileName)
{
  // read the image
  vtkSmartPointer<vtkNIFTIReader> reader =
    vtkSmartPointer<vtkNIFTIReader>::New();

  reader->SetFileName(fileName);
  reader->Update();

  double spacing[3];
  reader->GetOutput()->GetSpacing(spacing);
  spacing[0] = fabs(spacing[0]);
  spacing[1] = fabs(spacing[1]);
  spacing[2] = fabs(spacing[2]);

  // flip the image rows into a DICOM-style ordering
  vtkSmartPointer<vtkImageReslice> flip =
    vtkSmartPointer<vtkImageReslice>::New();

  flip->SetInputConnection(reader->GetOutputPort());
  flip->SetResliceAxesDirectionCosines(
    -1,0,0, 0,-1,0, 0,0,1);
  flip->SetOutputSpacing(spacing);
  flip->Update();

  //vtkImageData *image = flip->GetOutput();
  vtkImageData *image = reader->GetOutput();

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());

  // get the SForm or QForm matrix if present
  static double nMatrix[16] =
    { 1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
  if (reader->GetSFormMatrix())
    {
    vtkMatrix4x4::DeepCopy(nMatrix, reader->GetSFormMatrix());
    }
  else if (reader->GetQFormMatrix())
    {
    vtkMatrix4x4::DeepCopy(nMatrix, reader->GetQFormMatrix());
    }

  // generate the matrix, but modify to use DICOM coords
  static double xyFlipMatrix[16] =
    { -1, 0, 0, 0,  0, -1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
  // correct for the flip that was done earlier
  vtkMatrix4x4::Multiply4x4(nMatrix, xyFlipMatrix, *matrix->Element);
  // do the left/right, up/down dicom-to-minc transformation
  vtkMatrix4x4::Multiply4x4(xyFlipMatrix, *matrix->Element, *matrix->Element);
  matrix->Modified();
}
#endif /* AIRS_USE_NIFTI */

void SetViewFromMatrix(
  vtkRenderer *renderer,
  vtkInteractorStyleImage *istyle,
  vtkMatrix4x4 *matrix)
{
  istyle->SetCurrentRenderer(renderer);

  // This view assumes the data uses the DICOM Patient Coordinate System.
  // It provides a right-is-left view of axial and coronal images
  double viewRight[4] = { 1.0, 0.0, 0.0, 0.0 };
  double viewUp[4] = { 0.0, -1.0, 0.0, 0.0 };

  matrix->MultiplyPoint(viewRight, viewRight);
  matrix->MultiplyPoint(viewUp, viewUp);

  istyle->SetImageOrientation(viewRight, viewUp);
}

};

void printUsage(const char *cmdname)
{
    cout << "Usage 1: " << cmdname << " --nodisplay -o output.nii source.nii target.nii"
         << endl;
    cout << "Usage 2: " << cmdname << " --nodisplay -o output.nii dicomdir1/ dicomdir2/"
         << endl;
}

int main (int argc, char *argv[])
{
  if (argc < 3)
    {
    printUsage(argv[0]);
    return EXIT_FAILURE;
    }

  // -------------------------------------------------------
  // the files
  int argi = 1;
  const char *xfmfile = NULL;
  const char *outputfile = NULL;
  const char *sourcefile;
  const char *targetfile;
  bool display = true;

  if (strcmp(argv[argi], "--nodisplay") == 0)
    {
    display = false;
    argi++;
    }
  if (strcmp(argv[argi], "-o") == 0)
    {
    if (argc <= argi + 1)
      {
      cerr << argv[0] << " : missing output file after -o\n" << endl;
      return EXIT_FAILURE;
      }
    // is the output an xfm file or an image file?
    xfmfile = argv[argi + 1];
    argi += 2;
    size_t m = strlen(xfmfile);
    if (m < 4 || strcmp(&xfmfile[m-4], ".xfm") != 0)
      {
      // it isn't an .xfm file, assume that it is an image file
      outputfile = xfmfile;
      xfmfile = NULL;
      }
    }

  if (argc <= argi + 1)
    {
    printUsage(argv[0]);
    return EXIT_FAILURE;
    }
  targetfile = argv[argi];
  sourcefile = argv[argi + 1];

  // -------------------------------------------------------
  // parameters for registration

  int interpolatorType = vtkImageRegistration::Linear;
  double transformTolerance = 0.1; // tolerance on transformation result
  int numberOfBins = 64; // for Mattes' mutual information
  double initialBlurFactor = 4.0;

  // -------------------------------------------------------
  // load the images

  int n = 0;

  // Read the source image
  vtkSmartPointer<vtkImageData> sourceImage =
    vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkMatrix4x4> sourceMatrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  n = strlen(sourcefile);
  if (n > 4 && strcmp(&sourcefile[n-4], ".mnc") == 0)
    {
    ReadMINCImage(sourceImage, sourceMatrix, sourcefile);
    }
#ifdef AIRS_USE_NIFTI
  else if ((n > 4 && strcmp(&sourcefile[n-4], ".nii") == 0) ||
           (n > 7 && strcmp(&sourcefile[n-7], ".nii.gz") == 0))
    {
    ReadNIFTIImage(sourceImage, sourceMatrix, sourcefile);
    }
#endif
  else
    {
    ReadDICOMImage(sourceImage, sourceMatrix, sourcefile);
    }

  // Read the target image
  vtkSmartPointer<vtkImageData> targetImage =
    vtkSmartPointer<vtkImageData>::New();
  vtkSmartPointer<vtkMatrix4x4> targetMatrix =
    vtkSmartPointer<vtkMatrix4x4>::New();
  n = strlen(targetfile);
  if (n > 4 && strcmp(&targetfile[n-4], ".mnc") == 0)
    {
    ReadMINCImage(targetImage, targetMatrix, targetfile);
    }
#ifdef AIRS_USE_NIFTI
  else if ((n > 4 && strcmp(&targetfile[n-4], ".nii") == 0) ||
           (n > 7 && strcmp(&targetfile[n-7], ".nii.gz") == 0))
    {
    ReadNIFTIImage(targetImage, targetMatrix, targetfile);
    }
#endif
  else
    {
    ReadDICOMImage(targetImage, targetMatrix, targetfile);
    }

  // -------------------------------------------------------
  // display the images

  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> interactor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  vtkSmartPointer<vtkInteractorStyleImage> istyle =
    vtkSmartPointer<vtkInteractorStyleImage>::New();

  istyle->SetInteractionModeToImageSlicing();
  interactor->SetInteractorStyle(istyle);
  renderWindow->SetInteractor(interactor);
  renderWindow->AddRenderer(renderer);

  vtkSmartPointer<vtkImageSlice> sourceActor =
    vtkSmartPointer<vtkImageSlice>::New();
  vtkSmartPointer<vtkImageResliceMapper> sourceMapper =
    vtkSmartPointer<vtkImageResliceMapper>::New();
  vtkSmartPointer<vtkImageProperty> sourceProperty =
    vtkSmartPointer<vtkImageProperty>::New();

  sourceMapper->SET_INPUT_DATA(sourceImage);
  sourceMapper->SliceAtFocalPointOn();
  sourceMapper->SliceFacesCameraOn();
  sourceMapper->ResampleToScreenPixelsOff();

  double sourceRange[2];
  vtkSmartPointer<vtkImageHistogramStatistics> autoRange =
    vtkSmartPointer<vtkImageHistogramStatistics>::New();
  autoRange->SET_INPUT_DATA(sourceImage);
  autoRange->Update();
  autoRange->GetAutoRange(sourceRange);

  sourceProperty->SetInterpolationTypeToLinear();
  sourceProperty->SetColorWindow((sourceRange[1]-sourceRange[0]));
  sourceProperty->SetColorLevel(0.5*(sourceRange[0]+sourceRange[1]));
  sourceProperty->CheckerboardOn();
  sourceProperty->SetCheckerboardSpacing(40,40);

  sourceActor->SetMapper(sourceMapper);
  sourceActor->SetProperty(sourceProperty);
  sourceActor->SetUserMatrix(sourceMatrix);

  vtkSmartPointer<vtkImageSlice> targetActor =
    vtkSmartPointer<vtkImageSlice>::New();
  vtkSmartPointer<vtkImageResliceMapper> targetMapper =
    vtkSmartPointer<vtkImageResliceMapper>::New();
  vtkSmartPointer<vtkImageProperty> targetProperty =
    vtkSmartPointer<vtkImageProperty>::New();

  targetMapper->SET_INPUT_DATA(targetImage);
  targetMapper->SliceAtFocalPointOn();
  targetMapper->SliceFacesCameraOn();
  targetMapper->ResampleToS