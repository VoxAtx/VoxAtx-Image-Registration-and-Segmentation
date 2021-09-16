
/*=========================================================================

Program:   Atamai Image Registration and Segmentation
Module:    register.cxx

   This software is distributed WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=========================================================================*/

// Image registration is done first on a blurred, low-resolution version of
// the image before being done on the full resolution image, and is also
// done first with no interpolation before being done with linear interpolation.
// This multi-stage approach increases the robustness and often the speed of
// the registration.

#include <vtkSmartPointer.h>

#include <vtkImageReslice.h>
#include <vtkImageResize.h>
#include <vtkImageBSplineCoefficients.h>
#include <vtkImageBSplineInterpolator.h>
#include <vtkImageSincInterpolator.h>
#include <vtkImageHistogramStatistics.h>
#include <vtkImageThreshold.h>
#include <vtkImageCast.h>
#include <vtkROIStencilSource.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkIdTypeArray.h>
#include <vtkStringArray.h>
#include <vtkMath.h>
#include <vtkCommand.h>
#include <vtkMultiThreader.h>

#include <vtkMINCImageReader.h>
#include <vtkMINCImageWriter.h>
#include <vtkMNITransformReader.h>
#include <vtkMNITransformWriter.h>
#include <vtkDICOMImageReader.h>

#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageSlice.h>
#include <vtkImageStack.h>
#include <vtkImageResliceMapper.h>
#include <vtkImageProperty.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkJPEGWriter.h>

#include <vtkTimerLog.h>
#include <vtkVersion.h>

#include <vtksys/SystemTools.hxx>

#include "AIRSConfig.h"
#include "vtkITKXFMReader.h"
#include "vtkITKXFMWriter.h"
#include "vtkImageRegistration.h"
#include "vtkLabelInterpolator.h"

// optional readers
#ifdef AIRS_USE_DICOM
#define AIRS_USE_NIFTI
#include <vtkNIFTIReader.h>
#include <vtkNIFTIWriter.h>
#include <vtkDICOMReader.h>
#include <vtkDICOMSorter.h>
#include <vtkDICOMMRGenerator.h>
#include <vtkDICOMCTGenerator.h>
#include <vtkDICOMWriter.h>
#include <vtkDICOMMetaData.h>
#include <vtkGlobFileNames.h>
#endif

#include <vector>
#include <string>

// A macro to assist VTK 5 backwards compatibility
#if VTK_MAJOR_VERSION >= 6
#define SET_INPUT_DATA SetInputData
#define SET_STENCIL_DATA SetStencilData
#else
#define SET_INPUT_DATA SetInput
#define SET_STENCIL_DATA SetStencil
#endif

// Check for vtkImageReslice::SetOutputScalarType()
#if VTK_MAJOR_VERSION > 6 || (VTK_MAJOR_VERSION == 6 && VTK_MINOR_VERSION >= 2)
#define VTK_RESLICE_HAS_OUTPUT_SCALAR_TYPE
#endif

#if VTK_MAJOR_VERSION > 7 || (VTK_MAJOR_VERSION == 7 && VTK_MINOR_VERSION >= 0)
#define USE_SMP_THREADED_IMAGE_ALGORITHM
#endif

// parallel processing
enum { MultiThread = 1, ThreadPool = 2 };

// coord systems
enum { NativeCoords, DICOMCoords, NIFTICoords };

// file types
enum { DICOMImage, NIFTIImage, MINCImage,
         LastImageType = MINCImage,
       MNITransform, ITKTransform, CSVTransform, TXTTransform,
         LastTransformType = TXTTransform };

// internal methods for reading images, these methods read the image
// into the specified data object and also provide a matrix for converting
// the data coordinates into patient coordinates.
namespace {

void ComputeRange(vtkImageData *image, double range[2], double fill[2]);

// set the interpolator as requested
void SetInterpolator(vtkImageReslice *reslice, int interpolator)
{
  vtkSmartPointer<vtkImageBSplineInterpolator> bsplineInterpolator =
    vtkSmartPointer<vtkImageBSplineInterpolator>::New();

  vtkSmartPointer<vtkImageSincInterpolator> sincInterpolator =
    vtkSmartPointer<vtkImageSincInterpolator>::New();
  sincInterpolator->SetWindowFunctionToBlackman();

  vtkSmartPointer<vtkLabelInterpolator> labelInterpolator =
    vtkSmartPointer<vtkLabelInterpolator>::New();

  switch (interpolator)
    {
    case vtkImageRegistration::Nearest:
      reslice->SetInterpolationModeToNearestNeighbor();
      break;
    case vtkImageRegistration::Linear:
      reslice->SetInterpolationModeToLinear();
      break;
    case vtkImageRegistration::Cubic:
      reslice->SetInterpolationModeToCubic();
      break;
    case vtkImageRegistration::BSpline:
      reslice->SetInterpolator(bsplineInterpolator);
      break;
    case vtkImageRegistration::Sinc:
      reslice->SetInterpolator(sincInterpolator);
      break;
    case vtkImageRegistration::ASinc:
      sincInterpolator->AntialiasingOn();
      reslice->SetInterpolator(sincInterpolator);
      break;
    case vtkImageRegistration::Label:
      reslice->SetInterpolator(labelInterpolator);
      break;
    }
}

// use file extension to guess file type
int GuessFileType(const char *filename)
{
  size_t n = strlen(filename);

  if (n > 4 && strcmp(&filename[n-4], ".txt") == 0)
    {
    int ftype = TXTTransform;
    ifstream infile(filename);
    if (infile.good())
      {
      char firstline[32];
      memset(firstline, '\0', 32);
      infile.getline(firstline, 32);
      if (strncmp(firstline, "#Insight Transform File V1.0", 28) == 0)
        {
        ftype = ITKTransform;
        }
      }
    infile.close();
    return ftype;
    }
  if (n > 4 && strcmp(&filename[n-4], ".xfm") == 0)
    {
    return MNITransform;
    }
  if (n > 4 && strcmp(&filename[n-4], ".tfm") == 0)
    {
    return ITKTransform;
    }
  if (n > 4 && strcmp(&filename[n-4], ".mat") == 0)
    {
    return TXTTransform;
    }
  if (n > 4 && strcmp(&filename[n-4], ".csv") == 0)
    {
    return CSVTransform;
    }

  if (n > 4 && strcmp(&filename[n-4], ".mnc") == 0)
    {
    return MINCImage;
    }
  if ((n > 4 && strcmp(&filename[n-4], ".nii") == 0) ||
      (n > 7 && strcmp(&filename[n-7], ".nii.gz") == 0) ||
      (n > 4 && strcmp(&filename[n-4], ".hdr") == 0) ||
      (n > 4 && strcmp(&filename[n-4], ".img") == 0) ||
      (n > 7 && strcmp(&filename[n-7], ".img.gz") == 0))
    {
    return NIFTIImage;
    }

  return DICOMImage;
}

#ifdef AIRS_USE_DICOM
vtkDICOMReader *ReadDICOMImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *directoryName,
  int coordSystem)
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
    vtkSmartPointer<vtkDICOMSorter> sorter =
      vtkSmartPointer<vtkDICOMSorter>::New();
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

  if (coordSystem == NIFTICoords)
    {
    reader->SetMemoryRowOrderToBottomUp();
    }
  else
    {
    reader->SetMemoryRowOrderToFileNative();
    }

  reader->UpdateInformation();
  if (reader->GetErrorCode())
    {
    exit(1);
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

  return reader;
}

void WriteDICOMImage(
  vtkImageReader2 *sourceReader, vtkImageReader2 *targetReader,
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *directoryName,
  int vtkNotUsed(coordSystem))
{
  if (vtksys::SystemTools::FileExists(directoryName))
    {
    if (!vtksys::SystemTools::FileIsDirectory(directoryName))
      {
      fprintf(stderr, "option -o must give a DICOM directory, not a file.\n");
      exit(1);
      }
    }
  else if (!vtksys::SystemTools::MakeDirectory(directoryName))
    {
    fprintf(stderr, "Cannot create directory: %s\n", directoryName);
    exit(1);
    }

  // get the meta data
  vtkDICOMReader *reader = vtkDICOMReader::SafeDownCast(targetReader);
  vtkDICOMReader *reader2 = vtkDICOMReader::SafeDownCast(sourceReader);

  vtkSmartPointer<vtkDICOMMetaData> meta =
    vtkSmartPointer<vtkDICOMMetaData>::New();

  if (reader)
    {
    // copy the bulk of the meta data from the target image
    meta->DeepCopy(reader->GetMetaData());
    meta->SetAttributeValue(DC::SeriesNumber,
      meta->GetAttributeValue(DC::SeriesNumber).AsUnsignedInt() + 1000);
    std::string seriesDescription =
      meta->GetAttributeValue(DC::SeriesDescription).AsString() + " REG";
    if (seriesDescription.size() < 64)
      {
      meta->SetAttributeValue(DC::SeriesDescription, seriesDescription);
      }
    }
  if (reader2)
    {
    // set the frame of reference from the source image
    meta->SetAttributeValue(DC::FrameOfReferenceUID,
      reader2->GetMetaData()->GetAttributeValue(
      DC::FrameOfReferenceUID));
    }

  // make the generator
  vtkSmartPointer<vtkDICOMMRGenerator> mrgenerator =
    vtkSmartPointer<vtkDICOMMRGenerator>::New();
  vtkSmartPointer<vtkDICOMCTGenerator> ctgenerator =
    vtkSmartPointer<vtkDICOMCTGenerator>::New();
  vtkDICOMGenerator *generator = 0;
  if (reader)
    {
    std::string SOPClass =
      meta->GetAttributeValue(DC::SOPClassUID).AsString();
    if (SOPClass == "1.2.840.10008.5.1.4.1.1.2" ||
        SOPClass == "1.2.840.10008.5.1.4.1.1.2.1" ||
        SOPClass == "1.2.840.10008.5.1.4.1.1.2.2")
      {
      generator = ctgenerator;
      }
    else if (SOPClass == "1.2.840.10008.5.1.4.1.1.4" ||
             SOPClass == "1.2.840.10008.5.1.4.1.1.4.1" ||
             SOPClass == "1.2.840.10008.5.1.4.1.1.4.4")
      {
      generator = mrgenerator;
      }
    }

  // prepare the writer to write the image
  vtkSmartPointer<vtkDICOMWriter> writer =
    vtkSmartPointer<vtkDICOMWriter>::New();
  if (generator)
    {
    writer->SetGenerator(generator);
    }
  writer->SetMetaData(meta);
  writer->SetFilePrefix(directoryName);
  writer->SetFilePattern("%s/IM-0001-%04.4d.dcm");
  writer->TimeAsVectorOn();
  if (reader)
    {
    if (reader->GetTimeDimension() > 1)
      {
      writer->SetTimeDimension(reader->GetTimeDimension());
      writer->SetTimeSpacing(reader->GetTimeSpacing());
      }
    if (reader->GetRescaleSlope() > 0)
      {
      writer->SetRescaleSlope(reader->GetRescaleSlope());
      writer->SetRescaleIntercept(reader->GetRescaleIntercept());
      }
    writer->SetMemoryRowOrder(reader->GetMemoryRowOrder());
    }
  writer->SET_INPUT_DATA(data);
  writer->SetPatientMatrix(matrix);
  writer->Write();
}

#else

vtkDICOMImageReader *ReadDICOMImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *directoryName,
  int coordSystem)
{
  // read the image
  vtkDICOMImageReader *reader = vtkDICOMImageReader::New();

  reader->SetDirectoryName(directoryName);
  reader->Update();
  if (reader->GetErrorCode())
    {
    exit(1);
    }

  vtkSmartPointer<vtkImageData> image = reader->GetOutput();

  if (coordSystem != NIFTICoords)
    {
    // the reader flips the image and reverses the ordering, so undo these
    vtkSmartPointer<vtkImageReslice> flip =
      vtkSmartPointer<vtkImageReslice>::New();

    flip->SetInputConnection(reader->GetOutputPort());
    flip->SetResliceAxesDirectionCosines(
      1,0,0, 0,-1,0, 0,0,-1);
    flip->Update();

    image = flip->GetOutput();
    }

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

  if (coordSystem == NIFTICoords)
    {
    double spacing[3], origin[3];
    int extent[6];
    image->GetSpacing(spacing);
    image->GetOrigin(origin);
    image->GetExtent(extent);
    // account fo the y and z flips
    double point[4];
    point[0] = origin[0] + spacing[0]*extent[0];
    point[1] = origin[1] + spacing[1]*extent[3];
    point[2] = origin[2] + spacing[2]*extent[5];
    point[3] = 1.0;
    matrix->MultiplyPoint(point, point);
    for (int j = 0; j < 3; j++)
      {
      matrix->Element[j][1] = -matrix->Element[j][1];
      matrix->Element[j][2] = -matrix->Element[j][2];
      matrix->Element[j][3] = point[j];
      }
    // do the DICOM to NIFTI coord conversion
    for (int k = 0; k < 4; k++)
      {
      matrix->Element[0][k] = -matrix->Element[0][k];
      matrix->Element[1][k] = -matrix->Element[1][k];
      }
    }

  matrix->Modified();

  return reader;
}
#endif

vtkMINCImageReader *ReadMINCImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *fileName,
  int coordSystem)
{
  // read the image
  vtkMINCImageReader *reader = vtkMINCImageReader::New();

  reader->SetFileName(fileName);
  reader->Update();
  if (reader->GetErrorCode())
    {
    exit(1);
    }

  vtkSmartPointer<vtkImageData> image = reader->GetOutput();

  if (coordSystem == DICOMCoords)
    {
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

    image = flip->GetOutput();
    }

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());

  if (coordSystem == DICOMCoords)
    {
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
  else
    {
    matrix->DeepCopy(reader->GetDirectionCosines());
    }

  return reader;
}

void WriteMINCImage(
  vtkImageReader2 *vtkNotUsed(sourceReader),
  vtkImageReader2 *vtkNotUsed(targetReader),
  vtkImageData *data, vtkMatrix4x4 *vtkNotUsed(matrix), const char *fileName,
  int vtkNotUsed(coordSystem))
{
  fprintf(stderr, "Writing MINC images is not supported yet, "
          "the output file will have incorrect information\n");
  vtkSmartPointer<vtkMINCImageWriter> writer =
    vtkSmartPointer<vtkMINCImageWriter>::New();
  writer->SetFileName(fileName);
  writer->SET_INPUT_DATA(data);
  // the input matrix must be converted
  //writer->SetDirectionCosines(matrix);
  writer->Write();
}

#ifdef AIRS_USE_NIFTI
vtkNIFTIReader *ReadNIFTIImage(
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *fileName,
  int coordSystem)
{
  // read the image
  vtkNIFTIReader *reader = vtkNIFTIReader::New();

  reader->SetFileName(fileName);
  reader->Update();
  if (reader->GetErrorCode())
    {
    exit(1);
    }

  vtkSmartPointer<vtkImageData> image = reader->GetOutput();

  if (coordSystem == DICOMCoords)
    {
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

    image = flip->GetOutput();
    }

  // get the data
  data->CopyStructure(image);
  data->GetPointData()->PassData(image->GetPointData());

  // get the SForm or QForm matrix if present
  static double nMatrix[16] =
    { 1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
  if (reader->GetQFormMatrix())
    {
    vtkMatrix4x4::DeepCopy(nMatrix, reader->GetQFormMatrix());
    }
  else if (reader->GetSFormMatrix())
    {
    vtkMatrix4x4::DeepCopy(nMatrix, reader->GetSFormMatrix());
    }

  if (coordSystem == DICOMCoords)
    {
    // generate the matrix, but modify to use DICOM coords
    static double xyFlipMatrix[16] =
      { -1, 0, 0, 0,  0, -1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1 };
    // correct for the flip that was done earlier
    vtkMatrix4x4::Multiply4x4(nMatrix, xyFlipMatrix, *matrix->Element);
    // do the left/right, up/down dicom-to-minc transformation
    vtkMatrix4x4::Multiply4x4(xyFlipMatrix, *matrix->Element, *matrix->Element);
    matrix->Modified();
    }
  else
    {
    matrix->DeepCopy(nMatrix);
    }

  return reader;
}

void WriteNIFTIImage(
  vtkImageReader2 *sourceReader, vtkImageReader2 *targetReader,
  vtkImageData *data, vtkMatrix4x4 *matrix, const char *fileName,
  int vtkNotUsed(coordSystem))
{
  vtkNIFTIReader *sreader = vtkNIFTIReader::SafeDownCast(sourceReader);
  vtkNIFTIReader *treader = vtkNIFTIReader::SafeDownCast(targetReader);

  vtkSmartPointer<vtkNIFTIWriter> writer =
    vtkSmartPointer<vtkNIFTIWriter>::New();
  if (treader)
    {
    writer->SetNIFTIHeader(treader->GetNIFTIHeader());
    if (treader->GetTimeDimension() > 1)
      {
      writer->SetTimeDimension(treader->GetTimeDimension());
      writer->SetTimeSpacing(treader->GetTimeSpacing());
      }
    }
  if (sreader)
    {
    if (sreader->GetQFac() < 0)
      {
      writer->SetQFac(-1.0);
      }
    }
  writer->SET_INPUT_DATA(data);
  writer->SetQFormMatrix(matrix);
  writer->SetSFormMatrix(matrix);
  writer->SetFileName(fileName);
  writer->Write();
}

#endif /* AIRS_USE_NIFTI */

vtkImageReader2 *ReadImage(
  vtkImageData *image, vtkMatrix4x4 *matrix, double vrange[2],
  const char *filename, int coordSystem, int interpolator)
{
  int t = GuessFileType(filename);
  vtkImageReader2 *reader = 0;

  if (t == MINCImage)
    {
    reader = ReadMINCImage(image, matrix, filename, coordSystem);
    }
  else if (t == NIFTIImage)
    {
#ifdef AIRS_USE_NIFTI
    reader = ReadNIFTIImage(image, matrix, filename, coordSystem);
#else
    fprintf(stderr, "NIFTI files are not supported.\n");
    exit(1);
#endif
    }
  else
    {
    reader = ReadDICOMImage(image, matrix, filename, coordSystem);
    }