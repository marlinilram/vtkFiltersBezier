#include "vtkBREPReader.h"

#include "vtkCachedStreamingDemandDrivenPipeline.h"
#include "vtkControlPointArray.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkNURBSPatchAdaptor.h"
#include "vtkObjectFactory.h"
#include "vtkPatchInterpolation.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"

#include <fstream>

vtkStandardNewMacro(vtkBREPReader);

vtkBREPReader::vtkBREPReader()
{
  this->FileName = NULL;

  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

vtkBREPReader::~vtkBREPReader()
{
  delete[] this->FileName;
  this->FileName = NULL;
}

// ----------------------------------------------------------------------
void vtkBREPReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "File Name: "
    << (this->FileName ? this->FileName : "(none)") << "\n";
}

// ----------------------------------------------------------------------
vtkUnstructuredGrid* vtkBREPReader::GetOutput()
{
  return this->GetOutput(0);
}

// ----------------------------------------------------------------------
vtkUnstructuredGrid* vtkBREPReader::GetOutput(int port)
{
  return vtkUnstructuredGrid::SafeDownCast(this->GetOutputDataObject(port));
}

// ----------------------------------------------------------------------
void vtkBREPReader::SetOutput(vtkDataObject* d)
{
  this->GetExecutive()->SetOutputData(0, d);
}

// ----------------------------------------------------------------------
vtkDataObject* vtkBREPReader::GetInput()
{
  return this->GetInput(0);
}

// ----------------------------------------------------------------------
vtkDataObject* vtkBREPReader::GetInput(int port)
{
  return this->GetExecutive()->GetInputData(port, 0);
}

// ----------------------------------------------------------------------
int vtkBREPReader::ProcessRequest(vtkInformation* request,
                                  vtkInformationVector** inputVector,
                                  vtkInformationVector* outputVector)
{
  // generate the data
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
    {
    return this->RequestData(request, inputVector, outputVector);
    }

  if(request->Has(vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
    {
    return this->RequestUpdateExtent(request, inputVector, outputVector);
    }

  // execute information
  if(request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
    {
    return this->RequestInformation(request, inputVector, outputVector);
    }

  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

// ----------------------------------------------------------------------
int vtkBREPReader::FillOutputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  return 1;
}

// ----------------------------------------------------------------------
int vtkBREPReader::FillInputPortInformation(
  int vtkNotUsed(port), vtkInformation* info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  return 1;
}

// ----------------------------------------------------------------------
int vtkBREPReader::RequestInformation(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* vtkNotUsed(outputVector))
{
  return 1;
}

// ----------------------------------------------------------------------
int vtkBREPReader::RequestUpdateExtent(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector,
  vtkInformationVector* vtkNotUsed(outputVector))
{
  int numInputPorts = this->GetNumberOfInputPorts();
  for (int i=0; i<numInputPorts; i++)
    {
    int numInputConnections = this->GetNumberOfInputConnections(i);
    for (int j=0; j<numInputConnections; j++)
      {
      vtkInformation* inputInfo = inputVector[i]->GetInformationObject(j);
      inputInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);
      }
    }
  return 1;
}

// ----------------------------------------------------------------------
void vtkBREPReader::SetInputData(vtkDataObject* input)
{
  this->SetInputData(0, input);
}

// ----------------------------------------------------------------------
void vtkBREPReader::SetInputData(int index, vtkDataObject* input)
{
  this->SetInputDataInternal(index, input);
}

// ----------------------------------------------------------------------
void vtkBREPReader::AddInputData(vtkDataObject* input)
{
  this->AddInputData(0, input);
}

// ----------------------------------------------------------------------
void vtkBREPReader::AddInputData(int index, vtkDataObject* input)
{
  this->AddInputDataInternal(index, input);
}

// ----------------------------------------------------------------------
int vtkBREPReader::RequestData(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** vtkNotUsed(inputVector),
  vtkInformationVector* outputVector)
{
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  vtkUnstructuredGrid* output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkSmartPointer<vtkPoints> bezierPoints = vtkSmartPointer<vtkPoints>::New();
  output->SetPoints(bezierPoints);
  output->Allocate();

  if (!this->FileName)
    {
    vtkErrorMacro(<< "A FileName must be specified.");
    return 0;
    }

  //FILE* in = fopen(this->FileName, "r");
  std::ifstream in(this->FileName);

  if (!in.is_open())
    {
    vtkErrorMacro(<< "File " << this->FileName << "not found");
    return 0;
    }

  vtkDebugMacro(<<"Reading file");

  std::string str_buffer;

  while (!in.eof())
  {
    in >> str_buffer;
    if (str_buffer == "Surfaces")
    {
      int num_surfaces;
      in >> num_surfaces;
      for (int i = 0; i < num_surfaces; ++i)
      {
        int surface_type;
        in >> surface_type;
        if (surface_type == 9)
        {
          //std::cout<<"A NURBS Patch!\n";
        }
        this->HandleSurface(in, surface_type, output);
        //std::cout<<"surface "<<i+1<<" processed.\n";
      }
    }
  }

  return 1;
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleSurface(std::ifstream& in, int shapeType, vtkDataObject* out)
{
  enum {
    PLANE = 1,
    CYLINDER = 2,
    CONE = 3,
    SPHERE = 4,
    TORUS = 5,
    LINEAREXTRUSION = 6,
    REVOLUTION = 7,
    BEZIER = 8,
    BSPLINE = 9,
    RECTANGULAR = 10,
    OFFSET = 11
  };

   switch (shapeType)
   {
   case PLANE:
     this->HandlePlane(in, out);
     break;
   case CYLINDER:
     this->HandleCylinder(in, out);
     break;
   case CONE:
     this->HandleCone(in, out);
     break;
   case SPHERE:
     this->HandleShpere(in, out);
     break;
   case TORUS:
     this->HandleTorus(in, out);
     break;
   case LINEAREXTRUSION:
     this->HandleLinearextrusion(in, out);
     break;
   case REVOLUTION:
     this->HandleRevolution(in, out);
     break;
   case BEZIER:
     this->HandleBezierSurface(in, out);
     break;
   case BSPLINE:
     this->HandleBsplineSurface(in, out);
     break;
   case RECTANGULAR:
     this->HandleRectangular(in, out);
     break;
   case OFFSET:
     this->HandleOffset(in, out);
     break;
   default:
     break;
   }
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleCurve(std::ifstream& in, int shapeType, vtkDataObject* out)
{
  enum {
    LINE = 1,
    CIRCLE = 2,
    ELLIPSE = 3,
    PARABOLA = 4,
    HYPERBOLA = 5,
    BEZIER = 6,
    BSPLINE = 7,
    TRIMMED = 8,
    OFFSET = 9
  };

  switch (shapeType)
  {
  case LINE:
    this->HandleLine(in, out);
    break;
  case CIRCLE:
    this->HandleCircle(in, out);
    break;
  case ELLIPSE:
    this->HandleEllipse(in, out);
    break;
  case PARABOLA:
    this->HandleParabola(in, out);
    break;
  case HYPERBOLA:
    this->HandleHyperbola(in, out);
    break;
  case BEZIER:
    this->HandleBezierCurve(in, out);
    break;
  case BSPLINE:
    this->HandleBsplineCurve(in, out);
    break;
  case TRIMMED:
    this->HandleTrimmed(in, out);
  case OFFSET:
    this->HandleOffsetCurve(in, out);
  default:
    break;
  }
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandlePlane(std::ifstream& in, vtkDataObject* vtkNotUsed(out))
{
  std::string out_buffer;
  for (int i = 0; i < 12; ++i)
  {
    in >> out_buffer;
  }
}

// ----------------------------------------------------------------------------
void vtkBREPReader::HandleCylinder(
  std::ifstream& in,
  vtkDataObject* vtkNotUsed(out))
{
  std::string out_buffer;
  for (int i = 0; i < 13; ++i)
  {
    in >> out_buffer;
  }
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleCone(
  std::ifstream& in,
  vtkDataObject* vtkNotUsed(out))
{
  std::string out_buffer;
  for (int i = 0; i < 14; ++i)
  {
    in >> out_buffer;
  }
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleShpere(
  std::ifstream& in,
  vtkDataObject* vtkNotUsed(out))
{
  std::string out_buffer;
  for (int i = 0; i < 13; ++i)
  {
    in >> out_buffer;
  }
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleTorus(
  std::ifstream& in,
  vtkDataObject* vtkNotUsed(out))
{
  std::string out_buffer;
  for (int i = 0; i < 14; ++i)
  {
    in >> out_buffer;
  }
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleLinearextrusion(
  std::ifstream& in,
  vtkDataObject* out)
{
  double Dv[3];
  in >> Dv[0] >> Dv[1] >> Dv[2];

  int curve_type;
  in >> curve_type;
  this->HandleCurve(in, curve_type, out);
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleRevolution(std::ifstream& in, vtkDataObject* out)
{
  double P[3];
  double D[3];
  in >> P[0] >> P[1] >> P[2] >> D[0] >> D[1] >> D[2];

  int curve_type;
  in >> curve_type;
  this->HandleCurve(in, curve_type, out);
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleBezierSurface(std::ifstream& in, vtkDataObject* out)
{
  int sample[2] = {20,20};
  int rational_flag[2] = {0,0};
  in >> rational_flag[0] >> rational_flag[1];

  int degree[2] = {1,1};
  in >> degree[0] >> degree[1];

  std::string out_buffer;
  vtkSmartPointer<vtkDoubleArray> ctrlPts = vtkSmartPointer<vtkDoubleArray>::New();
  ctrlPts->SetNumberOfComponents(4);
  ctrlPts->SetNumberOfTuples((degree[0] + 1) * (degree[1] + 1));
  bool rational = ((rational_flag[0] + rational_flag[1]) > 0 ? true : false);

  for (int i = 0; i < degree[0] + 1; ++i)
  {
    for (int j = 0; j < degree[1] + 1; ++j)
    {
      double pts[4] = {1.0,1.0,1.0,1.0};
      in >> pts[0] >> pts[1] >> pts[2];
      if (rational)
      {
        in >> pts[3];
        pts[0] *= pts[3];
        pts[1] *= pts[3];
        pts[2] *= pts[3];
      }
      ctrlPts->SetTuple(i + j * (degree[0] + 1), pts);
    }
  }

  vtkSmartPointer<vtkPatchInterpolation> patchInterp = vtkSmartPointer<vtkPatchInterpolation>::New();
  patchInterp->GenerateShape(
    vtkUnstructuredGrid::SafeDownCast(out),
    sample, 2, ctrlPts, degree);
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleBsplineSurface(std::ifstream& in, vtkDataObject* out)
{
  int u_r_flag, v_r_flag;
  int u_period, v_period;
  in >> u_r_flag >> v_r_flag;
  in >> u_period >> v_period;

  int u_degree, v_degree;
  int u_ctrlPt_num, v_ctrlPt_num;
  int u_knot_len, v_knot_len;
  in >> u_degree >> v_degree;
  in >> u_ctrlPt_num >> v_ctrlPt_num;
  in >> u_knot_len >> v_knot_len;

  vtkStructuredGrid* nurbs = vtkStructuredGrid::New();

  vtkSmartPointer<vtkDoubleArray> ctrlPts = vtkSmartPointer<vtkDoubleArray>::New();
  ctrlPts->SetNumberOfComponents(4);
  ctrlPts->SetNumberOfTuples(u_ctrlPt_num * v_ctrlPt_num);
  ctrlPts->SetName("ControlPolygon");
  bool rational = ((u_r_flag + v_r_flag) > 0 ? true : false);

  for (int i = 0; i < u_ctrlPt_num; ++i)
  {
    for (int j = 0; j < v_ctrlPt_num; ++j)
    {
      double pts[4] = {1.0,1.0,1.0,1.0};
      in >> pts[0] >> pts[1] >> pts[2];
      if (rational)
      {
        in >> pts[3];
        pts[0] *= pts[3];
        pts[1] *= pts[3];
        pts[2] *= pts[3];
      }
      ctrlPts->SetTuple(i + j * u_ctrlPt_num, pts);
    }
  }

  vtkSmartPointer<vtkControlPointArray<double> > ctrlPtsInSitu =
    vtkSmartPointer<vtkControlPointArray<double> >::New();
  ctrlPtsInSitu->InitializeArray(ctrlPts);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetData(ctrlPtsInSitu);

  nurbs->GetPointData()->AddArray(ctrlPts);
  nurbs->SetDimensions(u_ctrlPt_num, v_ctrlPt_num, 1);
  nurbs->SetPoints(points);

  std::vector<double> knotsLen(3, 0);
  std::vector<double> knotsArr;

  // load u knot vector
  for (int i = 0; i < u_knot_len; ++i)
  {
    double knots;
    int multi;
    in >> knots >> multi;
    knotsArr.insert(knotsArr.end(), multi, knots);
  }

  knotsLen[0] = knotsArr.size();

  for (int i = 0; i < v_knot_len; ++i)
  {
    double knots;
    int multi;
    in >> knots >> multi;
    knotsArr.insert(knotsArr.end(), multi, knots);
  }

  knotsLen[1] = knotsArr.size() - knotsLen[0];

  vtkSmartPointer<vtkDoubleArray> array1 = vtkSmartPointer<vtkDoubleArray>::New();
  array1->SetNumberOfComponents((int)knotsArr.size());
  array1->SetName("Knots");
  array1->InsertNextTuple(&knotsArr[0]);

  vtkSmartPointer<vtkIntArray> array2 = vtkSmartPointer<vtkIntArray>::New();
  array2->SetNumberOfComponents(3);
  array2->SetName("KnotsLen");
  array2->InsertNextTuple(&knotsLen[0]);

  // put control point, knots array and knot length to vtkStructuredGrid
  vtkSmartPointer<vtkFieldData> fieldData = vtkSmartPointer<vtkFieldData>::New();
  fieldData->AddArray(array1);
  fieldData->AddArray(array2);

  nurbs->SetFieldData(fieldData);

  vtkSmartPointer<vtkNURBSPatchAdaptor> nurbsAdaptor = vtkSmartPointer<vtkNURBSPatchAdaptor>::New();
  nurbsAdaptor->SetControlPointData(nurbs);
  nurbsAdaptor->GetPatchShape(vtkUnstructuredGrid::SafeDownCast(out));
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleRectangular(std::ifstream& in, vtkDataObject* out)
{
  double umin, umax, vmin, vmax;
  in >> umin >> umax >> vmin >> vmax;

  int surface_type;
  in >> surface_type;
  this->HandleSurface(in, surface_type, out);
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleOffset(std::ifstream& in, vtkDataObject* out)
{
  double d;
  in >> d;

  int surface_type;
  in >> surface_type;
  this->HandleSurface(in, surface_type, out);
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleLine(
  std::ifstream& in,
  vtkDataObject* vtkNotUsed(out))
{
  std::string out_buffer;
  for (int i = 0; i < 6; ++i)
  {
    in >> out_buffer;
  }
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleCircle(
  std::ifstream& in,
  vtkDataObject* vtkNotUsed(out))
{
  std::string out_buffer;
  for (int i = 0; i < 13; ++i)
  {
    in >> out_buffer;
  }
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleEllipse(
  std::ifstream& in,
  vtkDataObject* vtkNotUsed(out))
{
  std::string out_buffer;
  for (int i = 0; i < 14; ++i)
  {
    in >> out_buffer;
  }
}
// ----------------------------------------------------------------------
void vtkBREPReader::HandleParabola(
  std::ifstream& in,
  vtkDataObject* vtkNotUsed(out))
{
  std::string out_buffer;
  for (int i = 0; i < 13; ++i)
  {
    in >> out_buffer;
  }
}
// ----------------------------------------------------------------------
void vtkBREPReader::HandleHyperbola(
  std::ifstream& in,
  vtkDataObject* vtkNotUsed(out))
{
  std::string out_buffer;
  for (int i = 0; i < 14; ++i)
  {
    in >> out_buffer;
  }
}
// ----------------------------------------------------------------------
void vtkBREPReader::HandleBezierCurve(
  std::ifstream& in,
  vtkDataObject* vtkNotUsed(out))
{
  int rational;
  in >> rational;

  int degree;
  in >> degree;

  for (int i = 0; i < degree + 1; ++i)
  {
    double pts[4] = {1.0,1.0,1.0,1.0};
    in >> pts[0] >> pts[1] >> pts[2];
    if (rational)
    {
      in >> pts[3];
      pts[0] *= pts[3];
      pts[1] *= pts[3];
      pts[2] *= pts[3];
    }
  }
}
// ----------------------------------------------------------------------
void vtkBREPReader::HandleBsplineCurve(
  std::ifstream& in,
  vtkDataObject* vtkNotUsed(out))
{
  int rational;
  int period;
  in >> rational >> period;

  int degree, ctrlPt_num, knot_len;
  in >> degree >> ctrlPt_num >> knot_len;

  for (int i = 0; i < ctrlPt_num; ++i)
  {
    double pts[4];
    in >> pts[0] >> pts[1] >> pts[2];
    if (rational)
    {
      in >> pts[3];
      pts[0] *= pts[3];
      pts[1] *= pts[3];
      pts[2] *= pts[3];
    }
  }

  for (int i = 0; i < knot_len; ++i)
  {
    double knots;
    int multi;
    in >> knots >> multi;
  }
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleTrimmed(std::ifstream& in, vtkDataObject* out)
{
  double u_min, u_max;
  in >> u_min >> u_max;

  int curve_type;
  in >> curve_type;
  this->HandleCurve(in, curve_type, out);
}

// ----------------------------------------------------------------------
void vtkBREPReader::HandleOffsetCurve(std::ifstream& in, vtkDataObject* out)
{
  double d;
  in >> d;

  double dir[3];
  in >> dir[0] >> dir[1] >> dir[2];

  int curve_type;
  in >> curve_type;
  this->HandleCurve(in, curve_type, out);
}