#include "vtkPatchInterpolation.h"
#include "vtkNURBSPatchAdaptor.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkIndent.h"
#include "vtkFieldData.h"
#include "vtkCellArray.h"
#include "vtkTriangleStrip.h"
#include "vtkPolyData.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStructuredGrid.h"
#include "vtkLineSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkPointData.h"
#include "vtkControlPointArray.h"
#include "vtkNew.h"

#include "vtkBREPReader.h"

#include "vtkMath.h"
#include "vtkTypeTraits.h"

#include <ctime>

bool testRationalLine()
{


  return false;
}

bool testRationalCircle()
{
  return false;
}

bool testGenerateShape()
{
  vtkSmartPointer<vtkUnstructuredGrid> shape = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  shape->SetPoints(points);
  shape->Allocate();
  std::cout<<points->GetNumberOfPoints()<<"\n";

  int* sample = nullptr;
  int* degree = nullptr;
  vtkSmartPointer<vtkDoubleArray> ctrl_pts = vtkSmartPointer<vtkDoubleArray>::New();
  ctrl_pts->SetNumberOfComponents(4);
  ctrl_pts->InsertNextTuple4(1,1,0,1);
  ctrl_pts->InsertNextTuple4(1,1,1,1);
  ctrl_pts->InsertNextTuple4(2,0,2,2);
  ctrl_pts->InsertNextTuple4(-1,1,0,1);
  ctrl_pts->InsertNextTuple4(-1,1,1,1);
  ctrl_pts->InsertNextTuple4(-2,0,2,2);
  int dim = 2;
  degree = new int[dim];
  degree[0] = 2;
  degree[1] = 2;
  sample = new int[dim];
  sample[0] = 20;
  sample[1] = 20;

  vtkSmartPointer<vtkDoubleArray> ctrl_pts2 = vtkSmartPointer<vtkDoubleArray>::New();
  ctrl_pts2->SetNumberOfComponents(4);
  ctrl_pts2->InsertNextTuple4(10,1,0,1);
  ctrl_pts2->InsertNextTuple4(10,1,1,1);
  ctrl_pts2->InsertNextTuple4(20,0,2,2);
  ctrl_pts2->InsertNextTuple4(5,1,0,1);
  ctrl_pts2->InsertNextTuple4(5,1,1,1);
  ctrl_pts2->InsertNextTuple4(10,0,2,2);

  vtkSmartPointer<vtkPatchInterpolation> patchInterpolation = vtkSmartPointer<vtkPatchInterpolation>::New();
  patchInterpolation->GenerateShape(shape, sample, dim, ctrl_pts, degree);
  patchInterpolation->GenerateShape(shape, sample, dim, ctrl_pts2, degree);

  // Create an actor and mapper
  vtkSmartPointer<vtkDataSetMapper> mapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
  mapper->SetInput(polydata);
#else
  mapper->SetInputData(shape);
#endif

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  // Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(actor);
  renderWindow->Render();
  renderWindowInteractor->Start();

  return true;
}

bool testPointsStorage()
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkDoubleArray> ctrl_pts = vtkSmartPointer<vtkDoubleArray>::New();
  ctrl_pts->SetNumberOfComponents(4);
  for (int j = 0; j < 5; ++j)
  {
    for (int i = 0; i < 5; ++i)
    {
      ctrl_pts->InsertNextTuple4(i,j,i*j,1);
    }
  }
  return true;
}

bool testNURBSEval()
{
  return true;
}

int TestPatchInterpolation(int argc, char* argv[])
{
  return testGenerateShape();

  //return testNURBSShape();

//  vtkSmartPointer<vtkPoints> points =
//    vtkSmartPointer<vtkPoints>::New();
//
//  //points->SetNumberOfPoints(4);
//  points->GetData()->InsertNextTuple3(0,0,0);
//  points->GetData()->InsertNextTuple3(0,1,0);
//  points->GetData()->InsertNextTuple3(1,0,0);
//  points->GetData()->InsertNextTuple3(1.5,1,0);
//  //points->InsertNextPoint(0,0,0);
//  //points->InsertNextPoint(0,1,0);
//  //points->InsertNextPoint(1,0,0);
//  //points->InsertNextPoint(1.5,1,0);
//
//  vtkSmartPointer<vtkTriangleStrip> triangleStrip =
//    vtkSmartPointer<vtkTriangleStrip>::New();
//  triangleStrip->GetPointIds()->SetNumberOfIds(4);
//  triangleStrip->GetPointIds()->SetId(0,0);
//  triangleStrip->GetPointIds()->SetId(1,1);
//  triangleStrip->GetPointIds()->SetId(2,2);
//  triangleStrip->GetPointIds()->SetId(3,3);
//
//  vtkSmartPointer<vtkCellArray> cells =
//    vtkSmartPointer<vtkCellArray>::New();
//  cells->InsertNextCell(triangleStrip);
//
//  vtkSmartPointer<vtkPolyData> polydata =
//    vtkSmartPointer<vtkPolyData>::New();
//  polydata->SetPoints(points);
//  polydata->SetStrips(cells);
//
//  // Create an actor and mapper
//  vtkSmartPointer<vtkDataSetMapper> mapper =
//    vtkSmartPointer<vtkDataSetMapper>::New();
//#if VTK_MAJOR_VERSION <= 5
//  mapper->SetInput(polydata);
//#else
//  mapper->SetInputData(polydata);
//#endif
//
//  vtkSmartPointer<vtkActor> actor =
//    vtkSmartPointer<vtkActor>::New();
//  actor->SetMapper(mapper);
//
//  // Create a renderer, render window, and interactor
//  vtkSmartPointer<vtkRenderer> renderer =
//    vtkSmartPointer<vtkRenderer>::New();
//  vtkSmartPointer<vtkRenderWindow> renderWindow =
//    vtkSmartPointer<vtkRenderWindow>::New();
//  renderWindow->AddRenderer(renderer);
//  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//    vtkSmartPointer<vtkRenderWindowInteractor>::New();
//  renderWindowInteractor->SetRenderWindow(renderWindow);
//
//  renderer->AddActor(actor);
//  renderWindow->Render();
//  renderWindowInteractor->Start();

  //vtkIndent indent;



  //// test insert knot
  //vtkSmartPointer<vtkDoubleArray> pointsIn = vtkSmartPointer<vtkDoubleArray>::New();
  //pointsIn->SetNumberOfComponents(4);
  //pointsIn->InsertNextTuple4(0,0,0,1);
  //pointsIn->InsertNextTuple4(1,0.8,0,1);
  //pointsIn->InsertNextTuple4(2,-1,0,1);
  //pointsIn->InsertNextTuple4(3,-1,0,1);
  //pointsIn->InsertNextTuple4(4,1,0,1);
  //pointsIn->InsertNextTuple4(3,2,0,1);
  //pointsIn->InsertNextTuple4(2,2,0,1);
  //pointsIn->InsertNextTuple4(1.5,1.5,0,1);
  //int knotsLen[3] = {12,0,0};
  //double knotsArr[12] = {0,0,0.5,1,1,2,3,4,5,5,5,5};
  //int ctrlPtsNum[3] = {8,1,1};

  //vtkSmartPointer<vtkDoubleArray> pointsOut = vtkSmartPointer<vtkDoubleArray>::New();
  //int knotsLenNew[3] = {0,0,0};
  //double knotsArrNew[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
  //int ctrlPtsNumNew[3] = {0,0,0};

  //int insertDim = 0;
  //double tNew = 0;

  //vtkSmartPointer<vtkNURBSPatchAdaptor> adaptor = vtkSmartPointer<vtkNURBSPatchAdaptor>::New();
  ////adaptor->InsertKnot(
  ////  pointsIn, knotsLen, knotsArr, ctrlPtsNum,
  ////  pointsOut, knotsLenNew, knotsArrNew, ctrlPtsNumNew,
  ////  tNew, insertDim);

  //std::cout<<"Knots Len New: ";
  //for (int i = 0; i < 3; ++i)
  //{
  //  std::cout<<knotsLenNew[i]<<"\t";
  //}
  //std::cout<<"\n";

  //std::cout<<"Knots Array New: \n";
  //int cnt = 0;
  //for (int i = 0; i < 3; ++i)
  //{
  //  std::cout<<"Knot Vector Dim "<<i<<": ";
  //  for (int j = 0; j < knotsLenNew[i]; ++j)
  //  {
  //    std::cout<<knotsArrNew[cnt]<<"\t";
  //    ++cnt;
  //  }
  //  std::cout<<"\n";
  //}
  //std::cout<<"\n";

  //std::cout<<"Control Points Number New: ";
  //for (int i = 0; i < 3; ++i)
  //{
  //  std::cout<<ctrlPtsNumNew[i]<<"\t";
  //}
  //std::cout<<"\n";

  //std::cout<<"Control Points New:\n";
  //for (int i = 0; i < pointsOut->GetNumberOfTuples(); ++i)
  //{
  //  std::cout<<"Points "<<i<<": ";
  //  double* curPt = pointsOut->GetTuple(i);
  //  for (int j = 0; j < 4; ++j)
  //  {
  //    std::cout<<curPt[j]<<"\t";
  //  }
  //  std::cout<<"\n";
  //}
  //std::cout<<"\n";
  //return 0;

  //// test vtkFieldData
  //vtkSmartPointer<vtkFieldData> fieldData = vtkSmartPointer<vtkFieldData>::New();

  //double testKnot[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
  //vtkSmartPointer<vtkDoubleArray> array1 = vtkSmartPointer<vtkDoubleArray>::New();
  //array1->SetNumberOfComponents(8);
  //array1->SetName("Knots");
  //array1->InsertNextTuple(testKnot);

  //double lenKnots[3] = {8,0,0};
  //vtkSmartPointer<vtkIntArray> array2 = vtkSmartPointer<vtkIntArray>::New();
  //array2->SetNumberOfComponents(3);
  //array2->SetName("KnotsLen");
  //array2->InsertNextTuple(lenKnots);

  //fieldData->AddArray(array1);
  //fieldData->AddArray(array2);

  //fieldData->GetArray("Knots")->PrintSelf(std::cout, indent);
  //return 0;

  //// test binary search
  ////double testKnot[8] = {0.0,1.0,3.0,3.0,4.0,5.0,6.0,7.0};
  ////std::cout<<"Find span id: "<<adaptor->FindSpan(-1.0, 8, testKnot)<<"\n";
  //return 0;

  //// declare
  //double* r = nullptr;
  //int* degree = nullptr;
  //vtkSmartPointer<vtkDoubleArray> ctrl_pts = vtkSmartPointer<vtkDoubleArray>::New();
  //ctrl_pts->SetNumberOfComponents(4);
  //vtkSmartPointer<vtkDoubleArray> output_array = vtkSmartPointer<vtkDoubleArray>::New();
  //output_array->SetNumberOfComponents(4);
  //vtkSmartPointer<vtkPatchInterpolation> test_class = vtkSmartPointer<vtkPatchInterpolation>::New();

  //// case 1: dim == 1
  //r = new double[1];
  //degree = new int[1];
  //// case 1.1: linear
  //r[0] = 0.2;
  //degree[0] = 1;

  //ctrl_pts->InsertNextTuple4(0, 0, 0, 1);
  //ctrl_pts->InsertNextTuple4(1, 1, 1, 1);


  //  // case 1.2: quadratic
  //degree[0] = 2;
  //ctrl_pts->InsertNextTuple4(0.3, 0.5, 0.7, 1);
  //  // case 1.3: cubic
  //degree[0] = 3;
  //ctrl_pts->InsertNextTuple4(0.7, 0.5, 0.3, 1);
  //delete r;
  //delete degree;
  //// case 2: dim == 2
  //  // case 2.1: linear
  //  // case 2.2: quadratic
  //  // case 2.3: cubic
  //// case 3: dim == 3
  //  // case 3.1: linear
  //  // case 3.2: quadratic
  //  // case 3.3: cubic

  //std::cout << "test hahaha\n";

  //return 0;
}
