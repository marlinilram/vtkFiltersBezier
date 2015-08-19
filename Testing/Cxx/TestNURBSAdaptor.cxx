#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkControlPointArray.h"
#include "vtkDataSetMapper.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkIntArray.h"
#include "vtkNew.h"
#include "vtkNURBSPatchAdaptor.h"
#include "vtkPointData.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"

#include "vtkTestUtilities.h"
#include "vtkRegressionTestImage.h"

#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"

int TestNURBSAdaptor(int argc, char* argv[])
{
  vtkNew<vtkStructuredGrid> nurbs;

  vtkNew<vtkFieldData> fieldData;

  double knotsLen[3] = {12,0,0};
  double knotsArr[12] = {
    0,0,0,0.25,0.25,0.5,0.5,0.75,0.75,1,1,1};

  vtkNew<vtkDoubleArray> array1;
  array1->SetNumberOfComponents(sizeof(knotsArr)/sizeof(*knotsArr));
  array1->SetName("Knots");
  array1->InsertNextTuple(knotsArr);

  vtkNew<vtkIntArray> array2;
  array2->SetNumberOfComponents(3);
  array2->SetName("KnotsLen");
  array2->InsertNextTuple(knotsLen);

  fieldData->AddArray(array1.GetPointer());
  fieldData->AddArray(array2.GetPointer());
  nurbs->SetFieldData(fieldData.GetPointer());

  vtkNew<vtkDoubleArray> ctrl_pts;
  ctrl_pts->SetNumberOfComponents(4);
  ctrl_pts->SetName("ControlPolygon");
  //for (int j = 0; j < 5; ++j)
  //{
  //  for (int i = 0; i < 5; ++i)
  //  {
  //      ctrl_pts->InsertNextTuple4(i,j,i*j,1);
  //  }
  //}
  //ctrl_pts->InsertNextTuple4(0,1,5,1);
  //ctrl_pts->InsertNextTuple4(0.5,1,4,1);
  //ctrl_pts->InsertNextTuple4(3,1,4,1);
  //ctrl_pts->InsertNextTuple4(4,1,3,1);
  //ctrl_pts->InsertNextTuple4(6,1,3,1);

  //ctrl_pts->InsertNextTuple4(0,2,5,1);
  //ctrl_pts->InsertNextTuple4(0.5,2,4,1);
  //ctrl_pts->InsertNextTuple4(3,2,4,1);
  //ctrl_pts->InsertNextTuple4(4,2,3,1);
  //ctrl_pts->InsertNextTuple4(6,2,3,1);

  //ctrl_pts->InsertNextTuple4(0,2.3,3,1);
  //ctrl_pts->InsertNextTuple4(0.5,2.3,2,1);
  //ctrl_pts->InsertNextTuple4(3,2.3,2,1);
  //ctrl_pts->InsertNextTuple4(4,2.3,0,1);
  //ctrl_pts->InsertNextTuple4(6,2.3,0,1);

  //ctrl_pts->InsertNextTuple4(0,3,3,1);
  //ctrl_pts->InsertNextTuple4(0.5,3,2,1);
  //ctrl_pts->InsertNextTuple4(3,3.1,2,1);
  //ctrl_pts->InsertNextTuple4(4,3.2,0,1);
  //ctrl_pts->InsertNextTuple4(6,3.2,0,1);

  //ctrl_pts->InsertNextTuple4(0,3.5,3,1);
  //ctrl_pts->InsertNextTuple4(0.5,3.5,2,1);
  //ctrl_pts->InsertNextTuple4(3,3.8,2,1);
  //ctrl_pts->InsertNextTuple4(4,3.8,0,1);
  //ctrl_pts->InsertNextTuple4(6,3.8,0,1);

  ctrl_pts->InsertNextTuple4(1,0,0,1);
  ctrl_pts->InsertNextTuple4(1/sqrt(2),1/sqrt(2),0,1/sqrt(2));
  ctrl_pts->InsertNextTuple4(0,1,0,1);
  ctrl_pts->InsertNextTuple4(-1/sqrt(2),1/sqrt(2),0,1/sqrt(2));
  ctrl_pts->InsertNextTuple4(-1,0,0,1);
  ctrl_pts->InsertNextTuple4(-1/sqrt(2),-1/sqrt(2),0,1/sqrt(2));
  ctrl_pts->InsertNextTuple4(0,-1,0,1);
  ctrl_pts->InsertNextTuple4(1/sqrt(2),-1/sqrt(2),0,1/sqrt(2));
  ctrl_pts->InsertNextTuple4(1,0,0,1);

  nurbs->GetPointData()->AddArray(ctrl_pts.GetPointer());
  nurbs->SetDimensions(9,1,1);

  vtkNew<vtkControlPointArray<double>> ctrlPts;
  ctrlPts->InitializeArray(ctrl_pts.GetPointer());

  vtkNew<vtkPoints> points;
  points->SetData(ctrlPts.GetPointer());
  nurbs->SetPoints(points.GetPointer());

  std::cout<<"number of points: "<<ctrlPts->GetNumberOfTuples()<<"\n";
  std::cout<<"number of components: "<<ctrlPts->GetNumberOfComponents()<<"\n";

  //for (int i = 0; i < ctrlPts->GetNumberOfTuples(); ++i)
  //{
  //  ctrlPts->GetTuple(i, pt);
  //  std::cout<<"Pt "<<i<<": "<<pt[0]<<"\t"<<pt[1]<<"\t"<<pt[2]<<"\n";
  //}

  vtkNew<vtkNURBSPatchAdaptor> nurbsAdaptor;
  nurbsAdaptor->SetControlPointData(nurbs.GetPointer());

  vtkNew<vtkUnstructuredGrid> bezierShape;
  vtkNew<vtkPoints> bezierPoints;
  bezierShape->SetPoints(bezierPoints.GetPointer());
  bezierShape->Allocate();

  nurbsAdaptor->GetPatchShape(bezierShape.GetPointer());

  // Create an actor and mapper
  vtkNew<vtkDataSetMapper> mapper;
  mapper->SetInputData(bezierShape.GetPointer());

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper.GetPointer());

  // Create a renderer, render window, and interactor
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer.GetPointer());
  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow.GetPointer());

  renderer->AddActor(actor.GetPointer());
  renderWindow->Render();
  renderWindowInteractor->Start();

  return true;
}
