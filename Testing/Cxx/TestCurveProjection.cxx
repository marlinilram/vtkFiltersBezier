/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestCurveProjection.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkControlPointArray.h"
#include "vtkDataSetMapper.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkIntArray.h"
#include "vtkLineSource.h"
#include "vtkNew.h"
#include "vtkNURBSPatchAdaptor.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkStructuredGrid.h"
#include "vtkSphereSource.h"
#include "vtkUnstructuredGrid.h"

#include "vtkTestUtilities.h"
#include "vtkRegressionTestImage.h"

#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"

int TestCurveProjection(int argc, char* argv[])
{
  vtkNew<vtkStructuredGrid> nurbs;

  vtkNew<vtkFieldData> fieldData;

  double knotsLen[3] = {10,0,0};
  double knotsArr[10] = {
    0,0,0,0.2, 0.4, 0.6, 0.8,1,1,1};

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

  ctrl_pts->InsertNextTuple4(-2,0.25,0,1);
  ctrl_pts->InsertNextTuple4(-1.75,1.5,0,1);
  ctrl_pts->InsertNextTuple4(1.5*2,1.75*2,0,2);
  ctrl_pts->InsertNextTuple4(0,0,0,1);
  ctrl_pts->InsertNextTuple4(-1.5*2,-1.75*2,0,2);
  ctrl_pts->InsertNextTuple4(1.75,-1.5,0,1);
  ctrl_pts->InsertNextTuple4(2,-0.25,0,1);

  nurbs->GetPointData()->AddArray(ctrl_pts.GetPointer());
  nurbs->SetDimensions(7,1,1);

  vtkNew<vtkControlPointArray<double>> ctrlPts;
  ctrlPts->InitializeArray(ctrl_pts.GetPointer());

  vtkNew<vtkPoints> points;
  points->SetData(ctrlPts.GetPointer());
  nurbs->SetPoints(points.GetPointer());

  //std::cout<<"number of points: "<<ctrlPts->GetNumberOfTuples()<<"\n";
  //std::cout<<"number of components: "<<ctrlPts->GetNumberOfComponents()<<"\n";

  vtkNew<vtkNURBSPatchAdaptor> nurbsAdaptor;
  nurbsAdaptor->SetControlPointData(nurbs.GetPointer());

  vtkNew<vtkUnstructuredGrid> bezierShape;
  vtkNew<vtkPoints> bezierPoints;
  bezierShape->SetPoints(bezierPoints.GetPointer());
  bezierShape->Allocate();

  nurbsAdaptor->GetPatchShape(bezierShape.GetPointer());

  // test projection
  std::vector<vtkSmartPointer<vtkActor>> actorsProjLine;
  double target[27] = {-2.1,-0.25,0,-1.5,0.5,0,-2,1.75,0,-0.25,2.0,0,0,0.25,0,1.6, 2.0,0,1.0,0,0,1.5,-1.25,0,2,0.25,0};
  for (int i = 0; i < 9; ++i)
  {
    double projParas[3] = {0,0,0};
    double projPt[3] = {0,0,0};
    nurbsAdaptor->PointInversion(&target[3*i], projParas);
    nurbsAdaptor->EvaluatePoint(projPt, projParas);
    vtkNew<vtkLineSource> proj_line;
    proj_line->SetPoint1(&target[3*i]);
    proj_line->SetPoint2(projPt);
    proj_line->Update();
    vtkNew<vtkPolyDataMapper> mapperProjLine;
    mapperProjLine->SetInputConnection(proj_line->GetOutputPort());
    vtkNew<vtkActor> actorProjLine;
    actorProjLine->SetMapper(mapperProjLine.GetPointer());
    actorsProjLine.push_back(actorProjLine.GetPointer());
  }

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
  for (int i = 0; i < actorsProjLine.size(); ++i)
  {
  renderer->AddActor(actorsProjLine[i]);
  }
  renderer->GetActiveCamera()->SetPosition(0,0,10);
  renderer->GetActiveCamera()->SetFocalPoint(0,0,0);
  renderer->GetActiveCamera()->SetViewUp(0,1,0);
  renderer->GetActiveCamera()->SetClippingRange(0.01, 1000);
  renderer->GetActiveCamera()->SetViewAngle(30);
  renderWindow->SetSize(300, 300);
  renderWindow->Render();

  int testStatus = vtkRegressionTestImage(renderWindow.GetPointer());
  if (testStatus == vtkRegressionTester::DO_INTERACTOR)
  {
    renderWindowInteractor->Start();
  }

  return (testStatus ? EXIT_SUCCESS : EXIT_FAILURE);
}
