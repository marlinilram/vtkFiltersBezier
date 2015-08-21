/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestSurfaceProjection.cxx

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
#include "vtkMath.h"
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

int TestSurfaceProjection(int argc, char* argv[])
{
  vtkNew<vtkStructuredGrid> nurbs;
  // Important: a nurbs data is stored in a vtkStructuredGrid
  // knot vector are stored in a 1-D array, so there should be
  // an additional array to store the length of knot vector for
  // each dimension.
  //
  // Knots array are stored with name "Knots".
  // Knots len are stored with name "KnoraLwn"
  // Both of them are stored in the vtkFieldData of vtkStructuredGrid.
  //
  // Control points are stored in a 4-component array with name
  // "ControlPolygon". It is held by vtkStructuredGrid::vtkPointData
  //
  // The vtkStructuredGrid::vtkPoints also hold a virtual copy of
  // the control points. Because vtkPoints only allow a vtkDataArray
  // with 3 components. We need to initialize a vtkControlPointArray
  // to hold the "ControlPolygon" array and then put the
  // vtkControlPointArray into the vtkStructuredGrid::vtkPoints.
  // The vtkControlPointArray is an in-situ vtkDataArray which interprets
  // the 4-component "ControlPolygon" array as a 3-component array.

  vtkNew<vtkFieldData> fieldData;

  double knotsLen[3] = {9,9,0};
  double knotsArr[18] = {
    0,0,0,0,0.5,1,1,1,1,0,0,0,0,0.5,1,1,1,1};

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

  ctrl_pts->InsertNextTuple4(-50.118762969970703, -50.593822479248047, 19.958602905273441, 1);
  ctrl_pts->InsertNextTuple4(-25.059381484985352, -50.593822479248047, 19.958602905273441,  1);
  ctrl_pts->InsertNextTuple4(0, -50.593822479248047, 19.958602905273441,  1);
  ctrl_pts->InsertNextTuple4(25.059381484985352, -50.593822479248047, 19.958602905273441,  1);
  ctrl_pts->InsertNextTuple4(50.118762969970703, -50.593822479248047, 19.958602905273441,  1);
  ctrl_pts->InsertNextTuple4(-50.118762969970703, -25.415678024291989, 19.958602905273441,  1);
  ctrl_pts->InsertNextTuple4(-25.059381484985352*5, -25.415676116943359*5, 0,  5);
  ctrl_pts->InsertNextTuple4(0, -25.415676116943359*5, 0,  5);
  ctrl_pts->InsertNextTuple4(25.059381484985352*5, -25.415676116943359*5, 0,  5);
  ctrl_pts->InsertNextTuple4(50.118762969970703, -25.415678024291989, 19.958602905273441,  1);
  ctrl_pts->InsertNextTuple4(-50.118762969970703, -0.23753062705691261, 19.958602905273441,  1);
  ctrl_pts->InsertNextTuple4(-25.059381484985352*5, -0.2375297546386719*5, 0,  5);
  ctrl_pts->InsertNextTuple4(0, -0.23753241037047701*10, 60.756065368652337*10,  10);
  ctrl_pts->InsertNextTuple4(25.059381484985352*5, -0.2375297546386719*5, 0,  5);
  ctrl_pts->InsertNextTuple4(50.118762969970703, -0.23753062705691261, 19.958602905273441,  1);
  ctrl_pts->InsertNextTuple4(-50.118762969970703, 24.940614700317379, 19.958602905273441,  1);
  ctrl_pts->InsertNextTuple4(-25.059381484985352*5, 24.940616607666019*5, 0,  5);
  ctrl_pts->InsertNextTuple4(0, 24.940616607666019*5, 0,  5);
  ctrl_pts->InsertNextTuple4(25.059381484985352*5, 24.940616607666019*5, 0,  5);
  ctrl_pts->InsertNextTuple4(50.118762969970703, 24.940614700317379, 19.958602905273441,  1);
  ctrl_pts->InsertNextTuple4(-50.118762969970703, 50.118762969970703, 19.958602905273441,  1);
  ctrl_pts->InsertNextTuple4(-25.059381484985352, 50.118762969970703, 19.958602905273441, 1);
  ctrl_pts->InsertNextTuple4(0, 50.118762969970703, 19.958602905273441,  1);
  ctrl_pts->InsertNextTuple4(25.059381484985352, 50.118762969970703, 19.958602905273441,  1);
  ctrl_pts->InsertNextTuple4(50.118762969970703, 50.118762969970703, 19.958602905273441,  1);

  nurbs->GetPointData()->AddArray(ctrl_pts.GetPointer());
  nurbs->SetDimensions(5,5,1);

  vtkNew<vtkControlPointArray<double>> ctrlPts;
  ctrlPts->InitializeArray(ctrl_pts.GetPointer());

  vtkNew<vtkPoints> points;
  points->SetData(ctrlPts.GetPointer());
  nurbs->SetPoints(points.GetPointer());

  vtkNew<vtkNURBSPatchAdaptor> nurbsAdaptor;
  nurbsAdaptor->SetControlPointData(nurbs.GetPointer());

  vtkNew<vtkUnstructuredGrid> bezierShape;
  vtkNew<vtkPoints> bezierPoints;
  bezierShape->SetPoints(bezierPoints.GetPointer());
  bezierShape->Allocate();

  nurbsAdaptor->GetPatchShape(bezierShape.GetPointer());

  // test projection
  std::vector<vtkSmartPointer<vtkActor>> actorsProjLine;
  double target[27] = {0,0,80,-15,30,30,0,35,55,-35,20,30,-10,-25,20,16, 20,30,40,-20,20,5,-25,20,20,25,10};
  for (int i = 0; i < 9; ++i)
    {
    double projParas[3] = {0,0,0};
    double projPt[3] = {0,0,0};
    double deriPt[6] = {0,0,0,0,0,0};
    double normPt[3] = {0,0,0};
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
    vtkNew<vtkSphereSource> proj_sphere;
    proj_sphere->SetCenter(&target[3*i]);
    proj_sphere->SetRadius(1.0);
    proj_sphere->Update();
    vtkNew<vtkPolyDataMapper> mapperProjSphere;
    mapperProjSphere->SetInputConnection(proj_sphere->GetOutputPort());
    vtkNew<vtkActor> actorProjSphere;
    actorProjSphere->SetMapper(mapperProjSphere.GetPointer());
    actorsProjLine.push_back(actorProjSphere.GetPointer());
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
  renderer->GetActiveCamera()->SetPosition(135,-135,100);
  renderer->GetActiveCamera()->SetFocalPoint(0,0,30);
  renderer->GetActiveCamera()->SetViewUp(0,0,1);
  renderer->GetActiveCamera()->SetClippingRange(0.01, 1000);
  renderer->GetActiveCamera()->SetViewAngle(30);
  renderWindow->SetSize(400, 400);
  renderWindow->Render();

  int testStatus = vtkRegressionTestImage(renderWindow.GetPointer());
  if (testStatus == vtkRegressionTester::DO_INTERACTOR)
    {
    renderWindowInteractor->Start();
    }

  return (testStatus ? EXIT_SUCCESS : EXIT_FAILURE);
}
