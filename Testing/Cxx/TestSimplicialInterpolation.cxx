/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestSimplicialInterpolation.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkDataSetMapper.h"
#include "vtkDoubleArray.h"
#include "vtkNew.h"
#include "vtkNURBSPatchAdaptor.h"
#include "vtkPatchInterpolation.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"

#include "vtkTestUtilities.h"
#include "vtkRegressionTestImage.h"

#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"

int TestSimplicialInterpolation(int argc, char* argv[])
{
  vtkSmartPointer<vtkUnstructuredGrid> shape = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  shape->SetPoints(points);
  shape->Allocate();

  int sample[2] = {20,20};

  vtkSmartPointer<vtkDoubleArray> ctrl_pts2 = vtkSmartPointer<vtkDoubleArray>::New();
  ctrl_pts2->SetNumberOfComponents(4);
  ctrl_pts2->InsertNextTuple4(-39,-10,0,1);
  ctrl_pts2->InsertNextTuple4(-13,-11,10,1);
  ctrl_pts2->InsertNextTuple4(15,-22,-12,1);
  ctrl_pts2->InsertNextTuple4(42,-13,0,1);
  ctrl_pts2->InsertNextTuple4(-25,2.2,0,1);
  ctrl_pts2->InsertNextTuple4(0,6,8,1);
  ctrl_pts2->InsertNextTuple4(26,2.2,0,1);
  ctrl_pts2->InsertNextTuple4(-20,13,0,1);
  ctrl_pts2->InsertNextTuple4(0,12,0,1);
  ctrl_pts2->InsertNextTuple4(-9,26,-7,1);

  vtkNew<vtkPatchInterpolation> patchInterpolation;
  patchInterpolation->GenerateSimplicialShape(shape, sample, 3, ctrl_pts2, 3);

  // Create an actor and mapper
  vtkNew<vtkDataSetMapper> mapper;
  mapper->SetInputData(shape);

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper.GetPointer());

  // Create a renderer, render window, and interactor
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer.GetPointer());
  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow.GetPointer());

  renderer->AddActor(actor.GetPointer());
  renderer->GetActiveCamera()->SetPosition(1.5,4.745,175.553);
  renderer->GetActiveCamera()->SetFocalPoint(1.5,4.745,-2.156);
  renderer->GetActiveCamera()->SetViewUp(0,1,0);
  renderer->GetActiveCamera()->SetClippingRange(166.292, 192.509);
  renderer->GetActiveCamera()->SetViewAngle(30);
  renderWindow->Render();

  int testStatus = vtkRegressionTestImage(renderWindow.GetPointer());
  if (testStatus == vtkRegressionTester::DO_INTERACTOR)
  {
    renderWindowInteractor->Start();
  }

  return (testStatus ? EXIT_SUCCESS : EXIT_FAILURE);
}
