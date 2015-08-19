/*=========================================================================

  Program:   Visualization Toolkit
  Module:    TestBREPReader.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkActor.h"
#include "vtkBREPReader.h"
#include "vtkCamera.h"
#include "vtkDataSetMapper.h"
#include "vtkNew.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkUnstructuredGrid.h"

#include "vtkTestUtilities.h"
#include "vtkRegressionTestImage.h"

#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"

int TestBREPReader(int argc, char* argv[])
{
  char* modelFile = vtkTestUtilities::ExpandDataFileName(argc, argv,
    "Data/Moto.brep");

  // Read file
  vtkNew<vtkBREPReader> brepReader;
  brepReader->SetFileName(modelFile);
  brepReader->Update();
  vtkUnstructuredGrid* outShape = brepReader->GetOutput();

  // Create an actor and mapper
  vtkNew<vtkDataSetMapper> mapper;
  mapper->SetInputData(outShape);

  vtkNew<vtkActor> actor;
  actor->SetMapper(mapper.GetPointer());

  // Create a renderer, render window, and interactor
  vtkNew<vtkRenderer> renderer;
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer.GetPointer());
  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow.GetPointer());

  renderer->AddActor(actor.GetPointer());
  renderer->GetActiveCamera()->SetPosition(-781.363, -5316.17, 263.986);
  renderer->GetActiveCamera()->SetFocalPoint(-781.363, 0, 263.986);
  renderer->GetActiveCamera()->SetViewUp(0,0,1);
  renderer->GetActiveCamera()->SetClippingRange(3000, 6000);
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