/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkBREPReader.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//.NAME vtkBREPReader ¨C reader for .brep file
//.SECTION Description
// vtkBREPReader can read and processed the bezier and nurbs patch in .brep
// files and return vtkUnstructuredGrid for visualization.
// Currently only a small part of shapes defined in brep file can be processed,
// the rest will be added in the future.
//
// The class should be initialized with the fileName and it will return a
// vtkUnstructuredGrid which contains the triangle mesh of the input file for
// visualization.

#ifndef vtkBREPReader_h
#define vtkBREPReader_h

#include "vtkAlgorithm.h"
#include "vtkFiltersBezierModule.h"

class vtkUnstructuredGrid;

class VTKFILTERSBEZIER_EXPORT vtkBREPReader : public vtkAlgorithm
{
public:
  static vtkBREPReader *New();
  vtkTypeMacro(vtkBREPReader, vtkAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify file name of data file (brep).
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Get the output data object for a port on this algorithm.
  vtkUnstructuredGrid* GetOutput();
  vtkUnstructuredGrid* GetOutput(int);
  virtual void SetOutput(vtkDataObject* d);

  // Description:
  // see vtkAlgorithm for details
  virtual int ProcessRequest(vtkInformation*,
                             vtkInformationVector**,
                             vtkInformationVector*);

  // this method is not recommended for use, but lots of old style filters
  // use it
  vtkDataObject* GetInput();
  vtkDataObject* GetInput(int port);

  // Description:
  // Assign a data object as input. Note that this method does not
  // establish a pipeline connection. Use SetInputConnection() to
  // setup a pipeline connection.
  void SetInputData(vtkDataObject*);
  void SetInputData(int, vtkDataObject*);

  // Description:
  // Assign a data object as input. Note that this method does not
  // establish a pipeline connection. Use AddInputConnection() to
  // setup a pipeline connection.
  void AddInputData(vtkDataObject*);
  void AddInputData(int, vtkDataObject*);

protected:
  vtkBREPReader();
  ~vtkBREPReader();

  // convenience method
  virtual int RequestInformation(vtkInformation* request,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector);

  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  virtual int RequestData(vtkInformation* request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);

  // Description:
  // This is called by the superclass.
  // This is the method you should override.
  virtual int RequestUpdateExtent(vtkInformation*,
                                  vtkInformationVector**,
                                  vtkInformationVector*);

  // see algorithm for more info
  virtual int FillOutputPortInformation(int port, vtkInformation* info);
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Description:
  // Process different type of shape in BREP files.
  // Only bezier and nurbs surface are implemented now.
  void HandleSurface(std::ifstream& in, int shapeType, vtkDataObject* out);
  void HandlePlane(std::ifstream& in, vtkDataObject* out);
  void HandleCylinder(std::ifstream& in, vtkDataObject* out);
  void HandleCone(std::ifstream& in, vtkDataObject* out);
  void HandleShpere(std::ifstream& in, vtkDataObject* out);
  void HandleTorus(std::ifstream& in, vtkDataObject* out);
  void HandleLinearextrusion(std::ifstream& in, vtkDataObject* out);
  void HandleRevolution(std::ifstream& in, vtkDataObject* out);
  void HandleBezierSurface(std::ifstream& in, vtkDataObject* out);
  void HandleBsplineSurface(std::ifstream& in, vtkDataObject* out);
  void HandleRectangular(std::ifstream& in, vtkDataObject* out);
  void HandleOffset(std::ifstream& in, vtkDataObject* out);
  void HandleCurve(std::ifstream& in, int shapeType, vtkDataObject* out);
  void HandleLine(std::ifstream& in, vtkDataObject* out);
  void HandleCircle(std::ifstream& in, vtkDataObject* out);
  void HandleEllipse(std::ifstream& in, vtkDataObject* out);
  void HandleParabola(std::ifstream& in, vtkDataObject* out);
  void HandleHyperbola(std::ifstream& in, vtkDataObject* out);
  void HandleBezierCurve(std::ifstream& in, vtkDataObject* out);
  void HandleBsplineCurve(std::ifstream& in, vtkDataObject* out);
  void HandleTrimmed(std::ifstream& in, vtkDataObject* out);
  void HandleOffsetCurve(std::ifstream& in, vtkDataObject* out);

  char* FileName;

private:
  vtkBREPReader(const vtkBREPReader&);
  void operator=(const vtkBREPReader&);
};

#endif