/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkBezierPatchAdaptor.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//.NAME vtkBezierPatchAdaptor ¨C Abstract class for spline patch adapter
//.SECTION Description
// vtkBezierPatchAdaptor defined interface for all spline patch adapter,
// including get degree, range of parameters and methods to iterate over
// patches.

#ifndef __vtkBezierPatchAdaptor_h
#define __vtkBezierPatchAdaptor_h

#include "vtkObject.h"
#include "vtkFiltersBezierModule.h"

class vtkDataObject;
class vtkPoints;

/**\adapter to iterate bezier patches
  *
  */

class VTKFILTERSBEZIER_EXPORT vtkBezierPatchAdaptor : public vtkObject
{
public:
  vtkTypeMacro(vtkBezierPatchAdaptor, vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void SetControlPointData(vtkDataObject* controlPointData);
  vtkGetObjectMacro(ControlPointData, vtkDataObject);

  // Description:
  // Methods to iterate over patches.
  virtual void Begin() = 0;
  virtual bool IsAtEnd() = 0;
  virtual bool GoToPatch(vtkIdType patchId) = 0; // random access iteration
  virtual bool Next() = 0;

  // Description:
  // Methods to access the current patch.
  virtual int GetPatchDimension() const = 0;
  virtual void GetPatchShape() const = 0;
    // returns VTK_TRIANGLE/VTK_QUAD when PatchDimension==2,
    // returns VTK_HEXAHEDRON/VTK_TETRA when PatchDimension==3
  virtual void GetPatchDegree(int* degreeOut) const = 0;
  virtual void GetPatchParameterRange(double* paramsOut) const = 0;
  virtual void GetPatchPoints(vtkPoints* pointsOut) const = 0;

  // Description:
  // Methods to evaluate points, derivatives and projection.
  virtual void PointInversion(double Pt[3], double* Paras);
  virtual void EvaluatePoint(double Pt[3], double* Paras);
  virtual void EvaluateDeriv(double* Pt, double* Paras, int d);

  // Some helper methods that could make Python wrappings usable:
  int GetPatchDegree(int dim) const;
  void GetPatchParameterRange(int dim, double paramRange[2]);

protected:
  vtkBezierPatchAdaptor();
  ~vtkBezierPatchAdaptor();
  vtkDataObject* ControlPointData;

private:
  vtkBezierPatchAdaptor(const vtkBezierPatchAdaptor&); // Not implemented.
  void operator = (const vtkBezierPatchAdaptor&); // Not implemented.
};
#endif