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

  // Methods to iterate over patches:
  virtual void Begin() = 0;
  virtual bool IsAtEnd() = 0;
  virtual bool GoToPatch(vtkIdType patchId) = 0; // random access iteration
  virtual bool Next() = 0;

  // Methods to access the current patch:
  virtual int GetPatchDimension() const = 0;
  virtual void GetPatchShape() const = 0;
    // returns VTK_TRIANGLE/VTK_QUAD when PatchDimension==2,
    // returns VTK_HEXAHEDRON/VTK_TETRA when PatchDimension==3
  virtual void GetPatchDegree(int* degreeOut) const = 0;
  virtual void GetPatchParameterRange(double* paramsOut) const = 0;
  virtual void GetPatchPoints(vtkPoints* pointsOut) const = 0;

  virtual void PointInversion(double Pt[3], double* Paras);
  virtual void EvaluatePoint(double Pt[3], double* Paras);
  virtual void EvaluateDeriv(double Pt[3], double* Paras, int d);

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