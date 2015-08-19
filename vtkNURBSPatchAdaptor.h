#ifndef vtkNURBSPatchAdaptor_h
#define vtkNURBSPatchAdaptor_h

#include "vtkBezierPatchAdaptor.h"
#include "vtkFiltersBezierModule.h"

#include <vector>

class vtkDataArray;
class vtkDoubleArray;
class vtkPoints;
class vtkUnstructuredGrid;

class VTKFILTERSBEZIER_EXPORT vtkNURBSPatchAdaptor : public vtkBezierPatchAdaptor
{
public:
  static vtkNURBSPatchAdaptor* New();

  vtkTypeMacro(vtkNURBSPatchAdaptor, vtkBezierPatchAdaptor);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual void SetControlPointData(vtkDataObject* controlPointData);
  vtkGetObjectMacro(ControlPointData, vtkDataObject);

  // Methods to iterate over patches:
  virtual void Begin();
  virtual bool IsAtEnd();
  virtual bool GoToPatch(vtkIdType patchId); // random access iteration
  virtual bool Next();

  // Methods to access the current patch:
  virtual int GetPatchDimension() const;
  virtual void GetPatchShape() const;
  // returns VTK_TRIANGLE/VTK_QUAD when PatchDimension==2,
  // returns VTK_HEXAHEDRON/VTK_TETRA when PatchDimension==3
  void GetPatchDegree(int* degreeOut) const;
  void GetPatchParameterRange(double* paramsOut) const;
  virtual void GetPatchPoints(vtkPoints* pointsOut) const;
  virtual void GetPatchShape(vtkUnstructuredGrid* outputSp);

  virtual void PointInversion(double Pt[3], double* Paras);
  virtual void EvaluatePoint(double Pt[3], double* Paras);
  virtual void EvaluateDeriv(double* Pt, double* Paras, int d);

protected:
  vtkNURBSPatchAdaptor();
  ~vtkNURBSPatchAdaptor();

  // dim in internal functions ranges from 0 - 2
  int FindSpan(double tNew, int lenKnot, double* knot);
  void FindSpanMult(
    double tNew, int lenKnot, double* knot,
    int& kSpan, int& sMult);
  void FindKnotStart(int* lenKnot, int dim, int& knotStart) const;
  void FetchBezierPatch(
    vtkUnstructuredGrid* shape,
    vtkDataArray* ctrlPts, int nurbsDim,
    int ctrlPtsNum[3], int ptAxisStart[3], int degree[3]);

  void InsertKnot(
    vtkDataArray* pointsIn, int* knotsLen, double* knotsArr, int* ctrlPtsNum,
    vtkDataArray* pointsOut, int* knotsLenNew, double* knotsArrNew, int* ctrlPtsNumNew,
    double tNew, int insertDim, int k = -1);
  void InsertKnotMulti(
    vtkDataArray* pointsIn, int* knotsLenOld, double* knotsArrOld, int* ctrlPtsNumOld,
    vtkDataArray* pointsOut, int* knotsLenNew, double* knotsArrNew, int* ctrlPtsNumNew,
    double tNew, int insertDim, int insertTimes, int kSpan = -1, int sMult = -1);
  void InsertKnotForEval(double* Paras, int& PtIdx);

  void EvaluateFirstOrderDeriv(
    double Dv[3],
    int degree, double U,
    double* Pts);
  void EvaluateSecondOrderDeriv(
    double Dv[3],
    int degree, double* U,
    double* Pts);
  void EvaluateSecondPartialDeriv(
    double Dv[3],
    int* degree, double* U,
    double* Pts);

  void EvaluateCurvDeriv(double* Pt, double* Paras, int d);
  void EvaluateCurveSecondDeriv(double* Dv, double* Paras);
  void EvaluateSurfaceDeriv(double* Pt, double* Paras, int d);
  void EvaluateSurfaceSecondDeriv(double* Dv, double* Paras);
  void EvaluateVolumeDeriv(double* Pt, double* Paras, int d);
  void EvaluateVolumeSecondDeriv(double* Dv, double* Paras);

  void PointInversionCurve(double Pt[3], double* Paras);
  void PointInversionSurface(double Pt[3], double* Paras);
  void PointInversionVolume(double Pt[3], double* Paras);

private:
  vtkNURBSPatchAdaptor(const vtkNURBSPatchAdaptor&);
  void operator = (const vtkNURBSPatchAdaptor&);

  vtkDoubleArray* PointsCache;
  std::vector<int> CtrlPtsNumCache;
  std::vector<int> KnotsLenCache;
  std::vector<double> KnotsArrCache;
};

#endif