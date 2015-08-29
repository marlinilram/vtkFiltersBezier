/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkNURBSPatchAdaptor.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//.NAME vtkNURBSPatchAdaptor ¨C adapter class for NURBS patch
//.SECTION Description
// vtkNURBSPatchAdaptor provides several methods to support computation and
// visualization of 1-d, 2-d and 3-d NURBS (Non-uniform Rational B-spline).
// It converts input NURBS patch to Bezier patch for visualization with
// methods to evaluate point, first order and second order derivatives of
// NURBS patches. The core algorithm for these methods is knot insertion.
// To calculate point inversion and projection, this class implements a Newton
// method to optimize the estimated parameters.
//
// This class provide all necessary methods for NURBS patch .
// It should be initialized by calling method SetControlPointData and it will
// add the reference count of the input vtkDataObject.
// The input vtkDataObject should be vtkStructuredGrid. How NURBS data is
// stored in vtkStructuredGrid can be seen in TestCurveProjection.cxx and
// TestSurfaceProjection.cxx
//
// All of the algorithms are highly based on:
// Piegl, Les; Tiller, Wayne (1997). The NURBS Book (2. ed.).
// Berlin: Springer. ISBN 3-540-61545-8.

#ifndef vtkNURBSPatchAdaptor_h
#define vtkNURBSPatchAdaptor_h

#include "vtkBezierPatchAdaptor.h"
#include "vtkFiltersBezierModule.h"

#include <vector>

class vtkDataArray;
class vtkDoubleArray;
class vtkPoints;
class vtkUnstructuredGrid;

class VTKFILTERSBEZIER_EXPORT vtkNURBSPatchAdaptor
  : public vtkBezierPatchAdaptor
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

  // Description:
  // Given a point, return the estimated parameter coordinate.
  virtual void PointInversion(double Pt[3], double* Paras);
  // Description:
  // Given a parameter coordinate, return point position.
  virtual void EvaluatePoint(double Pt[3], double* Paras);
  // Description:
  // Given a parameter coordinate and the order of derivative d,
  // return derivatives of the input order d.
  // Notice the Pt should allocate memory first and it depends
  // on both d and the dimension of the NURBS patch.
  // For example, for 2-D patch and first order derivatives
  // it will be [1, 0] and [0, 1] which means two partial
  // derivatives for each dimension; for 2-D patch and second
  // order derivatives, it will be [2, 0], [1, 1] and [0,2].
  virtual void EvaluateDeriv(double* Pt, double* Paras, int d);

protected:
  vtkNURBSPatchAdaptor();
  ~vtkNURBSPatchAdaptor();

  // dim in internal functions ranges from 0 - 2

  // Description:
  // Find span of which the knot value we want to insert lies in
  int FindSpan(double tNew, int lenKnot, double* knot);
  // Description:
  // Find span of which the new knot value lies in and the multiplicity
  // of this value.
  void FindSpanMult(
    double tNew, int lenKnot, double* knot,
    int& kSpan, int& sMult);
  // Description:
  // Find the start idx of indicate dim in the knotArr. Because
  // we store all the knot vector in a 1-D array successively, we need
  // to locate the start idx before we process the current dimension.
  void FindKnotStart(int* lenKnot, int dim, int& knotStart) const;
  // Description:
  // To visualize NURBS patch, we need to convert it into Bezier patch.
  // This method fetch the Bezier patch shape given the id of start corner
  // control points ptAxisStart and information of the NURBS patch (The
  // NURBS patch has been inserted knots to ensure every knot will have
  // the multiplicity degree times). The returned triangle mesh will stored
  // in vtkUnstructuredGrid.
  void FetchBezierPatch(
    vtkUnstructuredGrid* shape,
    vtkDataArray* ctrlPts, int nurbsDim,
    int ctrlPtsNum[3], int ptAxisStart[3], int degree[3]);

  // Description:
  // This method insert knot value tNew into the current knot vector.
  // Old information of NURBS patch is given by pointsIn, knotsLen, knotsArr
  // and ctrlPtsNum. New information of NURBS patch is returned by pointsOut,
  // knotsLenNew, knotsArrNew and ctrlPtsNumNew
  // k means the span id. It can be provided or can be computed inside the
  // function.
  void InsertKnot(
    vtkDataArray* pointsIn, int* knotsLen, double* knotsArr, int* ctrlPtsNum,
    vtkDataArray* pointsOut, int* knotsLenNew, double* knotsArrNew, int* ctrlPtsNumNew,
    double tNew, int insertDim, int k = -1);
  // Description:
  // This method insert knot value multiple times into the current knot vector.
  void InsertKnotMulti(
    vtkDataArray* pointsIn,
    std::vector<int>& knotsLenOld,
    std::vector<double>& knotsArrOld,
    std::vector<int>& ctrlPtsNumOld,
    vtkDataArray* pointsOut,
    std::vector<int>& knotsLenNew,
    std::vector<double>& knotsArrNew,
    std::vector<int>& ctrlPtsNumNew,
    double tNew, int insertDim, int insertTimes, int kSpan = -1, int sMult = -1);
  // Description:
  // This method is a helper function to compute points or derivatives given
  // the parameter coordinates. It is done by inserting given parameters as new
  // knot degree times. This will split the NURBS patch at the given parameters
  // and points or derivatives can be computed by end points formulation.
  // More details are provided in:
  // Piegl, Les; Tiller, Wayne (1997). The NURBS Book (2. ed.).
  void InsertKnotForEval(double* Paras, int& PtIdx);

  // Description:
  // Fundamental routine to compute first order derivative in end point
  // formulation.
  void EvaluateFirstOrderDeriv(
    double Dv[3],
    int degree, double U,
    double* Pts);
  // Description:
  // Fundamental routine to compute second order derivative in end point
  // formulation.
  void EvaluateSecondOrderDeriv(
    double Dv[3],
    int degree, double* U,
    double* Pts);
  // Description:
  // Fundamental routine to compute partial derivative in end point
  // formulation.
  void EvaluateSecondPartialDeriv(
    double Dv[3],
    int* degree, double* U,
    double* Pts);

  // Description:
  // Return derivative of 1-d NURBS curve given the order d.
  void EvaluateCurvDeriv(double* Pt, double* Paras, int d);
  // Description:
  // Return highest to second derivatives given the parameters.
  // Because we need zero and first order derivative to compute
  // second order derivative, to avoid call knotInsert redundant
  // times, we return all of them. The Dv need allocate memory
  // for 3 * 3 double.
  void EvaluateCurveSecondDeriv(double* Dv, double* Paras);
  // Description:
  // Return derivative of 2-d NURBS surface given the order d.
  void EvaluateSurfaceDeriv(double* Pt, double* Paras, int d);
  // Description:
  // Return highest to second derivatives given the parameters.
  // Because we need zero and first order derivative to compute
  // second order derivative, to avoid call knotInsert redundant
  // times, we return all of them. The Dv need allocate memory
  // for 3 * 6 double. [0, 0], [1, 0], [0, 1], [2, 0], [1, 1],
  // and [0, 2].
  void EvaluateSurfaceSecondDeriv(double* Dv, double* Paras);
  // Description:
  // Return derivative of 3-d NURBS volume given the order d.
  void EvaluateVolumeDeriv(double* Pt, double* Paras, int d);
  // Description:
  // Return highest to second derivatives given the parameters.
  // Because we need zero and first order derivative to compute
  // second order derivative, to avoid call knotInsert redundant
  // times, we return all of them. The Dv need allocate memory
  // for 3 * 10 double. [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
  // [2, 0, 0], [1, 1, 0], [0, 2, 0], [1, 0, 1], [0, 1, 1], [0, 0, 2].
  // and [0, 2].
  void EvaluateVolumeSecondDeriv(double* Dv, double* Paras);

  // Description:
  // Return estimated parameter given point position for curve.
  void PointInversionCurve(double Pt[3], double* Paras);
  // Description:
  // Return estimated parameter given point position for surface.
  void PointInversionSurface(double Pt[3], double* Paras);
  // Description:
  // Return estimated parameter given point position for volume.
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