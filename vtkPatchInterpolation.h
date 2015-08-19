/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPatchInterpolation.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
//.NAME vtkPatchInterpolation ¨C Point interpolation for bezier patch
//.SECTION Description
// vtkPatchInterpolation provides methods to compute point interpolation for
// rectangular/simplicial 1-d, 2-d and 3-d bezier patch; several methods to
// generate control points for conic curves (ellipse, parabola and hyperbola);
// it also provides a the method for tessellation to visualize the bezier patch.
// Longer description of class here.


#ifndef __vtkPatchInterpolation_h
#define __vtkPatchInterpolation_h

#include "vtkObject.h"
#include "vtkFiltersBezierModule.h" // For EXPORT macro

class vtkDataArray;
class vtkVector3d;
class vtkUnstructuredGrid;

/**\brief Provide methods to interpolate and evaluate Bezier patches.
  *
  */
class VTKFILTERSBEZIER_EXPORT vtkPatchInterpolation : public vtkObject
{
public:
  static vtkPatchInterpolation* New();
  vtkTypeMacro(vtkPatchInterpolation,vtkObject);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  void InterpolateOnPatch(
    vtkDataArray* outputPt,
    int dim, vtkDataArray* P_i, int* degree, double* r);
  void InterpolateOnSimplicialPatch(
    vtkDataArray* outputPt,
    vtkDataArray* P_i, int degree, double* r, int simplex);
  void GenerateShape(
    vtkUnstructuredGrid* outputSp,
    int* samples, int dim, vtkDataArray* ctrlPts, int* degree);
  void GenerateEllipseCtrlPt(
    vtkDataArray* P_i,
    double x_axis_rad, double y_axis_rad, int quadrant);
  void GenerateEllipseCtrlPt(
    vtkDataArray* P_i,
    vtkVector3d center, vtkVector3d majorAxis, vtkVector3d minorAxis,
    int quadrant);
  void GenerateHyperbolaCtrlPt(
    vtkDataArray* P_i,
    double semi_major_axis, double eccentricity);

protected:
  vtkPatchInterpolation();
  virtual ~vtkPatchInterpolation();

  template<typename Iterator>
  void EvalCoord(
    Iterator data_iter, const Iterator data_end,
    const int& dim, const int* const degree, const double* const r,
    double* coord);
  template<typename T_i, typename T_d>
  void EvalBernstein(const T_i& n, const T_i& v, const T_d& x, T_d& b);
  template<typename Iterator>
  void EvalTriangleCoord(
    Iterator data_iter, const Iterator data_end,
    const int degree, const double r[3],
    double* coord);
  template<typename Iterator>
  void EvalTetrahedronCoord(
    Iterator data_iter, const Iterator data_end,
    const int degree, const double r[4],
    double* coord);

private:
  vtkPatchInterpolation(const vtkPatchInterpolation&); // Not implemented.
  void operator = (const vtkPatchInterpolation&); // Not implemented.
};

#endif // __vtkPatchInterpolation_h
