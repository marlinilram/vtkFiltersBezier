/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPatchInterpolation.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkPatchInterpolation.h"

#include "vtkCellType.h"
#include "vtkDataArray.h"
#include "vtkDataArrayIteratorMacro.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVector.h"

#include <math.h>
#include <vector>

vtkStandardNewMacro(vtkPatchInterpolation);

vtkPatchInterpolation::vtkPatchInterpolation()
{

}

vtkPatchInterpolation::~vtkPatchInterpolation()
{

}

void vtkPatchInterpolation::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

/**\brief Interpolate \a r on a patch defined by \a dim, \a P_i, and \a degree.
  *
  * Given parametric coordinates \a r, append the world
  * coordinates for the given patch to the \a outputPt array.
  * The control point array \a P_i is ordered as a regular
  * grid of dimension \a dim whose sides have \a degree[i]
  * points along each dimension i.
  */
void vtkPatchInterpolation::InterpolateOnPatch(
  vtkDataArray* outputPt,
  int dim, vtkDataArray* P_i, int* degree, double* r)
{
  double output_pt[4] = {0.0, 0.0, 0.0, 0.0};

  if (P_i->GetNumberOfComponents() != 4)
    {
    std::cout
      << "Number of components: " << P_i->GetNumberOfComponents() << "\n";
    }

  switch (P_i->GetDataType())
    {
    vtkDataArrayIteratorMacro(
      P_i, this->EvalCoord(vtkDABegin, vtkDAEnd, dim, degree, r, output_pt));
    }

#if DEBUG
  std::cout
    << "Point Coordinates: " <<output_pt[0]
  << "\t" << output_pt[1]
  << "\t" << output_pt[2]
  << "\t" << output_pt[3]
  << "\n";
#endif

  switch (outputPt->GetNumberOfComponents())
    {
  case 3:
    outputPt->InsertNextTuple3(
      output_pt[0]/output_pt[3],
      output_pt[1]/output_pt[3],
      output_pt[2]/output_pt[3]);
    break;
  case 4:
    outputPt->InsertNextTuple4(
      output_pt[0], output_pt[1], output_pt[2], output_pt[3]);
    break;
  default:
    vtkErrorMacro(
      "Number of components should be 3 or 4, not "
      << outputPt->GetNumberOfComponents() << ".");
    break;
    }
}

void vtkPatchInterpolation::InterpolateOnSimplicialPatch(
  vtkDataArray* outputPt,
  vtkDataArray* P_i, int degree, double* r, int simplex)
{
  double output_pt[4] = {0.0, 0.0, 0.0, 0.0};

  if (P_i->GetNumberOfComponents() != 4)
    {
    std::cout
      << "Number of components: " << P_i->GetNumberOfComponents() << "\n";
    }

  if (simplex == 3)
    {
    if (P_i->GetNumberOfTuples() != (degree + 1) * (degree + 2) / 2)
      {
      vtkErrorMacro(
        "Number of tuples should be "
        << (degree + 1) * (degree + 2) / 2
        << ", not "
        << outputPt->GetNumberOfTuples() << ".");
      return;
      }

    switch (P_i->GetDataType())
      {
      vtkDataArrayIteratorMacro(
        P_i, this->EvalTriangleCoord(vtkDABegin, vtkDAEnd, degree, r, output_pt));
      }
    }
  else if (simplex == 4)
    {
    if (P_i->GetNumberOfTuples() != 
      degree + 1 + (degree + 1) * degree / 2 + degree * (degree + 1) * (degree + 2) / 6)
      {
      vtkErrorMacro(
        "Number of tuples should be "
        << degree + 1 + (degree + 1) * degree / 2 + degree * (degree + 1) * (degree + 2) / 6
        << ", not "
        << outputPt->GetNumberOfTuples() << ".");
      return;
      }

    switch (P_i->GetDataType())
      {
      vtkDataArrayIteratorMacro(
        P_i, this->EvalTetrahedronCoord(vtkDABegin, vtkDAEnd, degree, r, output_pt));
      }
  }
  else
    {
    vtkErrorMacro(
      "Only triangle (simplex = 3) and tetrahedron (simplex = 4) acceptable,"
      << " Simplex higher than 4 is not supported.");
    return;
    }

#if DEBUG
  std::cout
    << "Point Coordinates: " <<output_pt[0]
  << "\t" << output_pt[1]
  << "\t" << output_pt[2]
  << "\t" << output_pt[3]
  << "\n";
#endif

  switch (outputPt->GetNumberOfComponents())
    {
  case 3:
    outputPt->InsertNextTuple3(
      output_pt[0]/output_pt[3],
      output_pt[1]/output_pt[3],
      output_pt[2]/output_pt[3]);
    break;
  case 4:
    outputPt->InsertNextTuple4(
      output_pt[0], output_pt[1], output_pt[2], output_pt[3]);
    break;
  default:
    vtkErrorMacro(
      "Number of components should be 3 or 4, not "
      << outputPt->GetNumberOfComponents() << ".");
    break;
    }
}

/**\brief Evaluate the contribution of a control point to a world coordinate.
  *
  */
template<typename Iterator>
void vtkPatchInterpolation::EvalCoord(
  Iterator data_iter, const Iterator data_end,
  const int& dim, const int* const degree, const double* const r,
  double* coord)
{
  // idx keeps the current k for binomial n in each dimension for example if idxs = {2, 1, 0}, it
  // means we need to evaluate the coefficient for control point P_{2, 1, 0}. Hence we need to
  // evaluate the following Bernstein polynomials B_{2, degree[0]}(r[0]), B_{1, degree[1]}(r[1]),
  // B_{0, degree[2]}(r[2]) We assume that control points P are stored as
  // - P_{0, 0, 0}, ..., P_{degree[0], 0, 0};
  // - P_{0, 1, 0}, ..., P_{degree[0], 1, 0};
  // - ...
  // - P_{0, degree[1], 0}, ..., P_{degree[0], degree[1], 0};
  // - ...
  // - ...
  // - P_{0, degree[1], degree[2]}, ..., P_{degree[0], degree[1], degree[2]};
  std::vector<int> idxs(dim, 0);
  double coef = 1.0;
  coord[0] = coord[1] = coord[2] = coord[3] = 0.0;

  while (data_iter != data_end)
    {
    coef = 1.0;
    double coef_per_dim = 0.0;
    for (size_t i = 0; i < dim; ++i)
      {
      this->EvalBernstein(degree[i], idxs[i], r[i], coef_per_dim);
      coef *= coef_per_dim;
      }

    coord[0] += (*data_iter) * coef;
    ++data_iter;
    coord[1] += (*data_iter) * coef;
    ++data_iter;
    coord[2] += (*data_iter) * coef;
    ++data_iter;
    coord[3] += (*data_iter) * coef;
    ++data_iter;

    for (int j = 0; j < dim; ++j)
      {
      ++idxs[j];
      if (idxs[j] < degree[j]+1)
        {
        break;
        }
      idxs[j] = 0;
      }
    }

  //for (; data_iter != data_end; ++data_iter)
  //{
  //    std::cout<<(*data_iter)<<"\n";
  //}
}

template<typename Iterator>
void vtkPatchInterpolation::EvalTriangleCoord(
  Iterator data_iter, const Iterator data_end,
  const int degree, const double r[3], double* coord)
{
  // i j k step style:
  // - i = 0, j = 0, k = n;
  // - i = 0, j = 1, k = n-1;
  // - ...
  // - i = 0, j = n, k = 0;
  // - i = 1, j = 0, k = n-1;
  // - ...
  // - i = 1, j = n-1, k = 0;
  double coef = 1.0;
  coord[0] = coord[1] = coord[2] = coord[3] = 0.0;

  for (int i = 0; i <= degree; ++i)
    {
    for (int j = 0; j <= degree - i; ++j)
      {
      int k = degree - i - j;
      coef =
        vtkMath::Factorial(degree) /
        (vtkMath::Factorial(i) *
        vtkMath::Factorial(j) *
        vtkMath::Factorial(k)) *
        (pow(r[0], i) * pow(r[1], j) * pow(r[2], k));

      coord[0] += (*data_iter) * coef;
      ++data_iter;
      coord[1] += (*data_iter) * coef;
      ++data_iter;
      coord[2] += (*data_iter) * coef;
      ++data_iter;
      coord[3] += (*data_iter) * coef;
      ++data_iter;
      }
    }

  if (data_iter != data_end)
    {
    vtkErrorMacro(
      "Errors: data of control points doesn't match with supposed dimension.");
    }
}

template<typename Iterator>
void vtkPatchInterpolation::EvalTetrahedronCoord(
  Iterator data_iter, const Iterator data_end,
  const int degree, const double r[4], double* coord)
{
  // i j k step style:
  // - i = 0, j = 0, k = n;
  // - i = 0, j = 1, k = n-1;
  // - ...
  // - i = 0, j = n, k = 0;
  // - i = 1, j = 0, k = n-1;
  // - ...
  // - i = 1, j = n-1, k = 0;
  double coef = 1.0;
  coord[0] = coord[1] = coord[2] = coord[3] = 0.0;

  for (int i = 0; i <= degree; ++i)
    {
    for (int j = 0; j <= degree - i; ++j)
      {
      for (int k = 0; k <= degree - i - j; ++k)
        {
        int l = degree - i - j - k;
        coef =
          vtkMath::Factorial(degree) /
          (vtkMath::Factorial(i) *
          vtkMath::Factorial(j) *
          vtkMath::Factorial(k)) *
          vtkMath::Factorial(l) *
          (pow(r[0], i) * pow(r[1], j) * pow(r[2], k) * pow(r[3], l));

        coord[0] += (*data_iter) * coef;
        ++data_iter;
        coord[1] += (*data_iter) * coef;
        ++data_iter;
        coord[2] += (*data_iter) * coef;
        ++data_iter;
        coord[3] += (*data_iter) * coef;
        ++data_iter;
        }
      }
    }

  if (data_iter != data_end)
    {
    vtkErrorMacro(
      "Errors: data of control points doesn't match with supposed dimension.");
    }
}

/**\brief Evaluate a Bernstein polynomial.
  *
  */
template<typename T_i, typename T_d>
void vtkPatchInterpolation::EvalBernstein(
  const T_i& n, const T_i& v, const T_d& x, T_d& b)
{
  b = (vtkMath::Factorial(n) / vtkMath::Factorial(v) / vtkMath::Factorial(n - v)) * pow(x, v) * pow(1 - x, n - v);
}

void vtkPatchInterpolation::GenerateShape(
  vtkUnstructuredGrid* outputSp,
  int* samples, int dim, vtkDataArray* ctrlPts, int* degree)
{
  //outputSp->InsertNextCell()
  vtkPoints* points = outputSp->GetPoints();
  if (dim == 1)
  {
    double stepX = 1 / double(samples[0]);
    double r[1] = {0};
    vtkIdType curNumPts = points->GetNumberOfPoints();
    for (int i = 0; i < samples[0] + 1; ++i)
    {
      r[0] = i * stepX;
      //std::cout<<"r[0]: "<<r[0]<<"\nr[1]: "<<r[1]<<"\n";
      this->InterpolateOnPatch(points->GetData(), dim, ctrlPts, degree, r);
    }

    std::vector<vtkIdType> ptIds((samples[0]+1), 0);
    for (int i = 0; i < samples[0] + 1; ++i)
    {
      ptIds[i] = curNumPts + i;
    }
    outputSp->InsertNextCell(VTK_POLY_LINE, ptIds.size(), &ptIds[0]);
  }
  else if (dim == 2)
  {
    double stepX = 1 / double(samples[0]);
    double stepY = 1 / double(samples[1]);
    double r[2] = {0.0,0.0};
    vtkIdType curNumPts = points->GetNumberOfPoints();
    for (int i = 0; i < samples[1] + 1; ++i)
    {
      r[1] = i * stepY;
      for (int j = 0; j < samples[0] + 1; ++j)
      {
        r[0] = j * stepX;
        //std::cout<<"r[0]: "<<r[0]<<"\nr[1]: "<<r[1]<<"\n";
        this->InterpolateOnPatch(points->GetData(), dim, ctrlPts, degree, r);
      }
    }

    std::vector<vtkIdType> ptIds(2*(samples[0]+1), 0);
    for (int i = 0; i < samples[1]; ++i)
    {
      for (int j = 0; j < samples[0] + 1; ++j)
      {
        ptIds[2*j + 0] = curNumPts + i * (samples[0] + 1) + j;
        ptIds[2*j + 1] = curNumPts + (i + 1) * (samples[0] + 1) + j;
      }
      outputSp->InsertNextCell(VTK_TRIANGLE_STRIP, ptIds.size(), &ptIds[0]);
    }
  }
  else if (dim == 3)
  {

  }
  else
  {
    vtkErrorMacro(
      "Number of dim should be 1,2 or 3, not "
      << dim << ".");
  }
}


/**\Generate control points for ellipse
  *
  * Generate control points for ellipse given
  * radius for x and y axis.
  * The center of this ellipse locates in origin.
  */
void vtkPatchInterpolation::GenerateEllipseCtrlPt(
  vtkDataArray* P_i,
  double x_axis_rad, double y_axis_rad, int quadrant)
{
  // rational function for ellipse is
  // - x = b * (1 - u .* u) ./ (1 + u .* u);
  // - y = a * 2 * u ./ (1 + u .* u); 
  // - b is y_axis-rad and a is x_axis_rad

  //vtkDoubleArray* P_i = vtkDoubleArray::New();


  P_i->SetNumberOfComponents(4);
  P_i->SetNumberOfTuples(3);
  switch (quadrant)
    {
  case 1:
    P_i->SetTuple4(0, x_axis_rad, 0, 0, 1);
    P_i->SetTuple4(1, x_axis_rad, y_axis_rad, 0, 1);
    P_i->SetTuple4(2, 0, y_axis_rad*2, 0, 2);
    break;
  case 2:
    P_i->SetTuple4(0, 0, y_axis_rad*2, 0, 2);
    P_i->SetTuple4(1, -x_axis_rad, y_axis_rad, 0, 1);
    P_i->SetTuple4(2, -x_axis_rad, 0, 0, 1);
    break;
  case 3:
    P_i->SetTuple4(0, -x_axis_rad, 0, 0, 1);
    P_i->SetTuple4(1, -x_axis_rad, -y_axis_rad, 0, 1);
    P_i->SetTuple4(2, 0, -y_axis_rad*2, 0, 2);
    break;
  case 4:
    P_i->SetTuple4(0, 0, -y_axis_rad*2, 0, 2);
    P_i->SetTuple4(1, x_axis_rad, -y_axis_rad, 0, 1);
    P_i->SetTuple4(2, x_axis_rad, 0, 0, 1);
    break;
  default:
    vtkErrorMacro(
      "Quadrant should be 1, 2, 3 or 4, not "
      << quadrant << ".");
    break;
    }
}

/**\Generate control points for ellipse
  * in 3-d space
  * Generate control points for ellipse
  * given major axis and minor axis.
  * The minor axis should be orthogonalized
  * with major axis.
  */
void vtkPatchInterpolation::GenerateEllipseCtrlPt(
  vtkDataArray* P_i,
  vtkVector3d center, vtkVector3d majorAxis, vtkVector3d minorAxis,
  int quadrant)
{
  double x_axis_rad = majorAxis.Norm();
  double y_axis_rad = minorAxis.Norm();

  // build new coordinate system
  double GS_multiplier = minorAxis.Dot(majorAxis) / x_axis_rad;
  double* minorAxis_ptr = minorAxis.GetData();
  double* majorAxis_ptr = majorAxis.GetData();
  vtkVector3d minorAxis_orth(
    minorAxis_ptr[0] - GS_multiplier * majorAxis_ptr[0],
    minorAxis_ptr[1] - GS_multiplier * majorAxis_ptr[1],
    minorAxis_ptr[2] - GS_multiplier * majorAxis_ptr[2]);
  vtkVector3d zAxis;
  vtkMath::Cross(majorAxis.GetData(), minorAxis_orth.GetData(), zAxis.GetData());
  majorAxis.Normalize();
  minorAxis_orth.Normalize();
  zAxis.Normalize();
  double coordTransMat[3][3] = {
    {majorAxis.GetX(), minorAxis.GetX(), zAxis.GetX()},
    {majorAxis.GetY(), minorAxis.GetY(), zAxis.GetY()},
    {majorAxis.GetZ(), minorAxis.GetZ(), zAxis.GetZ()}};

  // transform control points into new coordinate system process 4 patches independently
  this->GenerateEllipseCtrlPt(P_i, x_axis_rad, y_axis_rad, quadrant);
  for (vtkIdType i_pt = 0; i_pt < P_i->GetNumberOfTuples(); ++i_pt)
    {
    double* i_tuple = P_i->GetTuple4(i_pt);
    double o_coord[3];
    vtkMath::Multiply3x3(coordTransMat, i_tuple, o_coord);
    P_i->SetTuple4(
      i_pt,
      o_coord[0] + center.GetX() * i_tuple[3],
      o_coord[1] + center.GetY() * i_tuple[3],
      o_coord[2] + center.GetZ() * i_tuple[3],
      i_tuple[3]);
    }
}

/**\Generate control points for hyperbola
  *
  * Generate control points for hyperbola given
  * semi major axis and eccentricity.
  * The center of this ellipse locates in origin.
  */
void vtkPatchInterpolation::GenerateHyperbolaCtrlPt(
  vtkDataArray* P_i,
  double semi_major_axis, double eccentricity)
{
  if (eccentricity < 1.0 + 1e-4)
    {
    vtkErrorMacro(
      "Eccentricity should larger than 1.0, not"
      << eccentricity << ".");
    }
  double a = semi_major_axis;
  double b = a * sqrt(eccentricity * eccentricity - 1);
  double w_0 = 0.75;
  double w_1 = 1.25;
  double w_2 = 0.75;

  P_i->SetNumberOfComponents(4);
  P_i->SetNumberOfTuples(3);

  P_i->SetTuple4(0, w_0 * 5 * a / 3, - w_0 * 4 * b / 3, 0, w_0);
  P_i->SetTuple4(1, w_1 * 3 * a / 5, 0, 0, w_1);
  P_i->SetTuple4(2, w_2 * 5 * a / 3, w_2 * 4 * b / 3, 0, w_2);
}