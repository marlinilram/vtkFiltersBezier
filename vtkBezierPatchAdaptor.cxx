#include "vtkBezierPatchAdaptor.h"

#include "vtkControlPointArray.h"
#include "vtkDataObject.h"
#include "vtkPoints.h"

vtkBezierPatchAdaptor::vtkBezierPatchAdaptor()
{
  this->ControlPointData = NULL;
}

vtkBezierPatchAdaptor::~vtkBezierPatchAdaptor()
{
  if (this->ControlPointData
   && this->ControlPointData->GetReferenceCount() == 1)
  {
    this->ControlPointData->Delete();
  }
}

void vtkBezierPatchAdaptor::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

// ----------------------------------------------------------------------------
void vtkBezierPatchAdaptor::SetControlPointData(vtkDataObject* controlPointData)
{
  this->ControlPointData = controlPointData;
  controlPointData->Register(this);
}

// ----------------------------------------------------------------------------
int vtkBezierPatchAdaptor::GetPatchDegree(int dim) const
{
  std::vector<int> degree(this->GetPatchDimension());
  this->GetPatchDegree(&degree[0]);

  return degree[dim];
}
// ----------------------------------------------------------------------------
void vtkBezierPatchAdaptor::GetPatchParameterRange(int dim, double paramRange[2])
{
  std::vector<int> range(2 * this->GetPatchDimension());
  this->GetPatchDegree(&range[0]);
  paramRange[0] = range[2 * dim];
  paramRange[1] = range[2 * dim + 1];
}

// ----------------------------------------------------------------------------
void vtkBezierPatchAdaptor::PointInversion(
  double vtkNotUsed(Pt)[3],
  double* vtkNotUsed(Paras))
{
  vtkWarningMacro("Not implemented yet\n");
}
// ----------------------------------------------------------------------------
void vtkBezierPatchAdaptor::EvaluatePoint(
  double vtkNotUsed(Pt)[3],
  double* vtkNotUsed(Paras))
{
  vtkWarningMacro("Not implemented yet\n");
}
// ----------------------------------------------------------------------------
void vtkBezierPatchAdaptor::EvaluateDeriv(
  double* vtkNotUsed(Pt),
  double* vtkNotUsed(Paras), int vtkNotUsed(d))
{
  vtkWarningMacro("Not implemented yet\n");
}