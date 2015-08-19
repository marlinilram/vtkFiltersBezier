#include "vtkNURBSPatchAdaptor.h"

#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkIntArray.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPatchInterpolation.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"

#include <iostream>
#include <limits>
#include <vector>

vtkStandardNewMacro(vtkNURBSPatchAdaptor);

vtkNURBSPatchAdaptor::vtkNURBSPatchAdaptor()
{
  this->ControlPointData = NULL;
  this->PointsCache = vtkDoubleArray::New();
  this->CtrlPtsNumCache.resize(3, 0);
  this->KnotsLenCache.resize(3, 0);
}

vtkNURBSPatchAdaptor::~vtkNURBSPatchAdaptor()
{
  if (this->ControlPointData && 
  this->ControlPointData->GetReferenceCount() == 1)
    {
    this->ControlPointData->Delete();
    }

  if (this->PointsCache)
    {
    this->PointsCache->Delete();
    }
}

void vtkNURBSPatchAdaptor::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

void vtkNURBSPatchAdaptor::SetControlPointData(vtkDataObject* controlPointData)
{
  this->ControlPointData = controlPointData;
  controlPointData->Register(this);
}

void vtkNURBSPatchAdaptor::Begin()
{

}

bool vtkNURBSPatchAdaptor::IsAtEnd()
{
  return false;
}

bool vtkNURBSPatchAdaptor::GoToPatch(vtkIdType patchId) // random access iteration
{
  patchId = 0;
  vtkWarningMacro("Not implmented yet, set patchId = 0.");
  return false;
}

bool vtkNURBSPatchAdaptor::Next()
{
  return false;
}

// Methods to access the current patch:
int vtkNURBSPatchAdaptor::GetPatchDimension() const
{
  vtkStructuredGrid* nurbs = vtkStructuredGrid::SafeDownCast(this->ControlPointData);
  std::vector<int> ctrlPtsNum(3, 0);
  nurbs->GetDimensions(&ctrlPtsNum[0]);
  if (ctrlPtsNum[2] > 1)
    {
    return 3;
    }
  else if (ctrlPtsNum[1] > 1)
    {
    return 2;
    }
  else if (ctrlPtsNum[0] > 1)
    {
    return 1;
    }
  else
    {
    std::cout << "Dimension of NURBS can not be zero.\n";
    return 0;
    }
}
void vtkNURBSPatchAdaptor::GetPatchShape() const
// returns VTK_TRIANGLE/VTK_QUAD when PatchDimension==2, returns VTK_HEXAHEDRON/VTK_TETRA when PatchDimension==3
{

}
void vtkNURBSPatchAdaptor::GetPatchDegree(int* degreeOut) const
{
  int nurbsDim = this->GetPatchDimension();
  vtkStructuredGrid* nurbs = vtkStructuredGrid::SafeDownCast(this->ControlPointData);
  int* knotsLen = static_cast<int*>(nurbs->GetFieldData()->GetArray("KnotsLen")->GetVoidPointer(0));

  std::vector<int> ctrlPtsNum(3, 0);
  nurbs->GetDimensions(&ctrlPtsNum[0]);
  for (int i = 0; i < nurbsDim; ++i)
    {
    degreeOut[i] = knotsLen[i] - ctrlPtsNum[i] - 1;
    }
}
void vtkNURBSPatchAdaptor::GetPatchParameterRange(double* paramsOut) const
{
  vtkStructuredGrid* nurbs = vtkStructuredGrid::SafeDownCast(this->ControlPointData);

  int* knotsLenPtr = static_cast<int*>(nurbs->GetFieldData()->GetArray("KnotsLen")->GetVoidPointer(0));
  double* knotsArrPtr = static_cast<double*>(nurbs->GetFieldData()->GetArray("Knots")->GetVoidPointer(0));
  std::vector<int> knotsLen(knotsLenPtr, knotsLenPtr + 3);
  std::vector<double> knotsArr(knotsArrPtr, knotsArrPtr + knotsLen[0] + knotsLen[1] + knotsLen[2]);

  std::vector<int> ctrlPtsNum(3, 0);
  nurbs->GetDimensions(&ctrlPtsNum[0]);

  int degree[3] = {0,0,0};
  this->GetPatchDegree(degree);
  int nurbsDim = this->GetPatchDimension();

  int knotStart = 0;
  for (int i = 0; i < nurbsDim; ++i)
    {
    this->FindKnotStart(&knotsLen[0], i, knotStart);
    paramsOut[2*i] = knotsArr[knotStart + degree[i]];
    paramsOut[2*i + 1] = knotsArr[knotStart + knotsLen[i] - 1 - degree[i]];
    }
}
void vtkNURBSPatchAdaptor::GetPatchPoints(vtkPoints* pointsOut) const
{
  vtkUnstructuredGrid* nurbs = vtkUnstructuredGrid::SafeDownCast(this->ControlPointData);
  pointsOut = nurbs->GetPoints();
}

// convert nurbs to bezier for rendering
void vtkNURBSPatchAdaptor::GetPatchShape(vtkUnstructuredGrid* outputSp)
{
  // Several instructions:
  // 1. knot vector representation
  // - a knot vector [1,2,3,4,5,6,7,8] for a n-degree curve means parameter ranges [t_{n-1}, t_{m-n+1}]
  // - for historical reasons, there will be a (meaningless) additional knot at the beginning and
  //   end of a knot vector e.g. [t_0,1,2,3,4,5,6,7,8,t_9]
  // 2. convert NURBS to bezier
  // - given a full knot vector [t_0,...,t_m] (with additional knots) find the parameter range which
  //   is [t_n, t_{m-n}]
  // - insert knot until each knot in the range has multiplicity n then we will get a series of
  //   bezier curves

  vtkStructuredGrid* nurbs = vtkStructuredGrid::SafeDownCast(this->ControlPointData);

  int* knotsLen = static_cast<int*>(nurbs->GetFieldData()->GetArray("KnotsLen")->GetVoidPointer(0));
  double* knotsArr = static_cast<double*>(nurbs->GetFieldData()->GetArray("Knots")->GetVoidPointer(0));

  std::vector<int> knotsLenOld(knotsLen, knotsLen + 3);
  std::vector<int> knotsLenNew = knotsLenOld;
  std::vector<double> knotsArrOld(knotsArr, knotsArr + knotsLen[0] + knotsLen[1] + knotsLen[2]);
  std::vector<double> knotsArrNew = knotsArrOld;
  std::vector<int> ctrlPtsNumOld(3, 0);
  nurbs->GetDimensions(&ctrlPtsNumOld[0]);
  std::vector<int> ctrlPtsNumNew = ctrlPtsNumOld;
  int nurbsDim = 1;
  if (ctrlPtsNumOld[0] == 1)
    std::cout << "Dimension of NURBS can not be zero.\n";
  else if (ctrlPtsNumOld[2] > 1)
    nurbsDim = 3;
  else if (ctrlPtsNumOld[1] > 1)
    nurbsDim = 2;

  vtkSmartPointer<vtkDoubleArray> pointsIn = vtkSmartPointer<vtkDoubleArray>::New();
  pointsIn->DeepCopy(nurbs->GetPointData()->GetArray("ControlPolygon"));
  vtkSmartPointer<vtkDoubleArray> pointsOut = vtkSmartPointer<vtkDoubleArray>::New();

  int knotStart = 0;
  bool tag = true;
  std::vector<double> tNews;
  for (int i = 0; i < nurbsDim; ++i)
    {
    // for each dimension

    // 1. find parameter range and new knots value tNews
    tNews.clear();
    int degree; double paraStart; double paraEnd;
    if (tag)
      {
      degree = knotsLenOld[i] - ctrlPtsNumOld[i] - 1;
      this->FindKnotStart(&knotsLenOld[0], i, knotStart);
      paraStart = knotsArrOld[knotStart + degree];
      paraEnd = knotsArrOld[knotStart + knotsLenOld[i] - 1 - degree];
      tNews.push_back(paraStart);
      for (int j = knotStart + degree + 1; j <= knotStart + knotsLenOld[i] - 1 - degree; ++j)
        {
        if (abs(knotsArrOld[j-1] - knotsArrOld[j]) > 1e-5)
          {
          tNews.push_back(knotsArrOld[j]);
          }
        }
      }
    else
      {
      degree = knotsLenNew[i] - ctrlPtsNumNew[i] - 1;
      this->FindKnotStart(&knotsLenNew[0], i, knotStart);
      paraStart = knotsArrNew[knotStart + degree];
      paraEnd = knotsArrNew[knotStart + knotsLenNew[i] - 1 - degree];
      tNews.push_back(paraStart);
      for (int j = knotStart + degree + 1; j <= knotStart + knotsLenNew[i] - 1 - degree; ++j)
        {
        if (abs(knotsArrNew[j-1] - knotsArrNew[j]) > 1e-5)
          {
          tNews.push_back(knotsArrNew[j]);
          }
        }
      }

    // insert new knots until all knot have multiplicity degree
    double tNew;
    for (int j = 0; j < tNews.size(); ++j)
      {
      tNew = tNews[j];
      int kSpan; int sMult;
      if (tag)
        {
        // use knotsLen[i] - 1 to ignore the last additional knot
        this->FindSpanMult(
          tNew, knotsLenOld[i]-1, &knotsArrOld[knotStart],
          kSpan, sMult);

        if (degree == sMult)
          continue;

        // insert knot tNew degree - sMult times
        knotsArrNew.resize(knotsArrOld.size() + degree - sMult);
        this->InsertKnotMulti(
          pointsIn, &knotsLenOld[0], &knotsArrOld[0], &ctrlPtsNumOld[0],
          pointsOut, &knotsLenNew[0], &knotsArrNew[0], &ctrlPtsNumNew[0],
          tNew, i, degree - sMult, kSpan, sMult);

        tag = false;
        }
      else
        {
        this->FindSpanMult(
          tNew, knotsLenNew[i]-1, &knotsArrNew[knotStart],
          kSpan, sMult);

        if (degree == sMult)
          continue;

        knotsArrOld.resize(knotsArrNew.size() + degree - sMult);
        this->InsertKnotMulti(
          pointsOut, &knotsLenNew[0], &knotsArrNew[0], &ctrlPtsNumNew[0],
          pointsIn, &knotsLenOld[0], &knotsArrOld[0], &ctrlPtsNumOld[0],
          tNew, i, degree - sMult, kSpan, sMult);

        tag = true;
        }
      }
    }

  if (tag) // final data is in xxxOld
    {
    pointsOut->DeepCopy(pointsIn);
    knotsLenNew = knotsLenOld;
    knotsArrNew = knotsArrOld;
    ctrlPtsNumNew = ctrlPtsNumOld;
    }

  //vtkSmartPointer<vtkPoints> outputSpPts = vtkSmartPointer<vtkPoints>::New();
  //outputSp->SetPoints(outputSpPts);
  //outputSp->Allocate();

  // outputSp should initialize its memory allocation outside of the function, for it would happen
  // that the whole shape consists of several nurbs surfaces

  int degree[3] = {0,0,0};
  int ptAxisStart[3] = {0,0,0};
  int patchAxisNum[3] = {1,1,1};
  knotStart = 0;
  for (int i = 0; i < nurbsDim; ++i)
    {
    degree[i] = knotsLenNew[i] - ctrlPtsNumNew[i] - 1;

    int firstPtIdx = this->FindSpan(
      knotsArrNew[knotStart + degree[i]],
      knotsLenNew[i]-1, &knotsArrNew[knotStart]) - degree[i];
    int lastPtIdx = this->FindSpan(
      knotsArrNew[knotStart + knotsLenNew[i] - 1 - degree[i]],
      knotsLenNew[i]-1, &knotsArrNew[knotStart]) - degree[i];

    ptAxisStart[i] = firstPtIdx;
    patchAxisNum[i] = (lastPtIdx - firstPtIdx) / degree[i];

    knotStart += knotsLenNew[i];
    }

  int ptAxisCurStart[3] = {ptAxisStart[0], ptAxisStart[1], ptAxisStart[2]};
  for (int k = 0; k < patchAxisNum[2]; ++k)
    {
    ptAxisCurStart[1] = ptAxisStart[1];
    for (int j = 0; j < patchAxisNum[1]; ++j)
      {
      ptAxisCurStart[0] = ptAxisStart[0];
      for (int i = 0; i < patchAxisNum[0]; ++i)
        {
        this->FetchBezierPatch(
          outputSp,
          pointsOut, nurbsDim,
          &ctrlPtsNumNew[0], ptAxisCurStart, degree);
        ptAxisCurStart[0] += degree[0];
        }
      ptAxisCurStart[1] += degree[1];
      }
    ptAxisCurStart[2] += degree[2];
    }
}

// Point Inversion related method
void vtkNURBSPatchAdaptor::PointInversion(double Pt[3], double* Paras)
{
  vtkStructuredGrid* nurbs = vtkStructuredGrid::SafeDownCast(this->ControlPointData);

  std::vector<int> ctrlPtsNum(3, 0);
  nurbs->GetDimensions(&ctrlPtsNum[0]);

  if (ctrlPtsNum[2] > 1)
    {
    this->PointInversionVolume(Pt, Paras);
    }
  else if (ctrlPtsNum[1] > 1)
    {
    this->PointInversionSurface(Pt, Paras);
    }
  else if (ctrlPtsNum[0] > 1)
    {
    this->PointInversionCurve(Pt, Paras);
    }
  else
    {
    std::cout << "Dimension of NURBS can not be zero.\n";
    }
}

// Point evaluation related method
void vtkNURBSPatchAdaptor::EvaluatePoint(double Pt[3], double* Paras)
{
  // we insert the given Paras to original knot vectors until the multiplicity reaches degree times
  // and get the corresponding control point it will be the Point on the patch

  int PtIdx = 1;

  this->InsertKnotForEval(Paras, PtIdx);

  // get the control point corresponding to the paras
  double* Ptw = this->PointsCache->GetTuple4(PtIdx - 1);
  Pt[0] = Ptw[0] / Ptw[3];
  Pt[1] = Ptw[1] / Ptw[3];
  Pt[2] = Ptw[2] / Ptw[3];
}
void vtkNURBSPatchAdaptor::EvaluateDeriv(double Pt[3], double* Paras, int d)
{
  if (d >= 3)
  {
    vtkErrorMacro(
      "Only allowed first and second derivatives now.\n");
  }

  int nurbsDim = this->GetPatchDimension();
  switch (nurbsDim)
  {
  case 1:
    this->EvaluateCurvDeriv(Pt, Paras, d);
    break;
  case 2:
    this->EvaluateSurfaceDeriv(Pt, Paras, d);
    break;
  case 3:
    this->EvaluateVolumeDeriv(Pt, Paras, d);
    break;
  default:
    break;
  }
}

// Internal computation method for find span, multiplicity and insertion
void vtkNURBSPatchAdaptor::FindKnotStart(int* lenKnot, int dim, int& knotStart) const
{
  knotStart = 0;
  for (int i = 0; i < dim; ++i)
    {
    knotStart += lenKnot[i];
    }
}
int vtkNURBSPatchAdaptor::FindSpan(double tNew, int lenKnot, double* knot)
{
  // Find where to insert the new knot.
  if (tNew < knot[0])
    {
    return -1;
    }
  if (tNew >= knot[lenKnot-1])
    {
    return lenKnot - 1;
    }

  int low = 0;
  int high = lenKnot - 1;
  int mid = (low + high) / 2;
  while (tNew < knot[mid] || tNew >= knot[mid+1])
    {
    if (tNew < knot[mid])
      {
      high = mid;
      }
    else
      {
      low = mid;
      }
    mid = (low + high) / 2;
    }
  return mid;
}
void vtkNURBSPatchAdaptor::FindSpanMult(
  double tNew, int lenKnot, double* knot,
  int& kSpan, int& sMult)
{
  kSpan = this->FindSpan(tNew, lenKnot, knot);

  sMult = 0;
  for (int j = kSpan; j >= 1 && abs(knot[j]-tNew)<=1e-5; --j)
    {
    // use j >= 1 to ignore the first additional knot
    ++sMult;
    }
}

void vtkNURBSPatchAdaptor::FetchBezierPatch(
  vtkUnstructuredGrid* shape,
  vtkDataArray* ctrlPts, int nurbsDim,
  int ctrlPtsNum[3], int ptAxisStart[3], int degree[3])
{
  vtkSmartPointer<vtkPatchInterpolation> patchInterpolation =
    vtkSmartPointer<vtkPatchInterpolation>::New();

  vtkSmartPointer<vtkDoubleArray> bezierCtrlPts
    = vtkSmartPointer<vtkDoubleArray>::New();
  bezierCtrlPts->SetNumberOfComponents(4);

  int sample[3] = {20,20,20};

  for (int k = ptAxisStart[2]; k <=  ptAxisStart[2] + degree[2]; ++k)
    {
    for (int j = ptAxisStart[1]; j <=  ptAxisStart[1] + degree[1]; ++j)
      {
      for (int i = ptAxisStart[0]; i <=  ptAxisStart[0] + degree[0]; ++i)
        {
        // point (i,j,k)
        vtkIdType ptIdx =
          i +
          j * ctrlPtsNum[0] +
          k * ctrlPtsNum[1] * ctrlPtsNum[2];
        double* pt = ctrlPts->GetTuple(ptIdx);
        bezierCtrlPts->InsertNextTuple4(
          pt[0], pt[1], pt[2], pt[3]);
        }
      }
    }

  patchInterpolation->GenerateShape(shape, sample, nurbsDim, bezierCtrlPts, degree);
}

void vtkNURBSPatchAdaptor::InsertKnot(
  vtkDataArray* pointsIn, int* knotsLen, double* knotsArr, int* ctrlPtsNum,
  vtkDataArray* pointsOut, int* knotsLenNew, double* knotsArrNew, int* ctrlPtsNumNew,
  double tNew, int insertDim, int k)
{
  if (insertDim >= 3 || insertDim < 0)
    {
    vtkErrorMacro(
      "Dimension should between 0-2, not"
      << insertDim << ".");
    }

  // find where knot vector start based on the dimension we want in insert for
  int knot_start = 0;
  for (int i = 0; i < insertDim; ++i)
    {
    knot_start += knotsLen[i];
    }

  // find the knot span where the new knot is
  if (k == -1)
    k = FindSpan(tNew, knotsLen[insertDim], &knotsArr[knot_start]);

  // copy old knots array to new knots array
  std::memcpy(knotsLenNew, knotsLen, 3*sizeof(knotsLen));
  knotsLenNew[insertDim] += 1;
  std::memcpy(knotsArrNew, knotsArr, (knot_start+k+1)*sizeof(knotsArr));
  knotsArrNew[knot_start+k+1] = tNew;
  std::memcpy(knotsArrNew+knot_start+k+2, knotsArr+knot_start+k+1,
    (knotsLen[0]+knotsLen[1]+knotsLen[2]-knot_start-k-1)*sizeof(knotsArr));

  // copy old control points to new control points
  std::memcpy(ctrlPtsNumNew, ctrlPtsNum, 3*sizeof(ctrlPtsNum));
  ctrlPtsNumNew[insertDim] += 1;
  vtkIdType totalCPoints = ctrlPtsNumNew[0]*ctrlPtsNumNew[1]*ctrlPtsNumNew[2];
  pointsOut->SetNumberOfComponents(4);
  pointsOut->SetNumberOfTuples(totalCPoints);

  int p = knotsLen[insertDim] - ctrlPtsNum[insertDim] - 1;
  int ptIter[3] = {0,0,0};
  for (; ptIter[2] < ctrlPtsNumNew[2]; ++ptIter[2])
    {
    for (; ptIter[1] < ctrlPtsNumNew[1]; ++ptIter[1])
      {
      for (; ptIter[0] < ctrlPtsNumNew[0]; ++ptIter[0])
        {
        // compute alpha
        double alpha;
        int i = ptIter[insertDim];
        vtkIdType ptIdxNew =
          ptIter[2]*ctrlPtsNumNew[1]*ctrlPtsNumNew[0] +
          ptIter[1]*ctrlPtsNumNew[0] +
          ptIter[0];

        if (i <= k - p)
          {
          vtkIdType ptIdxOld =
            ptIter[2]*ctrlPtsNum[1]*ctrlPtsNum[0] +
            ptIter[1]*ctrlPtsNum[0] +
            ptIter[0];

          double* ptOld = pointsIn->GetTuple(ptIdxOld);
          pointsOut->SetTuple4(
            ptIdxNew,
            ptOld[0], ptOld[1], ptOld[2], ptOld[3]);
          }
        else if (i <= k && i >= k - p + 1)
          {
          alpha = (tNew - knotsArr[knot_start + i])
            / (knotsArr[knot_start + i + p] - knotsArr[knot_start + i]);

          std::cout<<"i:\t"<<i<<"\talpha:\t"<<alpha<<"\n";

          vtkIdType ptIdxOld =
            ptIter[2]*ctrlPtsNum[1]*ctrlPtsNum[0] +
            ptIter[1]*ctrlPtsNum[0] +
            ptIter[0];
          --ptIter[insertDim];
          vtkIdType ptIdxOldPrev =
            ptIter[2]*ctrlPtsNum[1]*ctrlPtsNum[0] +
            ptIter[1]*ctrlPtsNum[0] +
            ptIter[0];
          ++ptIter[insertDim];

          double ptOld[4];
          pointsIn->GetTuple(ptIdxOld, ptOld);
          double ptOldPrev[4];
          pointsIn->GetTuple(ptIdxOldPrev, ptOldPrev);
          pointsOut->SetTuple4(
            ptIdxNew,
            alpha * ptOld[0] + (1 - alpha) * ptOldPrev[0],
            alpha * ptOld[1] + (1 - alpha) * ptOldPrev[1],
            alpha * ptOld[2] + (1 - alpha) * ptOldPrev[2],
            alpha * ptOld[3] + (1 - alpha) * ptOldPrev[3]);

          std::cout<<"Old Pt "<<ptIdxOld<<":\t"<<ptOld[0]<<"\t"<<ptOld[1]<<"\t"<<ptOld[2]<<"\t"<<ptOld[3]<<"\n";
          std::cout<<"Old Pt Prev"<<ptIdxOldPrev<<":\t"<<ptOldPrev[0]<<"\t"<<ptOldPrev[1]<<"\t"<<ptOldPrev[2]<<"\t"<<ptOldPrev[3]<<"\n";
          }
        else
          {
          --ptIter[insertDim];
          vtkIdType ptIdxOld =
            ptIter[2]*ctrlPtsNum[1]*ctrlPtsNum[0] +
            ptIter[1]*ctrlPtsNum[0] +
            ptIter[0];
          ++ptIter[insertDim];

          double* ptOld = pointsIn->GetTuple(ptIdxOld);
          pointsOut->SetTuple4(
            ptIdxNew,
            ptOld[0], ptOld[1], ptOld[2], ptOld[3]);
          }
        }
      }
    }
}
void vtkNURBSPatchAdaptor::InsertKnotMulti(
  vtkDataArray* pointsIn, int* knotsLenOld, double* knotsArrOld, int* ctrlPtsNumOld,
  vtkDataArray* pointsOut, int* knotsLenNew, double* knotsArrNew, int* ctrlPtsNumNew,
  double tNew, int insertDim, int insertTimes, int kSpan, int sMult)
{
  if (insertDim >= 3 || insertDim < 0)
    {
    vtkErrorMacro(
      "Dimension should between 0-2, not"
      << insertDim << ".");
    }

  // find where knot vector start based on the dimension we want in insert for
  int knot_start = 0;
  for (int i = 0; i < insertDim; ++i)
    {
    knot_start += knotsLenOld[i];
    }

  // find the knot span where the new knot is
  if (kSpan == -1)
    kSpan = FindSpan(tNew, knotsLenOld[insertDim], &knotsArrOld[knot_start]);
  if (sMult == -1)
    {
    sMult = 0; // multiplicity of tNew
    for (int i = kSpan; i >= 0 && abs(knotsArrOld[knot_start+i]-tNew)<=1e-5; --i)
      {
      ++sMult;
      }
    }
  //std::cout<<"kSpan: "<<kSpan<<"\tsMult: "<<sMult<<"\n";

  // copy old knots array to new knots array
  std::memcpy(knotsLenNew, knotsLenOld, 3*sizeof(knotsLenOld));
  knotsLenNew[insertDim] += insertTimes;
  std::memcpy(knotsArrNew, knotsArrOld, (knot_start+kSpan+1)*sizeof(knotsArrOld));
  for (int i = 1; i <= insertTimes; ++i)
    {
    knotsArrNew[knot_start+kSpan+i] = tNew;
    }
  std::memcpy(knotsArrNew+knot_start+kSpan+insertTimes+1, knotsArrOld+knot_start+kSpan+1,
    (knotsLenOld[0]+knotsLenOld[1]+knotsLenOld[2]-knot_start-kSpan-1)*sizeof(knotsArrOld));

  // copy old control points to new control points
  std::memcpy(ctrlPtsNumNew, ctrlPtsNumOld, 3*sizeof(ctrlPtsNumOld));
  ctrlPtsNumNew[insertDim] += insertTimes;
  vtkIdType totalCPoints = ctrlPtsNumNew[0]*ctrlPtsNumNew[1]*ctrlPtsNumNew[2];
  pointsOut->SetNumberOfComponents(4);
  pointsOut->SetNumberOfTuples(totalCPoints);


  // To iterate the axis based on which dimension we want to insert knots we define step length to
  // help read data in different way For example, if insertDim = 1 which is the Y axis we hope to
  // iterate the point data like this: first from (0,0,0) to (0,maxIdY,0) then (1,0,0) to
  // (1,maxIdY,0) ... until (maxIdX,0,0) to (maxIdX,maxIdY,0) then axis Z move one step (0,0,1) to
  // (maxIdX,maxIdY,1) So here majorStep is for axis X for it moves faster minorStep is for axis Z
  // for it moves slower insertStep is for axis Y for which we insert knots to make the code concise
  // use int ptStep[3] to store the variables ptStep[0] stores insertStep ptStep[1] stores majorStep
  // ptStep[2] stores minorStep
  int ptStepNew[3]; int ptMaxNew[3];
  int ptStepOld[3]; int ptMaxOld[3];
  switch (insertDim)
    {
  case 0:
    // insert from axis X
    ptStepNew[0] = 1; ptStepNew[1] = ctrlPtsNumNew[0]; ptStepNew[2] = ctrlPtsNumNew[0]*ctrlPtsNumNew[1];
    ptStepOld[0] = 1; ptStepOld[1] = ctrlPtsNumOld[0]; ptStepOld[2] = ctrlPtsNumOld[0]*ctrlPtsNumOld[1];

    ptMaxNew[0] = ctrlPtsNumNew[0]; ptMaxNew[1] = ctrlPtsNumNew[1]; ptMaxNew[2] = ctrlPtsNumNew[2];
    ptMaxOld[0] = ctrlPtsNumOld[0]; ptMaxNew[1] = ctrlPtsNumOld[1]; ptMaxNew[2] = ctrlPtsNumOld[2];
    break;
  case 1:
    // insert from axis Y
    ptStepNew[0] = ctrlPtsNumNew[0]; ptStepNew[1] = 1; ptStepNew[2] = ctrlPtsNumNew[0]*ctrlPtsNumNew[1];
    ptStepOld[0] = ctrlPtsNumOld[0]; ptStepOld[1] = 1; ptStepOld[2] = ctrlPtsNumOld[0]*ctrlPtsNumOld[1];

    ptMaxNew[0] = ctrlPtsNumNew[1]; ptMaxNew[1] = ctrlPtsNumNew[0]; ptMaxNew[2] = ctrlPtsNumNew[2];
    ptMaxOld[0] = ctrlPtsNumOld[1]; ptMaxNew[1] = ctrlPtsNumOld[0]; ptMaxNew[2] = ctrlPtsNumOld[2];
    break;
  case 2:
    // insert from axis Z
    ptStepNew[0] = ctrlPtsNumNew[0]*ctrlPtsNumNew[1]; ptStepNew[1] = 1; ptStepNew[2] = ctrlPtsNumNew[0];
    ptStepOld[0] = ctrlPtsNumOld[0]*ctrlPtsNumOld[1]; ptStepOld[1] = 1; ptStepOld[2] = ctrlPtsNumOld[0];

    ptMaxNew[0] = ctrlPtsNumNew[2]; ptMaxNew[1] = ctrlPtsNumNew[0]; ptMaxNew[2] = ctrlPtsNumNew[1];
    ptMaxOld[0] = ctrlPtsNumOld[2]; ptMaxNew[1] = ctrlPtsNumOld[0]; ptMaxNew[2] = ctrlPtsNumOld[1];
    break;
  default:
    break;
    }
  //std::cout<<"ptMaxNew[1]: "<<ptMaxNew[1]<<"\nptMaxNew[2]: "<<ptMaxNew[2]<<"\n";
  int p = knotsLenOld[insertDim] - ctrlPtsNumOld[insertDim] - 1;
  std::vector<std::vector<double>> Rw(p+1, std::vector<double>(4,0));
  std::vector<std::vector<double>> alpha(p-sMult, std::vector<double>(insertTimes+1, 0));
  //std::vector<std::vector<double>> alfa(ptMaxOld[0], std::vector<double>(ptMaxNew[0], 0));
  // implement oslo algorithm here
  // compute alfa
  //for (int i = 0; i < )
  //{
  //}


  // store alpha
  for (int j = 1; j <= insertTimes; ++j)
    {
    int L = kSpan - p + j;
    for (int i = 0; i <= p - j - sMult; ++i)
      {
      alpha[i][j] =
        (tNew - knotsArrOld[knot_start + L + i]) /
        (knotsArrOld[knot_start + kSpan + i + 1] - knotsArrOld[knot_start + L + i]);
      }
    }
  // compute new control points
  for (int i = 0; i < ptMaxNew[2]; ++i)
    {
    for (int j = 0; j < ptMaxNew[1]; ++j)
      {
      for (int k = 0; k <= kSpan - p; ++k)
        {
        vtkIdType ptIdxNew = k*ptStepNew[0] + j*ptStepNew[1] + i*ptStepNew[2];
        vtkIdType ptIdxOld = k*ptStepOld[0] + j*ptStepOld[1] + i*ptStepOld[2];
        double* ptOld = pointsIn->GetTuple(ptIdxOld);
        pointsOut->SetTuple4(
          ptIdxNew,
          ptOld[0], ptOld[1], ptOld[2], ptOld[3]);
        }
      //std::cout<<"test1\n";
      for (int k = kSpan - sMult; k < ptMaxOld[0]; ++k)
        {
        vtkIdType ptIdxNew = (k+insertTimes)*ptStepNew[0] + j*ptStepNew[1] + i*ptStepNew[2];
        vtkIdType ptIdxOld = k*ptStepOld[0] + j*ptStepOld[1] + i*ptStepOld[2];
        double* ptOld = pointsIn->GetTuple(ptIdxOld);
        pointsOut->SetTuple4(
          ptIdxNew,
          ptOld[0], ptOld[1], ptOld[2], ptOld[3]);
        }
      //std::cout<<"test2\n";
      for (int k = 0; k <= p - sMult; ++k)
        {
        vtkIdType ptIdxOld = (kSpan - p + k)*ptStepOld[0] + j*ptStepOld[1] + i*ptStepOld[2];
        double* ptOld = pointsIn->GetTuple(ptIdxOld);
        Rw[k][0] = ptOld[0]; Rw[k][1] = ptOld[1]; Rw[k][2] = ptOld[2]; Rw[k][3] = ptOld[3];
        }

      int L = 0;
      for (int kk = 1; kk <= insertTimes; ++kk)
        {
        L = kSpan - p + kk;
        for (int k = 0; k <= p - kk - sMult; ++k)
          {
          //double alpha =
          //  (tNew - knotsArrOld[knot_start + L + k]) /
          //  (knotsArrOld[knot_start + k + kSpan + 1] - knotsArrOld[knot_start + L + k]);
          Rw[k][0] = alpha[k][kk]*Rw[k + 1][0] + (1.0 - alpha[k][kk])*Rw[k][0];
          Rw[k][1] = alpha[k][kk]*Rw[k + 1][1] + (1.0 - alpha[k][kk])*Rw[k][1];
          Rw[k][2] = alpha[k][kk]*Rw[k + 1][2] + (1.0 - alpha[k][kk])*Rw[k][2];
          Rw[k][3] = alpha[k][kk]*Rw[k + 1][3] + (1.0 - alpha[k][kk])*Rw[k][3];
          }
        vtkIdType ptIdxNew = L*ptStepNew[0] + j*ptStepNew[1] + i*ptStepNew[2];
        pointsOut->SetTuple4(ptIdxNew, Rw[0][0], Rw[0][1], Rw[0][2], Rw[0][3]);
        ptIdxNew = (kSpan + insertTimes - kk - sMult)*ptStepNew[0] + j*ptStepNew[1] + i*ptStepNew[2];
        pointsOut->SetTuple4(
          ptIdxNew,
          Rw[p - kk - sMult][0],
          Rw[p - kk - sMult][1],
          Rw[p - kk - sMult][2],
          Rw[p - kk - sMult][3]);
        }
      //std::cout<<"test3\n";
      for (int k = L + 1; k < kSpan - sMult; ++k)
        {
        vtkIdType ptIdxNew = k*ptStepNew[0] + j*ptStepNew[1] + i*ptStepNew[2];
        pointsOut->SetTuple4(
          ptIdxNew,
          Rw[k - L][0],
          Rw[k - L][1],
          Rw[k - L][2],
          Rw[k - L][3]);
        }
      //std::cout<<"test4\n";
      //std::cout<<"i: "<<i<<"\nj: "<<j<<"\n";
      }
    }
}
void vtkNURBSPatchAdaptor::InsertKnotForEval(double* Paras, int& PtIdx)
{
  vtkStructuredGrid* nurbs = vtkStructuredGrid::SafeDownCast(this->ControlPointData);

  int* knotsLenPtr = static_cast<int*>(nurbs->GetFieldData()->GetArray("KnotsLen")->GetVoidPointer(0));
  double* knotsArrPtr = static_cast<double*>(nurbs->GetFieldData()->GetArray("Knots")->GetVoidPointer(0));
  std::vector<int> knotsLen(knotsLenPtr, knotsLenPtr + 3);
  std::vector<double> knotsArr(knotsArrPtr, knotsArrPtr + knotsLen[0] + knotsLen[1] + knotsLen[2]);

  std::vector<int> ctrlPtsNum(3, 0);
  nurbs->GetDimensions(&ctrlPtsNum[0]);

  vtkSmartPointer<vtkDoubleArray> pointsIn = vtkSmartPointer<vtkDoubleArray>::New();
  pointsIn->DeepCopy(nurbs->GetPointData()->GetArray("ControlPolygon"));

  // 1. get dimension and degree for each dimension of the nurbs patch
  int nurbDim = this->GetPatchDimension();
  int degree[3] = {0, 0, 0};
  this->GetPatchDegree(degree);

  // 2. insert given paras to each dimension multiple times
  PtIdx = 1;
  int kSpan, sMult;
  int knotStart;

  for (int i = 0; i < nurbDim; ++i)
    {
    this->FindKnotStart(&knotsLen[0], i, knotStart);
    this->FindSpanMult(
      Paras[i], knotsLen[i]-1, &knotsArr[knotStart],
      kSpan, sMult);
    PtIdx *= (kSpan - sMult + 1);

    this->KnotsArrCache.resize(
      knotsLen[0] + knotsLen[1] + knotsLen[2]
    + degree[i] - sMult);
    this->InsertKnotMulti(
      pointsIn, &knotsLen[0], &knotsArr[0], &ctrlPtsNum[0],
      this->PointsCache,
      &this->KnotsLenCache[0],
      &this->KnotsArrCache[0],
      &this->CtrlPtsNumCache[0],
      Paras[i], i, degree[i] - sMult, kSpan, sMult);

    pointsIn->DeepCopy(this->PointsCache);
    knotsLen = this->KnotsLenCache;
    knotsArr = this->KnotsArrCache;
    ctrlPtsNum = this->CtrlPtsNumCache;
    }
}

void vtkNURBSPatchAdaptor::EvaluateFirstOrderDeriv(
  double Dv[3],
  int degree, double U, double* Pts)
{
  // first order derivative needs U_p_1 and Pw_0, Pw_1, so Pts has 8 elements
  double coef =
    degree / U * Pts[7] / Pts[3];
  Dv[0] = coef * (Pts[4] / Pts[7] - Pts[0] / Pts[3]);
  Dv[1] = coef * (Pts[5] / Pts[7] - Pts[1] / Pts[3]);
  Dv[2] = coef * (Pts[6] / Pts[7] - Pts[2] / Pts[3]);
}
void vtkNURBSPatchAdaptor::EvaluateSecondOrderDeriv(
  double Dv[3],
  int degree, double* U, double* Pts)
{
  // second order derivative needs U_p_1, U_p_2 and Pw_0, Pw_1, Pw_2, so Pts has 12 elements

  double coef[3];
  coef[0] =
    degree * (degree - 1) / U[0] / U[0];
  coef[1] =
    degree * (degree - 1) / U[0]
  * (U[0] + U[1]) / (U[0] * U[1]);
  coef[2] =
    degree * (degree - 1) / U[0] / U[1];
  double A_2_0[3]; // A_2_0 means A''(0)
  A_2_0[0] =
    coef[0] * Pts[0]
  + coef[1] * Pts[4]
  + coef[2] * Pts[8];
  A_2_0[1] =
    coef[0] * Pts[1]
  + coef[1] * Pts[5]
  + coef[2] * Pts[9];
  A_2_0[2] =
    coef[0] * Pts[2]
  + coef[1] * Pts[6]
  + coef[2] * Pts[10];

  // w_1_0 means w'(0)
  double w_1_0 =
    degree / U[0] * (Pts[7] - Pts[3]);
  // w_2_0 means w''(0)
  double w_2_0 =
    coef[0] * Pts[3]
  + coef[1] * Pts[7]
  + coef[2] * Pts[11];

  double Dv_first[3];
  this->EvaluateFirstOrderDeriv(Dv_first, degree, U[0], Pts);

  //double d_zero_order[3];
  double Dv_zero[3];
  Dv_zero[0] = Pts[0] / Pts[3];
  Dv_zero[1] = Pts[1] / Pts[3];
  Dv_zero[2] = Pts[2] / Pts[3];

  Dv[0] =
    (
    A_2_0[0]
  - 2 * w_1_0 * Dv_first[0]
  - w_2_0 * Dv_zero[0])
    / Pts[3];

  Dv[1] =
    (
    A_2_0[1]
  - 2 * w_1_0 * Dv_first[1]
  - w_2_0 * Dv_zero[1])
    / Pts[3];

  Dv[2] =
    (
    A_2_0[2]
  - 2 * w_1_0 * Dv_first[2]
  - w_2_0 * Dv_zero[2])
    / Pts[3];
}
void vtkNURBSPatchAdaptor::EvaluateSecondPartialDeriv(
  double Dv[3],
  int* degree, double* U, double* Pts)
{
  // for second partial derivative needs P_0_0, P_0_1, P_1_0 and P_1_1
  double coef[3];
  coef[0] =
    degree[0] * degree[1] / Pts[3] / U[0] / U[1] * Pts[15];
  coef[1] =
    degree[0] * degree[1] / Pts[3] / U[0] / U[1]
  * Pts[7] * Pts[11] / Pts[3];
  coef[2] =
    degree[0] * degree[1] / Pts[3] / U[0] / U[1]
  * (2 * Pts[7] * Pts[11] / Pts[3] - Pts[15]);

  Dv[0] =
    coef[0] * Pts[12] / Pts[15]
  - coef[1] * (Pts[4] / Pts[7] + Pts[8] / Pts[11])
    + coef[2] * Pts[0] / Pts[3];

  Dv[1] =
    coef[0] * Pts[13] / Pts[15]
  - coef[1] * (Pts[5] / Pts[7] + Pts[9] / Pts[11])
    + coef[2] * Pts[1] / Pts[3];

  Dv[2] =
    coef[0] * Pts[14] / Pts[15]
  - coef[1] * (Pts[6] / Pts[7] + Pts[10] / Pts[11])
    + coef[2] * Pts[2] / Pts[3];
}

void vtkNURBSPatchAdaptor::EvaluateCurvDeriv(double* Pt, double* Paras, int d)
{
  std::vector<double> Dv(9, 0);
  this->EvaluateCurveSecondDeriv(&Dv[0], Paras);
  // Dv stores zero to second order derivatives successively
  Pt[0] = Dv[3 * d + 0];
  Pt[1] = Dv[3 * d + 1];
  Pt[2] = Dv[3 * d + 2];
}
void vtkNURBSPatchAdaptor::EvaluateCurveSecondDeriv(double* Dv, double* Paras)
{
  // currently only implement first order and second order derivatives
  // Dv should have space for 9 elements
  int PtIdx = 1;

  this->InsertKnotForEval(Paras, PtIdx);

  std::vector<int> degree(3, 0);
  this->GetPatchDegree(&degree[0]);

  int knotStart;
  int kSpan, sMult;
  this->FindKnotStart(&this->KnotsLenCache[0], 0, knotStart);
  this->FindSpanMult(
    Paras[0],
    this->KnotsLenCache[0] - 1,
    &this->KnotsArrCache[knotStart],
    kSpan, sMult);

  std::vector<double> U(2, 0.0);
  U[0] = this->KnotsArrCache[knotStart + kSpan + 1];
  U[1] = this->KnotsArrCache[knotStart + kSpan + 2];

  std::vector<double> Pts(12, 0.0);
  this->PointsCache->GetTuple(PtIdx - 1, &Pts[0]);
  this->PointsCache->GetTuple(PtIdx, &Pts[4]);
  this->PointsCache->GetTuple(PtIdx + 1, &Pts[8]);

  // store zero to d order derivatives to Dv Dv should have space for 3 * (d + 1)
  Dv[0] = Pts[0] / Pts[3];
  Dv[1] = Pts[1] / Pts[3];
  Dv[2] = Pts[2] / Pts[3];

  this->EvaluateFirstOrderDeriv(&Dv[3], degree[0], U[0], &Pts[0]);
  this->EvaluateSecondOrderDeriv(&Dv[6], degree[0], &U[0], &Pts[0]);
}
void vtkNURBSPatchAdaptor::EvaluateSurfaceDeriv(double* Pt, double* Paras, int d)
{
  std::vector<double> Dv(18, 0);
  this->EvaluateSurfaceSecondDeriv(&Dv[0], Paras);
  // Todo: return indicated derivative for order d
  Pt = Pt;
  d = d;
}
void vtkNURBSPatchAdaptor::EvaluateSurfaceSecondDeriv(double* Dv, double* Paras)
{
  // currently only implement first order and second order derivatives the Dv will contains 3 * 6 doubles
  // - P_1 Dv_zero namely the point position
  // - P_2 - P_3 First and second order derivatives for parametric space u Dv_1_u Dv_2_u
  // - P_4 - P_5 First and second order derivatives for parametric space v Dv_1_v Dv_2_v
  // - P_6 Partial derivatives for parametric space u and v Dv_1_u_v
  int PtIdx = 1;

  this->InsertKnotForEval(Paras, PtIdx);

  std::vector<int> degree(3, 0);
  this->GetPatchDegree(&degree[0]);

  std::vector<int> knotStart(2, 0);
  std::vector<int> kSpan(2, 0);
  std::vector<int> sMult(2, 0);
  std::vector<double> U(2, 0.0);
  std::vector<double> Pts(16, 0.0);
  int axisStep[2] = {1, this->CtrlPtsNumCache[0]};
  // we need to retrieve control points along different axis

  // point pos
  this->PointsCache->GetTuple(PtIdx - 1, &Pts[0]);
  Dv[0] = Pts[0] / Pts[3];
  Dv[1] = Pts[1] / Pts[3];
  Dv[2] = Pts[2] / Pts[3];

  // first order and second order derivative
  for (int i = 0; i < 2; ++i)
  {
    this->FindKnotStart(&this->KnotsLenCache[0], i, knotStart[i]);
    this->FindSpanMult(
      Paras[i],
      this->KnotsLenCache[i] - 1,
      &this->KnotsArrCache[knotStart[i]],
      kSpan[i], sMult[i]);

    U[0] = this->KnotsArrCache[knotStart[i] + kSpan[i] + 1];
    U[1] = this->KnotsArrCache[knotStart[i] + kSpan[i] + 2];

    this->PointsCache->GetTuple(PtIdx - 1 + axisStep[i], &Pts[4]);
    this->PointsCache->GetTuple(PtIdx - 1 + 2 * axisStep[i], &Pts[8]);

    this->EvaluateFirstOrderDeriv(&Dv[3 + i * 6], degree[i], U[0], &Pts[0]);
    this->EvaluateSecondOrderDeriv(&Dv[6 + i * 6], degree[i], &U[0], &Pts[0]);
  }

  // partial derivative
  U[0] = this->KnotsArrCache[knotStart[0] + kSpan[0] + 1];
  U[1] = this->KnotsArrCache[knotStart[1] + kSpan[1] + 1];

  this->PointsCache->GetTuple(PtIdx - 1 + axisStep[0], &Pts[4]);
  this->PointsCache->GetTuple(PtIdx - 1 + axisStep[1], &Pts[8]);
  this->PointsCache->GetTuple(PtIdx - 1 + axisStep[0] + axisStep[1], &Pts[12]);

  this->EvaluateSecondPartialDeriv(&Dv[15], &degree[0], &U[0], &Pts[0]);
}
void vtkNURBSPatchAdaptor::EvaluateVolumeDeriv(double* Pt, double* Paras, int d)
{
  std::vector<double> Dv(30, 0);
  this->EvaluateSurfaceSecondDeriv(&Dv[0], Paras);
  // Todo: return indicated derivative for order d
  Pt = Pt;
  d = d;
}
void vtkNURBSPatchAdaptor::EvaluateVolumeSecondDeriv(double* Dv, double* Paras)
{
  // currently only implement first order and second order derivatives the Dv will contains 3 * 10 doubles
  // - P_1 Dv_zero namely the point position
  // - P_2 - P_3 First and second order derivatives for parametric space u Dv_1_u Dv_2_u
  // - P_4 - P_5 First and second order derivatives for parametric space v Dv_1_v Dv_2_v
  // - P_6 - P_7 First and second order derivatives for parametric space w Dv_1_w Dv_2_w
  // - P_8 - P_10 Partial derivatives Dv_1_u_v, Dv_1_u_w, Dv_1_v_w
  int PtIdx = 1;

  this->InsertKnotForEval(Paras, PtIdx);

  std::vector<int> degree(3, 0);
  this->GetPatchDegree(&degree[0]);

  std::vector<int> knotStart(3, 0);
  std::vector<int> kSpan(3, 0);
  std::vector<int> sMult(3, 0);
  std::vector<double> U(2, 0.0);
  std::vector<double> Pts(16, 0.0);
  int axisStep[3] = {
    1,
    this->CtrlPtsNumCache[0],
    this->CtrlPtsNumCache[0] * this->CtrlPtsNumCache[1]};
  // we need to retrieve control points along different axis

  // point pos
  this->PointsCache->GetTuple(PtIdx - 1, &Pts[0]);
  Dv[0] = Pts[0] / Pts[3];
  Dv[1] = Pts[1] / Pts[3];
  Dv[2] = Pts[2] / Pts[3];

  // first order and second order derivative
  for (int i = 0; i < 3; ++i)
  {
    this->FindKnotStart(&this->KnotsLenCache[0], i, knotStart[i]);
    this->FindSpanMult(
      Paras[i],
      this->KnotsLenCache[i] - 1,
      &this->KnotsArrCache[knotStart[i]],
      kSpan[i], sMult[i]);

    U[0] = this->KnotsArrCache[knotStart[i] + kSpan[i] + 1];
    U[1] = this->KnotsArrCache[knotStart[i] + kSpan[i] + 2];

    this->PointsCache->GetTuple(PtIdx - 1 + axisStep[i], &Pts[4]);
    this->PointsCache->GetTuple(PtIdx - 1 + 2 * axisStep[i], &Pts[8]);

    this->EvaluateFirstOrderDeriv(&Dv[3 + i * 6], degree[i], U[0], &Pts[0]);
    this->EvaluateSecondOrderDeriv(&Dv[6 + i * 6], degree[i], &U[0], &Pts[0]);
  }

  // partial derivative
  int dim_comb[6] = {0,1,0,2,1,2};
  for (int i = 0; i < 3; ++i)
  {
    int id_0 = dim_comb[2 * i + 0];
    int id_1 = dim_comb[2 * i + 1];

    U[0] = this->KnotsArrCache[knotStart[id_0] + kSpan[id_0] + 1];
    U[1] = this->KnotsArrCache[knotStart[id_1] + kSpan[id_1] + 1];

    this->PointsCache->GetTuple(PtIdx - 1 + axisStep[id_0], &Pts[4]);
    this->PointsCache->GetTuple(PtIdx - 1 + axisStep[id_1], &Pts[8]);
    this->PointsCache->GetTuple(PtIdx - 1 + axisStep[id_0] + axisStep[id_1], &Pts[12]);

    int partial_degree[2] = {degree[id_0], degree[id_1]};
    this->EvaluateSecondPartialDeriv(&Dv[21 + 3 * i], &partial_degree[0], &U[0], &Pts[0]);
  }
}

void vtkNURBSPatchAdaptor::PointInversionCurve(double Pt[3], double* Paras)
{
  // 1. find a start u (in each span)
  // 2. insert u n-degree times
  // 3. compute derivatives at u

  // To find a start u, we compute n equally spaced parameter values on each span

  // compute distance for all candidate parameters to find the start u_0
  std::vector<double> parasRange(2, 0);
  this->GetPatchParameterRange(&parasRange[0]);
  parasRange[0] += 0.01;
  parasRange[1] -= 0.01;
  int n_para = 25;
  double step = (parasRange[1] - parasRange[0]) / n_para;
  double init_para = parasRange[0];
  std::vector<double> Dv(9, 0);
  double min_dist = std::numeric_limits<double>::max();

  for (int i = 0; i <= n_para; ++i)
  {
    double start_para = parasRange[0] + i * step;
    this->EvaluatePoint(&Dv[0], &start_para);

    double cur_dist = vtkMath::Distance2BetweenPoints(Pt, &Dv[0]);
    if (cur_dist < min_dist)
    {
      min_dist = cur_dist;
      init_para = start_para;
    }
  }

  // start from the init_para
  std::vector<double> diff_C_P(3,0);
  double ep_1 = 1e-4;
  double ep_2 = 1e-3;
  for (;;)
  {
    this->EvaluateCurveSecondDeriv(&Dv[0], &init_para);
    diff_C_P[0] = Dv[0] - Pt[0];
    diff_C_P[1] = Dv[1] - Pt[1];
    diff_C_P[2] = Dv[2] - Pt[2];
    double zero_cos =
      abs(vtkMath::Dot(&Dv[3], &diff_C_P[0]))
      / (sqrt(vtkMath::Dot(&Dv[3], &Dv[3])
      * sqrt(vtkMath::Dot(&diff_C_P[0], &diff_C_P[0]))));
    double pt_coincidence =
      sqrt(vtkMath::Dot(&diff_C_P[0], &diff_C_P[0]));
    if (pt_coincidence <= ep_1 && zero_cos <= ep_2)
    {
      break;
    }

    double para_step =
      vtkMath::Dot(&Dv[3], &diff_C_P[0])
      / (vtkMath::Dot(&Dv[6], &diff_C_P[0]) + vtkMath::Dot(&Dv[3], &Dv[3]));

    // assume it a non-closed curve
    if (init_para - para_step < parasRange[0])
    {
      para_step = init_para - parasRange[0];
    }
    if (init_para - para_step > parasRange[1])
    {
      para_step = init_para - parasRange[1];
    }

    init_para -= para_step;

    double para_change = abs(para_step) * sqrt(vtkMath::Dot(&Dv[3], &Dv[3]));
    if (pt_coincidence <= ep_1 || zero_cos <= ep_2 || para_change <= ep_1)
    {
      break;
    }
  }

  Paras[0] = init_para;
  Paras[1] = 0;
  Paras[2] = 0;
}
void vtkNURBSPatchAdaptor::PointInversionSurface(double Pt[3], double* Paras)
{
  std::vector<double> parasRange(4, 0);
  this->GetPatchParameterRange(&parasRange[0]);
  parasRange[0] += 0.01;
  parasRange[1] -= 0.01;
  parasRange[2] += 0.01;
  parasRange[3] -= 0.01;

  // find initial u and v
  int n_para = 25;
  std::vector<double> step(2, 0);
  std::vector<double> init_para(2, 0);
  for (int i = 0; i < 2; ++i)
  {
    step[i] = (parasRange[2 * i + 1] - parasRange[2 * i + 0]) / n_para;
    init_para[i] = parasRange[2 * i];
  }
  std::vector<double> Dv(18, 0);
  std::vector<double> start_para(2, 0);
  double min_dist = std::numeric_limits<double>::max();

  for (int j = 0; j <= n_para; ++j)
  {
    for (int i = 0; i <= n_para; ++i)
    {
      start_para[0] = parasRange[0] + i * step[0];
      start_para[1] = parasRange[2] + j * step[1];
      this->EvaluatePoint(&Dv[0], &start_para[0]);

      double cur_dist = vtkMath::Distance2BetweenPoints(Pt, &Dv[0]);
      if (cur_dist < min_dist)
      {
        min_dist = cur_dist;
        init_para = start_para;
      }
    }
  }

  // start from the init_para
  std::vector<double> diff_C_P(3,0);
  double ep_1 = 1e-4;
  double ep_2 = 1e-3;
  for (;;)
  {
    this->EvaluateSurfaceSecondDeriv(&Dv[0], &init_para[0]);
    diff_C_P[0] = Dv[0] - Pt[0];
    diff_C_P[1] = Dv[1] - Pt[1];
    diff_C_P[2] = Dv[2] - Pt[2];

    // check conditions
    // 1. point coincidence
    double pt_coincidence =
      sqrt(vtkMath::Dot(&diff_C_P[0], &diff_C_P[0]));
    // 2. zero cosine
    std::vector<double> zero_cos(2, 0);
    for (int i = 0; i < 2; ++i)
    {
      zero_cos[i] =
        abs(vtkMath::Dot(&Dv[3 + i * 6], &diff_C_P[0]))
        / (sqrt(vtkMath::Dot(&Dv[3 + i * 6], &Dv[3])
        * sqrt(vtkMath::Dot(&diff_C_P[0], &diff_C_P[0]))));
    }
    if (pt_coincidence <= ep_1 && zero_cos[0] <= ep_2 && zero_cos[1] <= ep_2)
    {
      break;
    }

    // compute new init_para
    double J[2][2];
    J[0][0] = vtkMath::Dot(&Dv[3], &Dv[3]) + vtkMath::Dot(&diff_C_P[0], &Dv[6]);
    J[0][1] = vtkMath::Dot(&Dv[3], &Dv[9]) + vtkMath::Dot(&diff_C_P[0], &Dv[15]);
    J[1][0] = J[0][1];
    J[1][1] = vtkMath::Dot(&Dv[9], &Dv[9]) + vtkMath::Dot(&diff_C_P[0], &Dv[12]);
    std::vector<double> delta(2, 0);
    delta[0] = - vtkMath::Dot(&diff_C_P[0], &Dv[3]);
    delta[1] = - vtkMath::Dot(&diff_C_P[0], &Dv[9]);
    double *J_ptr[2];
    J_ptr[0] = J[0];
    J_ptr[1] = J[1];
    if (!vtkMath::SolveLinearSystem(J_ptr, &delta[0], 2))
    {
      std::cout<<"Error in Solving linear system. Return the latest state.\n";
      Paras[0] = init_para[0];
      Paras[1] = init_para[1];
      Paras[2] = 0;

      return;
    }

    // assume non-closed
    for (int i = 0; i < 2; ++i)
    {
      if (init_para[i] + delta[i] < parasRange[i * 2 + 0])
      {
        delta[i] = parasRange[i * 2 + 0] - init_para[i];
      }
      if (init_para[i] + delta[i] > parasRange[i * 2 + 1])
      {
        delta[i] = parasRange[i * 2 + 1] - init_para[i];
      }
    }
    init_para[0] += delta[0];
    init_para[1] += delta[1];

    double para_change = abs(
      delta[0] * sqrt(vtkMath::Dot(&Dv[3], &Dv[3]))
      + delta[1] * sqrt(vtkMath::Dot(&Dv[9], &Dv[9])));
    if (pt_coincidence <= ep_1
      || (zero_cos[0] <= ep_2 && zero_cos[1] <= ep_2)
      || para_change <= ep_1)
    {
      break;
    }
  }

  Paras[0] = init_para[0];
  Paras[1] = init_para[1];
  Paras[2] = 0;
}
void vtkNURBSPatchAdaptor::PointInversionVolume(double Pt[3], double* Paras)
{
  std::vector<double> parasRange(6, 0);
  this->GetPatchParameterRange(&parasRange[0]);
  parasRange[0] += 0.01;
  parasRange[1] -= 0.01;
  parasRange[2] += 0.01;
  parasRange[3] -= 0.01;
  parasRange[4] += 0.01;
  parasRange[5] -= 0.01;

  // find initial u and v
  int n_para = 25;
  std::vector<double> step(3, 0);
  std::vector<double> start_para(3, 0);
  for (int i = 0; i < 3; ++i)
  {
    step[i] = (parasRange[2 * i + 1] - parasRange[2 * i + 0]) / n_para;
    start_para[i] = parasRange[2 * i];
  }
  std::vector<double> init_para = start_para;
  std::vector<double> Dv(30, 0);
  double min_dist = std::numeric_limits<double>::max();

  for (int k = 0; k <= n_para; ++k)
  {
    for (int j = 0; j <= n_para; ++j)
    {
      for (int i = 0; i <= n_para; ++i)
      {
        start_para[0] = parasRange[0] + i * step[0];
        start_para[1] = parasRange[2] + j * step[1];
        start_para[2] = parasRange[4] + k * step[2];
        this->EvaluatePoint(&Dv[0], &start_para[0]);

        double cur_dist = vtkMath::Distance2BetweenPoints(Pt, &Dv[0]);
        if (cur_dist < min_dist)
        {
          min_dist = cur_dist;
          init_para = start_para;
        }
      }
    }
  }

  // start from the init_para
  std::vector<double> diff_C_P(3,0);
  double ep_1 = 1e-4;
  double ep_2 = 1e-3;
  for (;;)
  {
    this->EvaluateVolumeSecondDeriv(&Dv[0], &init_para[0]);
    diff_C_P[0] = Dv[0] - Pt[0];
    diff_C_P[1] = Dv[1] - Pt[1];
    diff_C_P[2] = Dv[2] - Pt[2];

    // check conditions
    // 1. point coincidence
    double pt_coincidence =
      sqrt(vtkMath::Dot(&diff_C_P[0], &diff_C_P[0]));
    // 2. zero cosine
    std::vector<double> zero_cos(3, 0);
    for (int i = 0; i < 3; ++i)
    {
      zero_cos[i] =
        abs(vtkMath::Dot(&Dv[3 + i * 6], &diff_C_P[0]))
        / (sqrt(vtkMath::Dot(&Dv[3 + i * 6], &Dv[3])
        * sqrt(vtkMath::Dot(&diff_C_P[0], &diff_C_P[0]))));
    }
    if (pt_coincidence <= ep_1
      && zero_cos[0] <= ep_2
      && zero_cos[1] <= ep_2
      && zero_cos[2] <= ep_2)
    {
      break;
    }

    // compute new init_para
    double J[3][3];
    J[0][0] = vtkMath::Dot(&Dv[3], &Dv[3]) + vtkMath::Dot(&diff_C_P[0], &Dv[6]);
    J[0][1] = vtkMath::Dot(&Dv[3], &Dv[9]) + vtkMath::Dot(&diff_C_P[0], &Dv[21]);
    J[0][2] = vtkMath::Dot(&Dv[3], &Dv[15]) + vtkMath::Dot(&diff_C_P[0], &Dv[24]);
    J[1][0] = J[0][1];
    J[1][1] = vtkMath::Dot(&Dv[9], &Dv[9]) + vtkMath::Dot(&diff_C_P[0], &Dv[12]);
    J[1][2] = vtkMath::Dot(&Dv[9], &Dv[15]) + vtkMath::Dot(&diff_C_P[0], &Dv[27]);
    J[2][0] = J[0][2];
    J[2][1] = J[1][2];
    J[2][2] = vtkMath::Dot(&Dv[15], &Dv[15]) + vtkMath::Dot(&diff_C_P[0], &Dv[18]);
    std::vector<double> delta(3, 0);
    delta[0] = - vtkMath::Dot(&diff_C_P[0], &Dv[3]);
    delta[1] = - vtkMath::Dot(&diff_C_P[0], &Dv[9]);
    delta[2] = - vtkMath::Dot(&diff_C_P[0], &Dv[15]);
    double *J_ptr[3];
    J_ptr[0] = J[0];
    J_ptr[1] = J[1];
    J_ptr[2] = J[2];
    if (!vtkMath::SolveLinearSystem(J_ptr, &delta[0], 3))
    {
      std::cout<<"Error in Solving linear system. Return the latest state.\n";
      Paras[0] = init_para[0];
      Paras[1] = init_para[1];
      Paras[2] = init_para[2];

      return;
    }

    // assume non-closed
    for (int i = 0; i < 3; ++i)
    {
      if (init_para[i] + delta[i] < parasRange[i * 2 + 0])
      {
        delta[i] = parasRange[i * 2 + 0] - init_para[i];
      }
      if (init_para[i] + delta[i] > parasRange[i * 2 + 1])
      {
        delta[i] = parasRange[i * 2 + 1] - init_para[i];
      }
    }
    init_para[0] += delta[0];
    init_para[1] += delta[1];
    init_para[2] += delta[2];

    double para_change = abs(
      delta[0] * sqrt(vtkMath::Dot(&Dv[3], &Dv[3]))
      + delta[1] * sqrt(vtkMath::Dot(&Dv[9], &Dv[9]))
      + delta[2] * sqrt(vtkMath::Dot(&Dv[15], &Dv[15])));
    if (pt_coincidence <= ep_1
      || (zero_cos[0] <= ep_2 && zero_cos[1] <= ep_2 && zero_cos[2] <= ep_2)
      || para_change <= ep_1)
    {
      break;
    }
  }

  Paras[0] = init_para[0];
  Paras[1] = init_para[1];
  Paras[2] = init_para[2];
}