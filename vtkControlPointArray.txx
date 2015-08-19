#include "vtkControlPointArray.h"

#include "vtkDataArray.h"
#include "vtkIdList.h"
#include "vtkObjectFactory.h"
#include "vtkVariant.h"
#include "vtkVariantCast.h"

template <class Scalar> vtkControlPointArray<Scalar>*
  vtkControlPointArray<Scalar>::New()
{
  VTK_STANDARD_NEW_BODY(vtkControlPointArray<Scalar>);
}

template <class Scalar> void vtkControlPointArray<Scalar>
::PrintSelf(ostream& os, vtkIndent indent)
{
  this->vtkControlPointArray<Scalar>::Superclass::PrintSelf(
    os, indent);
}

template <class Scalar> void vtkControlPointArray<Scalar>
::InitializeArray(vtkDataArray* controlPoints)
{
  this->MaxId = controlPoints->GetNumberOfTuples() * 3 - 1;
  this->Size = this->MaxId + 1;
  this->NumberOfComponents = 3;
  this->ControlPoints = controlPoints;
}

template <class Scalar> void vtkControlPointArray<Scalar>
::Initialize()
{
  this->MaxId = -1;
  this->Size = 0;
  this->NumberOfComponents = 3;
}

template <class Scalar> void vtkControlPointArray<Scalar>
::GetTuples(vtkIdList *ptIds, vtkAbstractArray *output)
{
  vtkDataArray* da = vtkDataArray::FastDownCast(output);
  if (!da)
    {
    vtkWarningMacro(<<"Input is not a vtkDataArray");
    return;
    }

  if (da->GetNumberOfComponents() != this->GetNumberOfComponents())
    {
     vtkWarningMacro(<<"Incorrect number of components in input array.");
     return;
    }

  const vtkIdType numPoints = ptIds->GetNumberOfIds();
  for (vtkIdType i = 0; i < numPoints; ++i)
    {
    da->SetTuple(i, this->GetTuple(ptIds->GetId(i)));
    }
}

template <class Scalar> void vtkControlPointArray<Scalar>
::GetTuples(vtkIdType p1, vtkIdType p2, vtkAbstractArray *output)
{
  vtkDataArray* da = vtkDataArray::FastDownCast(output);
  if (!da)
    {
    vtkWarningMacro(<<"Input is not a vtkDataArray");
    return;
    }

  if (da->GetNumberOfComponents() != this->GetNumberOfComponents())
    {
    vtkWarningMacro(<<"Incorrect number of components in input array.");
    return;
    }

  for (vtkIdType daTupleId = 0; p1 <= p2; ++p1)
    {
    da->SetTuple(daTupleId++, this->GetTuple(p1));
    }
}

template <class Scalar> void vtkControlPointArray<Scalar>
::Squeeze()
{

}

template <class Scalar> vtkArrayIterator*
vtkControlPointArray<Scalar>::NewIterator()
{
  vtkErrorMacro(<<"Not implemented.");
  return NULL;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType vtkControlPointArray<Scalar>
::LookupValue(vtkVariant value)
{
  bool valid = true;
  Scalar val = vtkVariantCast<Scalar>(value, &valid);
  if (valid)
    {
    // return this->Lookup(val, 0);
    }

  // to suppress unreference warnings;
  val = val;
  return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::LookupValue(vtkVariant value, vtkIdList *ids)
{
  /*
  bool valid = true;
  Scalar val = vtkVariantCast<Scalar>(value, &valid);
  ids->Reset();
  if (valid)
    {
    vtkIdType index = 0;
    while ((index = this->Lookup(val, index)) >= 0)
      {
      ids->InsertNextId(index);
      ++index;
      }
    }
  */

  // to suppress unreference warnings;
  ids = ids;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkVariant vtkControlPointArray<Scalar>
::GetVariantValue(vtkIdType idx)
{
  return vtkVariant(this->GetValueReference(idx));
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::ClearLookup()
{
  // no-op, no fast lookup implemented.
}

//------------------------------------------------------------------------------
template <class Scalar> double* vtkControlPointArray<Scalar>
::GetTuple(vtkIdType i)
{
  this->GetTuple(i, this->TempDoubleArray);
  return this->TempDoubleArray;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::GetTuple(vtkIdType i, double *tuple)
{
  double* tupleH = this->ControlPoints->GetTuple(i);
  tuple[0] = tupleH[0] / tupleH[3];
  tuple[1] = tupleH[1] / tupleH[3];
  tuple[2] = tupleH[2] / tupleH[3];
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType vtkControlPointArray<Scalar>
::LookupTypedValue(Scalar value)
{
//  return this->Lookup(value, 0);
  // to suppress unreference warnings;
  value = value;
  return 0;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::LookupTypedValue(Scalar value, vtkIdList *ids)
{
  /*
  ids->Reset();
  vtkIdType index = 0;
  while ((index = this->Lookup(value, index)) >= 0)
    {
    ids->InsertNextId(index);
    ++index;
    }
  */

  // to suppress unreference warnings;
  ids = ids;
  value = value;
}

//------------------------------------------------------------------------------
template <class Scalar> Scalar vtkControlPointArray<Scalar>
::GetValue(vtkIdType idx)
{
  return this->GetValueReference(idx);
}

//------------------------------------------------------------------------------
template <class Scalar> Scalar& vtkControlPointArray<Scalar>
::GetValueReference(vtkIdType idx)
{
  //assert(this->ControlPoints);

  const vtkIdType tuple = idx / 3;
  const vtkIdType comp = idx % 3;

  // return reference here doesn't make real difference
  double* tupleH = this->ControlPoints->GetTuple(tuple);
  this->TempDoubleArray[comp] = tupleH[comp] / tupleH[3];

  return this->TempDoubleArray[comp];
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::GetTupleValue(vtkIdType tupleId, Scalar *tuple)
{
  double* tupleH = this->ControlPoints->GetTuple(tupleId);
  tuple[0] = tupleH[0] / tupleH[3];
  tuple[1] = tupleH[1] / tupleH[3];
  tuple[2] = tupleH[2] / tupleH[3];
}

//------------------------------------------------------------------------------
template <class Scalar> int vtkControlPointArray<Scalar>
::Allocate(vtkIdType, vtkIdType)
{
  vtkErrorMacro("Read only container.")
    return 0;
}

//------------------------------------------------------------------------------
template <class Scalar> int vtkControlPointArray<Scalar>
::Resize(vtkIdType)
{
  vtkErrorMacro("Read only container.")
    return 0;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::SetNumberOfTuples(vtkIdType)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::SetTuple(vtkIdType, vtkIdType, vtkAbstractArray *)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::SetTuple(vtkIdType, const float *)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::SetTuple(vtkIdType, const double *)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::InsertTuple(vtkIdType, vtkIdType, vtkAbstractArray *)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::InsertTuple(vtkIdType, const float *)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::InsertTuple(vtkIdType, const double *)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::InsertTuples(vtkIdList *, vtkIdList *, vtkAbstractArray *)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::InsertTuples(vtkIdType, vtkIdType, vtkIdType, vtkAbstractArray *)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType vtkControlPointArray<Scalar>
::InsertNextTuple(vtkIdType, vtkAbstractArray *)
{
  vtkErrorMacro("Read only container.")
    return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType vtkControlPointArray<Scalar>
::InsertNextTuple(const float *)
{

  vtkErrorMacro("Read only container.")
    return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType vtkControlPointArray<Scalar>
::InsertNextTuple(const double *)
{
  vtkErrorMacro("Read only container.")
    return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::DeepCopy(vtkAbstractArray *)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::DeepCopy(vtkDataArray *)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::InterpolateTuple(vtkIdType, vtkIdList *, vtkAbstractArray *, double *)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::InterpolateTuple(vtkIdType, vtkIdType, vtkAbstractArray*, vtkIdType,
                   vtkAbstractArray*, double)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::SetVariantValue(vtkIdType, vtkVariant)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::RemoveTuple(vtkIdType)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::RemoveFirstTuple()
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::RemoveLastTuple()
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::SetTupleValue(vtkIdType, const Scalar*)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::InsertTupleValue(vtkIdType, const Scalar*)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType vtkControlPointArray<Scalar>
::InsertNextTupleValue(const Scalar *)
{
  vtkErrorMacro("Read only container.")
    return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::SetValue(vtkIdType, Scalar)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkIdType vtkControlPointArray<Scalar>
::InsertNextValue(Scalar)
{
  vtkErrorMacro("Read only container.")
    return -1;
}

//------------------------------------------------------------------------------
template <class Scalar> void vtkControlPointArray<Scalar>
::InsertValue(vtkIdType, Scalar)
{
  vtkErrorMacro("Read only container.")
    return;
}

//------------------------------------------------------------------------------
template <class Scalar> vtkControlPointArray<Scalar>
::vtkControlPointArray()
  : ControlPoints(0)
{
}

//------------------------------------------------------------------------------
template <class Scalar> vtkControlPointArray<Scalar>
::~vtkControlPointArray()
{
}