/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkFluentWriter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#include "vtkFluentWriter.h"

#include "vtkAppendFilter.h"
#include "vtkCellType.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCompositeDataSet.h"
#include "vtkDataObject.h"
#include "vtkDataObjectTreeIterator.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkIdTypeArray.h"
#include "vtkMath.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtksys/SystemTools.hxx"
#include "vtkThreshold.h"
#include "vtkUnstructuredGrid.h"

#include <sstream>

vtkObjectFactoryNewMacro(vtkFluentWriter);

namespace
{
unsigned int GetNumberOfDigits(unsigned int i)
{
  if (i < 10)
  {
    return 1;
  }
  return GetNumberOfDigits(i/10)+1;
}
}

//----------------------------------------------------------------------------
vtkFluentWriter::vtkFluentWriter()
{
  this->FileName = nullptr;
  this->BinaryFile = 0;

  this->caseFd = nullptr;
  this->dataFd = nullptr;

  this->GhostLevel = 0;
  this->WriteAllTimeSteps = 0;
  this->TransientGeometry = 0;

  this->NumberOfTimeSteps = 0;
  this->CurrentTimeIndex = 0;

  this->EntityOffset = 2;
  this->FluentCellArrays = true;

  this->BlockIdArrayName = "BlockIdScalars";

  this->UnstructuredOutput = vtkSmartPointer<vtkUnstructuredGrid>::New();

  this->LoadVariableNames();
}

//----------------------------------------------------------------------------
vtkFluentWriter::~vtkFluentWriter()
{
  this->SetFileName(nullptr);
}

//----------------------------------------------------------------------------
void vtkFluentWriter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "FileName " << (this->FileName ? this->FileName : "(none)") << endl;
  os << indent << "GhostLevel " << this->GhostLevel << endl;
  os << indent << "WriteAllTimeSteps " << this->WriteAllTimeSteps << endl;
  os << indent << "TransientGeometry" << this->TransientGeometry << endl;
  os << indent << "BlockIdArrayName " << 
    (!this->BlockIdArrayName.empty() ? this->BlockIdArrayName : "(none)") << endl;
}

//----------------------------------------------------------------------------
int vtkFluentWriter::ProcessRequest(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector)
{
  if (request->Has(vtkDemandDrivenPipeline::REQUEST_INFORMATION()))
  {
    return this->RequestInformation(request, inputVector, outputVector);
  }
  else if (request->Has(
      vtkStreamingDemandDrivenPipeline::REQUEST_UPDATE_EXTENT()))
  {
    return this->RequestUpdateExtent(request, inputVector, outputVector);
  }
  // generate the data
  else if (request->Has(vtkDemandDrivenPipeline::REQUEST_DATA()))
  {
    return this->RequestData(request, inputVector, outputVector);
  }

  return this->Superclass::ProcessRequest(request, inputVector, outputVector);
}

//----------------------------------------------------------------------------
int vtkFluentWriter::RequestInformation(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector,
  vtkInformationVector* vtkNotUsed(outputVector))
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  if (inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
  {
    this->NumberOfTimeSteps =
      inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  }
  else
  {
    this->NumberOfTimeSteps = 0;
  }

  return 1;
}

//----------------------------------------------------------------------------
int vtkFluentWriter::RequestUpdateExtent(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector,
  vtkInformationVector* vtkNotUsed(outputVector))
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  if (this->WriteAllTimeSteps &&
      inInfo->Has(vtkStreamingDemandDrivenPipeline::TIME_STEPS()))
  {
    double* timeSteps =
      inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
    double timeReq = timeSteps[this->CurrentTimeIndex];
    inputVector[0]->GetInformationObject(0)->Set
      (vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP(), timeReq);
  }
  return 1;
}

//----------------------------------------------------------------------------
int vtkFluentWriter::FillInputPortInformation(
  int vtkNotUsed(port),
  vtkInformation* info)
{
  info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),
    "vtkCompositeDataSet");
  return 1;
}

//----------------------------------------------------------------------------
int vtkFluentWriter::RequestData(
  vtkInformation* request,
  vtkInformationVector** inputVector,
  vtkInformationVector* vtkNotUsed(outputVector))
{
  if (!this->FileName)
  {
    return 1;
  }

  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  this->OriginalInput = vtkDataObject::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  // is this the first request
  if (this->CurrentTimeIndex == 0 && this->WriteAllTimeSteps)
  {
    // Tell the pipeline to start looping.
    request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 1);
  }

  this->WriteData();

  ++this->CurrentTimeIndex;
  if (this->CurrentTimeIndex >= this->NumberOfTimeSteps)
  {
    this->CloseFluentFile();
    this->CurrentTimeIndex = 0;
    if (this->WriteAllTimeSteps)
    {
      // Tell the pipeline to stop looping.
      request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 0);
    }
  }

  this->CloseFluentFile();

  int localContinue = request->Get(
    vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING());
  if (this->GlobalContinueExecuting(localContinue) != localContinue)
  {
    // Some other node decided to stop the execution.
    assert(localContinue == 1);
    request->Set(vtkStreamingDemandDrivenPipeline::CONTINUE_EXECUTING(), 0);
  }
  return 1;
}

//----------------------------------------------------------------------------
int vtkFluentWriter::GlobalContinueExecuting(int localContinueExecution)
{
  return localContinueExecution;
}

//----------------------------------------------------------------------------
void vtkFluentWriter::WriteData()
{
  this->BaseName = vtksys::SystemTools::GetFilenamePath(this->FileName) + "/" +
    vtksys::SystemTools::GetFilenameWithoutLastExtension(this->FileName);

  int group = 0;
  this->FlattenedInput.clear();
  this->FlattenedNames.clear();

  if (!this->FlattenHierarchy(this->OriginalInput, "", group))
  {
    vtkErrorMacro(
      "vtkFluentWriter::WriteData Unable to flatten hierarchy");
    return;
  }

  this->RemoveGhostCells();

  this->CheckParametes();

  this->MergeBlocks();

  // move check parameters up here and then if there's a change, new file.
  if (this->WriteAllTimeSteps && this->CurrentTimeIndex != 0)
  {
    if (!this->WriteNextTimeStep())
    {
      vtkErrorMacro("vtkFluentWriter::WriteData results");
    }
    return;
  }
  else
  {
    // Close out the old file, if we have one
    if (this->CurrentTimeIndex > 0)
    {
      this->CloseFluentFile();
    }

    this->OpenFluentFile("cas");
    this->WriteGeometry();

    if (this->FluentCellArrays)
    {
      this->OpenFluentFile("dat");
      this->WriteDataFile();
    }
  }
}

//----------------------------------------------------------------------------
int vtkFluentWriter::FlattenHierarchy(vtkDataObject* input, const char *name, int& group)
{
  if (input->IsA("vtkMultiBlockDataSet"))
  {
    vtkMultiBlockDataSet* castObj = vtkMultiBlockDataSet::SafeDownCast(input);
    vtkSmartPointer<vtkDataObjectTreeIterator> iter;
    iter.TakeReference(castObj->NewTreeIterator());
    iter->VisitOnlyLeavesOff();
    iter->TraverseSubTreeOff();
    iter->SkipEmptyNodesOff();
    for (iter->InitTraversal();
         !iter->IsDoneWithTraversal();
         iter->GoToNextItem())
    {
      name = iter->GetCurrentMetaData()->Get(vtkCompositeDataSet::NAME());
      if (name != nullptr && strstr(name, "Sets") != nullptr)
      {
        continue;
      }
      if (name == nullptr)
      {
        // avoid null references in std::string
        name = "";
      }
      if (iter->GetCurrentDataObject() && !this->FlattenHierarchy(iter->GetCurrentDataObject(), name, group))
      {
        return 0;
      }
    }
  }
  else if (input->IsA("vtkCompositeDataSet"))
  {
    vtkCompositeDataSet* castObj = vtkCompositeDataSet::SafeDownCast(input);
    vtkSmartPointer<vtkCompositeDataIterator> iter;
    iter.TakeReference(castObj->NewIterator());
    vtkDataObjectTreeIterator* treeIter = vtkDataObjectTreeIterator::SafeDownCast(castObj);
    if (treeIter)
    {
      treeIter->VisitOnlyLeavesOff();
      treeIter->TraverseSubTreeOff();
      treeIter->SkipEmptyNodesOff();
    }
    for (iter->InitTraversal();
         !iter->IsDoneWithTraversal();
         iter->GoToNextItem())
    {
      if(iter->GetCurrentDataObject() && !this->FlattenHierarchy(iter->GetCurrentDataObject(), name, group))
      {
        return 0;
      }
    }
  }
  else if (input->IsA("vtkDataSet"))
  {
    vtkSmartPointer<vtkUnstructuredGrid> output =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
    if (input->IsA("vtkUnstructuredGrid"))
    {
      output->ShallowCopy(input);
    }
    else
    {
      this->ToUnstructuredGrid(input, output);
    }

    vtkIdType numCells = output->GetNumberOfCells();
    vtkUnsignedCharArray* cArray = vtkUnsignedCharArray::New();
    cArray->SetNumberOfTuples(numCells);

    for (int i = 0; i < numCells; i++)
    {
      cArray->SetValue(i, group);
    }
    cArray->SetName(this->BlockIdArrayName.c_str());
    output->GetCellData()->AddArray(cArray);
    cArray->Delete();

    this->FlattenedInput.push_back(output);
    if (strcmp(name, "") == 0)
    {
      // Setting an arbitrary name for datasets that have not been assigned one.
      char str[512];
      sprintf(str, "block_%d", group);
      name = str;
    }
    this->FlattenedNames.push_back(name);
    ++group;
  }
  else
  {
    vtkErrorMacro(<< "Incorrect class type " << input->GetClassName() << " on input");
    return 0;
  }
  return 1;
}

//----------------------------------------------------------------------------
FILE* vtkFluentWriter::OpenFile(const char* name)
{
  FILE* fd = fopen(name, "wb");

  if (fd == nullptr)
  {
    vtkErrorMacro("Error opening " << name
      << ": " << strerror(errno));
    return nullptr;
  }
  return fd;
}

//----------------------------------------------------------------------------
void vtkFluentWriter::OpenFluentFile(char* extension)
{
  std::string ext = extension;
  std::ostringstream myFileName;
  myFileName << this->BaseName;

  if (this->NumberOfProcesses == 1)
  {
    if (this->WriteAllTimeSteps == false)
    {
      myFileName << "." << ext;
    }
    else
    {
      myFileName << std::setfill('0') << std::setw(6) <<
        this->CurrentTimeIndex << std::setw(0) << "." << ext;
    }
  }
  else
  {
    if (this->WriteAllTimeSteps == false || this->CurrentTimeIndex == 0)
    {
      myFileName << ".";
    }
    else
    {
      myFileName << "." << std::setfill('0') << std::setw(6) <<
        this->CurrentTimeIndex << std::setw(0) << ".";
    }
    unsigned int numDigits = GetNumberOfDigits(
      static_cast<unsigned int>(this->NumberOfProcesses-1));
    myFileName << this->NumberOfProcesses << "." << std::setfill('0')
               << std::setw(numDigits) << this->MyRank << "." << ext;
  }

  if (ext == "cas")
  {
    this->caseFd = this->OpenFile(myFileName.str().c_str());
  }
  else if (ext == "dat")
  {
    this->dataFd = this->OpenFile(myFileName.str().c_str());
  }
  else
  {
    vtkErrorMacro("Unkonwn fluent extension: " << ext);
    return;
  }
}

//----------------------------------------------------------------------------
void vtkFluentWriter::CloseFluentFile()
{
  if (this->caseFd != nullptr)
  {
    fclose(this->caseFd);
    this->caseFd = nullptr;
  }

  if (this->dataFd != nullptr)
  {
    fclose(this->dataFd);
    this->dataFd = nullptr;
  }
}

//----------------------------------------------------------------------------
void vtkFluentWriter::LoadVariableNames()
{
  this->VariableMap["PRESSURE"] = 1;
  this->VariableMap["MOMENTUM"] = 2;
  this->VariableMap["TEMPERATURE"] = 3;
  this->VariableMap["ENTHALPY"] = 4;
  this->VariableMap["X_VELOCITY"] = 111;
  this->VariableMap["Y_VELOCITY"] = 112;
  this->VariableMap["Z_VELOCITY"] = 113;
}

//----------------------------------------------------------------------------
void vtkFluentWriter::RemoveGhostCells()
{
  for (size_t i = 0; i < this->FlattenedInput.size(); ++i)
  {
    vtkUnsignedCharArray *da = this->FlattenedInput[i]->GetCellGhostArray();

    if (da)
    {
      vtkThreshold *t = vtkThreshold::New();
      t->SetInputData(this->FlattenedInput[i]);
      t->ThresholdByLower(0);
      t->SetInputArrayToProcess(
        0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, vtkDataSetAttributes::GhostArrayName());

      t->Update();

      this->FlattenedInput[i] = vtkSmartPointer<vtkUnstructuredGrid>(t->GetOutput());
      t->Delete();

      this->FlattenedInput[i]->GetCellData()->RemoveArray(vtkDataSetAttributes::GhostArrayName());
      this->FlattenedInput[i]->GetPointData()->RemoveArray(vtkDataSetAttributes::GhostArrayName());

      this->GhostLevel = 1;
    }
    else
    {
      this->GhostLevel = 0;
    }
  }
}

//----------------------------------------------------------------------------
int vtkFluentWriter::CheckParametersInternal(int _NumberOfProcesses, int _MyRank)
{
  if (!this->FileName)
  {
    vtkErrorMacro("No filename specified.");
    return 0;
  }

  this->NumberOfProcesses = _NumberOfProcesses;
  this->MyRank = _MyRank;

  this->ConstructBlockInfo();
  this->CheckFluentCellArrays();

  return 1;
}

//----------------------------------------------------------------------------
int vtkFluentWriter::CheckParametes()
{
  return this->CheckParametersInternal(1, 0);
}

//----------------------------------------------------------------------------
void vtkFluentWriter::CheckFluentCellArrays()
{
  std::list<std::string> names0, names;
  for (size_t i = 0; i < this->FlattenedInput.size(); ++i)
  {
    vtkCellData *cellData = this->FlattenedInput[i]->GetCellData();
    int NumCellArrays = cellData->GetNumberOfArrays();

    if (NumCellArrays == 1) // there is a BlockIdScalars array only.
    {
      this->FluentCellArrays = false;
      return;
    }

    names.clear();
    int count = 0;
    for (vtkIdType j = 0; j < NumCellArrays; ++j)
    {
      auto iter = this->VariableMap.find(cellData->GetArray(j)->GetName());
      if (iter != end(this->VariableMap))
      {
        names.push_back(cellData->GetArray(j)->GetName());
        ++count;
      }
    }

    if (count == 0)
    {
      this->FluentCellArrays = false;
      return;
    }

    names.sort();
    if (i == 0)
    {
      names0.assign(names.begin(), names.end());
    }
    else
    {
      if (!(this->FluentCellArrays = names0.size() == names.size() &&
        std::equal(names0.cbegin(), names0.cend(), names.cbegin())))
      {
        this->FluentCellArrays = false;
      }
    }
    if (this->FluentCellArrays == false)
    {
      return;
    }
  }

  this->VariableNames.assign(names0.begin(), names0.end()); 
}

//----------------------------------------------------------------------------
void vtkFluentWriter::ConstructBlockInfo()
{
  this->BoundaryFaceBlockIds.clear();
  this->VolumeBlockIds.clear();

  for (size_t i = 0; i < this->FlattenedInput.size(); ++i)
  {
    vtkIdType numCells = this->FlattenedInput[i]->GetNumberOfCells();
    this->NumBlockElments.push_back(numCells);

    // It is assumed that the dimension of all cells is same through a dataset
    vtkCell* cell = this->FlattenedInput[i]->GetCell(0);
    if (cell->GetCellDimension() == 2)
    {
      this->BoundaryFaceBlockIds.push_back(i);
    }
    else if (cell->GetCellDimension() == 3)
    {
      this->VolumeBlockIds.push_back(i);
    }
  }
}

//----------------------------------------------------------------------------
void vtkFluentWriter::MergeBlocks()
{
  vtkAppendFilter* appender = vtkAppendFilter::New();
  for (size_t i = 0; i < this->FlattenedInput.size(); ++i)
  {
    appender->AddInputData(this->FlattenedInput[i]);
  }
  appender->MergePointsOn();
  appender->Update();

  this->UnstructuredOutput->ShallowCopy(appender->GetOutput());
  //this->UnstructuredOutput->GetFieldData()->PassData(cd->GetFieldData());

  if (this->VolumeBlockIds.size() == 0)
  {
    vtkErrorMacro("Volumemric elments dose not exist.")
    return;
  }

  // for a 3d mesh, all faces on a volumetric element should be attached to the
  // face of another volumetric element or to a surface element. 
  if (this->BoundaryFaceBlockIds.size() == 0)
  {
    int group = this->FlattenedInput.size();

    vtkDataSetSurfaceFilter *surface = vtkDataSetSurfaceFilter::New();
    surface->SetInputData(this->UnstructuredOutput);
    surface->Update();

    vtkSmartPointer<vtkUnstructuredGrid> output =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

    this->ToUnstructuredGrid(surface->GetOutputDataObject(0), output);

    vtkIdType numCells = output->GetNumberOfCells();
    vtkUnsignedCharArray* cArray = vtkUnsignedCharArray::New();
    cArray->SetNumberOfTuples(numCells);

    for (int i = 0; i < numCells; i++)
    {
      cArray->SetValue(i, group);
    }
    cArray->SetName(this->BlockIdArrayName.c_str());
    output->GetCellData()->AddArray(cArray);
    cArray->Delete();

    this->FlattenedInput.push_back(output);
    this->FlattenedNames.push_back("boundary");
    this->BoundaryFaceBlockIds.push_back(group);
    this->NumBlockElments.push_back(numCells);

    appender->AddInputData(output);
    appender->Modified();
    appender->Update();
    this->UnstructuredOutput->ShallowCopy(appender->GetOutput());

    vtkIdType numCells2 = UnstructuredOutput->GetNumberOfCells();

    surface->Delete();
  }
  appender->Delete();

  this->UnstructuredOutput->BuildLinks();
}

//----------------------------------------------------------------------------
void vtkFluentWriter::ToUnstructuredGrid(vtkDataObject *input, vtkUnstructuredGrid *output)
{
  vtkDataSet* castObj = vtkDataSet::SafeDownCast(input);

  output->GetFieldData()->ShallowCopy(castObj->GetFieldData());
  output->GetPointData()->ShallowCopy(castObj->GetPointData());
  output->GetCellData()->ShallowCopy(castObj->GetCellData());

  vtkIdType numPoints = castObj->GetNumberOfPoints();
  vtkSmartPointer<vtkPoints> outPoints = vtkSmartPointer<vtkPoints>::New();
  outPoints->SetNumberOfPoints(numPoints);
  for (vtkIdType i = 0; i < numPoints; i++)
  {
    outPoints->SetPoint(i, castObj->GetPoint(i));
  }
  output->SetPoints(outPoints);

  int numCells = castObj->GetNumberOfCells();
  output->Allocate(numCells);

  vtkIdList* ptIds = vtkIdList::New();
  for (int i = 0; i < numCells; i++)
  {
    castObj->GetCellPoints(i, ptIds);
    output->InsertNextCell(castObj->GetCellType(i), ptIds);
  }
  ptIds->Delete();
}

//----------------------------------------------------------------------------
int vtkFluentWriter::GetCellTypeNumber(int type)
{
  int typeNum = -1;

  switch (type)
  {
    case VTK_TRIANGLE:
      typeNum = 1;
      break;
    case VTK_QUAD:
      typeNum = 3;
      break;
    case VTK_TETRA:
      typeNum = 2;
      break;
    case VTK_HEXAHEDRON:
    case VTK_VOXEL:
      typeNum = 4;
      break;
    case VTK_PYRAMID:
      typeNum = 5;
      break;
    case VTK_WEDGE:
      typeNum = 6;
      break;
    case VTK_POLYGON:
    case VTK_POLYHEDRON:
      typeNum = 7;
      break;
    default:
      typeNum = -1;
      vtkErrorMacro("Unknown cell type is found.");
      break;
  }

  return typeNum;
}

//----------------------------------------------------------------------------
int vtkFluentWriter::WriteGeometry()
{
  int rc = 1;

  char str[512];

  this->WriteStringToFile("(0 \"VTK to Fluent File\")\n", this->caseFd);
  this->WriteStringToFile("(0 \"Dimension:\")\n", this->caseFd);
  this->WriteStringToFile("(2 3)\n\n", this->caseFd);

  int entityId = this->EntityOffset;

  vtkUnsignedCharArray* boundaryDataArray = vtkUnsignedCharArray::SafeDownCast(
    this->UnstructuredOutput->GetCellData()->GetArray(this->BlockIdArrayName.c_str()));

  if (!boundaryDataArray)
  {
    vtkErrorMacro("BlockIdScalars array does not exist.");
    return 0;
  }

  this->NumOfPoints = this->UnstructuredOutput->GetNumberOfPoints();
  vtkIdType numOfCells = this->UnstructuredOutput->GetNumberOfCells();
  vtkIdType numOfFaces = 0;
  vtkIdType numOfBoundaryFaces = 0;
  this->NumOfVolumes = 0;
  vtkIdList* volumeCellIdArray = vtkIdList::New();
  for (vtkIdType cId = 0; cId < numOfCells; ++cId)
  {
    vtkCell* cell = this->UnstructuredOutput->GetCell(cId);
    if (cell->GetCellDimension() == 2)
    {
      ++numOfBoundaryFaces;
      continue;
    }

    if (cell->GetCellDimension() == 3)
    {
      volumeCellIdArray->InsertNextId(cId);
      numOfFaces += cell->GetNumberOfFaces();
      this->NumOfVolumes += 1;
    }
  }

  vtkIdType numOfInteriorFaces = numOfFaces - numOfBoundaryFaces;
  numOfInteriorFaces /= 2;
  this->NumOfFaces = numOfInteriorFaces + numOfBoundaryFaces;

  sprintf(str, "(10 (0 1 %x 0))\n\n", this->NumOfPoints);
  this->WriteStringToFile(str, this->caseFd);

  sprintf(str, "(12 (0 1 %x 0))\n", this->NumOfVolumes);
  this->WriteStringToFile(str, this->caseFd);

  int volOffset = 1;
  for (size_t i = 0; i < this->VolumeBlockIds.size(); ++i)
  {
    int blockId = this->VolumeBlockIds[i];

    if (this->BinaryFile)
    {
      sprintf(str, "(3012 (%x %x %x 1 0)(", entityId, volOffset, volOffset + NumBlockElments[blockId] - 1);
    }
    else
    {
      sprintf(str, "(12 (%x %x %x 1 0)(\n", entityId, volOffset, volOffset+NumBlockElments[blockId]-1);
    }
    this->WriteStringToFile(str, this->caseFd);

    if (this->BinaryFile)
    {
      int *cellSet = new int[NumBlockElments[blockId]];
      for (int i = 0, id = 0; i < this->NumOfVolumes; ++i)
      {
        vtkIdType cellId = volumeCellIdArray->GetId(i);

        if (boundaryDataArray->GetTuple1(cellId) != blockId)
        {
          continue;
        }
        cellSet[id++] = this->GetCellTypeNumber(this->UnstructuredOutput->GetCell(cellId)->GetCellType());
      }
      this->WriteIntToFile(cellSet, NumBlockElments[blockId], this->caseFd);
      delete[] cellSet;
    }
    else
    {
      int colCounts = 0;
      for (int i = 0; i < this->NumOfVolumes; ++i)
      {
        vtkIdType cellId = volumeCellIdArray->GetId(i);

        if (boundaryDataArray->GetTuple1(cellId) != blockId)
        {
          continue;
        }

        int typeNum = this->GetCellTypeNumber(this->UnstructuredOutput->GetCell(cellId)->GetCellType());
        sprintf(str, " %d", typeNum);
        this->WriteStringToFile(str, this->caseFd);

        if (++colCounts % 10 == 0)
        {
          this->WriteStringToFile("\n", this->caseFd);
        }
      }
    }
    if (this->BinaryFile)
    {
      this->WriteStringToFile(")\nEnd of Binary Section   3012)\n\n", this->caseFd);
    }
    else
    {
      this->WriteStringToFile("))\n\n", this->caseFd);
    }
    volOffset += this->NumBlockElments[blockId];
    ++entityId;
  }

  this->WriteStringToFile("(0 \"Node Section\")\n", this->caseFd);
  if (this->BinaryFile)
  {
    sprintf(str, "(3010 (1 1 %x 1 3) (", this->NumOfPoints);
  }
  else
  {
    sprintf(str, "(10 (1 1 %x 1 3) (\n", this->NumOfPoints);
  }
  this->WriteStringToFile(str, this->caseFd);

  vtkPoints *pts = this->UnstructuredOutput->GetPoints();
  if (pts)
  {
    vtkDataArray *da = pts->GetData();
    if (this->BinaryFile)
    {
      double *coords = new double[this->NumOfPoints*3];
      for (int j = 0; j < this->NumOfPoints; ++j)
      {
        int id = 3*j;
        coords[id+0] = da->GetComponent(j, 0);
        coords[id+1] = da->GetComponent(j, 1);
        coords[id+2] = da->GetComponent(j, 2);
      }
      this->WriteDoubleToFile(coords, this->NumOfPoints *3, this->caseFd);
      this->WriteStringToFile(")\nEnd of Binary Section   3010)\n\n", this->caseFd);
      delete[] coords;
    }
    else
    {
      for (int j = 0; j < this->NumOfPoints; ++j)
      {
        sprintf(str, "  %17.10e  %17.10e  %17.10e\n",
            da->GetComponent(j, 0),
            da->GetComponent(j, 1),
            da->GetComponent(j, 2));
        this->WriteStringToFile(str, this->caseFd);
      }
      this->WriteStringToFile("))\n\n", this->caseFd);
    }
  }

  vtkIdList* volumeCellIdMap = vtkIdList::New();
  for (int i = 0; i < numOfCells; ++i)
  {
    volumeCellIdMap->InsertNextId(-1);
  }

  for (int i = 0; i < this->NumOfVolumes; ++i)
  {
    volumeCellIdMap->SetId(volumeCellIdArray->GetId(i), i);
  }  

  this->WriteStringToFile("(0 \"Faces Section\")\n", this->caseFd);
  sprintf(str, "(13 (0 1 %x 0))\n\n", this->NumOfFaces);
  this->WriteStringToFile(str, this->caseFd);

  int faceOffset = 1;
  vtkIdList* neighborCellIds = vtkIdList::New();
  for (size_t i = 0; i < this->BoundaryFaceBlockIds.size(); ++i)
  {
    int blockId = this->BoundaryFaceBlockIds[i];

    if (this->BinaryFile)
    {
      sprintf(str, "(3013 (%x %x %x 3 0)(", entityId, faceOffset, faceOffset + NumBlockElments[blockId] - 1);
      this->WriteStringToFile(str, this->caseFd);
    }
    else
    {
      sprintf(str, "(13 (%x %x %x 3 0)(\n", entityId, faceOffset, faceOffset+NumBlockElments[blockId]-1);
      this->WriteStringToFile(str, this->caseFd);
    }

    for (vtkIdType cellId = 0; cellId < numOfCells; ++cellId)
    {
      if (boundaryDataArray->GetTuple1(cellId) != blockId)
      {
        continue;
      }

      vtkIdList* facePointIds = this->UnstructuredOutput->GetCell(cellId)->GetPointIds();
      this->UnstructuredOutput->GetCellNeighbors(cellId, facePointIds, neighborCellIds);
      vtkIdType neighborCellId = neighborCellIds->GetId(0);
      vtkIdType numIds = facePointIds->GetNumberOfIds();

      if (this->BinaryFile)
      {
        int *faceSet = new int[3+numIds];
        faceSet[0] = numIds;
        for (vtkIdType pntId = 0; pntId < numIds; ++pntId)
        {
          faceSet[pntId+1] = (int)facePointIds->GetId(pntId) + 1;
        }
        faceSet[numIds+1] = 0;
        faceSet[numIds+2] = (int)volumeCellIdMap->GetId(neighborCellId) + 1;
        this->WriteIntToFile(faceSet, 3+numIds, this->caseFd);
        delete[] faceSet;
      }
      else
      {
        sprintf(str, " %d", numIds);
        this->WriteStringToFile(str, this->caseFd);
        for (vtkIdType pntId = 0; pntId < numIds; ++pntId)
        {
          sprintf(str, " %x", (int)facePointIds->GetId(pntId)+1);
          this->WriteStringToFile(str, this->caseFd);
        }
        sprintf(str, " 0 %x\n", (int)volumeCellIdMap->GetId(neighborCellId)+1);
        this->WriteStringToFile(str, this->caseFd);
      }
    }
    if (this->BinaryFile)
    {
      this->WriteStringToFile(")\nEnd of Binary Section   3013)\n\n", this->caseFd);
    }
    else
    {
      this->WriteStringToFile("))\n\n", this->caseFd);
    }

    faceOffset += this->NumBlockElments[blockId];
    ++entityId;
  }

  if (this->BinaryFile)
  {
    sprintf(str, "(3013 (%x %x %x 2 0)(", (int)entityId, faceOffset, faceOffset + numOfInteriorFaces - 1);
  }
  else
  {
    sprintf(str, "(13 (%x %x %x 2 0)(\n", (int)entityId, faceOffset, faceOffset+numOfInteriorFaces-1);
  }
  this->WriteStringToFile(str, this->caseFd);

  for (size_t i = 0; i < this->VolumeBlockIds.size(); ++i)
  {
    int blockId = this->VolumeBlockIds[i];

    for (int i = 0; i < this->NumOfVolumes; ++i)
    {
      vtkIdType cellId = volumeCellIdArray->GetId(i);

      if (boundaryDataArray->GetTuple1(cellId) != blockId)
      {
        continue;
      }

      vtkCell *cell = this->UnstructuredOutput->GetCell(cellId);
      vtkIdType numOfFaces = cell->GetNumberOfFaces();
      for (vtkIdType fId = 0; fId < numOfFaces; ++fId)
      {
        vtkCell* face = cell->GetFace(fId);
        vtkIdList* facePointIds = face->GetPointIds();
        int numIds = facePointIds->GetNumberOfIds();

        if (cell->GetCellType() == VTK_VOXEL)
        {
          vtkIdType id2 = facePointIds->GetId(2);
          vtkIdType id3 = facePointIds->GetId(3);
          facePointIds->SetId(2, id3);
          facePointIds->SetId(3, id2);
        }

        this->UnstructuredOutput->GetCellNeighbors(cellId, facePointIds, neighborCellIds);
        int numNeighbors = neighborCellIds->GetNumberOfIds();
        if (numNeighbors != 1 ||
            volumeCellIdMap->GetId(neighborCellIds->GetId(0)) == -1)
        {
          continue;
        }
        if (neighborCellIds->GetId(0) < cellId)
        {
          continue;
        }

        if (this->BinaryFile)
        {
          int *faceSet = new int[3 + numIds];
          faceSet[0] = numIds;
          for (vtkIdType pntId = 0; pntId < numIds; ++pntId)
          {
            faceSet[pntId + 1] = (int)facePointIds->GetId(pntId) + 1;
          }
          faceSet[numIds + 1] = (int)volumeCellIdMap->GetId(neighborCellIds->GetId(0)) + 1;
          faceSet[numIds + 2] = (int)volumeCellIdMap->GetId(cellId) + 1;
          this->WriteIntToFile(faceSet, 3 + numIds, this->caseFd);
          delete[] faceSet;
        }
        else
        {
          sprintf(str, " %d", numIds);
          this->WriteStringToFile(str, this->caseFd);
          for (int pntId = 0; pntId < numIds; ++pntId)
          {
            sprintf(str, " %x", (int)facePointIds->GetId(pntId) + 1);
            this->WriteStringToFile(str, this->caseFd);
          }
          sprintf(str, " %x %x\n", (int)volumeCellIdMap->GetId(neighborCellIds->GetId(0)) + 1, (int)volumeCellIdMap->GetId(cellId) + 1);
          this->WriteStringToFile(str, this->caseFd);
        }
      }
    }
    if (this->BinaryFile)
    {
      this->WriteStringToFile(")\nEnd of Binary Section   3013)\n\n", this->caseFd);
    }
    else
    {
      this->WriteStringToFile("))\n\n", this->caseFd);
    }
    faceOffset += numOfInteriorFaces;
    ++entityId;
  }
  neighborCellIds->Delete();

  this->WriteStringToFile("(0 \"Zone Sections\")\n", this->caseFd);

  entityId = this->EntityOffset;

  for (size_t i = 0; i < this->VolumeBlockIds.size(); ++i)
  {
    int blockId = this->VolumeBlockIds[i];
    std::string name = this->FlattenedNames[blockId];    
    sprintf(str, "(45 (%x fluid %s)())\n", entityId, name.c_str());
    this->WriteStringToFile(str, this->caseFd);
    ++entityId;
  }

  for (size_t i = 0; i < this->BoundaryFaceBlockIds.size(); ++i)
  {
    int blockId = this->BoundaryFaceBlockIds[i];
    std::string name = this->FlattenedNames[blockId]; 
    sprintf(str, "(45 (%x wall %s)())\n", entityId, name.c_str());
    this->WriteStringToFile(str, this->caseFd);
    ++entityId;
  }

  for (size_t i = 0; i < this->VolumeBlockIds.size(); ++i)
  {
    sprintf(str, "(45 (%x interior %s)())\n", entityId, "default-interior");
    this->WriteStringToFile(str, this->caseFd);
    ++entityId;
  }

  volumeCellIdArray->Delete();
  volumeCellIdMap->Delete();

  return rc;
}

//----------------------------------------------------------------------------
int vtkFluentWriter::WriteDataFile()
{
  int rc = 1;
  char str[512];

  vtkUnsignedCharArray* boundaryDataArray = vtkUnsignedCharArray::SafeDownCast(
    this->UnstructuredOutput->GetCellData()->GetArray(this->BlockIdArrayName.c_str()));

  if (!boundaryDataArray)
  {
    vtkErrorMacro("BlockIdScalars array does not exist.");
    return 0;
  }

  vtkIdType numOfCells = this->UnstructuredOutput->GetNumberOfCells();
  vtkIdList* volumeCellIdArray = vtkIdList::New();
  for (vtkIdType cId = 0; cId < numOfCells; ++cId)
  {
    vtkCell* cell = this->UnstructuredOutput->GetCell(cId);
    if (cell->GetCellDimension() == 3)
    {
      volumeCellIdArray->InsertNextId(cId);
    }
  }

  sprintf(str, "(33 (%d %d %d))\n\n", this->NumOfVolumes, this->NumOfFaces, this->NumOfPoints);
  this->WriteStringToFile(str, this->dataFd);

  for (std::string name : this->VariableNames)
  {
    vtkDataArray *da = this->UnstructuredOutput->GetCellData()->GetArray(name.c_str());
    if (!da)
    {
      vtkErrorMacro(<< "cell data: " << name << " is invalid.");
      return 0;
    }

    int volOffset = 1;
    int entityId = this->EntityOffset;

    for (size_t i = 0; i < this->VolumeBlockIds.size(); ++i)
    {
      int blockId = this->VolumeBlockIds[i];

      if (this->BinaryFile)
      {
        sprintf(str, "(3300 (%d %d 1 0 0 %d %d)(", this->VariableMap[name], entityId,
          volOffset, volOffset + NumBlockElments[blockId] - 1);
      }
      else
      {
        sprintf(str, "(300 (%d %d 1 0 0 %d %d)(\n", this->VariableMap[name], entityId,
          volOffset, volOffset + NumBlockElments[blockId] - 1);
      }
      this->WriteStringToFile(str, this->dataFd);

      if (this->BinaryFile)
      {
        double *cellSet = new double[NumBlockElments[blockId]];
        for (int i = 0, id = 0; i < this->NumOfVolumes; ++i)
        {
          vtkIdType cellId = volumeCellIdArray->GetId(i);

          if (boundaryDataArray->GetTuple1(cellId) != blockId)
          {
            continue;
          }
          cellSet[id++] = da->GetTuple1(cellId);
        }
        this->WriteDoubleToFile(cellSet, NumBlockElments[blockId], this->dataFd);
        delete[] cellSet;
      }
      else
      {
        int colCounts = 0;
        for (int i = 0; i < this->NumOfVolumes; ++i)
        {
          vtkIdType cellId = volumeCellIdArray->GetId(i);

          if (boundaryDataArray->GetTuple1(cellId) != blockId)
          {
            continue;
          }

          sprintf(str, "%17.10e\n", da->GetTuple1(cellId));
          this->WriteStringToFile(str, this->dataFd);
        }
      }
      if (this->BinaryFile)
      {
        this->WriteStringToFile(")\nEnd of Binary Section   3300)\n\n", this->dataFd);
      }
      else
      {
        this->WriteStringToFile("))\n\n", this->dataFd);
      }
      volOffset += this->NumBlockElments[blockId];
      ++entityId;
    }
  }
  
  volumeCellIdArray->Delete();

  return rc;
}

//----------------------------------------------------------------------------
int vtkFluentWriter::WriteNextTimeStep()
{
  int rc = 1;

  float tsv = 0.0;
  if (this->GetInput()->GetInformation()->Has(vtkDataObject::DATA_TIME_STEP()))
  {
    tsv = this->GetInput()->GetInformation()->Has(vtkDataObject::DATA_TIME_STEP());
  }

  // Close out the old file, if we have one
  if (this->CurrentTimeIndex > 0)
  {
    this->CloseFluentFile();
  }

  this->CheckParametes();

  if (this->ShouldWriteGeometry())
  {
    this->OpenFluentFile("cas");
    this->WriteGeometry();
  }

  if (this->FluentCellArrays)
  {
    this->OpenFluentFile("dat");
    this->WriteDataFile();
  }

  return rc;
}

//----------------------------------------------------------------------------
void vtkFluentWriter::WriteStringToFile(const char* cstring, FILE* file)
{
  fwrite(cstring, sizeof(char), std::min(strlen(cstring), static_cast<size_t>(512)), file);
}

//----------------------------------------------------------------------------
void vtkFluentWriter::WriteDoubleToFile(const double* d, const int num, FILE* file)
{
  fwrite(d, sizeof(double), num, file);
}

//----------------------------------------------------------------------------
void vtkFluentWriter::WriteIntToFile(const int* i, const int num, FILE* file)
{
  fwrite(i, sizeof(int), num, file);
}

//----------------------------------------------------------------------------
bool vtkFluentWriter::ShouldWriteGeometry()
{
  return (this->TransientGeometry || (this->CurrentTimeIndex == 0));
}
