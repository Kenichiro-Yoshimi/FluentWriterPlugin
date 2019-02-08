/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkFluentWriter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
 * @class   vtkFluentWriter
 * @brief   Write Fluent files for 3D.
 *
 *     This is a vtkWriter that writes it's vtkUnstructuredGrid
 *     input out to a Fluent file as a 3D mesh. Since Fluent format treats with
 *     volumetric elements in 3D, polygons which are not attached to the face
 *     of a volumetric element are not supported.
*/

#ifndef vtkFluentWriter_h
#define vtkFluentWriter_h

#include "vtkSmartPointer.h"
#include "vtkWriter.h"

#include <list>
#include <map>
#include <vector>

class vtkUnstructuredGrid;

class vtkFluentWriter : public vtkWriter
{
public:
  static vtkFluentWriter *New();
  vtkTypeMacro(vtkFluentWriter,vtkWriter);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  //@{
  /**
   * Name for the output file.
   */
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
  //@}

  //@{
  /**
   * Is the file to be written in binary format (as opposed to ascii).
   */
  vtkSetMacro(BinaryFile, int);
  vtkGetMacro(BinaryFile, int);
  vtkBooleanMacro(BinaryFile, int);
  //@}

  //@{
  /**
   * We never write out ghost cells.  This variable is here to satisfy
   * the behavior of ParaView on invoking a parallel writer.
   */
  vtkSetMacro(GhostLevel, int);
  vtkGetMacro(GhostLevel, int);
  //@}

  //@{
  /**
   * When WriteAllTimeSteps is turned ON, the writer is executed once for
   * each timestep available from the reader.
   */
  vtkSetMacro(WriteAllTimeSteps, int);
  vtkGetMacro(WriteAllTimeSteps, int);
  vtkBooleanMacro(WriteAllTimeSteps, int);
  //@}

  //@{
  /**
   * Specify whetehr the geometry changes each timestep
   * if false, geometry is only written at timestep 0.
   */
  vtkSetMacro(TransientGeometry, int);
  vtkGetMacro(TransientGeometry, int);
  vtkBooleanMacro(TransientGeometry, int);
  //@}

  //
  //  Structures
  //
  struct Face;

protected:
  vtkFluentWriter();
  ~vtkFluentWriter();

  char *FileName;
  int BinaryFile;
  std::string BaseName;
  FILE* caseFd;
  FILE* dataFd;

  int NumberOfProcesses;
  int MyRank;

  int GhostLevel;
  int WriteAllTimeSteps;
  int TransientGeometry;

  int NumberOfTimeSteps;
  int CurrentTimeIndex;

  bool FluentCellArrays;

  int NumOfPoints;
  int NumOfVolumes;
  int NumOfFaces;

  int EntityOffset;

  vtkSmartPointer<vtkUnstructuredGrid> UnstructuredOutput;
  vtkDataObject *OriginalInput;
  std::vector<vtkSmartPointer<vtkUnstructuredGrid>> FlattenedInput;

  std::vector<std::string> FlattenedNames;
  std::vector<int> BoundaryFaceBlockIds;
  std::vector<int> VolumeBlockIds;
  std::vector<int> NumBlockElments;
  std::map<std::string, int> VariableMap;
  std::list<std::string> VariableNames;
  std::string BlockIdArrayName;
  std::map<int, std::vector<Face>> Faces;

  int ProcessRequest(vtkInformation* request,
                     vtkInformationVector** inputVector,
                     vtkInformationVector* outputVector) override;

  int RequestInformation(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVecotr);

  virtual int RequestUpdateExtent(vtkInformation* request,
                                  vtkInformationVector** inputVector,
                                  vtkInformationVector* outputVector);

  int FillInputPortInformation(int port, vtkInformation* info) override;

  int RequestData(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector) override;

  void WriteData() override;

  int FlattenHierarchy(vtkDataObject* input, const char *name, int& group);
  FILE* OpenFile(const char* name);

  void OpenFluentFile(char* ext);
  void CloseFluentFile();

  virtual void LoadVariableNames();
  int GetCellTypeNumber(int type);
  void RemoveGhostCells();
  int CheckParametes();
  int CheckParametersInternal(int NumberOfProcesses, int MyRank);
  void CheckFluentCellArrays();

  void ConstructBlockInfo();
  void MergeBlocks();
  void ToUnstructuredGrid(vtkDataObject *input, vtkUnstructuredGrid *output);

  /**
   * If writing in parallel multiple time steps exchange after each time step
   * if we should continue the execution. Pass local continueExecution as a
   * parameter and return the global continueExecution.
   */
  virtual int GlobalContinueExecuting(int localContinueExecution);

  virtual int WriteGeometry();
  virtual int WriteDataFile();
  virtual int WriteNextTimeStep();

  virtual void WriteStringToFile(const char* string, FILE* file);
  virtual void WriteDoubleToFile(const double* d, const int num, FILE* file);
  virtual void WriteIntToFile(const int* i, const int num, FILE* file);

  virtual bool ShouldWriteGeometry();

private:
  vtkFluentWriter(const vtkFluentWriter&) = delete;
  void operator=(const vtkFluentWriter&) = delete;
};

#endif
