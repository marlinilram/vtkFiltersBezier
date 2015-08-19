set(Module_SRCS
  vtkPatchInterpolation.cxx
  vtkBezierPatchAdaptor.cxx
  vtkNURBSPatchAdaptor.cxx
  vtkBREPReader.cxx
  )

set_source_files_properties(
  vtkBezierPatchAdaptor
  ABSTRACT
  )
  
vtk_module_library(vtkFiltersBezier ${Module_SRCS})
