vtk_add_test_cxx(${vtk-module}CxxTests tests
  TestBREPReader.cxx
  TestNURBSAdaptor.cxx
  TestPatchInterpolation.cxx
  )
vtk_test_cxx_executable(${vtk-module}CxxTests tests)
