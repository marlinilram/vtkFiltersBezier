vtk_add_test_cxx(${vtk-module}CxxTests tests
  TestBREPReader.cxx
  TestCurveProjection.cxx
  TestSurfaceProjection.cxx
  TestSimplicialInterpolation.cxx
  )
vtk_test_cxx_executable(${vtk-module}CxxTests tests)
