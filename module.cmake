vtk_module(vtkFiltersBezier
  GROUPS
    StandAlone
  DEPENDS
    vtkFiltersCore
  TEST_DEPENDS
    vtkIOXML
    vtkRendering${VTK_RENDERING_BACKEND}
    vtkTestingRendering
    vtkInteractionStyle
  KIT
    vtkFilters
  )
