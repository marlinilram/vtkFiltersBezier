# vtkFiltersBezier

## These codes are for VTK/GSoC2015 "spline curve- and surface-based visualizations":


Abstract
This project aims at adding a few extensions to the VTK library to support spline curve- and surface-based visualizations. A spline is a numeric function that is piecewise-defined by polynomial functions and a B¨¦zier curve is a special case which only have one segment and is evaluated by Bernstein basis functions. In this project, we will focus on rational B¨¦zier patches which are parametric surfaces generated from the Cartesian product of two B¨¦zier curves, since isogeometric spline simulations of many different types can generate them. Supporting this widely-used parametric representation will give VTK the ability to visualize meshes and simulations processed from CAD models.