import sys
import vtk
from vtk.test import Testing

class PatchInterpolation(Testing.vtkTest):
    """Test interpolation of values on various-dimensional patches."""

    def verifyRationalCurve(self, ctrlPts, f, prange, ns, name):
        """Evaluate points along a curve.

        ctrlPts  is a list of 4-tuple control points defining a rational
                 Bezier curve (the 4-th coordinate is the weight).
        f        is an implicit function that takes a point in and returns
                 zero if it lies on the curve.
        prange   is list containing minimum and maximum parameter values.
        ns       is the number of points to sample between the min/max.
        """
        import vtk
        npts = len(ctrlPts)
        degree = [npts - 1,0,0]
        self.assertGreater(degree[0], 0,
            'Need at least 2 control points, got {n}'.format(n=npts))

        cpt = vtk.vtkDoubleArray()
        cpt.SetNumberOfComponents(len(ctrlPts[0]))
        cpt.SetNumberOfTuples(npts)
        [cpt.SetTuple(i, ctrlPts[i]) for i in range(npts)]
        erp = vtk.vtkPatchInterpolation()
        pts = vtk.vtkPoints()
        pd = vtk.vtkPolyData()
        ln = vtk.vtkCellArray()
        pd.SetPoints(pts)
        pd.SetLines(ln)
        params = [0., 0., 0.]
        delta = (prange[1] - prange[0]) / (ns - 1.)
        for r in [prange[0] + delta * x for x in range(ns)]:
            params[0] = r
            erp.InterpolateOnPatch(pts.GetData(), 1, cpt, degree, params)

            # Test that f(xyz) = 0 to within 6 decimal places. We use 6 places
            # because vtkPoints stores values as 'float' by default and floats
            # are at places only accurate to 7 or 8 decimal digits.
            #xyz = pts.GetPoint(pts.GetNumberOfPoints() - 1)
            #self.assertAlmostEqual(f(xyz), 0., 6,
            #    'Point does not satisfy implicit function; f{p} = {f} != 0'.format(p=xyz, f=f(xyz)))

        # For debugging, uncomment the lines below. They will
        # write a VTK polydata file containing the output points.
        ##   Insert line segments connecting the points and
        ##   write out the resulting polydata to a file:
        [ln.InsertNextCell(2, [i, i+1]) for i in range(ns - 1)]
        #wri = vtk.vtkPolyDataWriter()
        #wri.SetInputDataObject(pd)
        #wri.SetFileName('bezier-{shape}.vtk'.format(shape=name))
        #wri.Write()
		
		# Create Instance for Rendering
        #
        #pdMapper = vtk.vtkPolyDataMapper()
        #pdMapper.SetInputData(pd)
        #pdActor = vtk.vtkActor()
        #pdActor.SetMapper(pdMapper)

        # Create rendering stuff
        #
        #ren = vtk.vtkRenderer()
        #renWin = vtk.vtkRenderWindow()
        #renWin.AddRenderer(ren)
        #ren.AddActor(pdActor)
        #renWin.SetSize(400, 150)
        #renWin.Render()
        #winToImg = vtk.vtkWindowToImageFilter()
        #winToImg.SetInput(renWin)
        #winToImg.Update()
        #imgWri = vtk.vtkPNGWriter()
        #imgWri.SetFileName('bezier-{shape}.png'.format(shape=name))
        #imgWri.SetInputConnection(winToImg.GetOutputPort())
        #imgWri.Write()
        
        #img_file = 'bezier-{shape}.png'.format(shape=name)
        #vtk.test.Testing.compareImage(renWin, vtk.test.Testing.getAbsImagePath(img_file), threshold=25)
        
        return pd

    def createRationalLine(self):
        return self.verifyRationalCurve(
            [[0,0,0,1], [1,1,1,1]],
            lambda x: x[0] - x[1],
            [0, 1], 3, 'line')

    def createRationalCircle(self):
        return self.verifyRationalCurve(
            [[1,0,0,1], [1,1,0,1], [0,2,0,2]],
            lambda x: 1. - x[0]**2 - x[1]**2,
            [0, 1], 21, 'circle')
            
    def createRationalHyperbola(self):
        import vtk
        erp = vtk.vtkPatchInterpolation()
        cpt = vtk.vtkDoubleArray()
        pts = vtk.vtkPoints()
        pd = vtk.vtkPolyData()
        ln = vtk.vtkCellArray()
        pd.SetPoints(pts)
        pd.SetLines(ln)
        prange = [0, 1]
        ns = 21;
        
        # generate control points of ellipse for one quadrant
        erp.GenerateHyperbolaCtrlPt(cpt, 1., 3.)
        
        npts = cpt.GetNumberOfTuples()
        degree = [npts - 1,0,0]
        self.assertGreater(degree[0], 0,
            'Need at least 2 control points, got {n}'.format(n=npts))

        # compute points coordinates
        params = [0., 0., 0.]
        delta = (prange[1] - prange[0]) / (ns - 1.)
        for r in [prange[0] + delta * x for x in range(ns)]:
            params[0] = r
            erp.InterpolateOnPatch(pts.GetData(), 1, cpt, degree, params)

        [ln.InsertNextCell(2, [i, i+1]) for i in range(ns - 1)]
        
        #wri = vtk.vtkPolyDataWriter()
        #wri.SetInputDataObject(pd)
        #wri.SetFileName('bezier-{shape}.vtk'.format(shape='hyperbola'))
        #wri.Write()
        
        return pd
            
    def createRationalEllipse(self):
        import vtk
        # declare vals
        erp = vtk.vtkPatchInterpolation()
        cpt = vtk.vtkDoubleArray()
        pts = vtk.vtkPoints()
        pd = vtk.vtkPolyData()
        ln = vtk.vtkCellArray()
        pd.SetPoints(pts)
        pd.SetLines(ln)
        prange = [0, 1]
        ns = 21
        center = vtk.vtkVector3d()
        center.SetX(0.0)
        center.SetY(0.0)
        center.SetZ(0.0)
        majorAxis = vtk.vtkVector3d()
        majorAxis.SetX(1.0)
        majorAxis.SetY(0.0)
        majorAxis.SetZ(0.0)
        minorAxis = vtk.vtkVector3d()
        minorAxis.SetX(0.0)
        minorAxis.SetY(1.414)
        minorAxis.SetZ(1.414)        
        
        for quadrant in range(1, 5):
            # generate control points of ellipse for one quadrant
            erp.GenerateEllipseCtrlPt(cpt, center, majorAxis, minorAxis, quadrant)
            
            npts = cpt.GetNumberOfTuples()
            degree = [npts - 1,0,0]
            self.assertGreater(degree[0], 0,
                'Need at least 2 control points, got {n}'.format(n=npts))

            # compute points coordinates
            params = [0., 0., 0.]
            delta = (prange[1] - prange[0]) / (ns - 1.)
            for r in [prange[0] + delta * x for x in range(ns)]:
                params[0] = r
                erp.InterpolateOnPatch(pts.GetData(), 1, cpt, degree, params)

            [ln.InsertNextCell(2, [(quadrant-1)*ns + i, (quadrant-1)*ns + i+1]) for i in range(ns - 1)]
        
        wri = vtk.vtkPolyDataWriter()
        wri.SetInputDataObject(pd)
        wri.SetFileName('bezier-{shape}.vtk'.format(shape='ellipse3d'))
        wri.Write()
        
        return pd
            
    def testRationalPatches(self):
        pdLine = self.createRationalLine()
        pdCircle = self.createRationalCircle()
        pdEllipse = self.createRationalEllipse()
        pdHyperbola = self.createRationalHyperbola()
        
        # Create Instance for Rendering
        #
        pdLineMapper = vtk.vtkPolyDataMapper()
        pdLineMapper.SetInputData(pdLine)
        pdLineActor = vtk.vtkActor()
        pdLineActor.SetMapper(pdLineMapper)
        pdLineActor.AddPosition(0.0, 0.0, 0.0)
        
        pdCircleMapper = vtk.vtkPolyDataMapper()
        pdCircleMapper.SetInputData(pdCircle)
        pdCircleActor = vtk.vtkActor()
        pdCircleActor.SetMapper(pdCircleMapper)
        pdCircleActor.AddPosition(1.0, 0.0, 0.0)
        
        pdEllipseMapper = vtk.vtkPolyDataMapper()
        pdEllipseMapper.SetInputData(pdEllipse)
        pdEllipseActor = vtk.vtkActor()
        pdEllipseActor.SetMapper(pdEllipseMapper)
        pdEllipseActor.AddPosition(4.0, 0.0, 0.0)
        
        pdHyperbolaMapper = vtk.vtkPolyDataMapper()
        pdHyperbolaMapper.SetInputData(pdHyperbola)
        pdHyperbolaActor = vtk.vtkActor()
        pdHyperbolaActor.SetMapper(pdHyperbolaMapper)
        pdHyperbolaActor.AddPosition(5.0, 0.0, 0.0)

        # Create rendering stuff
        #
        ren = vtk.vtkRenderer()
        ren.AddActor(pdLineActor)
        ren.AddActor(pdCircleActor)
        ren.AddActor(pdEllipseActor)
        ren.AddActor(pdHyperbolaActor)
        
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)        
        renWin.SetSize(400, 150)
        renWin.Render()
        
        #winToImg = vtk.vtkWindowToImageFilter()
        #winToImg.SetInput(renWin)
        #winToImg.Update()
        #imgWri = vtk.vtkPNGWriter()
        #imgWri.SetFileName('PatchInterpolation.png')
        #imgWri.SetInputConnection(winToImg.GetOutputPort())
        #imgWri.Write()
        
        img_file = "PatchInterpolation.png"
        vtk.test.Testing.compareImage(renWin, vtk.test.Testing.getAbsImagePath(img_file), threshold=25)
        vtk.test.Testing.interact()

if __name__ == '__main__':
    Testing.main([(PatchInterpolation,'test'),])
