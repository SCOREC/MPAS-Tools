from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
import jigsawpy
from jigsawpy.savejig import savejig

from mpas_tools.logging import check_call


def jigsaw_driver(cellWidth, x, y, on_sphere=True, earth_radius=6371.0e3,
                  geom_points=None, geom_edges=None, geom_bounds=None,
                  preserve_geometry=False, write_vtk=False, logger=None):
    """
    A function for building a jigsaw mesh

    Parameters
    ----------
    cellWidth : ndarray
        The size of each cell in the resulting mesh as a function of space

    x, y : ndarray
        The x and y coordinates of each point in the cellWidth array (lon and
        lat for spherical mesh)

    on_sphere : logical, optional
        Whether this mesh is spherical or planar

    earth_radius : float, optional
        Earth radius in meters

    geom_points : ndarray, optional
        list of point coordinates for geometric model vertices

    geom_edges : ndarray, optional
        list of edges between points in geom_points that define geometric
        model edges

    geom_bounds : ndarray, optional
        list of geometric model entities that define bounding geometric features

    preserve_geometry : logical, optional
        by default Jigsaw does not ensure that all model entities have a mesh
        entitiy associated with them (i.e., a mesh vertex at the location of
        each model vertex).  Setting this option to `True` will ensure that all
        model entities are represented in the mesh via Jigsaw's `mesh_top1`
        option.

    write_vtk : logical, optional
        write the triangular mesh to a VTK file for vizualization in ParaView

    logger : logging.Logger, optional
        A logger for the output if not stdout
    """
    # Authors
    # -------
    # Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    print("----foo")

    # setup files for JIGSAW
    opts = jigsawpy.jigsaw_jig_t()
    opts.geom_file = 'mesh.msh'
    opts.jcfg_file = 'mesh.jig'
    opts.mesh_file = 'mesh-MESH.msh'
    opts.hfun_file = 'mesh-HFUN.msh'

    # save HFUN data to file
    hmat = jigsawpy.jigsaw_msh_t()
    if on_sphere:
       hmat.mshID = 'ELLIPSOID-GRID'
       hmat.xgrid = numpy.radians(x)
       hmat.ygrid = numpy.radians(y)
    else:
       hmat.mshID = 'EUCLIDEAN-GRID'
       hmat.xgrid = x
       hmat.ygrid = y
    hmat.value = cellWidth
    jigsawpy.savemsh(opts.hfun_file, hmat)

    # define JIGSAW geometry
    geom = jigsawpy.jigsaw_msh_t()
    if on_sphere:
       geom.mshID = 'ELLIPSOID-MESH'
       geom.radii = earth_radius*1e-3*numpy.ones(3, float)
    else:
       geom.mshID = 'EUCLIDEAN-MESH'
       geom.vert2 = geom_points
       geom.edge2 = geom_edges
       geom.bound = geom_bounds
    jigsawpy.savemsh(opts.geom_file, geom)

    # build mesh via JIGSAW!
    opts.hfun_scal = 'absolute'
    opts.hfun_hmax = float("inf")
    opts.hfun_hmin = 0.0
    opts.mesh_dims = +2  # 2-dim. simplexes
    opts.optm_qlim = 0.9375
    opts.verbosity = +1
    opts.mesh_top1 = preserve_geometry

    savejig(opts.jcfg_file, opts)
    check_call(['jigsaw', opts.jcfg_file], logger=logger)
  
    if write_vtk:
       mesh = jigsawpy.jigsaw_msh_t()
       jigsawpy.loadmsh(opts.mesh_file, mesh)
       jigsawpy.savevtk("mesh.vtk", mesh)
