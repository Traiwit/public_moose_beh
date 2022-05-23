# -*- coding: utf-8 -*-
"""Base classes for other topo classes.
"""


############## Topography classes ########################################
#{Topology Classes

class TopoBase(object):
    """Base class for topo objects such as L{PointCloud}, L{StructuredGrid},
    (unstructured FE-) L{Mesh<bae.field_02.topo.mesh.Mesh>}.
    """

    def getSubTopo(self, item):
        """Mandatory service method selecting from the space dimension of a
        FieldsCollection. For FieldsCollection.__getitem__()

        @param item: a serial index or slice or the like
        """
        raise NotImplementedError(
            "Derived class of type %s must supply a getSubTopo method!"
            % type(self))

    def getPoints(self):
        raise NotImplementedError(
            "Derived class of type %s must supply a getPoints method!"
            % type(self))

    def __len__(self):
        raise NotImplementedError(
            "Derived class of type %s must supply a __len__ method!"
            % type(self))

    def getBoundingBox(self):
        raise NotImplementedError(
            "Derived class of type %s must supply a getBoundingBox method!"
            % type(self))

#} end of Topology Classes
