# -*- coding: utf-8 -*-
"""Topology classes for unstructured meshes

A "mesh" is a collection of nodes and elements. There are meshes with all the
same type of elements, they are of type L{HomoMesh} or one of its subclasses
like L{TriMesh} or L{TetMesh}. And there are meshes with mixed types of
elements of type L{HeteroMesh}.

A mesh object is not directly suitable as topo object of a
L{FieldsCollection}-object. Instead topo objects are of type L{MeshNodes},
L{MeshElems} or L{MeshGaussPts}. Those topo objects all have a mesh attribute
of one of the mesh types, i.e. L{HomoMesh} or L{HeteroMesh}.

Element type information is stored with the mesh object as an object of type
L{ElementType}. An ElementType-object holds all relevant (and known/determined)
attributes of a certain element type like shape, nodesPerElem, midSideNodes,
gaussPtElCoords, Abaqus name.


More specific aspects:
 - Nodes and elements are identified through their index in mesh.nodeCoords
   and mesh.elNodes.
 - Node and element numbers/labels are not stored in the mesh object. They can
   be stored in a corresponding field. (Hence storing a mesh with node and
   element labels would require at least two FieldsCollection objects
   --one with a MeshNodes topo for the node lables and one with a MeshElems
   topo for the element labels.
 - A mesh as well as a topo is an immutable object. Adding and removing nodes
   or elements is possible only through creating a new mesh object.


Programming notes and style guide:
 - consider mesh_01.Mesh names
"""
