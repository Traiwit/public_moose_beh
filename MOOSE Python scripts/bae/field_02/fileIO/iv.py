"""Service module for writing OpenInventor ".iv" files.
WARNING: Don't use this directly. The interface of this module is subject to
changes. Instead use 
"""

todo = """

- for writing binary files:
  . didn't find a verbatim description of the format, so only way is to check
    source code of implemnetation, e.g. coin3D.
  . ivcat from Inventor Tools (https://sourceforge.net/p/inventor-tools/)
    converts ascii to binary.
  . Seems like an object of class SoOutput gets data from someone who appears
    to be unaware of the write-mode binary or ascii.
    Here is the code that apparently writes the iv file for ivcat.
    https://sourceforge.net/p/inventor-tools/code/HEAD/tree/trunk/ivcat/ivcat.cpp#l311
---snip---
    // write stuff
    SoOutput out;
    out.setBinary(binary);

    OPEN_OUTPUT_FILE(&out, outputfile, &print_usage);

    if (verbose)
        fprintf(stderr, "Writing...\n");

    SoWriteAction writer(&out);

    for (int i = 0; i < nodeList.getLength(); i++) {
        writer.apply(nodeList[i]);
    }

    CLOSE_OUTPUT_FILE(&out, outputfile);
---snap---

    So maybe the structure of ascii and binary files is actually the same
    only the actual "formatting" differs. And this formatting seems to be
    done in SoOutput and its various write-methods:
    https://github.com/coin3d/coin/blob/master/src/io/SoOutput.cpp

    Maybe this is not true but then I don't know yet how the SoAction.apply()
    mechanism works. Apparently (from the code above) it's a SoWriteAction
    object that does the actual work through its apply method which is
    inherited from SoAction:

    https://github.com/coin3d/coin/blob/master/include/Inventor/actions/SoAction.h
    https://github.com/coin3d/coin/blob/master/include/Inventor/actions/SoWriteAction.h
    https://github.com/coin3d/coin/blob/master/src/actions/SoWriteAction.cpp
    https://github.com/coin3d/coin/blob/master/src/actions/SoAction.cpp

    Check out the SoAction.apply method. This should use SoWriteAction-methods
    or attributes...
"""

class IvWriter(object):
    """For writing files in OpenInventor (.iv) format. This is suitable for
    import as shapes into voxler.

    Usage:
     >>> writer = IvWriter(open("output.iv", "w"))
     >>> writer.writeTriMeshToIv(mesh1, colour=(0.8,0.2,0.2))
     >>> writer.writeTriMeshToIv(mesh2, colour=(0.0,0.8,0.8))
     >>> del writer
    """

    def __init__(self, output):
        """@param output: open file like object
        """
        self.output = output
        self.output.write("#Inventor V2.0 ascii\n")
        self.output.write("Separator {\n")

    def close(self):
        self.output.write("}\n")

    def __del__(self):
        self.close()

    def writeTriMeshToIv(self, mesh, colour=None, transparency=0):
        """Export the triangle surface to the iv file.

        @param mesh: a L{TriMesh} object
        @param colour: RGB-tuple, three floats between 0 and 1
        @param transparency: float between 0 and 1
        """

        self.output.write("  Separator {\n")

        # write material
        if colour is not None or transparency>0:
            if colour is None:
                colour = (1,1,1)
            else:
                colour = tuple(colour)
            self.output.write("    Material {\n")
            self.output.write("      ambientColor 0 0 0\n")
            self.output.write("      diffuseColor %g %g %g\n" % colour)
            self.output.write("      specularColor 0 0 0\n")
            self.output.write("      emissiveColor 0 0 0\n")
            self.output.write("      shininess 1\n")
            self.output.write("      transparency %g\n" % transparency)
            self.output.write("    }\n")
            self.output.write("    ShapeHints {\n")
            self.output.write("      vertexOrdering COUNTERCLOCKWISE\n")
            self.output.write("      shapeType UNKNOWN_SHAPE_TYPE\n")
            self.output.write("      faceType CONVEX\n")
            self.output.write("      creaseAngle 0\n")
            self.output.write("    }\n")

        # write point / node data
        self.output.write("    Coordinate3 {\n    point [\n")
        for coords in mesh.nodeCoords:
            self.output.write("  %g %g %g,\n" % tuple(coords))
        self.output.write("      ]\n  }\n")

        # write face / triangle data
        self.output.write("    IndexedFaceSet {\n    coordIndex [\n")
        for pts in mesh.elNodes:
            self.output.write("  %d,%d,%d,-1,\n" % tuple(pts))
        self.output.write("      ]\n  }\n")

        # fini
        self.output.write("  }\n")
