"""
Module that provides means to run and manipulate Voxler.

Usage
=====
 >>> from bae.voxler_01 import BEVoxler
 >>> vxl = BEVoxler()
 >>> vxl.createModule("Annotation")
 'Annotation' module created.
 >>> vxl.showModule("Annotation", False)
 'Annotation' module created.

Known bugs/problems
===================
initialising the BEVoxler by using BEVoxler("name.voxb") does not open
the desired file but manipulates an already open voxler or opens a new voxler.
Usually, an insatance of a random Voxler file (usually the first one that has 
been openned) is craeted. If this Voxler instance is closed, a new Voxler file 
will be openned if an instance of BEVoxler class is called again.
"""

__version__ = "1.0"

_version_history_ = r"""
1.0 BG: basic implementation of core functions. Almost all the needed
        functionality of Voxler can be built based on the core functions
        which exist in this version.
        Enjoy! :)
1.1 BG: added command_sequence function. This function sometimes works better
        especially when there are many commands to execute, e.g. modifyGridder,
        however, with some changes command_sequence can be avoided.
        added new functions: importCsv, modifyMath, exportModule, doGridder,
                             modifyGridModuleFltSlip
"""


_todo_ = """
The focus should be to implement a method to reduce usage of time consuming
functions in Voxler, i.e. trying to process the data before importing them into
Voxler would be more efficient. For example, running a gridder module in voxler
takes a long time.
Or some specific functionality needs to be implemented in a more user friendly
platform, like Paraview.

- setModuleOffset functions can be grouped into one function with additional
  input parameter as the moduleType. This can be done in Version2.
ideas:
======

"""

import win32com
import win32com.client
import os

###############################################################################
### ------------ Classes from Golden Software
###############################################################################
class ApiError(Exception):
    """
    Class to capture the Api Errors.
    """
    def __init__(self, msg):
        super(ApiError, self).__init__()
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


class Command(object):
    """
    This class encapsulates the error checking and ensures that exceptions are
    thrown on failure.
    
    @Note:the core functionality of "command" which is used in BEVoxler class 
    is defined here.
    """

    def __init__(self, api, cmd):
        self.api = api
        self.cmd = cmd

        # Catch all exceptions, get all bad return codes.
        cmd_constructed = False
        try:
            cmd_constructed = self.api.Construct(self.cmd)
        except:
            pass

        if not cmd_constructed:
            raise ApiError("Failed to construct command '{0}'".format(self.cmd))

    def option(self, key, value):
        opt_valid = False
        try:
            opt_valid = self.api.Option(key, value)
        except:
            pass

        if not opt_valid:
            raise ApiError("Failed to set option '{0}':'{1}'".format(key, value))

    def do(self):
        done = False
        try:
            done = self.api.Do()
        except:
            pass

        if not done:
            raise ApiError("Failed to execute command '{0}'".format(self.cmd))
# --- No need to have "doOnce", "do" does the same!
#    def doOnce(self):
#        done = False
#        try:
#            done = self.api.DoOnce()
#        except:
#            pass
#
#        if not done:
#            raise ApiError("Failed to execute command '{0}'".format(self.cmd))


###############################################################################
### ------------ BEVoxler class
###############################################################################

class BEVoxler(object):
    """
    This class is the main class to manipulate Voxler.
    
     >>> from vxl_core import BEVoxler
     >>> vxl = BEVoxler()
     >>> vxl.createModule("Annotaion")
     'Annotation' module created.
    
    """
    def __init__(self, vxbFilename=None, clearAllModules=False):
    
        try:
            self.app = win32com.client.GetObject(Class="Voxler.Application")
        except:
            self.app = win32com.client.Dispatch("Voxler.Application")

        self.app.Visible = 1
        self.api = self.app.CommandApi
        
        ### todo: do something differnet depending on vxbFilename
        #cmd = Command(self.api, "New")
        #cmd.do()
        self.vxbFilename = vxbFilename
        if False:
            defaultOpts = []
            if vxbFilename is not None:
                defaultOpts.append( ("Path", vxbFilename) )
            
            self.command(
                cmd_name="New", 
                opts=dict(defaultOpts),
                )

        if clearAllModules:
            self.command(
                cmd_name="DeleteAllModules", 
                opts=dict(),
                )
        
    def command(self, cmd_name, opts):
        """ 
        Executes a single command, and sets all of the options.
        
        @Note: Basically there is no need to use this function. Additional 
        functionality can be added to BEVoxler by using this function.
        
        The following example creates an Axes module that connects to the 
        ScatterPlot:
         >>> voxler.command('CreateModule', {'Type':'Axes', 
                            'AutoConnect':'True', 
                            'SourceModule':'ScatterPlot'})
        """
        cmd = Command(self.api, cmd_name)

        for k, v in opts.items():
            cmd.option(k, v)

        cmd.do()

    def command_sequence(self, cmd_name, opts, opts_seq):
        """Executes a command + opts, once for every option in opts_seq.

        Note: the options are not set in order.  The idea for this function is to
        set a whole bunch of options for a module or whatever, and not have to
        keep calling Construct, Option, Do.

        Example:
            Modifies the Axes Module, and for each of the key/values in seqOpts
            executes a new ModifyModule command on the "Axes" module, setting each
            of the options individually.

            voxler.command(
                # This command is created and this option is set
                "ModifyModule", {"Module":"Axes"},

                # For each key/value pair, the above command is executed.
                {
                    "AxesXAxisShow":"True",
                    "AxesXAxisColor":"light green",
                    "AxesXAxisBeg":"20",
                    "AxesXAxisEnd":"110",
                    "AxesXAxisCross1":"0",
                    "AxesXAxisCross2":"50",
                    "AxesXAxisPlane":"45",
                    "AxesXAxisTitle":"X Axis Title",
                    "AxesXAxisShowLabels":"True",
                    "AxesXAxisFlipLabelsX":"True",
                    "AxesXAxisFlipLabelsY":"True"
                }
            )
        """

        for cmdOptKey, cmdOptVal in opts_seq.items():

            # Keep things clean and make the command each time.
            cmd = Command(self.api, cmd_name)

            # Set the initial options common to all commands in this sequence.
            for k, v in opts.items():
                cmd.option(k, v)

            # Set the option for this key/value pair.
            cmd.option(cmdOptKey, cmdOptVal)

            cmd.do()

# --- no need to have "commandOnce", "command" does the same - similar to 
# "do" & "doOnce" commands
#    def commandOnce(self, cmd_name, opts):
#        """ Executes a single command, and sets all of the options.
#
#        Example:
#
#        # Creates an Axes module that connects to the ScatterPlot
#        voxler.command('CreateModule', {'Type':'Axes', 'AutoConnect':'True', 'SourceModule':'ScatterPlot'})
#        """
#        cmd = Command(self.api, cmd_name)
#
#        for k, v in opts.items():
#            cmd.option(k, v)
#
#        cmd.doOnce()

# --- following functions have no particular use in voxler automation
#    def path(self):
#        """Returns the path that Voxler is running from."""
#        return self.app.Path
#    
#    def name(self):
#        return self.app.Name
#
#    def fullname(self):
#        return self.app.FullName
#
#    def version(self):
#        return self.app.Version
#    
#    # one should not ever use this save function. One cannot undo after using
#    # python-voxler automation, thus, once a file is saved it is changed 
#    # eternally!
#    def save(self):
#        self.app.Save()
#    
#    def quit(self):
#        """Quits Voxler, obviously."""
#        self.app.Quit()
    
    
##########################################
#----------------------------------------- Common use functions
########################################## 
    def createModule(self,
                   moduleType,
                   ):
        """
        @param moduleType: STRING, one of the deined module types in Voxler
        """
        self.command(cmd_name="CreateModule",
                         opts=dict([
                                 ("Type", moduleType),
                                 ]),
                        )
        print("'%s' module created." % moduleType)
        return
    
    def deleteModule(self,
                   moduleName,
                   ):
        """
        @param moduleType: STRING
        """
        self.command(cmd_name="DeleteModule",
                         opts=dict([
                                 ("Module", moduleName),
                                 ]),
                        )
        print("Module '%s' deleted." % moduleName)
        return
    
    def connectModules(self,
                   sourceModule,
                   targetModule,
                   targetPort = 0,
                   ):
        """
        @param sourceModule: STRING
        @param targetModule: STRING
        @param targetPort: INT, if does not support multiple ports, use default
        i.e. 0.
        """
        self.command(cmd_name="ConnectModules",
                         opts=dict([
                                 ("SourceModule", sourceModule),
                                 ("TargetModule",targetModule),
                                 ("TargetPort",targetPort),
                                 ]),
                        )
        print("Module '%s' connected to '%s' via port %d."
              % (sourceModule, targetModule, int(targetPort)))
        return
    
    def connectMultiModules(self, 
                            modulePairs, 
                            targetPort = None,
                            ):
        """
        This function connects pairs of modules defined in a dictionary.

        @Note: use the following to generate modulePairs from raw strings
        without frame numbers
         >>>for m in modulePairs.keys():
         >>>    modulePairs[m % timeStamp] = modulePairs.pop(m)

        @param modulePairs: DICTIONARY name of the modules.
        For example:
         >>> modulePairs = {
         >>>          "VTK_B01"      : "JNC_B01",
         >>>          "VTK_B02"      : "JNC_B02", }

        @param TargetPort: LIST, If None, set to a list of zeros.
        """
        if targetPort is None:
            targetPort = [0] * len(modulePairs.keys())
        
        for i,module0 in enumerate(modulePairs.keys() ):
            module1 = modulePairs[module0]
            target = targetPort[i]
            
            self.command(cmd_name="ConnectModules",
                             opts=dict([
                                     ("SourceModule", module0),
                                     ("TargetModule", module1),
                                     ("TargetPort", target),
                                     ]),
                            )
            print("Module '%s' connected to '%s' via port %d." 
                  %(module0, module1, target))
        return
    
    def renameModules(self,
                      oldName,
                      newName,
                      ):
        """
        @param oldName: STR, the old name!
        @param newName: STR, the new name!
        """
        self.command(cmd_name="RenameModule",
                         opts=dict([
                                 ("Module", oldName),
                                 ("Name", newName),
                                 ("GuiEnabled", "False"),
                                 ]),
                        )
        print("Module name '%s' changed to '%s'." %(oldName, newName))
        return
    
    def changeAnnotationText(self,
                         annotationModuleName,
                         text,
                         ):
        """
        @param annotationModuleName: STR, name of the annotation module
        @param text: STR, the text to be shown. This text can be generated
        in a loop in a driver script.
        """
        self.command(cmd_name="ModifyModule",
                         opts=dict([
                                 ("Module", annotationModuleName),
                                 ("AnnotationText", text),
                                 ]),
                        )
        print("Annotation text in module '%s' changed to '%s'." 
              %(annotationModuleName,text))
        return
    
    def showModule(self,
                   moduleName,
                   show = True,
                   ):
        """
        This function returns a warning message if the module is on and the 
        operator tries to turn it on and vice versa. 
        
        @param moduleName: STR
        @param show: BOOLEAN, True = on, False = off
        """
        try:
            self.command(cmd_name="ShowModule",
                             opts=dict([
                                     ("Module",moduleName),
                                     ("Show",show),
                                     ]),
                            )
            print("Module '%s' show variable set to '%s'." 
                  %(moduleName,str(show)))
        except ApiError, exc:
            if str(exc).startswith('"Failed to execute command'):
                print("*WARNING: module '%s' show option is already set to'%s'."
                      %(moduleName, str(show)))
            else:
                raise
        return
    
    def changeIsosurfaceValue(self,
                              moduleName,
                              isoValue,
                              ):
        """
        @param moduleName: STRING
        @param isoValue: FLOAT or INT 
        """
        self.command(cmd_name="ModifyModule",
                         opts=dict([
                                 ("Module",moduleName),
                                 ("IsosurfaceIsovalue",float(isoValue)),
                                 ]),
                        )
        print("Isosurface value for module '%s' changed to '%f'." 
              %(moduleName,float(isoValue)))
        return
    
    def changeIsosurfaceValueMulti(self,
                                   isovals,
                                   ):
        """
        Uses changeIsosurfaceValue function.

        @param isovals: DICT, pairs of isosurface module name and value.
        Pairs dhould be {STRING : FLOAT or INT}
         >>> isovals = {"LOW_15E3_D01" : 15000,
         >>>            "MED_250E3_D01" : 150000,
         >>>            "HIGH_500E3_D01" : 750000,
         >>>            }
        """
        for i in range(len(isovals.keys())):
            name = isovals.keys()[i]
            val = float(isovals[name])
            self.changeIsosurfaceValue(name, val)
        return
    
    def changeDrawStyle(self,
                        moduleName,
                        drawStyle = 0,
                        ):
        """
        @param moduleName: STRING, only modules which support drawStyle
        @param drawStyle: INT, as follow
                            drawStyle = 0 --> As is
                            drawStyle = 1 --> Shaded
                            drawStyle = 2 --> Lines
                            drawStyle = 3 --> Points
        """
        self.command(cmd_name="ModifyModule",
                         opts=dict([
                                 ("Module",moduleName),
                                 ("GeomSrcDrawStyle",drawStyle),
                                 ]),
                        )
        print("DrawStyle in module %s changed to %d." 
              %(moduleName,int(drawStyle)))    
        return
    
    def moveModule(self,
                   moduleName,
                   positionX,
                   positionY,
                   ):
        """
        @param moduleName: STR
        @param positionX: INT, position in pixel
        @param positionY:INT, position in pixel
        """
        
        self.command(cmd_name="MoveModule",
                         opts=dict([
                                 ("Module", moduleName),
                                 ("UseOriginalPosition", "False"),
                                 ("XPosition", str(positionX)),
                                 ("YPosition", str(positionY)),
                                 ]),
                        )
        print("Module %s moved to (X = %d, Y = %d)" 
              %(moduleName, positionX, positionY))
        return
    
    def scatterPlotOps(self,
                       moduleName,
                       component,
                       ):
        """
        This function can be expanded by adding different parameters of
        scatterPlot module to it. At the moment, it can only::
        1. change input component

        @param moduleName: STR
        @param component: INT, the input component.
        """
        self.command(cmd_name="ModifyModule",
                         opts=dict([
                                 ("Module", moduleName),
                                 ("ScatterPlotColorComponent", component)
                                 ]),
                        )
        print("scatterPlot module '%s', colour component changed to component: %d." 
              % (moduleName,component))
        return
    
    def setObliqueOffset(self,
                   moduleName,
                   offset,
                   ):
        """
        @param moduleName: STR
        @param Offset: FLOAT, is the offset value from centre (cutting plane).
        """
        self.command(cmd_name="ModifyModule",
                         opts=dict([
                                 ("Module", moduleName),
                                 ("ObliqueSlicePlaneOffset", float(offset)),
                                 ]),
                        )
        print("'%s' module offset set to %d." % (moduleName, float(offset)) )
        return
    
    def setContourOffset(self,
                         moduleName,
                         offset,
                         ):
        """
        @param moduleName: STR
        @param Offset: FLOAT, is the offset value from centre (cutting plane).
        """
        self.command(cmd_name="ModifyModule",
                         opts=dict([
                                 ("Module", moduleName),
                                 ("ContourPlaneOffset", float(offset)),
                                 ]),
                        )
        print("'%s' module offset set to %d." % (moduleName, float(offset)) )
        return
    
    def setClipPlaneOffset(self,
                           moduleName,
                           distance,
                           ):
        """
        Changes the clippingPlane module offset. This should not be mixed up
        with setObliqueOffset. The clippingPlane and ObliqueImage modules
        work completely separately.

        @param moduleName: STR
        @param distance: the distance value from centre (FLOAT)
        """
        self.command(cmd_name="ModifyModule",
                     opts=dict([
                         ("Module", moduleName),
                         ("ClipPlaneDistance", float(distance)),
                     ]),
                 )
        print("'%s' module clip plane distance set to %d."
              % (moduleName, float(distance)) )
        return
    
    def modifyMath(self,
                   moduleName,
                   compDict = {1 : "A1"}
                   ):
        """
        Changes the output components of the math modules.
        
        @param moduleName: STR, name of the Math module
        @param compDict: DIC, a dictionary of component number vs expression.
        
        @Note: Voxler automation can take maximum of 20 components. If more
            than 20 components are given as input, an expection will be raised.
        
        To do:
            Geomtery of the math can be changed as well. This function can be 
            modified to do other modifications to math modules.
        """
        compNum = len(compDict.keys())
        
        if compNum > 20:
            raise Exception("Math module can take maximum of 20 components in"
                            "Voxler automation, input has got %d components" 
                            %compNum)
        
        self.command(cmd_name="ModifyModule",
                     opts=dict([
                             ("Module", moduleName),
                             ("LatticeMathOutComps", compNum),
                             ])
                    )
        for comp in compDict.keys():
            option = "LatticeMathExprComp" + str(comp)
            self.command(cmd_name="ModifyModule",
                     opts=dict([
                             ("Module", moduleName),
                             (option, compDict[comp]),
                             ])
                    )
        print("module %s has been modified successfully." %moduleName)
        return

##########################################
#----------------------------------------- gridding functions
########################################## 

    def doGridder(self,
                  moduleName
                  ):
        """
        Does the same as "Begin Gridding" in gridder module.
        
        @param moduleName: STR, name of the Math module
        """
        print("Gridding %s..." %moduleName)
        self.command(cmd_name="ModifyModule",
                         opts=dict([
                                 ("Module", moduleName),
                                 ("Gridder3DoIt", True),
                                 ]),
                        )
        return
    
    def modifyGridModuleFltSlip(self,
                                moduleName,
                                inComponent,
                                gridGeom,
                                gridRes,
                                gridSearchR = 7,
                                gridSearchMax = 50,
                                gridSearchMin = 1,
                                doIt = True,
                                ):
        """
        Changes the parameters in gridder modules. The parameters to modify 
        here are specific to faultSlip RER. That is why this function is 
        called ...FltSlip.
        
        @param moduleName: STR, name of the gridder module
        @param inComponent: INT, component number of the input (starts from 0)
        @param gridGeom: LIST of Floats, with format of [X_min, X_max,
                                                         Y_min, Y_max,
                                                         Z_min, Z_max]
        @param gridSearchR: INT, the radius search in simple search.
        @param gridSearchMax: INT, the max count in simple search.
        @param gridSearchMin: INT, the min count in simple search.
        @param doIt: Boolean, if True grids the data. Does the same as "Begin 
            Gridding" in gridder module.
            
        To do:
            This function can only perform the matrics-maximum gridding.
            Modification can be done to be flexible or other functions to be 
            written which do the other types of gridding.
        """
        print("Setting gridder options for module '%s'." %moduleName)
        self.command_sequence(cmd_name="ModifyModule",
                             opts=dict([
                                     ("Module", moduleName),
                                     ]),
                         opts_seq=dict([
                                 ("Gridder3Component", str(inComponent)),
                                 ("Gridder3Method", "0"),
                                 ("Gridder3Metric", "4"),
                                 ("Gridder3XMin", str(gridGeom[0])),
                                 ("Gridder3XMax", str(gridGeom[1])),
                                 ("Gridder3YMin", str(gridGeom[2])),
                                 ("Gridder3YMax", str(gridGeom[3])),
                                 ("Gridder3ZMin", str(gridGeom[4])),
                                 ("Gridder3ZMax", str(gridGeom[5])),
                                 ("Gridder3XSpacing", str(gridRes[0])),
                                 ("Gridder3YSpacing", str(gridRes[1])),
                                 ("Gridder3ZSpacing", str(gridRes[2])),
                                 ("Gridder3Search", "1"), 
                                 ("Gridder3SearchRadius", str(gridSearchR)),
                                 ("Gridder3SearchMax", str(gridSearchMax)),
                                 ("Gridder3SearchMin", str(gridSearchMin)),
                                 ]),
                        )
        if doIt:
            print("Gridding %s..." %moduleName)
            self.command(cmd_name="ModifyModule",
                             opts=dict([
                                     ("Module", moduleName),
                                     ("Gridder3DoIt", doIt),
                                     ]),
                            )
        return

##########################################
#----------------------------------------- Export utilities
########################################## 

#=========================================  screenshot
    def screenshot(self, 
            fileName,
            exportPath = None,
            multiShots = False,
            ):
        """
        This function can generate "jpg" and "png" output.
        
        @Note: fileType is automatically taken from the extension of the 
        fileName input. If no extension is given, extension is automatically 
        set to "jpg".
        
        @param fileName: STR, file name of the screenshot to be saved.
        @param exportPath: STR, default is the same route as the py script
        @param multiShots: BOOLEAN, default set to False. An extra 
        functionality to be able to capture multiple screenshots.
        It adds "_Img%02d" to the end of defined name.
        """
        
        fileType = fileName.split(".")[-1]
        if fileType.upper() != "JPG" and  fileType.upper() != "PNG":
            print("Invalid file extension (fileName = '%s') for screenshot. FileType is set to 'jpg'" %fileType)
            fileType = 'jpg'
            fileName = fileName.split(".")[0] + "." + fileType
        
        if exportPath is None:
           exportPath = os.path.abspath(os.getcwd()).replace(os.sep,"\\")
           print("exportPath input is None. exportPath is set to: %s" % exportPath)
        
        if not os.path.isdir(exportPath):
            print("Export path does not exist! Export was UNSUCCESSFULL!")
            print("***NOTE: if special characters exist (e.g. %s) use r in front of string name." %r"\v")
            print("         for example use r'%s' instead of '%s'." %(r"\voxler",r"\voxler"))
            return
        
        if not multiShots:
            # grab a screenshot and save it to the defined path
            self.command(cmd_name="Export",
                             opts=dict([
                                     ("GuiEnabled", "False"),
                                     ("ClearOptions", "False"),
                                     ("Filter", fileType),
                                     ("Height", "1051"),
                                     ("Width", "796"),
                                     ("Path",exportPath + "\\" + fileName),
                                     ("PersistOptions", "True"),
                                     ("ModuleID","0"),
                                     ]),
                            )
            print("Screenshot saved to: \n%s" % (exportPath + "\\" + fileName))
        
        if multiShots:
            cnt = 1
            while True:
                key = raw_input("Choose a suitable view angle and enter 'C' to screenshot or any other keys to pass: ").upper()
                if key == "C":
                    nameExt = "_Img%02d." %cnt
                    self.command(cmd_name="Export",
                                     opts=dict([
                                             ("GuiEnabled", "False"),
                                             ("ClearOptions", "False"),
                                             ("Filter", fileType),
                                             ("Height", "1051"),
                                             ("Width", "796"),
                                             ("Path",exportPath + "\\" + fileName.split(".")[0] + nameExt + fileName.split(".")[-1]),
                                             ("PersistOptions", "True"),
                                             ("ModuleID","0"),
                                             ]),
                            )
                    print("Screenshot saved to: \n%s" % (exportPath + "\\" + fileName.split(".")[0] + nameExt + fileName.split(".")[-1]))
                    cnt += 1
                elif key != "C":
                    break
                else:
                    print("Uknown input!")                
        return

#========================================= exportScreenAsIV
    
    def exportScreenAsIV(self,
                 fileName,
                 exportPath = None,
                 ):
        """
        Exports everything shown in the SCREEN as IV, if they can be exported
        as IV format, e.g. Annotations cannot be exported as IV but 
        obliqueImages and isosurfaces can be.
        
        @Warning: results in huge IV files.
        
        @param fileName: STR, The name of the file. If extension is not given,
        ".iv" will be added to the fileName.                         
        @param exportPath: STRING variable. 
        """
        if exportPath is None:
           exportPath = os.path.abspath(os.getcwd()).replace(os.sep,"\\")
           print("exportPath input is None. exportPath is set to: \n %s" % exportPath)
        
        if not os.path.isdir(exportPath):
            print("Export path does not exist! Export was UNSUCCESSFULL!")
            print("***NOTE: if special characters exist (e.g. %s) use r in front of string name." %r"\v")
            print("         for example use r'%s' instead of '%s'." %(r"\voxler",r"\voxler"))
            return
        
        if fileName.split(".")[-1].lower() != "iv":
            fileName = fileName + ".iv"
        
        self.command(cmd_name="Export",
                     opts=dict([
                             ("GuiEnabled", "False"),
                             ("ClearOptions", "True"),
                             ("Filter", "iv"),
                             ("Path", exportPath + "\\" + fileName),
                             ("PersistOptions", "True"),
                             ("ModuleId", "0"),
                             ])
                    )   
        print("Export file saved to: \n %s" % (exportPath + "\\" + fileName))
        return
    
    
#========================================= exportScreenAsDXF
    
    def exportModuleAsDXF(self,
                         moduleName,
                         fileName = None,
                         exportPath = None,
                         ):
        """
        This function only works if the module is exportable to dxf format.
        Usually used for isosurface modules.
        
        @param moduleName: STR, the name of the module to save as dxf.  
        @param fileName: STR, the name of the file including extension. 
        If extension is not given, ".dxf" will be added to the fileName.
        @param exportPath: STRING variable. 
        """
        
        if fileName is None:
            fileName = moduleName + ".dxf"
        
        if fileName.split(".")[-1].lower() != "dxf":
            fileName = fileName + ".dxf"
        
        if exportPath is None:
           exportPath = os.path.abspath(os.getcwd()).replace(os.sep,"\\")
           print("exportPath input is None. exportPath is set to: \n %s" % exportPath)
        
        if not os.path.isdir(exportPath):
            print("Export path does not exist! Export was UNSUCCESSFULL!")
            print("***NOTE: if special characters exist (e.g. %s) use r in front of string name." %r"\v")
            print("         for example use r'%s' instead of '%s'." %(r"\voxler",r"\voxler"))
            return
        
        self.command(cmd_name="Export",
                     opts=dict([
                             ("GuiEnabled", "False"),
                             ("ClearOptions", "False"),
                             ("Filter", "dxf"),
                             ("Module", moduleName),
                             ("ModuleId", "0"),
                             ("Options","FileCompatibility=21"),
                             ("Options", "FormatASCII=1"),
                             ("Options", "MaxBitmapSizeInMB=10"),
                             ("Options","AllColorsSame=0"),
                             ("Options","AllCStylesSame=0"),
                             ("Options","AllWidthsSame=0"),
                             ("Options", "AllTextToAreas=0"),
                             ("Options","FillSolidAreas=0"),
                             ("Options", "UseSpatialInfo=0"),
                             ("Options","ColorMapping=0"),
                             ("Options", "RenderMarkers=0"),
                             ("Path", exportPath + "\\" + fileName),
                             ("PersistOptions", "True"),
                             ])
                    )   
        print("Export file saved to: \n %s" % (exportPath + "\\" + fileName))
        return
    
    
#========================================= general export module (workes with
#========================================= VDAT and DAT format)
    
    def exportModule(self,
                     moduleName,
                     fileType = "vdat",
                     fileName = None,
                     exportPath = None,
                     ):
        """
        This function is a combination of exportDat, exportVdat and exportVtk.
        Other file formats may work.
        
        @param moduleName: STR, name of the Math module
        @param fileType: STR, fileType without "." (lower or upper case)
        @param fileName: STR, Optional. Name of the output. If None, module 
            name will be used.
        @param exportPath: STR, Optional. Output path. If None, current folder
            will be used.
        
        @Note: VDAT is a voxler data file format which can only be used in 
            Voxler. It is an efficient format for Voxler.
        @Note: as of 2018JUL10, there is a bug in VTK export from Voxler.It 
            replaces the NAN values by a random large number instead of zeros.
        @Note: DAT file format is an inefficient and large format. Often 
            results in large files. Can be viewed in standard text editors.
        """
        if fileName is None:
            fileName = moduleName + "." + fileType
        
        if exportPath is None:
           exportPath = os.path.abspath(os.getcwd()).replace(os.sep,"\\")
           print("exportPath input is None. exportPath is set to: \n %s" 
                 % exportPath)
        
        if not os.path.isdir(exportPath):
            print("Export path does not exist! Export was UNSUCCESSFULL!")
            print("***NOTE: if special characters exist (e.g. %s) use r in "
                  "front of string name." %r"\v")
            print("         for example use r'%s' instead of '%s'." 
                  %(r"\voxler",r"\voxler"))
            return
        
        # different print options based on the fileType 
        if fileType.upper() == "DAT":
            # DAT file needs a delimiter to be specified.
            self.command(cmd_name="Export",
                         opts=dict([
                                 ("GuiEnabled", "False"),
                                 ("ClearOptions", "False"),
                                 ("Filter", fileType.lower() ),
                                 ("Delimiter", "comma"),
                                 ("Path", exportPath + "\\" + fileName),
                                 ("PersistOptions", "True"),
                                 ("Module", moduleName),
                                 ("ModuleId", "0"),
                                 ])
                        )
        elif fileType.upper() == "VTK" or fileType.upper() == "VDAT":
            
            if fileType.upper() == "VTK":
                print("\n***WARNING***\n")
                print("vtk export from Voxler does not work properly!")
            
            self.command(cmd_name="Export",
                     opts=dict([
                             ("GuiEnabled", "False"),
                             ("ClearOptions", "False"),
                             ("Filter", fileType.lower()),
                             ("Path", exportPath + "\\" + fileName),
                             ("PersistOptions", "True"),
                             ("Module", moduleName),
                             ("ModuleId", "0"),
                             ])
                        )
        else:
            print("\n***WARNING***\n")
            print("expected dat, vdat or vtk as fileType but got '%s'" 
                  %fileType)
            print("The requested file type has not been trialled and may cause"
                  " an error.")
            self.command(cmd_name="Export",
                     opts=dict([
                             ("GuiEnabled", "False"),
                             ("ClearOptions", "False"),
                             ("Filter", fileType.lower()),
                             ("Path", exportPath + "\\" + fileName),
                             ("PersistOptions", "True"),
                             ("Module", moduleName),
                             ("ModuleId", "0"),
                             ])
                        )
        print("Export file saved to: \n %s" % (exportPath + "\\" + fileName))
        return


##########################################
#----------------------------------------- Import utilities
########################################## 

#========================================= general module input-All file types

    def importModule(self,
                     filePath,
                     fileName,
                     ):
        """
        Imports single modules to random location (Voxler's default location).
        This function can be used along with "moveModule" function to import
        multiple files and put them in the location one needs.

        @Note: fileType is automatically extracted from fileName

        @param filePath: STR, The path of the file.
        @param fileName: STR, The name of the file.
        """
        fileType = fileName.split(".")[-1]
        filePath += "/"
        
        if not os.path.isfile(filePath + fileName):
            print("Input path does not exist! Import was UNSUCCESSFULL!")
            print("***NOTE: if special characters exist (e.g. %s) use r in front of string name." %r"\v")
            print("         for example use r'%s' instead of '%s'." %(r"\voxler",r"\voxler"))
            return
        
        else:
            self.command(cmd_name="Import",
                             opts=dict([
                                     ("AutoConnect", "False"),
                                     ("ClearOptions", "True"),
                                     ("ClearPath", "True"),
                                     ("Filter", fileType),
                                     ("GuiEnabled", "False"),
                                     ("PersistOptions", "False"),
                                     ("ProgressEnabled", "False"),
                                     ("UndoRedoEnabled", "False"),
                                     ("Path", filePath + fileName),
                                     ]),
                            )
        return

#=========================================  import CSV

    def importCsv(self,
                  inputFile,
                  compVsColDict = None,
                  overWriteWarning = False,
                  ):
        """
        Imports a CSV file and returns the module name if successful.
        
        @param filePath: STR. path + file name.
        
        @param compVsColDict: DIC, Optional. If None, default values will be 
            used. column component vs column header count. For exampl:
                sample: x, y, z, col1, col2, col3,...,col50
                compVsColDict = {1 : 4, 2 : 5, ..., 50 : 54}
                Tip: to create such input use {i : i+3 for i in range(1,50)}
                     if a flag column exists use i+4
            WARNING: Voxler takes ages to perform this action as it updates
                     the database after each and every component assignment! 
                     Only for large CSVs. It usually imports and assigns 
                     the correct columns to XYZ and other components. So 
                     hopefully, in most cases, this is not needed.
        
        @param overWriteWarning: Boolean, if False a warning will be shown and
            the user is asked to continue or not. If True no user interaction
            will occur. This is particularly useful if this function needs to
            be used as a part of a larger script which we don't want user 
            interaction.
            
        @NOTE: fileType is automatically extracted from fileName and checked if 
            is CSV or not. If not, an exception will be raised.
        @NOTE: if the len(compVsColDict.keys()) is larger than the number of 
            components - 3 (it is assumed that 3 columns of XYZ exist in the 
            input data) Voxler returns an error. This can be checked before
            executing this function by means of pandas lib. For example:
                >>> import pandas as pd
                >>> df = pd.read_csv(inFile)
                >>> threshold = len(df.columns) - 3
            
        """
        fileType = inputFile.split(".")[-1]
        if fileType.upper() != "CSV":
            raise Exception("Expected 'CSV' but got '%s'" %fileType)
        
        if not os.path.isfile(inputFile):
            raise Exception("Input file '%s' does not exist." %inputFile)
        
        
        self.command(cmd_name="Import",
                         opts=dict([
                                 ("AutoConnect", "False"),
                                 ("ClearOptions", "True"),
                                 ("ClearPath", "True"),
                                 ("Filter", fileType),
                                 ("GuiEnabled", "False"),
                                 ("PersistOptions", "False"),
                                 ("ProgressEnabled", "False"),
                                 ("UndoRedoEnabled", "False"),
                                 ("Path", inputFile),
                                 ]),
                        )
        moduleName = inputFile.split("\\")[-1]
        print("File '%s' imported successfully." %inputFile)
        
        if compVsColDict is not None:
            
            if not overWriteWarning:
                # print warning and ask user to continue or not
                print("\n***** WARNING *****\n"
                      "Changing input components will take a long time"
                      " especially if CSV is a large file.")
                inp = raw_input("Would you like to continue (Y or N)?").upper()
                if inp != "Y":
                    print("Termination request received - not changing Voxler "
                          "default input parameters.")
                    return
            
            # set the number of components
            self.command(cmd_name="ModifyModule",
                         opts=dict([
                                 ("Module", moduleName),
                                 ("PtColCompCount", len(compVsColDict.keys())),
                                 ]),
                        )
            
            # create the input string for modifying imported csv
            importOptions = ["PtColCompCount=" + str(len(compVsColDict.keys()))
                             ] 
            for comp in compVsColDict.keys():
                optName = "ColComponent-" + str(comp)
                optVal = str(compVsColDict[comp])
                importOptions.append(optName + "=" + optVal)
            # executing the input string using command_sequence
            self.command_sequence(
                cmd_name="ModifyModule",
                opts=dict([
                          ("Module", moduleName),
                          ]),
                opts_seq=dict([ tuple(opt.split("=")) for opt in importOptions 
                               ])
                )
            print("Components in module '%s' options have been modified"
                  "successfully." %moduleName)
        return moduleName