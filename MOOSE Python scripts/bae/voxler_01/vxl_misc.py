"""
A more complicated implementation of BEVoxler which makes use of functions
defined in BEVoxler library (inherits from BEVoxler class).

@Warning: This class contains functions which have not been implemented
completely. Proceed with extreme care!
"""

_todo_ = """ 
Test and implement some important and useful functions.

ideas:
======

"""


from .vxl_core import BEVoxler
import os
from Tkinter import Tk
from tkFileDialog import askopenfilenames

class BEVoxler_misc(BEVoxler):
    
    def __init__(self):
        """
        @warning: This class contains functions which have not been implemented
        completely. Proceed with extreme care!
        """
        print("\n**************************************")
        print("---------------WARNING----------------")
        print("YOU ARE USING NOT IMPLEMENTED METHODs")        
        print("**************************************\n")
    
    
#========================================= View manipulation

    def changeView(self,
                   view,
                   ):
        """
        Set the ShowWindow option to Front, Back, Left, Right, Top, Bottom, 
        or Default
        
        @param view: STR, one of the defined views.
        """
        self.command(cmd_name="ViewDirection",
                         opts=dict([
                                 ("ViewDirection", view),
                                 ]),
                        )
        print("view changed to '%s'." % view)
        return
    
    
    def changeViewSize(self,
                   X, Y
                   ):
        """
        Changes the Voxler viewer window size
        
        @param X: INT, pixels in X direction
        @param Y: INT, pixels in Y direction
        """
        self.command(cmd_name="ViewSize",
                         opts=dict([
                                 ("XSize", X),
                                 ("YSize", Y),
                                 ]),
                        )
        print("view size changed to [%d, %d]." % (X,Y) )
        return
    
    def fitToWindow(self,
                   ):
        """
        Fits the view to window to include everything in the viewer window.
        """
        self.command(cmd_name="ViewFittoWindow",
                         opts=dict([
                                 ]),
                        )
        print("View fit to window.")
        return

#========================================= screenshotManual

    def screenshotManual(self, 
            fileName,
            size = [2548,1220],
            exportPath = None,
            multiShots = False,
            ):
        """
        This function can generate "jpg" and "png" output.
        NOTE: fileType is automatically taken from the extension of the fileName input.
              If no extension is given, the function will ask user to correct that.
        
        @param fileName: STRING variable. file name of the screenshot to be saved.
                         EXAMPLE: use a template like "P2017029_R012_F%03d.jpg"
        
        @param size: LIST variable.Default set to [2548,1220] (32" wide screen)
                     Width and Height of the screenshot (bottom right point)
                     The PIL library Image function is used to crop to the defined size.
        
        @param exportPath: STRING variable. default is the same route as the py script
        
        @param multiShots: BOOLEAN variable. default set to False. An extra 
                           functionality to be able to capture multiple screenshots.
                           It adds "_Img%02d" to the end of defined name.
        """
        from PIL import Image
        
        fileType = fileName.split(".")[-1]
        if fileType.upper() != "JPG" and  fileType.upper() != "PNG":
            print("Invalid file extension (fileName = '%s') for screenshot. FileType is set to 'jpg'" %fileType)
            fileType = 'jpg'
            fileName = fileName.split(".")[0] + "." + fileType
        
        # check if size is not defined or defined wrongly
        try:
            size
        except NameError:
            print("JPG size is not defined, default size is chosen as [w = 2548, h = 1220]")
            size = [2548,1220]
        if len(size) < 2:
            print("JPG size is not defined properly, default size is chosen as [w = 2548, h = 1220]")
            size = [2548,1220]
        
        if exportPath is None:
           exportPath = os.path.abspath(os.getcwd()).replace(os.sep,"\\")
           print("exportPath input is None. exportPath is set to: %s" % exportPath)
        
        if not os.path.isdir(exportPath):
            print("Export path does not exist! Export was UNSUCCESSFULL!")
            print("***NOTE: if special characters exist (e.g. %s) use r in front of string name." %r"\v")
            print("         for example use r'%s' instead of '%s'." %(r"\voxler",r"\voxler"))
            return
        
        if not multiShots:
            self.command(cmd_name="Export",
                             opts=dict([
                                     ("GuiEnabled", "False"),
                                     ("ClearOptions", "True"),
                                     ("Filter", fileType),
                                     #("Height", size[0]),
                                     #("Width", size[1]),
                                     ("KeepAspect","0"),
                                     ("Path",exportPath + "\\" + fileName),
                                     ("PersistOptions", "True"),
                                     ("ModuleID","0"),
                                     ]),
                            )
            img = Image.open(exportPath + "\\" + fileName)
            crp_img = img.crop((0, 0, size[0], size[1]))
            print("Image dimension is set to %s." % str(size))
            crp_img.save(exportPath + "\\" + fileName)
            print("Screenshot saved to: \n%s" % (exportPath + "\\" + fileName))
        
        if multiShots:
            cnt = 1
            while True:
                key = raw_input("Choose a suitable view angle and enter 'C' to screenshot or any other keys to return: ").upper()
                if key == "C":
                    nameExt = "_Img%02d." %cnt
                    self.command(cmd_name="Export",
                                     opts=dict([
                                             ("GuiEnabled", "False"),
                                             ("ClearOptions", "False"),
                                             ("Filter", fileType),
                                             #("Height", size[0]),
                                             #("Width", size[1]),
                                             ("Path",exportPath + "\\" + fileName.split(".")[0] + nameExt + fileName.split(".")[-1]),
                                             ("PersistOptions", "True"),
                                             ("ModuleID","0"),
                                             ]),
                            )
                    print("Image dimension is set to %s." % str(size))
                    img = Image.open(exportPath + "\\" + fileName.split(".")[0] + nameExt + fileName.split(".")[-1])
                    crp_img = img.crop((0, 0, size[0], size[1]))
                    crp_img.save(exportPath + "\\" + fileName.split(".")[0] + nameExt + fileName.split(".")[-1])
                    print("Screenshot saved to: \n%s" % (exportPath + "\\" + fileName.split(".")[0] + nameExt + fileName.split(".")[-1]))
                    cnt += 1
                elif key != "C":
                    break
                else:
                    print("Uknown input!")                
        return


    
##########################################
#----------------------------------------- Import utilities
########################################## 



#=========================================  importIVFiles

    def importIVFiles(self, 
            listOfFiles,
            modulePosFirst = None,
            modulePosOffset = (0,33),
            moduleShow = True,
            ):
        """
        @param listOfFiles: list of IV filenames (fullpath). An Exception is
          raised if a filename within the list does not exist.
          
        @param modulePosFirst: optional argument, defaults to None. 
          If None, then automatic placement of the module within the 
          NetworkManager. Can be a tuple of NetworkManager coords (x,y).

        @param modulePosOffset: optional argument, defaults to (0,25).
          gives the offset between two imported modules within the 
          NetworkManager. Only used if modulePosFirst is not None.
          
        @param moduleShow: optional argument, defaults to True. Can be either
          a boolean value, or a list of boolean values.
        """
        # function to return the name of imported modules
        listOfModules = list()
        
        if moduleShow is None:
            moduleShow = True
        if isinstance(moduleShow, bool):
            flag = bool(moduleShow)
            moduleShow = [ flag for i in range(len(listOfFiles)) ]
        elif isinstance(moduleShow, list):
            assert(len(moduleShow)==len(listOfFiles))
            assert(all([ isinstance(flag, bool) for flag in moduleShow ])==True)
            
        for i,(fn,show) in enumerate(zip(listOfFiles,moduleShow)):
    
            ### make sure fn exists to avoid non-meaningful ApiError exception
            if not os.path.exists(fn):
                raise Exception("try to import non-existing file '%s'. FAIL!"%fn)
                
            fn = fn.replace(os.sep,"/")
            moduleName = fn.rsplit("/",1)[-1]
            ### AF comment: it would be very nice if the voxler_lib.Voxler.command
            ### would return the assigned moduleName or moduleId
            ### Here, we can only hope the moduleName remains to be the filename
            self.command(
            #self.command(
                cmd_name="Import",
                opts=dict([
                    ("AutoConnect", "False"),
                    ("ClearOptions", "True"),
                    ("ClearPath", "True"),
                    ("Filter", "iv"),
                    ("GuiEnabled", "False"),
                    ("PersistOptions", "False"),
                    ("ProgressEnabled", "False"),
                    ("UndoRedoEnabled", "False"),
                    ("Path", str(fn)),
                    #("DefaultPosition", "False"),
                    #("XPosition", "%d"%(5)),
                    #("YPosition", "%d"%(70+25*i)),
                    ]),
                )
            if modulePosFirst is None:
                modulePosFirst = (0,33)
            
            self.command(
                cmd_name="MoveModule",
                opts=dict([
                    ("Module", moduleName),
                    ("UseOriginalPosition", "False"),
                    ("XPosition", str(modulePosFirst[0]+i*modulePosOffset[0])),
                    ("YPosition", str(modulePosFirst[1]+i*modulePosOffset[1])),
                    ]),
                )
            if show==False:
                #self.command(
                self.command(
                    cmd_name="ShowModule",
                    opts=dict([
                        ("Module", moduleName),
                        ("Show", "False"),
                        ]),
                    )
            print("IV file '%s' imported successfully." % moduleName)
            listOfModules.append(moduleName)
        return(listOfModules)

#=========================================  importIVFilesGUI
    def importIVFilesGUI(self, 
                         modulePosFirst = None,
                         modulePosOffset = (0,33),
                         moduleShow = True,
                         ):
        """
        DESCRIPTION
        """
        Tk().withdraw()
        listOfFiles = list(askopenfilenames(title="Choose IV files (press 'Ctrl' to choose multiple files)") )
        listOfFiles.sort()
        
        listOfModules = self.importIVFiles(listOfFiles,
                                           modulePosFirst,
                                           modulePosOffset,
                                           moduleShow)
        return(listOfFiles, listOfModules)
    
#=========================================  importVTKFiles
    
    def importVTKFiles(self, 
            listOfFiles,
            modulePosFirst = (0,38),
            modulePosOffset = (0,33),
            modulePosAuto = True,
            #moduleShow = True,
            ):
        """
        @param listOfFiles: LIST of VTK filenames (fullpath). An Exception is
          raised if a filename within the list does not exist.
          
        @param modulePosFirst: optional argument, defaults to (0,38). 
          If None, then (0,38). Can be a tuple of NetworkManager coords (x,y).
        
        @param modulePosOffset: optional argument, defaults to (0,33).
          gives the offset between two imported modules within the 
          NetworkManager. Only used if modulePosFirst is not None.
        
        @param modulePosAuto = optional argument, defaults to None. 
          If False then it follows the auto positioning defined by the previous 
          parameter, i.e. "modulePosOffset".
          True means to position each vtk module depending on its frame
          number, i.e. frame number 30 will sit in front of of frame 30 of IVs.

        """
        

        if modulePosAuto is None:
           modulePosAuto = True 
        
        if modulePosFirst is None:
            modulePosFirst = (0,38)
        #if moduleShow is None:
        #    moduleShow = True
        #if isinstance(moduleShow, bool):
        #    flag = bool(moduleShow)
        #    moduleShow = [ flag for i in range(len(listOfFiles)) ]
        #elif isinstance(moduleShow, list):
        #    assert(len(moduleShow)==len(listOfFiles))
        #    assert(all([ isinstance(flag, bool) for flag in moduleShow ])==True)
            
        for i,fn in enumerate(listOfFiles):
    
            ### make sure fn exists to avoid non-meaningful ApiError exception
            if not os.path.exists(fn):
                raise Exception("try to import non-existing file '%s'. FAIL!"%fn)
                
            fn = fn.replace(os.sep,"/")
            moduleName = fn.rsplit("/",1)[-1]
            #print(moduleName)
            self.command(cmd_name="Import",
                             opts=dict([
                                     ("AutoConnect", "False"),
                                     ("ClearOptions", "True"),
                                     ("ClearPath", "True"),
                                     ("Filter", "vtk"),
                                     ("GuiEnabled", "False"),
                                     ("PersistOptions", "False"),
                                     ("ProgressEnabled", "False"),
                                     ("UndoRedoEnabled", "False"),
                                     ("Path", str(fn)),
                                     #("DefaultPosition", "False"),
                                     #("XPosition", "%d"%(5)),
                                     #("YPosition", "%d"%(70+25*i)),
                                     ]),
                            )
            print("VTK file '%s' imported successfully." % moduleName)
    
            if modulePosAuto is False:
                self.command(cmd_name="MoveModule",
                                 opts=dict([
                                         ("Module", moduleName),
                                         ("UseOriginalPosition", "False"),
                                         ("XPosition", str(modulePosFirst[0]+i*modulePosOffset[0])),                            
                                         ("YPosition", str(modulePosFirst[1]+i*modulePosOffset[1])),
                                         ]),
                                )
            elif modulePosAuto is True:
                frame = int(moduleName [-7:-4])
                self.command(cmd_name="MoveModule",
                                 opts=dict([
                                         ("Module", moduleName),
                                         ("UseOriginalPosition", "False"),
                                         ("XPosition", str(modulePosFirst[0]+i*modulePosOffset[0])),
                                         ("YPosition", str(modulePosFirst[1]+(frame-1)*modulePosOffset[1])),
                                         ]),
                                )
                               
        return
    
##########################################
#----------------------------------------- preview IV
########################################## 

    def previewIV(self,
                  frameNumbers,
                  ivFiletmpl,
                  frameLabel = 'frame_label',
                  drawStyleToLines = True,
                  #screenCapture = False,
                  ):
        """
        @param frameNumbers: LIST of frame numbers.
        
        @param ivFiletmpl: STRING, e.g. "FRA2017_G02_Q03_ALL_F%03d.iv"

        @param frameLabel: STRING, optional parameter set to "frame_label" as
                           default. It is the name of the Annotation module.
                           Annotation text is derived from IV name automatically. 
        
        @param drawStyleToLines: BOOLEAN, True to set to lines and False to set to AsIs.
                                 This is only for demonstration purpose.
        
        
        """
        # asking the user to make sure IVs are off
        raw_input("**** IMPORTANT **** \n Make sure all the IVs are turned OFF, then press ENTER to continue...")
            
        
        #turn first IVs on
        cnt = 0
        frame0 = frameNumbers[cnt]
        frame1 = frameNumbers[cnt+1]
        print('Frame F%03d is being shown' % frame1)
        # turn IV_0 on
        self.showModule(ivFiletmpl %frame0, True)
        if drawStyleToLines:
            self.changeDrawStyle(ivFiletmpl %frame0, 0)
        # turn IV_1 on
        self.showModule(ivFiletmpl %frame1, True)
        # change IV_1 to lines
        if drawStyleToLines:
            self.changeDrawStyle(ivFiletmpl %frame1, 2)
        # change the frame_label
        self.changeAnnotationText(frameLabel, "_".join(ivFiletmpl.split('_')[0:3]) + '_F%03d' %frame1)
        
        while True:
            
            # ask user for an input
            print("------------------------")
            key = raw_input("press F for Forward, B for Backward, T for Terminate, C to screen Capture the current frame (current frame is F%03d): " %frame1).upper()
            
            if key == "C":
                
                raw_input("Select a suitable view point and press ENTER to screen capture.")
                self.screenshot("_".join(ivFiletmpl.split('_')[0:3]) + '_F%03d.jpg' %frame1)
            
            
            elif key == "F":
                                
                cnt += 1
                # the following condition is to prevent out of list errors
                if cnt < len(frameNumbers)-1:
                    print('Frame F%03d is being shown' % (frame1+1))
                    # turn previous IVs off
                    self.showModule(ivFiletmpl %frame0, False)
                    if drawStyleToLines:
                        self.changeDrawStyle(ivFiletmpl %frame0, 0)
                    self.showModule(ivFiletmpl %frame1, False)
                    if drawStyleToLines:
                        self.changeDrawStyle(ivFiletmpl %frame1, 0)
                    
                    frame0 = frameNumbers[cnt]
                    frame1 = frameNumbers[cnt+1]
                    # turn IV_0 on
                    self.showModule(ivFiletmpl %frame0, True)
                    # turn IV_1 on
                    self.showModule(ivFiletmpl %frame1, True)
                    # change IV_1 to lines
                    if drawStyleToLines:
                        self.changeDrawStyle(ivFiletmpl %frame1, 2)
                    # change the frame_label
                    self.changeAnnotationText(frameLabel, "_".join(ivFiletmpl.split('_')[0:3]) + '_F%03d' %frame1)
                    
                    
                else:
                    print("To go forward F%03d needs to be shown, which is out of the defined range. " % (frameNumbers[cnt]+1))
                    print("The defined range is from F%03d to F%03d." %(frameNumbers[0],frameNumbers[-1]))
                    cnt -= 1
                    
            elif key == "B":
                                
                cnt -= 1
                # the following condition is to prevent out of list errors
                if cnt >= 0:
                    print('Frame F%03d is being shown' % frame0)
                    # turn previous IVs off
                    self.showModule(ivFiletmpl %frame0, False)
                    if drawStyleToLines:
                        self.changeDrawStyle(ivFiletmpl %frame0, 0)
                    self.showModule(ivFiletmpl %frame1, False)
                    if drawStyleToLines:
                        self.changeDrawStyle(ivFiletmpl %frame1, 0)
                    
                    frame0 = frameNumbers[cnt]
                    frame1 = frameNumbers[cnt+1]
                    # turn IV_0 on
                    self.showModule(ivFiletmpl %frame0, True)
                    # turn IV_1 on
                    self.showModule(ivFiletmpl %frame1, True)
                    # change IV_1 to lines
                    if drawStyleToLines:
                        self.changeDrawStyle(ivFiletmpl %frame1, 2)
                    # change the frame_label
                    self.changeAnnotationText(frameLabel, "_".join(ivFiletmpl.split('_')[0:3]) + '_F%03d' %frame1)
                    
                else:
                    print("To go back F%03d needs to be shown, which is out of the defined range. " % (frameNumbers[cnt+1]-1))
                    print("The defined range is from F%03d to F%03d." %(frameNumbers[0],frameNumbers[-1]))
                    cnt += 1
        
            elif key == "T":
                print("Termination request has been received. Terminating...")
                self.showModule(ivFiletmpl %frame0, False)
                if drawStyleToLines:
                    self.changeDrawStyle(ivFiletmpl %frame0, 0)
                self.showModule(ivFiletmpl %frame1, False)
                if drawStyleToLines:
                    self.changeDrawStyle(ivFiletmpl %frame1, 0)
                print("\npreviewIV has been terminated successfully.")
                return
            else:
                print("Unknown input!")
                
        return
    
##########################################
#----------------------------------------- PgUp PgDn functionality
########################################## 

#========================================= switchFrame

    def switchFrame(self,
                   frameNumbers,
                   ivFiletmpl,
                   modulePairs,
                   labelOpt = ("frame_label","F%03d"),
                   screenCapture = False,
                   ):
        """
        @param frameNumbers: LIST, list of frames to go through
        
        @param ivFiletmpl: STRING, e.g. "FRA2017_G02_Q03_ALL_F%03d.iv"

        @param modulePairs: DICTIONARY, pairs of module names template and destination
                            modulePairs ={"MgB1_F%03d_F%03d"      : "JNC_RMRER_B01",
                                          "Merge_F%03d_F%03d" : "JNC_faultslip", 
                                          .
                                          .
                                          .
                                          }
                            NOTE: modulePairs have to follow same frame numbering.
                            There are often multiple VTK boxes with similar naming.
        
        @param labelOpt: TUPLE of two STRINGS, optional parameter set to 
                         ('frame_label','RER_F%03d-F%03d') as default. 
                         The two STRINGS are the Annotation module name and annotation text template.
        
        @param screencapture: BOOLEAN, default to False. True to capture screen.
                              Note that user will be asked to enter other screencapture function parameters.
        
        """
                # asking the user to make sure IVs are off
        raw_input("**** IMPORTANT **** \nMake sure all the IVs are turned OFF, then press ENTER to continue...")
        
        # asking the user for screenshot options if screenCapture is set to True
        if screenCapture:
            ImgPath = raw_input("Enter a path to save screenshots or press N to choose the current folder: ")
            if ImgPath == "N":
                ImgPath = None
            MultiShots = raw_input("Would you like to take multiple screenshots for each frame? (Y or N) ").upper()
            if MultiShots == "Y":
                MultiShots = True
            else:
                MultiShots = False
        
        print("================================================================================================")
        print("INPUT PARAMTERS ARE:")
        print("  frameNumbers : %s" %str(frameNumbers))
        #print("  frameNumberOffset : %d" %frameNumberOffset)
        if screenCapture:
            print("Screenshot options are:")
            print("  Image name template = %s" %labelOpt[1])
            print("  Path = %s" %str(ImgPath))
            print("  Multiple screenshots = %s" %MultiShots)
        print("================================================================================================\n")
        
        #turn first IVs on
        cnt = 0
        frame = frameNumbers[cnt]

        print("Frame F%03d is being shown..." % (frameNumbers[cnt]))
        
        print("\nTurning ON current IV frame...")
        self.showModule(ivFiletmpl % frame, True)
        self.changeDrawStyle(ivFiletmpl % frame, 0)
        
        print("\nConnecting modules...")
        # connecting module pairs
        modulePairs_temp = dict()
        for m in modulePairs.keys():
            modulePairs_temp.update({m  %frame : modulePairs[m]})
        self.connectMultiModules(modulePairs_temp)
        
        print #skipping a line in IPython console
        self.changeAnnotationText(labelOpt[0], labelOpt[1] %frame)
                
        # screenshot here if needed
        if screenCapture:
            print("\nCapturing the screen...")
            self.screenshot(labelOpt[1] %frame, ImgPath, MultiShots)
        
        while True:
            
            print("------------------------")
            # ask user for an input
            key = raw_input("Enter F for forward, B for backward, S to skip to a defined frame, T to terminate (current frame is F%03d): " %frame).upper()
            if key == "F":
                                
                cnt += 1
                # the following condition is to prevent out of list errors
                if cnt <= len(frameNumbers)-1:
                    print("Frame F%03d is being shown" %frameNumbers[cnt])
                    
                    frame0 = frame # this is to turn current IV on before turning previous IV off
                    frame = frameNumbers[cnt]
                    # turn current IVs on
                    print("\nTurning ON current IV frame...")
                    self.showModule(ivFiletmpl % frame, True)
                    self.changeDrawStyle(ivFiletmpl % frame, 0)
                    # turn previous IVs off
                    print("\nTurning OFF previous IV frame (F%03d)." %frame0)
                    self.showModule(ivFiletmpl %frame0, False)
                    self.changeDrawStyle(ivFiletmpl %frame0, 0)
                    
                    
                    print("\nConnecting modules...")
                    # connecting module pairs
                    modulePairs_temp = dict()
                    for m in modulePairs.keys():
                        modulePairs_temp.update({m  %frame : modulePairs[m]})
                    self.connectMultiModules(modulePairs_temp)
                    
                    print #skipping a line in IPython console
                    self.changeAnnotationText(labelOpt[0], labelOpt[1] %frame)
                    
                    # screenshot here if needed
                    if screenCapture:
                        print("\nCapturing the screen...")
                        self.screenshot(labelOpt[1] %frame, ImgPath, MultiShots)
                    
                    
                else:
                    print("Last frame reached (F%03d). Go backward, skip or terminate by pressing B, S or T " % (frameNumbers[cnt-1]))
                    print("The requested frameNumbers are %s." %str(frameNumbers))
                    cnt -= 1
            
            
            elif key == "B":
                                
                cnt -= 1
                # the following condition is to prevent out of list errors
                if cnt >= 0:
                    print("Frame F%03d is being shown" %frameNumbers[cnt])
                    
                    # turn previous IVs off
                    print("\nTurning OFF previous IV frame (F%03d)." %frame)
                    self.showModule(ivFiletmpl %frame, False)
                    self.changeDrawStyle(ivFiletmpl %frame, 0)
                    
                    frame = frameNumbers[cnt]
                    
                    print("\nTurning ON current IV frame...")
                    self.showModule(ivFiletmpl % frame, True)
                    self.changeDrawStyle(ivFiletmpl % frame, 0)
                    
                    print("\nConnecting modules...")
                    # connecting module pairs
                    modulePairs_temp = dict()
                    for m in modulePairs.keys():
                        modulePairs_temp.update({m  %frame : modulePairs[m]})
                    self.connectMultiModules(modulePairs_temp)
                    
                    print #skipping a line in IPython console
                    self.changeAnnotationText(labelOpt[0], labelOpt[1] %frame)
                    
                    # screenshot here if needed
                    if screenCapture:
                        print("\nCapturing the screen...")
                        self.screenshot(labelOpt[1] %frame, ImgPath, MultiShots)
                        
                else:
                    print("First frame reached (F%03d). Go forward, skip or terminate by pressing F, S or T " % (frameNumbers[cnt+1]))
                    print("The requested frameNumbers are %s." %str(frameNumbers))
                    cnt += 1
            
            
            elif key == 'S':
                
                
                # check if the new entry exists in the defined list
                while True:
                    try:
                        frameNew = input('Choose a frame number (format = integer) or press 0 to terminate: ')
                        if frameNew == 0:
                            print('Termination request has been received.')
                            return
                        if frameNew not in frameNumbers:
                            print("*** ERROR *** \nFrame F%03d is not in the input list." %frameNew)
                            print("Defined frameNumbers are: %s" %str(frameNumbers))
                        else:
                            break
                    except (NameError,SyntaxError):
                        print("Make sure to enter an integer as the frame number.")
                    
                    #except SyntaxError:
                    #    print("Make sure to enter an integer as the frame number.")
                        
                
                # assigning the index of new entry to cnt
                cnt = frameNumbers.index(frameNew)
                
                print("Frame F%03d is being shown" %frameNumbers[cnt])
                    
                # turn previous IVs off
                print("\nTurning OFF previous IV frame (F%03d)." %frame)
                self.showModule(ivFiletmpl %frame, False)
                self.changeDrawStyle(ivFiletmpl %frame, 0)
                
                frame = frameNumbers[cnt]
                
                print("\nTurning ON current IV frame...")
                self.showModule(ivFiletmpl % frame, True)
                self.changeDrawStyle(ivFiletmpl % frame, 0)
                
                print("\nConnecting modules...")
                # connecting module pairs
                modulePairs_temp = dict()
                for m in modulePairs.keys():
                    modulePairs_temp.update({m  %frame : modulePairs[m]})
                self.connectMultiModules(modulePairs_temp)
                
                print #skipping a line in IPython console
                self.changeAnnotationText(labelOpt[0], labelOpt[1] %frame)
                
                # screenshot here if needed
                if screenCapture:
                    print("\nCapturing the screen...")
                    self.screenshot(labelOpt[1] %frame, ImgPath, MultiShots)
                
                
            # manual termination (BG: do we need to disconnect modules?)
            elif key == 'T':
                print('Termination request has been received.')
                return
            
            # wrong keyword
            elif key != 'F' and key != 'B' and key != 'S' and key != 'T':
                print('Unknown keyword! Enter another keyword.')
                
            # automatic termination
            else:
                print('Automatic termination option is being executed!')
                return
        return
        
        

##########################################
#----------------------------------------- RER TOOLS
########################################## 

#========================================= switchFrameRER

    def switchFrameRER(self,
                       frameNumbers,
                       frameNumberOffset,
                       ivFiletmpl,
                       modulePairs,
                       isoVals,
                       labelOpt = ("frame_label","RER_F%03d-F%03d"),
                       screenCapture = False
                       ):
        """
        NOTE: in RER investigation, IV for (frame0-1) should be displayed as shaded
              and IV for frame1 should be displayed as lines.
              
        @param frameNumbers: LIST variable, a list of THE FIRST frame numbers for a desired period.
        
        @param frameNumberOffset: INT variable, the offset number for a desired period.
        
        @param ivFiletmpl: STRING variable, e.g. "FRA2017_G02_Q03_ALL_F%03d.iv"
        
        @param modulePairs: DICTIONARY variable, pairs of module names template and destination
                            modulePairs ={"MgB1_F%03d_F%03d"      : "JNC_RMRER_B01",
                                          "Merge_F%03d_F%03d" : "JNC_faultslip", 
                                          .
                                          .
                                          .
                                          }
                            NOTE: modulePairs have to follow same frame numbering.
                            
        @param isovals: DICTIONARY, pairs of isosurface module name and value
                        rm_vals = {"LOW_15E3_D01" : 15000,
                                   "MED_250E3_D01" : 150000,
                                   "HIGH_500E3_D01" : 750000,
                                   .
                                   .
                                   .
                                   }
        
        @param labelOpt: TUPLE of two STRINGS, optional parameter set to 
                         ('frame_label','RER_F%03d-F%03d') as default. 
                         The two parameters are the Annotation module name and annotation text template.
        
        @param screencapture: BOOLEAN, default to False. True to capture screen.
                              Note that user will be asked to enter other screencapture function parameters.
        
        """
        # asking the user to make sure IVs are off
        raw_input("**** IMPORTANT **** \nMake sure all the IVs are turned OFF, then press ENTER to continue...")
        
        # asking the user for screenshot options if screenCapture is set to True
        if screenCapture:
            ImgPath = raw_input("Enter a path to save screenshots or press N to choose the current folder: ").upper()
            if ImgPath == "N":
                ImgPath = None
            MultiShots = raw_input("Would you like to take multiple screenshots for each frame? (Y or N) ").upper()
            if MultiShots == "Y":
                MultiShots = True
            else:
                MultiShots = False
        
        
        print("================================================================================================")
        print("INPUT PARAMTERS ARE:")
        print("  frameNumbers : %s" %str(frameNumbers))
        print("  frameNumberOffset : %d" %frameNumberOffset)
        if screenCapture:
            print("Screenshot options are:")
            print("  Image name template = %s" %labelOpt[1])
            print("  Path = %s" %str(ImgPath))
            print("  Multiple screenshots = %s" %MultiShots)
        print("================================================================================================\n")
        
        #turn first IVs on
        cnt = 0
        frame0 = frameNumbers[cnt]
        frame1 = frame0 + frameNumberOffset
        print("Frame F%03d-F%03d is being shown..." % (frameNumbers[cnt], frameNumbers[cnt] + frameNumberOffset))
        
        print("\nTurning on IV frames F%03d and F%03d." %((frame0-1),frame1))
        self.showModule(ivFiletmpl % (frame0-1), True)
        self.changeDrawStyle(ivFiletmpl % (frame0-1), 0)
        self.showModule(ivFiletmpl % frame1, True)
        self.changeDrawStyle(ivFiletmpl % frame1, 2)
        
        print("\nConnecting modules...")
        # connecting module pairs
        modulePairs_temp = dict()
        for m in modulePairs.keys():
            modulePairs_temp.update({m  % (frame0,frame1):modulePairs[m]})
        self.connectMultiModules(modulePairs_temp)
        
        print #skipping a line in IPython console
        self.changeAnnotationText(labelOpt[0], labelOpt[1] %(frame0,frame1))
        
        print("\nRe-setting isovalues...")
        self.changeIsosurfaceValue2(isoVals)
        
        # screenshot here if needed
        if screenCapture:
            print("\nCapturing the screen...")
            self.screenshot(labelOpt[1] %(frame0,frame1), ImgPath, MultiShots)
                    
        while True:
            
            # ask user for an input
            print("------------------------")
            key = raw_input("press F for forward, B for Backward, S to skip to a frame, T for terminate (current frame is F%03d-F%03d): " %(frame0,frame1)).upper()
            if key == "F":
                                
                cnt += 1
                # the following condition is to prevent out of list errors
                if cnt <= len(frameNumbers)-1:
                    print("Frame F%03d-F%03d is being shown" % (frameNumbers[cnt], frameNumbers[cnt] + frameNumberOffset))
                    
                    # turn previous IVs off
                    print("\nTurning OFF previous IV frames F%03d and F%03d." %((frame0-1),frame1))
                    self.showModule(ivFiletmpl % (frame0-1), False)
                    self.changeDrawStyle(ivFiletmpl % (frame0-1), 0)
                    self.showModule(ivFiletmpl % frame1, False)
                    self.changeDrawStyle(ivFiletmpl % frame1, 0)
                    
                    frame0 = frameNumbers[cnt]
                    frame1 = frame0 + frameNumberOffset
                    
                    print("\nTurning on IV frames F%03d and F%03d." %((frame0-1),frame1))
                    self.showModule(ivFiletmpl % (frame0-1), True)
                    self.changeDrawStyle(ivFiletmpl % (frame0-1), 0)
                    self.showModule(ivFiletmpl % frame1, True)
                    self.changeDrawStyle(ivFiletmpl % frame1, 2)
                    
                    print("\nConnecting modules...")
                    # connecting module pairs
                    modulePairs_temp = dict()
                    for m in modulePairs.keys():
                        modulePairs_temp.update({m  % (frame0,frame1):modulePairs[m]})
                    self.connectMultiModules(modulePairs_temp)
                    
                    print #skipping a line in IPython console
                    self.changeAnnotationText(labelOpt[0], labelOpt[1] %(frame0,frame1))
                    
                    print("\nRe-setting isovalues...")
                    self.changeIsosurfaceValue2(isoVals)
                    
                    # screenshot here if needed
                    if screenCapture:
                        print("\nCapturing the screen...")
                        self.screenshot(labelOpt[1] %(frame0,frame1), ImgPath, MultiShots)
                    
                    
                else:
                    print("Last frame reached (F%03d). Go backward or terminate by pressing B or T " % (frameNumbers[cnt-1]))
                    print("The requested frameNumbers are %s." %str(frameNumbers))
                    cnt -= 1
                    
                    
            elif key == "B":
                                
                cnt -= 1
                # the following condition is to prevent out of list errors
                if cnt >= 0:
                    print("Frame F%03d-F%03d is being shown" % (frameNumbers[cnt], frameNumbers[cnt] + frameNumberOffset))
                    
                    # turn previous IVs off
                    print("\nTurning OFF previous IV frames F%03d and F%03d." %((frame0-1),frame1))
                    self.showModule(ivFiletmpl % (frame0-1), False)
                    self.changeDrawStyle(ivFiletmpl % (frame0-1), 0)
                    self.showModule(ivFiletmpl % frame1, False)
                    self.changeDrawStyle(ivFiletmpl % frame1, 0)
                    
                    frame0 = frameNumbers[cnt]
                    frame1 = frame0 + frameNumberOffset
                    
                    print("\nTurning on IV frames F%03d and F%03d." %((frame0-1),frame1))
                    self.showModule(ivFiletmpl % (frame0-1), True)
                    self.changeDrawStyle(ivFiletmpl % (frame0-1), 0)
                    self.showModule(ivFiletmpl % frame1, True)
                    self.changeDrawStyle(ivFiletmpl % frame1, 2)
                    
                    print("\nConnecting modules...")
                    # connecting module pairs
                    modulePairs_temp = dict()
                    for m in modulePairs.keys():
                        modulePairs_temp.update({m  % (frame0,frame1):modulePairs[m]})
                    self.connectMultiModules(modulePairs_temp)
                    
                    print #skipping a line in IPython console
                    self.changeAnnotationText(labelOpt[0], labelOpt[1] %(frame0,frame1))
                    
                    print("\nRe-setting isovalues...")
                    self.changeIsosurfaceValue2(isoVals)
                    
                    # screenshot here if needed
                    if screenCapture:
                        print("\nCapturing the screen...")
                        self.screenshot(labelOpt[1] %(frame0,frame1), ImgPath, MultiShots)
                        
                else:
                    print("First frame reached (F%03d). Go forward or terminate by pressing F or T " % (frameNumbers[cnt+1]))
                    print("The requested frameNumbers are %s." %str(frameNumbers))
                    cnt += 1
                    
                    
            elif key == "S":
                
                # check if the new entry exists in the defined list
                while True:
                    try:
                        frameNew = input('Choose a *FIRST* frame number (format = integer) or press 0 to terminate: ')
                        if frameNew == 0:
                            print('Termination request has been received.')
                            return
                        if frameNew not in frameNumbers:
                            print("*** ERROR *** \nFrame F%03d is not in the input list." %frameNew)
                            print("Defined frameNumbers are: %s" %str(frameNumbers))
                        else:
                            break
                    except (NameError,SyntaxError):
                        print("Make sure to enter an integer as the frame number.")
                
                # user has to input the start frame number of the pairs
                cnt = frameNumbers.index(frameNew)
                
                print("Frame F%03d-F%03d is being shown" % (frameNumbers[cnt], frameNumbers[cnt] + frameNumberOffset))
                
                # turn previous IVs off
                print("\nTurning OFF previous IV frames F%03d and F%03d." %((frame0-1),frame1))
                self.showModule(ivFiletmpl % (frame0-1), False)
                self.changeDrawStyle(ivFiletmpl % (frame0-1), 0)
                self.showModule(ivFiletmpl % frame1, False)
                self.changeDrawStyle(ivFiletmpl % frame1, 0)
                
                frame0 = frameNumbers[cnt]
                frame1 = frame0 + frameNumberOffset
                
                print("\nTurning on IV frames F%03d and F%03d." %((frame0-1),frame1))
                self.showModule(ivFiletmpl % (frame0-1), True)
                self.changeDrawStyle(ivFiletmpl % (frame0-1), 0)
                self.showModule(ivFiletmpl % frame1, True)
                self.changeDrawStyle(ivFiletmpl % frame1, 2)
                
                print("\nConnecting modules...")
                # connecting module pairs
                modulePairs_temp = dict()
                for m in modulePairs.keys():
                    modulePairs_temp.update({m  % (frame0,frame1):modulePairs[m]})
                self.connectMultiModules(modulePairs_temp)
                
                print #skipping a line in IPython console
                self.changeAnnotationText(labelOpt[0], labelOpt[1] %(frame0,frame1))
                
                print("\nRe-setting isovalues...")
                self.changeIsosurfaceValue2(isoVals)
                
                # screenshot here if needed
                if screenCapture:
                    print("\nCapturing the screen...")
                    self.screenshot(labelOpt[1] %(frame0,frame1), ImgPath, MultiShots)
                    
                    
            elif key == "T":
                print("Termination request has been received.")
                return
            else:
                print("Unknown input!")
                
        return


#========================================= switchFrameRER

    def switchFrameRERLoop(self,
                           frameNumbers,
                           frameNumberOffset,
                           ivFiletmpl,
                           modulePairs,
                           isoVals,
                           labelOpt = ("frame_label","RER_F%03d-F%03d"),
                           screenCapture = True,
                           ):
        """
        This function is essentially very similar to switchFrameRER but it 
        doesn'take input from the user and loop over all the frames to take
        screenshots
        """
        # asking the user to make sure IVs are off
        raw_input("**** IMPORTANT **** \nMake sure all the IVs are turned OFF, then press ENTER to continue...")
        
        # asking the user for screenshot options if screenCapture is set to True
        if screenCapture:
            ImgPath = raw_input("Enter a path to save screenshots or press N to choose the current folder: ").upper()
            if ImgPath == "N":
                ImgPath = None
            MultiShots = raw_input("Would you like to take multiple screenshots for each frame? (Y or N) ").upper()
            if MultiShots == "Y":
                MultiShots = True
            else:
                MultiShots = False
        
        
        print("================================================================================================")
        print("INPUT PARAMTERS ARE:")
        print("  frameNumbers : %s" %str(frameNumbers))
        print("  frameNumberOffset : %d" %frameNumberOffset)
        if screenCapture:
            print("Screenshot options are:")
            print("  Image name template = %s" %labelOpt[1])
            print("  Path = %s" %str(ImgPath))
            print("  Multiple screenshots = %s" %MultiShots)
        print("================================================================================================\n")
        
        #turn first IVs on
        cnt = 0
        frame0 = frameNumbers[cnt]
        frame1 = frame0 + frameNumberOffset
                    
        while True:
            
            print("------------------------")
            
            # the following condition is to prevent out of list errors
            if cnt <= len(frameNumbers)-1:
                print("Frame F%03d-F%03d is being shown" % (frameNumbers[cnt], frameNumbers[cnt] + frameNumberOffset))
                
                if cnt != 0:
                    # turn previous IVs off
                    print("\nTurning OFF previous IV frames F%03d and F%03d." %((frame0-1),frame1))
                    self.showModule(ivFiletmpl % (frame0-1), False)
                    self.changeDrawStyle(ivFiletmpl % (frame0-1), 0)
                    self.showModule(ivFiletmpl % frame1, False)
                    self.changeDrawStyle(ivFiletmpl % frame1, 0)
                
                frame0 = frameNumbers[cnt]
                frame1 = frame0 + frameNumberOffset
                
                print("\nTurning on IV frames F%03d and F%03d." %((frame0-1),frame1))
                self.showModule(ivFiletmpl % (frame0-1), True)
                self.changeDrawStyle(ivFiletmpl % (frame0-1), 0)
                self.showModule(ivFiletmpl % frame1, True)
                self.changeDrawStyle(ivFiletmpl % frame1, 2)
                
                print("\nConnecting modules...")
                # connecting module pairs
                modulePairs_temp = dict()
                for m in modulePairs.keys():
                    modulePairs_temp.update({m  % (frame0,frame1):modulePairs[m]})
                self.connectMultiModules(modulePairs_temp)
                
                print #skipping a line in IPython console
                self.changeAnnotationText(labelOpt[0], labelOpt[1] %(frame0,frame1))
                
                print("\nRe-setting isovalues...")
                self.changeIsosurfaceValue2(isoVals)
                
                # screenshot here if needed
                if screenCapture:
                    print("\nCapturing the screen...")
                    self.screenshot(labelOpt[1] %(frame0,frame1), ImgPath, MultiShots)
                
                cnt += 1
                
            else:
                break
                
        return


#========================================= RERconnect

    def RERconnect(self,
                   framePair,
                   ivFiletmpl,
                   modulePairs,
                   isoVals,
                   labelOpt = ("frame_label","RER_F%03d-F%03d"),
                   screenCapture = False,
                   ):
        """
        NOTE: in RER investigation, IV for (frame0-1) should be displayed as shaded
              and IV for frame1 should be displayed as lines.
              
        @param framePair: TUPLE variable, the two target frame numbers (int)
                          for desired period, e.g. (10 , 12)
        
        @param ivFiletmpl: STRING variable, e.g. "FRA2017_G02_Q03_ALL_F%03d.iv"
        
        @param modulePairs: DICTIONARY variable, pairs of module names template and destination
                            modulePairs ={"MgB1_F%03d_F%03d"      : "JNC_RMRER_B01",
                                          "Merge_F%03d_F%03d" : "JNC_faultslip", 
                                          .
                                          .
                                          .
                                          }
                            NOTE: modulePairs have to follow same frame numbering.
                            
        @param isovals: DICTIONARY, pairs of isosurface module name and value
                        rm_vals = {"LOW_15E3_D01" : 15000,
                                   "MED_250E3_D01" : 150000,
                                   "HIGH_500E3_D01" : 750000,
                                   .
                                   .
                                   .
                                   }
        
        @param labelOpt: TUPLE of two STRINGS, optional parameter set to 
                           ('frame_label','RER_F%03d-F%03d') as default. 
                           The two parameters are the Annotation module name and annotation text template.
        
        @param screencapture: BOOLEAN, default to False. True to capture screen.
                              Note that user will be asked to enter other screencapture function parameters.
        
        """
        # asking the user to make sure IVs are off
        raw_input("**** IMPORTANT **** \nMake sure all the IVs are turned OFF, then press ENTER to continue...")
        
        # asking the user for screenshot options if screenCapture is set to True
        if screenCapture:
            ImgPath = raw_input("Enter a path to save screenshots or press N to choose the current folder: ").upper()
            if ImgPath == "N":
                ImgPath = None
            MultiShots = raw_input("Would you like to take multiple screenshots for each frame? (Y or N) ").upper()
            if MultiShots == "Y":
                MultiShots = True
            else:
                MultiShots = False
        
        
        print("================================================================================================")

        
        frame0 = framePair[0]
        frame1 = framePair[1]
        
        print("Frame F%03d is being shown..." % frame1)
        
        print("\nTurning on IV frames F%03d and F%03d." %(frame0,frame1))
        self.showModule(ivFiletmpl %(frame0-1), True)
        self.changeDrawStyle(ivFiletmpl %(frame0-1), 0)
        self.showModule(ivFiletmpl % frame1, True)
        self.changeDrawStyle(ivFiletmpl % frame1, 2)
        
        print("\nConnecting modules...")
        # connecting module pairs
        modulePairs_temp = dict()
        for m in modulePairs.keys():
            modulePairs_temp.update({m  % (frame0,frame1):modulePairs[m]})
        self.connectMultiModules(modulePairs_temp)
        
        print #skipping a line in IPython console
        self.changeAnnotationText(labelOpt[0], labelOpt[1] %(frame0,frame1))
        
        print("\nRe-setting isovalues...")
        self.changeIsosurfaceValue2(isoVals)
        
        # screenshot here if needed
        if screenCapture:
            print("\nCapturing the screen...")
            self.screenshot(labelOpt[1] %(frame0,frame1), ImgPath, MultiShots)
        
                
        return
    
    
    
    
















##########################################
#----------------------------------------- SUPERSEDED - DON'T USE!!!
##########################################    
##########################################
#----------------------------------------- SUPERSEDED - DON'T USE!!!
##########################################    
##########################################
#----------------------------------------- SUPERSEDED - DON'T USE!!!
##########################################    
##########################################
#----------------------------------------- SUPERSEDED - DON'T USE!!!
##########################################    
##########################################
#----------------------------------------- SUPERSEDED - DON'T USE!!!
##########################################    
##########################################
#----------------------------------------- SUPERSEDED - DON'T USE!!!
##########################################    
  
    def SUPERSEDED_resetIsovalues(self,
                       rmIsovals,
                       fslipIsovals,
                       numOfBoxes,
                       moduleShow = None,
                       ):
        """
        @param rmIsovals: LIST of Rockmass isovalues, which has to have this
        format:
                [(first domain values:) #LOW, #MED, #HIGH,
                (second domain values:) #LOW, #MED, #HIGH,
                ...
                ]
                
        @param fslipIsovals: LIST of faultslip isovalues, which has to have 
        this format:
                [LOW, MED, HIGH]
                
        @param numOfBoxes: number of VTK boxes. NOTE that the isovalues for 
        different domains in different boxes have to be the same.
        
        @param moduleShow: String, defualt is 'True', set to 'False' to turn  
        modules off. Once set to False, the script will not change the isovalues.
        
        ###############################################
        ### COMPULSORY NAMING CONVENTION IN VOXLER: ###
        ###############################################
        1- fault slip isosurface modules, three isosurface modules with names:
            "LOW_###", "MED_###" and "HIGH_###" where the number is the isovalue
            e.g. "LOW_75000" is the isosurface module with isovalue of 75000
        
        2- Rockmass isosurface modules. There are three isosurface modules for
            each domain in each VTK box. So the number of isosurface modules 
            will be 3 x numDomains x numBoxes. Three isosurface modules for 
            each domain follow this naming convention:
                "LOW_###_B##", "MED_###_B##" and "HIGH_###_B##" where ### is 
                the isovalue and B## is the box number and domain number
                e.g. "MED_325000_B12" is the isosurface module with isovalue of
                325000 for second domain of box 1 (B12)
        """
        if moduleShow is None:
            moduleShow = 'True'
        
        numberOfDomains = len(rmIsovals)/3 # there are always three modules per domain (low, med, high)
        isoPrefix = ['LOW_','MED_','HIGH_']
        
        # rockmass isovalues
        for i in range(0,numOfBoxes):
            cnt = 0
            for ii in range(0,numberOfDomains):
                for iii,p in enumerate(isoPrefix):
                    name = p + str(rmIsovals[cnt]) + '_B' + str(i+1) + str(ii+1)
                    val = str(rmIsovals[cnt])
                    if moduleShow == 'True':
                        self.command(cmd_name="ShowModule",
                                        opts=dict([
                                                ("Module",name),
                                                ("Show",moduleShow),
                                                ]),
                                        )
                        self.command(cmd_name="ModifyModule",
                                        opts=dict([
                                                ("Module",name),
                                                ("IsosurfaceIsovalue",val),
                                                ]),
                                        )
                    elif moduleShow == 'False':
                        self.command(cmd_name="ShowModule",
                                        opts=dict([
                                                ("Module",name),
                                                ("Show",moduleShow),
                                                ]),
                                        )            
                    cnt += 1
        # faultslip isovalues
        for i,p in enumerate(isoPrefix):
            name = p + str(fslipIsovals[i])
            val = str(fslipIsovals[i])
            if moduleShow == 'True':
                self.command(cmd_name="ShowModule",
                                opts=dict([
                                        ("Module",name),
                                        ("Show",moduleShow),
                                        ]),
                                )
                self.command(cmd_name="ModifyModule",
                                opts=dict([
                                        ("Module",name),
                                        ("IsosurfaceIsovalue",val),
                                        ]),
                                )
            elif moduleShow == 'False':
                self.command(cmd_name="ShowModule",
                                 opts=dict([
                                         ("Module",name),
                                         ("Show",moduleShow),
                                         ]),
                                )
        
        
    def SUPERSEDED_RERconnect(self,
                   frameNumber,
                   offsetNumber,
                   ivFiletmpl,
                   vtkFiletmpl,
                   rmIsovals,
                   fslipIsovals,
                   vtkNumberOfBoxes = None,
                   vtkJunctiontmpl = None,
                   labeltmpl = None,
                   fslipModuletmpl = None,
                   fslipJunctiontmpl = None,
                   ):
        """
        @param frameNumber: integer, which is the frame number which we want
        to observe the results
        
        @param offsetNumber: integer, RER frame number offset
        
        @param ivFiletmpl: iv file template, e.g. "PH2017042_G01_Q02_excavUG_F%03d.iv"
        
        @param vtkFiletmpl: VTK file template. 
        IMPORTANT NOTE:
        1-  if there are more than two VTK boxes two variables are needed in 
            the VTK file template as box number and the frame number,
            e.g. "PH2017042_R06_Q02_B%02d_H3_F2%03d.vtk"
        2-  otherwise, only one variable is needed in the VTK file template as 
            the frame number, e.g. "PH2017042_R06_Q02_H3_F2%03d.vtk"
        
        @param rmIsovals: LIST of Rockmass isovalues, which has to have this
        format:
                [(first domain values:) #LOW, #MED, #HIGH,
                (second domain values:) #LOW, #MED, #HIGH,
                ...
                ]
        
        @param fslipIsovals: LIST of faultslip isovalues, which has to have 
        this format:
                [LOW, MED, HIGH]
        
        @param vtkNumberOfBoxes: integer, number of VTK boxes. Default = 2
        IMPORTANT NOTE: if there is only one box, the vtkJunctiontmpl will be
        changed inside the script.
        
        @param vtkJunctiontmpl: string with 1 variable (BOX number). 
        Default = "JUNCTION_RM_B%d"
        Junction module to connect to RM isosurfaces. 
        IMPORTANT NOTE: 
        1-  This input can have any format similar to "JUNCTION_RM_B%d",which 
            has the VTK box number as a variable in its name.
        2-  If only one box exists it will be overwritten to "JUNCTION_RM".
        
        @param labeltmpl: string with 1 variable. Deafualt = "Frame_%03d", 
        varible is the frameNumber.
                        
        @param fslipModuletmpl: string with 1 variable. Deafult is "F%03d_Merge"
        variable is frameNumber.
        IMPORTANT NOTE: can be any name with the format similar to "F%03d_Merge"
        
        @param fslipJunctiontmpl: string. Default is "JUNCTION_COH" 
        
        ###############################################
        ### COMPULSORY NAMING CONVENTION IN VOXLER: ###
        ###############################################
        1- JUNCTIONS:
            junction to connect RM RER values for >1 box: "JUNCTION_RM_B%d":
                box number is the variable.
            junction to connect RM RER values for 1 box: "JUNCTION_RM"
            junction to connect faultslip RER values: "JUNCTION_COH"
        
        2- faultslip Merge (or Grid modules if merge does not exist):
            'F%03d_Merge', where frame number is the variable.
        
        3- fault slip isosurface modules, three isosurface modules with names:
            "LOW_###", "MED_###" and "HIGH_###" where the number is the isovalue
            e.g. "LOW_75000" is the isosurface module with isovalue of 75000
        
        4- Rockmass isosurface modules. There are three isosurface modules for
            each domain in each VTK box. So the number of isosurface modules 
            will be 3 x numDomains x numBoxes. Three isosurface modules for 
            each domain follow this naming convention:
                "LOW_###_B##", "MED_###_B##" and "HIGH_###_B##" where ### is 
                the isovalue and B## is the box number and domain number
                e.g. "MED_325000_B12" is the isosurface module with isovalue of
                325000 for second domain of box 1 (B12)
        """
        
        # set default values to None variables
        if vtkNumberOfBoxes is None:
            vtkNumberOfBoxes = 2
        if vtkJunctiontmpl is None:
            vtkJunctiontmpl = "JUNCTION_RM_B%d"
        if labeltmpl is None:
            labeltmpl = "Frame_%03d"
        if fslipModuletmpl is None:
            fslipModuletmpl = "F%03d_Merge"
        if fslipJunctiontmpl is None:
            fslipJunctiontmpl = "JUNCTION_COH"        
        
        ### turn IVs on
        self.command(cmd_name="ShowModule",
                         opts=dict([
                                 ("Module",ivFiletmpl %frameNumber),
                                 ("Show","True"),
                                 ]),
                        )
        self.command(cmd_name="ShowModule",
                         opts=dict([
                                 ("Module",ivFiletmpl %(frameNumber-offsetNumber)),
                                 ("Show","True"),
                                 ]),
                        )
        ### change drawstyle of the second IV
        self.command(cmd_name="ModifyModule",
                        opts=dict([
                                ("Module",ivFiletmpl %frameNumber),
                                ("GeomSrcDrawStyle","2"),
                                ]),
                        )
        ### rename "frame_label" module (an annotation module)
        self.command(cmd_name="ModifyModule",
                        opts=dict([
                                ("Module","frame_label"),
                                ("AnnotationText",labeltmpl %frameNumber),
                                ]),
                        )
        ### connect modules
        # vtk modules
        if vtkNumberOfBoxes > 1:
            for i in range(0,vtkNumberOfBoxes):
                #vtkModule = vtkFiletmpl %(i+1,frameNumber)
                #vtkJunction = vtkJunctiontmpl %(i+1)
                self.command(cmd_name="ConnectModules",
                                 opts=dict([
                                         ("SourceModule",vtkFiletmpl %(i+1,frameNumber)),
                                         ("TargetModule",vtkJunctiontmpl %(i+1)),
                                         ]),
                                )
        else:
            #vtkModule = vtkFiletmpl %frameNumber
            #vtkJunction = "JUNCTION_RM" # CONSTANT variable
            vtkJunctiontmpl = "JUNCTION_RM"
            self.command(cmd_name="ConnectModules",
                             opts=dict([
                                     ("SourceModule",vtkFiletmpl %frameNumber),
                                     ("TargetModule",vtkJunctiontmpl),
                                     ]),
                            )
        # for Fault Slip connect Merge or any other name with similar format (name format: 'F%03d_Merge') to fslipJunctiontmpl
        self.command(cmd_name="ConnectModules",
                        opts=dict([
                                ("SourceModule", fslipModuletmpl %frameNumber),
                                ("TargetModule",fslipJunctiontmpl),
                                ]),
                        )
        ### re-set isovalues
        # rockmass isovalues
        numberOfDomains = len(rmIsovals)/3
        isoPrefix = ['LOW_','MED_','HIGH_']
        for i in range(0,vtkNumberOfBoxes):
            cnt = 0
            for ii in range(0,numberOfDomains):
                for iii,p in enumerate(isoPrefix):
                    name = p + str(rmIsovals[cnt]) + '_B' + str(i+1) + str(ii+1)
                    val = str(rmIsovals[cnt])
                    self.command(cmd_name="ModifyModule",
                                     opts=dict([
                                             ("Module",name),
                                             ("IsosurfaceIsovalue",val),
                                             ]),
                                    )
                    cnt += 1
        # faultslip isovalues
        for i,p in enumerate(isoPrefix):
            name = p + str(fslipIsovals[i])
            val = str(fslipIsovals[i])
            self.command(cmd_name="ModifyModule",
                             opts=dict([
                                     ("Module",name),
                                     ("IsosurfaceIsovalue",val),
                                     ]),
                            )



##########################################
#----------------------------------------- WORK IN PROGRESS
##########################################    

    def INCOMPLETE_importCSVFiles(self, 
            listOfFiles,
            importOptions = [],
            modulePosFirst = None,
            modulePosOffset = (0,25),
            moduleShow = True,
            ):
        """
        @param listOfFiles: list of CSV filenames (fullpath). An Exception is
          raised if a filename within the list does not exist.
          
        @param importOptions: optional argument, defaults to empty list. list 
          of options (as strings). For full description of CSV-file options
          see the Voxler-Automation Help (Construct Types > Import Options >
          CSV Import Options).
          
        @param modulePosFirst: optional argument, defaults to None. 
          If None, then automatic placement of the module within the 
          NetworkManager. Can be a tuple of NetworkManager coords (x,y).

        @param modulePosOffset: optional argument, defaults to (0,25).
          gives the offset between two imported modules within the 
          NetworkManager. Only used if modulePosFirst is not None.
          
        @param moduleShow: optional argument, defaults to True. Can be either
          a boolean value, or a list of boolean values.
        """

        if moduleShow is None:
            moduleShow = True
        if isinstance(moduleShow, bool):
            flag = bool(moduleShow)
            moduleShow = [ flag for i in range(len(listOfFiles)) ]
        elif isinstance(moduleShow, list):
            assert(len(moduleShow)==len(listOfFiles))
            assert(all([ isinstance(flag, bool) for flag in moduleShow ])==True)

        for i,(fn,show) in enumerate(zip(listOfFiles,moduleShow)):
    
            ### make sure fn exists to avoid non-meaningful ApiError exception
            if not os.path.exists(fn):
                raise Exception("try to import non-existing file '%s'. FAIL!"%fn)

            fn = fn.replace(os.sep,"/")
            moduleName = fn.rsplit("/",1)[-1]
            ### AF comment: it would be very nice if the voxler_lib.Voxler.command
            ### would return the assigned moduleName or moduleId
            ### Here, we can only hope the moduleName remains to be the filename

            ### create a csv-module via Import on the NetworkManager
            defaultOpts = [
                ("AutoConnect", "False"),
                ("ClearOptions", "True"),
                ("ClearPath", "True"),
                ("Filter", "csv"),
                ("GuiEnabled", "False"),
                ("PersistOptions", "False"),
                ("ProgressEnabled", "False"),
                ("UndoRedoEnabled", "False"),
                ("Path", str(fn)),
                #("DefaultPosition", "False")',
                #("XPosition", str(5+250)),
                #("YPosition", str(70+25*i)),
                ]
            if importOptions:
                defaultOpts.append(
                        ("Options", "; ".join(importOptions))
                    )
                #defaultOpts.extend([
                #    ("Option", opt) for opt in importOptions
                #    ])
            #self.command(
            self.command(
                cmd_name="Import",
                opts=dict(defaultOpts),
                )
            ### the following is wrong, this is more for sth like 'ModifyModule' with several api.do()'s
            #self.command_sequence(
            #    cmd_name="Import",
            #    opts=dict(defaultOpts),
            #    opts_seq=dict([ tuple(opt.split("=")) for opt in importOptions ])
            #    )
            
                
            ### if specified, re-position the created csv-module on the NetworkManager
            if modulePosFirst is not None:
                self.command(
                    cmd_name="MoveModule",
                    opts=dict([
                        ("Module", moduleName),
                        ("UseOriginalPosition", "False"),
                        ("XPosition", str(modulePosFirst[0]+i*modulePosOffset[0])),
                        ("YPosition", str(modulePosFirst[1]+i*modulePosOffset[1])),
                        ]),
                    )
                    
            ### if specified, turn the visibility of the created csv-module off
            if show==False:
                #self.command(
                self.command(
                    cmd_name="ShowModule",
                    opts=dict([
                        ("Module", moduleName),
                        ("Show", "False"),
                        ]),
                    )

        return




##########################################
##########################################
##########################################

