"""Library to create series of image compositions.

Usage
=====

Example: Script to compose series of multiple images with captions
 >>> import os, sys
 >>> from bae.compose_pngs_01 import PngImageComposer
 >>>
 >>> ### define frameList(s)
 >>> frameList = range(0,12+1,1)
 >>> frameIds = [ "F%03d"%frm for frm in frameList ]
 >>>
 >>> ### create a PngImageComposer Instance
 >>> w,h = 450,450
 >>> sep = 5
 >>> outStack = PngImageComposer(
 >>>     size=(2*w+1*sep,h),
 >>>     mode="RGBA",
 >>>     #bgcolor=None,
 >>>     bgcolor="white",
 >>>     )
 >>> outStack.setDefaultFont(fontFGColor="white", fontSize=18)
 >>> font1 = outStack.getDefaultFont()
 >>> font2 = dict(
 >>>     fontStyle="arial.ttf",
 >>>     fontSize=18,
 >>>     fontFGColor="black",
 >>>     fontBGColor="white",
 >>>     )
 >>>
 >>> ### positioning items
 >>> x = 0
 >>> y = 0
 >>> cropbox=(210,120,210+789,120+789)
 >>> outStack.insertImageStack(
 >>>     insertPos=(x,y),
 >>>     imgFilenames=[(
 >>>         "../R13H/PNGs_ss_zeroPWP/ViewE_VSectW250/PWP/"
 >>>         "P2016030_R13H_ss_zeroPWP_ViewE_VSectW250_PWP_2e+06_0_%s.png"
 >>>         )%frmName
 >>>         for frmName in frameIds ],
 >>>     imgCropBox=cropbox,
 >>>     insertSize=(w,h),
 >>>     )
 >>> outStack.insertTextStatic(
 >>>     insertPos=(x+10,y+20),
 >>>     text="R13H_ss_zeroPWP",
 >>>     fontDict=font1,
 >>>     )
 >>> outStack.insertTextStack(
 >>>     insertPos=(x+10,y+h-30),
 >>>     textStack=frameIds,
 >>>     fontDict=font1,
 >>>     )
 >>>
 >>> x += w+sep
 >>> outStack.insertImageStack(
 >>>     insertPos=(x,y),
 >>>     imgFilenames=[(
 >>>         "../R13H/PNGs_cn_ks1e-12/ViewE_VSectW250/PWP/"
 >>>         "P2016030_R13H_cn_ks1e-12_ViewE_VSectW250_PWP_2e+06_0_%s.png"
 >>>         )%frmName
 >>>         for frmName in frameIds ],
 >>>     imgCropBox=cropbox,
 >>>     insertSize=(w,h),
 >>>     )
 >>> outStack.insertTextStatic(
 >>>     text="R13H_cn_ks1e-12",
 >>>     insertPos=(x+10,y+20),
 >>>     fontDict=font1,
 >>>     )
 >>>
 >>> ### generate composed PNGs
 >>> outDir = sys.argv[0].replace(".py","")
 >>> outFilenames = [
 >>>     os.path.join(outDir,"DDP2016_comparePWP_EAST%02d_%s.png"%(i+1,frmId))
 >>>     for i,frmId in enumerate(frameIds)
 >>>     ]
 >>> outStack.generatePNGs(outFilenames)
"""

__version__ = "1.2"

_version_history_ = """
Versions:
=========

1.0 AF new
1.1 GP added getPngFileList
1.2 GP added PngImageComposerChequer
"""

import os
from PIL import Image
from PIL import ImageDraw, ImageFont

from bae.log_01 import msg


class _ImageStack(object):
    """
    Holds the information about a stack of images, ie., the list
    of image filenames, and the optional crop and/or resize information
    """

    def __init__(self, insertPos, imgFilenames, imgCropBox, imgReSize):
        """
        @param imgFilenames: list of image file paths. There can be None
           items then this item will be missing in the corresponding image
           of the resulting series.

        @param imgReSize: (width, height)-tuple. The image will be scaled to
           this size.
           If height (second item of the tuple) is None then the image will
           be scale to the specified width preserving the aspect ratio.
           If width (first item of the tuple) is None then the image will be
           scaled to the specified height preserving the aspect ratio.
        """
        self.insertPos = insertPos
        self.imgFilenames = imgFilenames
        self.imgCropBox = imgCropBox
        self.imgReSize = imgReSize
        #
        self.stackLength = len(imgFilenames)
        self.bgcolor = None

    def getItem(self, cnt):
        try:
            fn = self.imgFilenames[cnt]
        except IndexError:
            raise IndexError(
                "Index %d for this stack not applicable. returning None.")

        if not fn:
            return (None, None)

        if not os.path.isfile(fn):
            raise IOError("Could not find file '%s'" % fn)

        # load the png image
        image = Image.open(fn)

        # transparancy for PNGs
        if fn.endswith(".png") and (image.mode=="RGB" or image.mode=="P"):
            image = image.convert("RGBA")

        # crop image (if requested)
        if self.imgCropBox:
            image = image.crop(self.imgCropBox)

        # resize image (if requested)
        if self.imgReSize:
            w, h = self.imgReSize
            if h is None:
                wpercent = w / float(image.size[0])
                h = int(float(image.size[1]) * wpercent)
            if w is None:
                hpercent = w / float(image.size[1])
                w = int(float(image.size[0]) * hpercent)
            image = image.resize((w,h), Image.ANTIALIAS)

        return (image, self.insertPos)


class _TextStack(object):
    """
    Holds the information about a stack of texts, ie., the list
    of text and the font
    """

    def __init__(self, insertPos, textStack, fontDict):
        """
        @param fontDict: a dictionary with:
           "fontstyle": e.g. "arial.ttf"
           "fontSize": an integer
           "fontBGColor": "white", "black" or so
           "fontFGColor": dito
        """
        self.insertPos = insertPos
        self.textStack = textStack
        self.fontDict = fontDict
        #
        self.stackLength = len(textStack)
        try:
            self.imgFont = ImageFont.truetype(fontDict["fontStyle"],
                                              fontDict["fontSize"])
        except IOError:
            msg("WARNING: could not find the specified font '%s' and using the"
                " PIL-default instead." % fontDict["fontStyle"])
            self.imgFont = ImageFont.load_default().font
        self.bgcolor = fontDict["fontBGColor"]
        self.fgcolor = fontDict["fontFGColor"]
        self.lineSep = 5

    def getItem(self, cnt):
        try:
            text = self.textStack[cnt]
        except IndexError:
            raise IndexError(
                "Index %d for this stack not applicable. returning None." % cnt)
            return None
        else:
            strLines = text.split("\n")
            cntLines = len(strLines)

            ### determine text box size (generating each line of text on white
            ### background)
            image = Image.new("RGB", (0,0), "white")
            draw = ImageDraw.Draw(image)
            imgSize = [0,0]
            ypos = [0,]
            for strLine in strLines:
                lineSize = draw.textsize(strLine, font=self.imgFont)
                ypos.append(lineSize[1]+self.lineSep)
                imgSize = map(max, imgSize, lineSize)
            imgSize[1] += int(self.lineSep*(cntLines-1))
            imgSize = tuple(imgSize)
            del image

            ### actually write text as image
            image = Image.new("RGBA", imgSize, self.bgcolor)
            draw = ImageDraw.Draw(image)
            for cntLine,strLine in enumerate(strLines):
                draw.text(
                    (0, ypos[cntLine]),
                    strLine,
                    fill=self.fgcolor,
                    font=self.imgFont,
                    )

        return (image, self.insertPos)


class _TextStatic(object):
    """
    Holds the information about a single text, ie., the text and the font
    """

    def __init__(self, insertPos, text, fontDict):
        self.insertPos = insertPos
        self.text = text
        self.fontDict = fontDict
        #
        self.stackLength = 1

        try:
            self.imgFont = ImageFont.truetype(fontDict["fontStyle"],
                                              fontDict["fontSize"])
        except IOError:
            msg("WARNING: could not find the specified font '%s' and using the"
                " PIL-default instead." % fontDict["fontStyle"])
            self.imgFont = ImageFont.load_default().font

        self.bgcolor = fontDict["fontBGColor"]
        self.fgcolor = fontDict["fontFGColor"]
        self.lineSep = 5
        self.image = None

    def getItem(self, cnt):
        if self.image is None:

            strLines = self.text.split("\n")
            cntLines = len(strLines)

            ### determine text box size (generating each line of text on white
            ### background)
            image = Image.new("RGB", (0,0), "white")
            draw = ImageDraw.Draw(image)
            imgSize = [0,0]
            ypos = [0,]
            for strLine in strLines:
                lineSize = draw.textsize(strLine, font=self.imgFont)
                ypos.append(lineSize[1]+self.lineSep)
                imgSize = map(max, imgSize, lineSize)
            imgSize[1] += int(self.lineSep*(cntLines-1))
            imgSize = tuple(imgSize)
            del image

            ### actually write text as image
            image = Image.new("RGBA", imgSize, self.bgcolor)
            draw = ImageDraw.Draw(image)
            for cntLine,strLine in enumerate(strLines):
                draw.text(
                    (0, ypos[cntLine]),
                    strLine,
                    fill=self.fgcolor,
                    font=self.imgFont,
                    )
            self.image = image

        return (self.image, self.insertPos)


class PngImageComposer(object):
    """
    Class defining and generating the stack of result (composed) PNGs
    """

    def __init__(self, size, mode="RGBA", bgcolor=None):
        """
        @param size: (width, height)-tuple in points, e.g. (1680, 830)
        @param mode: e.g. "RGBA", "RGB"
        @param bgcolor: e.g. "white"
        """
        self.resultSize = size
        self.resultMode = mode
        self.resultBGColor = bgcolor
        self.maxStackLength = None
        self.insertedComponents = []
        self.defaultFont = dict(
            fontStyle="arial.ttf",
            fontSize=32,
            fontFGColor="black",
            fontBGColor=None,
            )
        return

    def setDefaultFont(self, fontStyle="arial.ttf", fontSize=32,
                       fontFGColor="black", fontBGColor=None):
        self.defaultFont = dict(
            fontStyle=fontStyle,
            fontSize=fontSize,
            fontFGColor=fontFGColor,
            fontBGColor=fontBGColor,
            )

    def getDefaultFont(self):
        return self.defaultFont

    def updateMaxStackLength(self):
        """Service function being called by the insert...() methods.
        Usually there is no need to call it explicitly.
        """
        maxStackLength = 0
        for imgComponent in self.insertedComponents:
            curStackLength = imgComponent.stackLength
            maxStackLength = max(maxStackLength,curStackLength)
        self.maxStackLength = maxStackLength

    def insertImageStack(self, insertPos=(0,0), imgFilenames=[],
                         imgCropBox=None, insertSize=None):
        """
        @param imgFilenames: list of image file paths. There can be None
           items then this item will be missing in the corresponding image
           of the resulting series.

        @param insertSize: (width, height)-tuple. The image will be scaled to
           this size.
           If height (second item of the tuple) is None then the image will be
           scaled to the specified width preserving the aspect ratio.
           If width (first item of the tuple) is None then the image will be
           scaled to the specified height preserving the aspect ratio.
        """
        imgComponent = _ImageStack(
            insertPos=insertPos,
            imgFilenames=imgFilenames,
            imgCropBox=imgCropBox,
            imgReSize=insertSize,
            )
        self.insertedComponents.append(imgComponent)
        self.updateMaxStackLength()

    def insertTextStack(self, insertPos=(0,0), textStack=[], fontDict=None):
        """
        @param textStack: List of texts for each output frame. There can be
           None items then this item will be missing in the corresponding
           frame / image of the resulting series.

        @param fontDict: optional dictionary with:
           "fontstyle": e.g. "arial.ttf"
           "fontSize": an integer
           "fontBGColor": "white", "black" or so
           "fontFGColor": dito
        """
        if fontDict is None:
            fontDict = self.getDefaultFont()

        imgComponent = _TextStack(
            insertPos=insertPos,
            textStack=textStack,
            fontDict=fontDict,
            )
        self.insertedComponents.append(imgComponent)
        self.updateMaxStackLength()

    def insertTextStatic(self, insertPos=(0,0), text="", fontDict=None):
        """
        @param fontDict: optional dictionary with:
           "fontstyle": e.g. "arial.ttf"
           "fontSize": an integer
           "fontBGColor": "white", "black" or so
           "fontFGColor": dito
        """
        if fontDict is None:
            fontDict = self.getDefaultFont()
        imgComponent = _TextStatic(
            insertPos=insertPos,
            text=text,
            fontDict=fontDict,
            )
        self.insertedComponents.append(imgComponent)
        self.updateMaxStackLength()

    def generatePNGs(self, outFilenames=None, test=None):
        """
        @param outFilenames: list of file names for the resulting composed
           images. If not given generic names will be generated:
           "imageComposerOutput%04d.png" % index

           This list determines the number of images that will be created. I.e.
           if outFileNames is shorter than the any of the content-stacks then
           excess items will be skipped.

        @param test: Optional index in outFilenames identifying the one and
           only page to be created for a test.
        """
        if outFilenames is None:
            outFilenames = [
                "imageComposerOutput%04d.png" % i
                for i in range(self.maxStackLength)
                ]
        if test is None:
            outList = [ (i,fn) for i,fn in enumerate(outFilenames) ]
        else:
            idx = int(test)
            if idx<0:
                idx += len(outFilenames)
            outList = [ (idx,outFilenames[idx]), ]

        ### loop over all output files
        for cnt,outFilename in outList:

            msg("  > composing '%s'..." % outFilename)

            ### initialize new result image
            result = Image.new(
                self.resultMode, self.resultSize, self.resultBGColor)

            ### loop over all items to be included
            for cntItem,insertItem in enumerate(self.insertedComponents):
                msg("    - inserting Item %d..." % cntItem)
                try:
                    img,pos = insertItem.getItem(cnt)
                except Exception as e:
                    msg("   WARNING: Item not inserted due to the following"
                        " error:\n   %s" % repr(e))
                    continue

                # skip None item in this image-stack
                if not img:
                    continue

                try:
                    bgcolor = insertItem.bgcolor
                except:
                    bgcolor = None
                if bgcolor is None:
                    result.paste(img, pos, img)
                    #result.paste(img, pos)
                else:
                    result.paste(img, pos)

            ### display (in test mode only)
            if test is not None:
                result.show()

            ### save
            tmp = outFilename.replace(os.sep,"/").rsplit("/",1)
            if len(tmp)==2:
                outDir = tmp[0]
                if not os.path.isdir(outDir):
                    os.makedirs(outDir)
            result.save(outFilename)
            msg("    saved\n")


class PngImageComposerChequer(PngImageComposer):
    """Combine images (and texts) in a chequerboard pattern.

     >>> from bae.compose_pngs_01 import \
     >>>     PngImageComposerChequer as Composer, getPngFileList
     >>>
     >>> frameList = [(2,i) for i in range(1, 25)]
     >>> composer = Composer(
     >>>     widths=[450,450,450], heights=[50,360],
     >>>     separator=5, bgcolor="white")
     >>> composer.insertTextStack(
     >>>     (0,0), ["Frame F%d%03d  option 1" % f for f in frameList])
     >>> composer.insertTextStatic((0,1), "option 2")
     >>> composer.insertTextStatic((0,2), "option 3")
     >>> composer.insertImageStack(
     >>>     (1,0), getPngFileList("PNG/opt1/", frameList)
     >>> composer.insertImageStack(
     >>>     (1,1), getPngFileList("PNG/opt2/", frameList)
     >>> composer.insertImageStack(
     >>>     (1,2), getPngFileList("PNG/opt3/", frameList)
     >>>
     >>> composer.generatePNGs(
     >>>     ["PNG/compare/Compare_opt123_F%d%03d" % f for f in frameList])

    Note: L{bae.avi_01} default image size is 1280x1024, ratio is 5:4.
    """
    def __init__(self, widths, heights, separator, mode="RGBA", bgcolor=None):
        """
        @param widths: list of widths of the subplots
        @param heights: list of heights of the subplots
        @param separator: distance between neighbouring images
        @param mode: e.g. "RGBA", "RGB"
        @param bgcolor: e.g. "white"
        """
        self.widths = widths
        self.heights = heights
        self.separator = separator
        totalwidth = sum(widths) + separator*(len(widths)-1)
        totalheight = sum(heights) + separator*(len(heights)-1)
        PngImageComposer.__init__(
            self, size=(totalwidth, totalheight),
            mode=mode, bgcolor=bgcolor)

    def getInsertPosFromIndex(self, locationIndex):
        """Compute the insertPos argument from an index tuple.
        Used by L{insertImageStack} and similar methods.
        """
        i,j = locationIndex  # abbreviation
        return (
            sum(self.widths[:i])+i*self.separator,
            sum(self.heights[:j])+j*self.separator)

    def insertImageStack(self, locationIndex, imgFilenames=[],
                         imgCropBox=None):
        """
        @param locationIndex: (0,0) is the top left subplot, (0,1) is its
           right neighbour...

        @param imgFilenames: list of image file paths. There can be None
           items then this item will be missing in the corresponding image
           of the resulting series.
        """
        i,j = locationIndex  # abbreviation
        PngImageComposer.insertImageStack(
            self,
            insertPos=self.getInsertPosFromIndex(locationIndex),
            imgFilenames=imgFilenames, imgCropBox=imgCropBox,
            insertSize=(self.widths[i],self.heights[j]))

    def insertTextStack(self, locationIndex, textStack=[], fontDict=None):
        """
        @param locationIndex: (0,0) is the top left subplot, (0,1) is its
           right neighbour...

        @param textStack: List of texts for each output frame. There can be
           None items then this item will be missing in the corresponding
           frame / image of the resulting series.

        @param fontDict: optional dictionary with:
           "fontstyle": e.g. "arial.ttf"
           "fontSize": an integer
           "fontBGColor": "white", "black" or so
           "fontFGColor": dito
        """
        i,j = locationIndex  # abbreviation
        PngImageComposer.insertTextStack(
            self,
            insertPos=self.getInsertPosFromIndex(locationIndex),
            textStack=textStack, fontDict=fontDict)

    def insertTextStatic(self, locationIndex, text, fontDict=None):
        """
        @param locationIndex: (0,0) is the top left subplot, (0,1) is its
           right neighbour...
        @param fontDict: optional dictionary with:
           "fontstyle": e.g. "arial.ttf"
           "fontSize": an integer
           "fontBGColor": "white", "black" or so
           "fontFGColor": dito
        """
        i,j = locationIndex  # abbreviation
        PngImageComposer.insertTextStatic(
            self,
            insertPos=self.getInsertPosFromIndex(locationIndex),
            text=text, fontDict=fontDict)


def getPngFileList(imagesDir, frameList, stepFrameToName="_F%d%03d"):
    """Return a list of png-file names in a given folder that correspond
    to the frames in frameList.

    For each frame in frameList the stepFrameToName expression or -function is
    evaluated to obtain a frame identifier string (e.g. (2,3) -> "F2003").
    A file is associated to a frame if this frame identifier string is part of
    the file name.

    The Nth entry in the resulting list is the file name of the image file
    corresponding to the Nth frame in frameList. If such a file could not be
    found then the corresponding item is None. The file names are complete
    paths, i.e. including the directory as specified by imagesDir.

    @param imagesDir: name (path) of the directory to search for image files
    @param frameList: list of (stepNb, frameNb)-tuples
    @param stepFrameToName: Defines how to build a frame identifier string
      (aka. frame name) from a (stepNb, frameNb)-tuple.
      Might be a function (stepNb, frameNb)-tuple --> frameName or a template
      string (like the default "_F%d%03d").
      The resulting string is used as search string in image file names to
      identify the image file corresponding to a particular frame.
    @returns: a list of paths to the image files.
    """

    # check argument stepFrameToName (preprocess if string)
    # if stepFrameToName is format string create a function
    if isinstance(stepFrameToName, basestring):
        stepFrameToNameStr = stepFrameToName
        def stepFrameToName(stepNb,frameNb):
            return stepFrameToNameStr % (stepNb, frameNb)

    # list of png files in imagesDir
    filesInDir = [x for x in os.listdir(imagesDir)
                  if x.lower().endswith(".png")]
    msg("Files in imagesDir %s:\n%s" % (imagesDir, filesInDir),
        debugLevel=10)

    fileList = list()
    for stepNb, frameNb in frameList:
        frameName = stepFrameToName(stepNb,frameNb)
        found = None
        for fn in filesInDir:
            if frameName in fn:
                found = fn
                break

        if found:
            filesInDir.remove(fn)
            path = os.path.join(imagesDir, fn)
            fileList.append(path)
            msg("For frame %s found file %s" % (frameName, path),
                debugLevel=10)
        else:
            fileList.append(None)
            msg("For frame %s did not find any file.", debugLevel=10)
    return fileList
