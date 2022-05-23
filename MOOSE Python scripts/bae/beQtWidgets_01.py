#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 08:37:13 2020

@author: tobias
"""
import sys
try:
    from PyQt4.QtGui import (QApplication,
                             QWidget, QCheckBox, QComboBox,
                             QMenu, QDoubleSpinBox, QSpinBox,
                             QInputDialog,QDialog,
                             QFormLayout, QLineEdit, QFileDialog,
                             QPushButton, QPlainTextEdit,
                             QVBoxLayout, QHBoxLayout, QSizePolicy,
                             QListWidget, QAbstractItemView,
                             QListWidgetItem,
                             QProgressBar, QGroupBox, QLabel,
                             QFrame)
    from PyQt4 import QtCore

except ImportError:
    from PyQt5.QtWidgets import (QApplication,
                                 QWidget, QCheckBox, QComboBox,
                                 QMenu, QDoubleSpinBox, QSpinBox,
                                 QInputDialog,QDialog,
                                 QFormLayout, QLineEdit, QFileDialog,
                                 QPushButton, QPlainTextEdit,
                                 QVBoxLayout, QHBoxLayout, QSizePolicy,
                                 QListWidget, QAbstractItemView,
                                 QListWidgetItem,
                                 QProgressBar, QGroupBox, QLabel,
                                 QFrame)

    from PyQt5 import QtCore

try:
    # raise ImportError
    from PyQt5.Qsci import (QsciScintilla, QsciLexerPython, QsciLexerFortran,
                            QsciLexerCPP)
except ImportError:
    QsciScintilla = None



import os, errno
import sip
from collections import OrderedDict

### support functions #########################################################
def clearLayout(layout):
    """Deletes a layout and all its items/subLayouts.
    Found here:
        U{https://stackoverflow.com/questions/22623151/python-how-to-unassign-layout-from-groupbox-in-pyqt}
    """
    while layout.count():
        item = layout.takeAt(0)
        widget = item.widget()
        if widget is not None:
            widget.deleteLater()
        else:
            clearLayout(item.layout())
    sip.delete(layout)


### Widgets ###################################################################
class RichQDoubleSpinBox(QDoubleSpinBox):
    '''Subclassed QDoubleSpinBox having a rightclick-popUp menu to set step'''
    def __init__(self, name='', value=1., upper=None, lower=None, step=1,
                 parent=None, decimals=None, suffix=None):
        super(RichQDoubleSpinBox, self).__init__(parent=parent)
        self.name = name

        self.setValue(value)
        if upper is not None:
            self.setMaximum(upper)
        if lower is not None:
            self.setMinimum(lower)
        if decimals is not None:
            self.setDecimals(decimals)
        if suffix is not None:
            self.setSuffix(suffix)
        self.setSingleStep(step)

        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.showMenu)


    def showMenu(self, event):
        '''Showes a simple popup menu and calls its action/callback on press'''
        menu = QMenu()
        setStepAction = menu.addAction("Set step")
        action = menu.exec_(self.mapToGlobal(event))
        if action == setStepAction:
            self.setStep()

    def setStep(self):
        '''Opens a simple dialog to enter the step of spinbox'''
        currentStep = self.singleStep()
        d, okPressed = QInputDialog.getDouble(self, "Set step",
                            "Value:", currentStep, 0, 1E32, 3)
        if okPressed:
            self.setSingleStep(d)


class RichQSpinBox(QSpinBox):
    '''Subclassed QSpinBox having a rightclick-popUp menu to set step'''
    def __init__(self, name='', value=1, upper=None, lower=None, step=1,
                 parent=None):
        super(RichQSpinBox, self).__init__(parent=parent)
        self.name = name
        if upper is not None:
            self.setMaximum(int(upper))
        if lower is not None:
            self.setMinimum(int(lower))
        self.setValue(int(value))
        self.setSingleStep(int(step))

        self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.showMenu)


    def showMenu(self, event):
        '''Showes a simple popup menu and calls its action/callback on press'''
        menu = QMenu()
        setStepAction = menu.addAction("Set step")
        action = menu.exec_(self.mapToGlobal(event))
        if action == setStepAction:
            self.setStep()

    def setStep(self):
        '''Opens a simple dialog to enter the step of spinbox'''
        currentStep = self.singleStep()
        d, okPressed = QInputDialog.getDouble(self, "Set step",
                            "Value:", currentStep, 0, int(1E32), 3)
        if okPressed:
            self.setSingleStep(d)


class Path(str):
    """A string-object enriched by some os fuctions
    """
    def exists(self):
        return os.path.exists(self)

    def dir(self, newDir=None, rel=False):
        """Set or gets the directory of a file"""
        if newDir is None:
            return os.path.dirname(self)
        else:
            if rel:
                return Path(os.path.join(self.dir(), newDir, self.file()))
            else:
                return Path(os.path.join(newDir, self.file()))

    def file(self):
        """Returns fileName without path"""
        return os.path.basename(self)

    def ext(self):
        """Returns extension of a file-path"""
        return os.path.splitext(self)[1]

    def createDir(self):
        """Creates directory (of file) if not exists"""
        try:
            os.makedirs(self.dir())
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise


class FileLineEdit(QWidget):
    """Extended FileLineEdit which will check the existance of a filePath.
    On double click it will open the fileOpenDialog.
    """
    _notExistCol = "color: rgb(224, 96, 96);"
    _existCol = "color: rgb(96, 224, 96);"
    _plainCol = "color: rgb(0, 0, 0);"

    def __init__(self, path='', initSearchPath='.', description='',
                 colorExists=True, filter='All files (*.*)',
                 evalInstant=True, folder=False):
        super(FileLineEdit, self).__init__()
        self.lineEdit = QLineEdit(parent=self)
        self.value = Path(path)
        self.description = description
        self.filter = filter
        self.selFolder = folder
        self.initSearchPath = initSearchPath
        self.colorExists = colorExists

        self.lineEdit.setText(self.value)
        if evalInstant:
            self.lineEdit.textChanged.connect(self.onTextChanged)
        else:
            self.lineEdit.editingFinished.connect(self.onTextChanged)

        self.button = QPushButton('Set',parent=self)
        self.button.clicked.connect(self.onSetPath)
        hLayout = QHBoxLayout()
        hLayout.addWidget(self.lineEdit)
        hLayout.addWidget(self.button)
        self.setLayout(hLayout)

    def onSetPath(self):
        old = self.lineEdit.text()
        if self.selFolder:
            fileName = QFileDialog.getExistingDirectory(self,
                                                        self.description,
                                                        self.initSearchPath,
                                            QFileDialog.ShowDirsOnly
                                             | QFileDialog.DontResolveSymlinks)
        else:
            fileName = QFileDialog.getOpenFileName(self,
                                                   self.description,
                                                   self.initSearchPath,
                                                   self.filter)[0]

        if fileName:
            self.lineEdit.setText(fileName)
        else:
            self.lineEdit.setText('')
            self.lineEdit.setText(old)

    def onTextChanged(self, newPath):
        """Callback on text changed - if file exists-->green"""
        self.value = Path(newPath)
        if self.colorExists:
            if not self.value.exists():
                self.lineEdit.setStyleSheet(self._notExistCol)
            else:
                self.lineEdit.setStyleSheet(self._existCol)


class PlainParameterWidget(QWidget):
    def __init__(self, parent=None, parameters={}, columns=1):
        """Creates a simple 2-column-layout with the first column showing
        the parameter key and the second column holding an appropieate
        input QT-input widget to show/mutate the parmeter value.
        Depending on the parameter type the following QTWidgets will be set
        and adjusted by optional keywords passed:
            - bool -> QCheckbox
            - float/int -> QSpinbox (aditional keywords: upper, lower, step)
            - string -> QTextBox or QCombobox if additional keyword
            - keys=['keyA','keyB',...] is passed
        """
        super(PlainParameterWidget, self).__init__(parent=parent)
        self.mainUI = parent
        self.parameters = OrderedDict(parameters)
        self.valueBoxes = {}
        self.formLayout = QFormLayout()
        self.updateLayout()
        self.setLayout(self.formLayout)

    def updateLayout(self):
        self.clearForm()
        for ii,(name, data) in enumerate(self.parameters.items()):
            value = data.get('value',0.) # default value
            keys = data.get('keys', [])

            ## value is bool --> checkbox
            if isinstance(value, bool):
                valueBox = QCheckBox(parent=self)
                valueBox.setChecked(value)
                valueBox.stateChanged.connect(
                    lambda _,name=name: self.onValueChanged(name))

            ## value is int or float --> (Rich)QSpinBox
            elif isinstance(value, (float, int)):
                upper = data.get('upper',None)
                lower = data.get('lower',None)
                step = data.get('step',1.)
                decimals = data.get('decimals', 3)
                if isinstance(value, float):
                    valueBox = RichQDoubleSpinBox(name, value,
                                                  upper, lower, step,
                                                  decimals=decimals,
                                                  parent=self)
                else:
                    valueBox = RichQSpinBox(name, value,
                                            upper, lower, step,
                                            parent=self)
                valueBox.valueChanged.connect(
                    lambda _,name=name: self.onValueChanged(name))

            ## value is a Path
            elif isinstance(value, Path):
                kwargs = dict((key,val) for key,val in data.iteritems()
                              if not key=='value')
                kwargs['path'] = value
                valueBox = FileLineEdit(**kwargs)
                valueBox.editingFinished.connect(
                        lambda name=name: self.onValueChanged(name))

            ## value is string
            elif isinstance(value, basestring):
                if keys: # --> comboBox
                    valueBox = QComboBox(parent=self)
                    valueBox.addItems(map(str, keys))
                    valueBox.setCurrentIndex(keys.index(value))
                    valueBox.currentIndexChanged.connect(
                            lambda ii, name=name: self.onValueChanged(name))
                else:
                    valueBox = QLineEdit(parent=self)
                    valueBox.setText(value)
                    valueBox.editingFinished.connect(
                            lambda name=name: self.onValueChanged(name))

            else:
                raise NotImplementedError('value type not supported yet')

            self.valueBoxes[name] = valueBox

            self.formLayout.addRow(name,valueBox)

    def clearForm(self):
        """Removes all parameter widgets from formlayout"""
        while self.formLayout.count():
            widget = self.formLayout.takeAt(0).widget()
            widget.deleteLater()


    def onValueChanged(self, name):
        """CallBack for changed parameters widgets which will update the
        self.parameters dictionary.
        """
        valueBox = self.valueBoxes[name]

        # value was float or int -- QSpinBox
        if isinstance(valueBox, RichQSpinBox):
            self.parameters[name]['value'] = valueBox.value()

        # value was bool -- QCheckBox
        elif isinstance(valueBox, QCheckBox):
            self.parameters[name]['value'] = bool(valueBox.isChecked())

        elif isinstance(valueBox, FileLineEdit):
            self.parameters[name]['value'] = str(valueBox.text())
            FileLineEdit.onTextChanged(valueBox, valueBox.text())

        # value was string (with keys) -- QComboBox
        elif isinstance(valueBox, QComboBox):
            self.parameters[name]['value'] = str(valueBox.currentText())

        # value was string (without keys) --
        elif isinstance(valueBox, QLineEdit):
            self.parameters[name]['value'] = str(valueBox.text())


class SimpleCodeEditor(QWidget):
    """Simple (code)Editor widget that allows
        - Syntax highlighting*
        - reasonable tab-formating*
        - saving and loading of codefiles
    * To enable these very usefull features you'll need to install
    U{QsciScintilla<https://qscintilla.com/#home>}.

        >>> pip install qscintilla --no-cache-dir

    or

        >>> conda install -c conda-forge qscintilla2

    The current text can be accessed via L{getText}-function. To modify the
    text programmatically use L{setText}.

    To add function buttons next to load/save buttons you can add
    these to the QHorzLayout menuLayout. E.g. you add this widget to a QDialog
    (self in following snippet):

        >>> self.editor = SimpleCodeEditor(parent=self, initText='blub')
        >>> menuLayout = self.editor.menuLayout
        >>> button = QPushButton('PrintBlob')
        >>> button.pressed.connect(lambda s: self.editor.setText('blob'))
        >>> menuLayout.addWidget(button)

    """



    extensions = {
            'python' : 'python script (*.py)',
            'fortran' : 'fortran source code (*.f)',
            'c'       : 'c source code (*.c)',
            'cpp'     : 'c++ source code (*.cpp)',
            'c++'     : 'c++ source code (*.cpp)',
            'abaqus'  : 'abaqus input (*.inp)',}


    def __init__(self, parent=None, initText=None, language='python'):
        super(SimpleCodeEditor, self).__init__(parent=parent)

        self.language = language

        # layout
        menuLayout = QHBoxLayout()
        loadButton = QPushButton('Load script')
        loadButton.pressed.connect(self.onLoadPressed)
        saveButton = QPushButton('Save script')
        saveButton.pressed.connect(self.onSavePressed)
        menuLayout.addWidget(loadButton)
        menuLayout.addWidget(saveButton)
        self.menuLayout = menuLayout

        self.layout = QVBoxLayout()
        sizePolicy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.setSizePolicy(sizePolicy)

        if QsciScintilla:
            # QScintilla editor setup
            # ------------------------
            self.__editor = QsciScintilla()
            self.__editor.setUtf8(True)  # Set encoding to UTF-8

            # TabSettings
            self.__editor.setTabWidth(4)
            self.__editor.setIndentationsUseTabs(False)
            self.__editor.setIndentationGuides(True)
            self.__editor.setAutoIndent(True)

            # Set lexer
            if self.language.lower() == 'python':
                self.__lexer = QsciLexerPython(self.__editor)
            elif self.language.lower() == 'fortran':
                self.__lexer = QsciLexerFortran(self.__editor)
            elif self.language.lower() in ['c','c++','cpp']:
                self.__lexer = QsciLexerCPP(self.__editor)

            self.__editor.setLexer(self.__lexer)

            # Brace matching: enable for a brace immediately before or after
            # the current position
            self.__editor.setBraceMatching(QsciScintilla.SloppyBraceMatch)

            self.__editor.SendScintilla(QsciScintilla.SCI_SETHSCROLLBAR, 0)

        else:
            self.__editor = QPlainTextEdit()

        # Add editor to layout
        self.layout.addLayout(menuLayout)
        self.layout.addWidget(self.__editor)
        self.setLayout(self.layout)

        if initText:
            self.setText(initText)

    def setText(self, text, append=True):
        """Sets the displayed text.
        @param text: text to add/set

        @param append: if True, text will be appended, if False, all will be
        replaced
        """
        if QsciScintilla:
            if append:
                self.__editor.append(text)
            else:
                self.__editor.setText(text)
        else:
            if append:
                self.__editor.appendPlainText(text)
            else:
                self.__editor.setPlainText(text)


    def getText(self):
        """Returns the current editor text"""
        return self.__editor.text()


    def onSavePressed(self):
        """CallBack for save button"""
        fileName = QFileDialog.getSaveFileName(self, 'Set FileName', './',
                                               self.extensions[self.language])
        if isinstance(fileName, tuple):
            fileName = fileName[0]
        if not fileName:
            return
        fid = open(fileName, 'w')
        fid.write(self.getText())


    def onLoadPressed(self):
        """CallBack for load button"""
        fileName = QFileDialog.getOpenFileName(self,
                                               'Select File to load',
                                               self.extensions[self.language])
        if isinstance(fileName, tuple):
            fileName = fileName[0]
        if not fileName:
            return
        fid = open(fileName, 'r')
        self.setText('\n'.join(fid.readlines()), append=False)


class DataQListWidgetItem(QListWidgetItem):
    def __init__(self, *args, **kwargs):
        self.itemData = kwargs.pop('itemData', {})
        super(DataQListWidgetItem,self).__init__(*args, **kwargs)


class DragDropQListWidget(QListWidget):
    """QListWidget which sets Drag/Drop ability via init-proceedure.
    The functions L{iterAllItems} and L{itemList} make the current content
    of widget easily accessable
    """
    def __init__(self, **kwargs):
        """
        @kwarg parent: qt-parentWidget (default is None)

        @kwarg items: list of strings wich will be set as listItems

        @kwarg dragMode: sets the drag behavior
            0 = no draggin allowed
            1 = allow dragging, dragged items remain in list (default)
            2 = allow dragging, dragged items will be removed from list

        @kwarg dropMode: sets the drop behavior
            0 = no dropping allowed
            1 = allow dropping and allow duplicated items
            2 = allow dropping but reject duplicated items

        @kwarg rejectItems: list of itemNames which will be rejected on
            dragAndDrop operation. Warning! These items can be set during
            __init__ or via addItem(s) but will removed on next drop operation

        @kwarg selMode: sets selection logic
            0 = only single item can be selected
            1 = multiple items can be selected
                (ctrl - to select arbitary items, shift to select all between)

        @kwarg itemData: list with same length as items, where each item will
            stored in listItem.itemData
            To use this data as a sorter, each item has to be a dictionary
            with a key named according to requested sorterKey

        @kwarg sorter: a list of sorters (string), "name" will sort items by
            item.text(), all other sorterKeys will use itemData for sorting

        Example:
            >>> myList = DragDropQListWidget(
            ...                      items=['item%d' % d for d in range(6)],
            ...                      itemData=[{'rand': random()}
            ...                                 for d in range(6)],
            ...                      sorter=['name', 'rand'])
        """

        parent = kwargs.get('parent', None)
        items  = kwargs.get('items', [])
        itemData = kwargs.get('itemData', [])
        rejectItems = kwargs.get('rejectItems', [])


        if itemData and not len(itemData) == len(items):
            raise ValueError('itemData must have same size as items')

        # drag/dropMode:
        #   0: not drag/dropAble
        #   1: if drag items remains, if dropitem will be appended
        #   2: if drag items will be removed, if drop unique(set) add item
        self.dragMode = kwargs.get('dragMode', 1)
        self.dropMode = kwargs.get('dropMode', 1)
        self.selMode  = kwargs.get('selMode', 1) #0 single, 1 multiple
        self.sorter = kwargs.get('sorter', ['name',])

        self.rejectItems = rejectItems
        self.allowDropsFrom = None

        super(DragDropQListWidget, self).__init__(parent=parent)
        self._parent = parent

        self.setupUI()
        if items:
            for ii,item in enumerate(items):
                item = DataQListWidgetItem(item)
                if itemData:
                    item.itemData = itemData[ii]
                self.addItem(item)

    def addItems(self, items, itemData=[]):
        for ii, item in enumerate(items):
            if itemData:
                data = itemData[ii]
            else:
                data = {}
            self.addItem(DataQListWidgetItem(item, itemData=data))

    def setupUI(self):
        if self.dragMode and self.dropMode:
            self.setDragDropMode(QAbstractItemView.DragDrop)
        elif self.dragMode == 1 and not self.dropMode:
            self.setDragDropMode(QAbstractItemView.DragOnly)
        elif not self.dragMode and self.dropMode:
            self.setDragDropMode(QAbstractItemView.DropOnly)

        if self.dragMode == 2:
            self.setDefaultDropAction(QtCore.Qt.MoveAction)

        if self.selMode:
            self.setSelectionMode(QListWidget.ExtendedSelection)

        if self.sorter:
            self.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
            self.customContextMenuRequested.connect(self.showContextMenu)

    def iterAllItems(self):
        """Returns itterator over all current widgetItems"""
        for ii in range(self.count()):
            yield self.item(ii)

    def itemList(self):
        """Returns a list of all current widgetStrings"""
        return [str(item.text()) for item in self.iterAllItems()]

    def dragEnterEvent(self, e):
        if not self.allowDropsFrom is None:
            if e.source() not in self.allowDropsFrom + [self,]:
                e.ignore()
            else:
                e.accept()
        else:
            e.accept()

    def dropEvent(self, e):
        # try to store itemData from eventSource (other DragDropQListWidget)
        try:
            itemData = dict((item.text(), item.itemData)
                            for item in e.source().iterAllItems())
        except AttributeError:
            itemData = dict((item.text(), {})
                            for item in e.source().iterAllItems())

        # now run default drop procedure and see which items are new
        QListWidget.dropEvent(self, e)

        newItemList, alreadyInList = [], []
        for item in self.iterAllItems():
            itemText = item.text()
            newItem = DataQListWidgetItem(itemText)

            if itemText in self.rejectItems:
                continue

            # skip duplicated item if dropmode does not allow duplicates
            if self.dropMode == 2 and itemText in alreadyInList:
                continue
            else:
                alreadyInList.append(itemText)

            # try to get itemData from source list
            try:
                newItem.itemData = itemData[item.text()]
            except KeyError:
                newItem.itemData = item.itemData
            newItemList.append(newItem)

        # rebuild listView
        self.clear()
        for newItem in newItemList:
            self.addItem(newItem)

    def showContextMenu(self, event):
        menu = QMenu(self)

        sorters = {}
        for sorter in self.sorter:
            for direction, verbal in zip([1,-1], ['ascending', 'descending']):
                action = menu.addAction("Sort by %s %s" % (sorter, verbal))
                sorters[(sorter, direction)] = action

        action = menu.exec_(self.mapToGlobal(event))

        try:
            idx = sorters.values().index(action)
        except ValueError:
            return

        sortKey, direction = sorters.keys()[idx]
        self.sortBy(sortKey, direction)

    def sortBy(self, sortKey, direction=1):
        unsortedItems = []
        if sortKey == 'name':
            sorter = [(item.text(), item)
                        for ii,item in enumerate(self.iterAllItems())]
        else:
            sorter = []
            for item in self.iterAllItems():
                try:
                    sorter.append((item.itemData[sortKey],item))
                except KeyError:
                    unsortedItems.append(item)

        newOrder = [item for _,item in sorted(sorter)[::direction]]

        if unsortedItems:
            print('Could not sort all items because of missing data')

        newItems = []
        for item in newOrder + unsortedItems:
            newItems.append(self.takeItem(self.row(item)))

        for item in newItems:
            self.addItem(item)


class ProgressWindow(QDialog):
    def __init__(self, parent=None, progresses=[{'title':'Progress',
                                                 'complete':1},],
                windowTitle='Progress'):
        super(ProgressWindow, self).__init__(parent=parent)
        self.setWindowTitle(windowTitle)
        self.mainLayout = QVBoxLayout()

        self.progresses = []
        for progress in progresses:
            group = QGroupBox(progress.get('title', 'Progress'))
            progressBar = QProgressBar()
            setattr(progressBar, 'text', QLabel())
            setattr(progressBar, 'complete', progress.get('complete', 1))
            self.progresses.append(progressBar)
            _tmpLayout = QVBoxLayout()
            _tmpLayout.addWidget(progressBar)
            _tmpLayout.addWidget(progressBar.text)
            group.setLayout(_tmpLayout)
            self.mainLayout.addWidget(group)

        self.setLayout(self.mainLayout)

    def update(self, barIndex, value, text=None):
        pBar = self.progresses[barIndex]
        pVal = int(100* value / float(pBar.complete))
        pBar.setValue(pVal)
        if not text is None:
            pBar.text.setText(text)



class QHLine(QFrame):
    def __init__(self):
        super(QHLine, self).__init__()
        self.setFrameShape(QFrame.HLine)
        self.setFrameShadow(QFrame.Sunken)


class QVLine(QFrame):
    def __init__(self):
        super(QVLine, self).__init__()
        self.setFrameShape(QFrame.VLine)
        self.setFrameShadow(QFrame.Sunken)


##### TestFunctions ##########################################################
def _testFileLineEdit():
    class TestForm(QDialog):
        def __init__(self):
            super(TestForm, self).__init__()
            vLayout = QVBoxLayout()
            self.fileSel = FileLineEdit(description='select a file',
                        filter='abaqus (*.inp) ;; csv (*.inp) ;; any (*.*)')
            self.folderSel = FileLineEdit(description='select a folder',
                                     folder=True)
            vLayout.addWidget(self.fileSel)
            vLayout.addWidget(self.folderSel)
            self.setLayout(vLayout)

    app =  QApplication(sys.argv)
    form = TestForm()
    form.show()                         # Show the form
    app.exec_()
    print form.fileSel.value
    print form.folderSel.value
def _testProgressWindow():
    app =  QApplication(sys.argv)
    form = ProgressWindow(windowTitle='TestProgress',
                          progresses=[{'title':'p1', 'complete':12345},
                                      {'title':'p2', 'complete':12.345},])
    form.show()                         # Show the form

    for ii in range(12345):
        form.update(0, ii)
    for ii in range(12345):
        form.update(1, ii/1., text= 'now %d of 123' % (ii+1))
    app.exec_()
    form.close()


def _testSimpleCodeEditor():
    class TestForm(QDialog):
        def __init__(self, parent=None, initText=None):
            super(TestForm, self).__init__(parent)
            self.editor = SimpleCodeEditor(parent=self, initText=initText)
            layout = QVBoxLayout()
            layout.addWidget(self.editor)
            layout.setContentsMargins(0,0,0,0)
            self.setLayout(layout)

    app =  QApplication(sys.argv)
    form = TestForm(initText= '#SimpleTest\ndef(input):\n\tprint "HalloWorld"')
    form.show()                         # Show the form
    app.exec_()


def _testPlainParameterWidget():
    parameters = {'myA'     : {},#dict(val=12,step=.5),
                  'myBool'  : dict(value=True),
                  'myFloat'     : dict(value=.5, step=.01, upper=1, lower=0),
                  'myInt'     : dict(value=-29, step=17, upper=100, lower=-30),
                  'myText'  : dict(value='Blub'),
                  'myStringSel' : dict(value='bla', keys=['foo', 'bla']),
                  'myPath'  : dict(value=Path(__file__),)
                  }

   # parameters = {}

    class TestForm(QDialog):
        def __init__(self, parent=None, parameters={}):
            super(TestForm, self).__init__(parent)
            layout = QVBoxLayout()
            self.PlainParameters = PlainParameterWidget(parent=self,
                                                        parameters=parameters)
            self.PlainParameters2 = PlainParameterWidget(parent=self,
                                                        parameters=parameters)
            hLayout = QHBoxLayout()
            hLayout.addWidget(self.PlainParameters)
            hLayout.addWidget(self.PlainParameters2)
            layout.addLayout(hLayout)
            self.b1 = QPushButton("add data")
            self.b1.clicked.connect(self.onAddData)
            layout.addWidget(self.b1)
            self.setLayout(layout)

        def onAddData(self):
            params = self.PlainParameters.parameters
            params['new_%d' % len(params)] = {'value':'bla'}
            self.PlainParameters.updateLayout()


    app =  QApplication(sys.argv)
    form = TestForm(parameters=parameters)
    form.show()                         # Show the form
    app.exec_()
    print form.PlainParameters.parameters


def _testDragDropQListWidget():
    from random import random
    class Test_DragDropQListWidget(QDialog):
        def __init__(self):
            super(Test_DragDropQListWidget, self).__init__()

            layout = QHBoxLayout()
            self.aa = DragDropQListWidget(parent=self,
                                     items=['item%d' % d for d in range(6)],
                                     itemData=[{'rand': random()}
                                                  for d in range(6)],
                                     dragMode=1,
                                     sorter=['name', 'rand'])

            self.bb = DragDropQListWidget(parent=self,
                                          sorter=['name', 'rand'],
                                          dragMode=2,
                                          rejectItems=['item3',])

            self.bb.allowDropsFrom = [self.aa,]
            self.bb.itemChanged.connect(lambda s: self.onItemChanged())
            layout.addWidget(self.aa)
            layout.addWidget(self.bb)

            self.setLayout(layout)

        def onItemChanged(self):
            print self.bb.itemList()
            print 'changed'

    app =  QApplication(sys.argv)
    form = Test_DragDropQListWidget()
    form.show()                         # Show the form
    app.exec_()
    print form.aa.itemList()
    print form.bb.itemList()


if __name__ == '__main__':
    print('No syntax errors')
    #_testSimplePythonEditor()
    #_testPlainParameterWidget()
    #_testFileLineEdit()
    #_testDragDropQListWidget()
    _testFileLineEdit()
    #_testProgressWindow()
