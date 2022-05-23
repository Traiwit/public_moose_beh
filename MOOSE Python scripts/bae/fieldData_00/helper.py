#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Some common fuctions
"""

import sys



import numpy as np

#check module is loaded by epydoc-build
#--> epyDocBuild=True to disable decorators
if 'epydoc' in sys.modules:
    epyDocBuild = True
else:
    epyDocBuild = False

#{ decorator functions
###GP: I'd prefer to replace the decorator by a verbatim call to
###GP: self._checkConsistency() in the function body
def _checkConsistencyBefore(func):
    """Calles self._checkConsistency before excecuting function"""
    if epyDocBuild:
        return func

    def func_wrapper(self, *args, **kwargs):
        self._checkConsistency()
        return func(self, *args, **kwargs)
    return func_wrapper


def _checkConsistencyAfter(func):
    """Calles self._checkConsistency after excecuting function"""
    if epyDocBuild:
        return func

    def func_wrapper(self, *args, **kwargs):
        ret = func(self, *args, **kwargs)
        self._checkConsistency()
        return ret
    return func_wrapper

#} #End decorator functions


def getBaseDataType(obj, numericStrict=False, stringStrict=False):
    if not numericStrict:
        atts = ['__add__', '__sub__', '__mul__', '__div__', '__pow__']
        #python3:
        #atts = ['__add__', '__sub__', '__mul__', '__truediv__', '__pow__']
        if all(hasattr(obj, attr) for attr in atts):
            return float

    if not stringStrict:
        if (isinstance(obj, str) or
            (isinstance(obj, np.ndarray) and obj.dtype.type is np.string_)):
            return str

    if isinstance(obj, np.ndarray):
        return obj.dtype.type

    else:
        return type(obj)

