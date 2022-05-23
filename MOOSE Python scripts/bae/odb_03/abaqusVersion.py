"""Module to determine the available Abaqus versions and things like that

Primarily for internal use.
"""
__version__ = "1.02"
_version_history_ = """
Versions:

1.0 GP extracted from odb_03/__init__.py
1.01 GP added the new abq2018-style
1.02 GP also accept abq614 as version
"""
from bae.log_01 import msg

def getAbqVersionStr(version=None):
    """Get Abaqus version string, something like "abq613".

    @param version: string as of the subdirectory for Abaqus like "6.13-2"
        For new style use something like "abq2018", ignoring the hot-fix
        postfix e.g. "hf4".
        If not specified then get it from what currently gets called with the
        shell command "abaqus".
    @returns: string identifying the version of the compiled extension module
        like "abq613" (Note that the minor number "-2" is ignored.)
    """

    if version is None:
        from subprocess import Popen, PIPE
        import os
        msg("getAbqVersionStr: determining path to abaqus exe", debugLevel=10)
        basepath = Popen(["which", "abaqus"], stdout=PIPE).communicate()[0]
        basepath = basepath.strip()  # line breaks seem to confuse realpath
        msg("getAbqVersionStr: path to abaqus exe: %s" % basepath,
            debugLevel=10)

        path = os.path.realpath(basepath)
        found = False
        version = 'dummy'
        while path and version and not found:
            path, version = os.path.split(path)
            found = path.endswith("abaqus") and version.startswith("6.")
            #### DEBUG
            # print(
            #     "Trying to find the base directory of the Abaqus"
            #     " installation. Current path %s, assuming I found it: %s"
            #     % (path, found))
        if not found:
            path = basepath
            while 1:
                version = os.path.basename(path)
                if version.startswith("abq"):
                    found = True
                    break
                if not os.path.islink(path):
                    found = False
                    break

                # for next iteration follow link one step
                newpath = os.readlink(path)
                if not os.path.isabs(newpath):
                    newpath = os.path.join(os.path.dirname(path), newpath)
                path = newpath

        if not found:
            raise ValueError(
                "Could not properly determine current Abaqus version.")
        msg("getAbqVersionStr: determined version %s" % version, debugLevel=10)

    # process 6.13-2 - style version string
    if "." in version:
        mainVersion = version.split("-",1)[0]
        abqVersionStr = "abq%s%s" % tuple(mainVersion.split("."))
    # if something like abq2018hf4 then strip off the "hf" part
    elif version.startswith("abq20"):
        abqVersionStr = version.split("hf")[0]
    elif version.startswith("abq61"):
        abqVersionStr = version[:6]
    else:
        abqVersionStr = version

    return abqVersionStr


def getAbqExe(version=None):
    """Get the path of the Abaqus executable, something like
    "/usr/abaqus/6.14-3/code/bin/abq6143".

    @param version: string as of the subdirectory for Abaqus like "6.13-2"
        For new style use something like "abq2018", ignoring the hot-fix
        postfix e.g. "hf4".
        If not specified then get it from what currently gets called with the
        shell command "abaqus".
    @returns: Abaqus executable -full path- of the specified version.
    """

    from subprocess import Popen, PIPE
    import os.path
    from os import environ, defpath, access, X_OK
    from glob import glob

    def which(exe):
        """try to find an executable using 'which'
        @returns: the full path of the executable or None
        """
        # Note: function takes variable 'version' from context
        path = None
        p = Popen(["which", exe], stdout=PIPE, stderr=PIPE)
        path, errmsg = p.communicate()
        if p.returncode!=0:
            path = None
            msg("Could not find an Abaqus executable for version %s:\n%s"
                % (version, errmsg), debugLevel=5)
        return path
    
    # try to find an executable...
    path = None  # not needed, just for peace of mind
    if version is None:
        # if version is not given try to find the executable "abaqus"
        path = which("abaqus")
    elif "." in version:
        # find Abaqus executable up to 6.14-3
        # according to the old convention: version="6.14-3" -> abq6143
        path = which("abq%s" % version.replace(".", "").replace("-", ""))
    else:
        # find an executable that is named exactly like version (e.g.
        # "abq2018hf6")
        path = which(version)
        if not path and "hf" in version:
            # try again with "hf" stripped off
            path = which(version.split("hf")[0])

    # If none of the above worked then assume a version starting with "abq".
    # Search the PATH environment variable to search for Abaqus executables
    # like abq2017, abq2018hf4, and abq6132.
    if not path and version.startswith("abq"):

        # path = directories that contain executables
        try:
            dirList = environ["PATH"]
        except KeyError:
            dirList = defpath

        # find matching executables
        for dirname in dirList.split(":"):
            pathList = [
                x for x in glob(os.path.join(dirname, version+"*"))
                if access(x, X_OK)]
            try:
                path = pathList[0]
            except IndexError:
                path = None
                continue
            else:
                break

    # diagnostic output
    if path:
        msg("Found executable for Abaqus version %s:\n%s"
            % (version, path), debugLevel=5)
    else:
        raise ValueError(
            "Don't know what to do with Abaqus version %s" % version)

    path = path.strip()  # some line breaks or so seem to confuse realpath
    path = os.path.realpath(path)
    return path
