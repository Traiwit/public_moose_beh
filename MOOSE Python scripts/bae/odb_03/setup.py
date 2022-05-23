"""Description for the odb_03.OdbReader_ext extension module
to be imported by the global setup.py script
"""

from distutils.core import Extension
import os
from subprocess import check_output
import shutil
from glob import glob
import inspect
from collections import defaultdict

#---------------------------------------------
#--- The result of this file are these lists:
ext_modules = list()
package_data = list()
extraInstallList = list()

#-----------------------------------------
#--- Abaqus version independent module(s)

ext_modules.append(
    Extension(
        name = "bae.odb_03.odbReader_ext",
        sources = [
            'bae/odb_03/fieldValueMap.cc',
            'bae/odb_03/interpolation.cc',
            'bae/odb_03/odbReader_ext.cc',],
    )
)

#--------------------------------------------------------
#--- build stand-alone Abaqus executable odbAccessServer

def _build_odbAccessServer(install_dir):

    moduleRelDir = "bae/odb_03"


    #--- find different Abaqus versions
    abqExeNames = check_output(
        'compgen -c abq', shell=True, executable='/bin/bash')
    abqExeNames = abqExeNames.split()
    versions = dict()
    for exeName in abqExeNames:
        if exeName.startswith("abq61"):
            mainVersion = exeName[:6]  # "abq614" for exeample
            try:
                subversion = int(exeName[6:])  # "abq6143" -> 3
            except ValueError:
                subversion = 0
        elif exeName.startswith("abq20"):
            mainVersion = exeName[:7]  # "abq2018" for exeample
            try:
                subversion = int(exeName.split("hf",1)[1])  # "abq2018hf11"->11
            except (IndexError, ValueError):
                subversion = 0
        else:
            print("WARNING: Found command %s that I can't interpret as Abaqus"
                  " command." % exeName)
            continue

        try:
            oldExe, oldSubversion = versions[mainVersion]
        except KeyError:
            # nothing stored so far for this main-version
            versions[mainVersion] = [exeName, subversion]
        else:
            if subversion>oldSubversion:
                # overwrite with later version
                versions[mainVersion] = [exeName, subversion]

    # finalize (get rid of numerical subversion, diagnostic output)
    versions = sorted(
        (mainVersion, exeName)
        for mainVersion, (exeName, subversion) in versions.iteritems())
    if versions:
        print("Preparing odb_03 extension module for: %s" % ", ".join(
            version[1] for version in versions))
    else:
        print("Found no Abaqus versions.")
        return


    #--- build odbAccessServer executables for the different Abaqus versions

    # scriptHomeDir: directory where this setup.py file lies
    # this script is being executed with execfile, that's why we need inspect
    # to determine its location. __file__ would point to the main setup-script
    # that started the execfile command...
    scriptHomeDir = inspect.getfile(inspect.currentframe())

    # directory of the odbAccessServer sources
    srcDir = os.path.join(os.path.dirname(scriptHomeDir), "odbAccessServer")

    # change the current working directory to the odbAccessServer sources
    oldWd = os.getcwd()
    install_dir = os.path.abspath(install_dir)
    os.chdir(srcDir)

    # determine all source files
    sourceFiles = glob("./*.cc") + glob("./*.h")
    if not sourceFiles:
        print "WARNING: No source files for odbAccessServer."
        os.chdir(oldWd)  # change back to old working directory
        return
    newestPre = max(os.path.getmtime(fn) for fn in sourceFiles)

    # now build executable for each available Abaqus version
    for versionStr, abqExeName in versions:

        targetName = "odbAccessServer_%s" % versionStr
        targetPath = os.path.join(install_dir, moduleRelDir, targetName)

        if (os.access(targetPath, os.F_OK)
                and os.path.getmtime(targetPath)>=newestPre):
            print "%s is already up-to-date" % targetPath
            continue

        else:
            print ("Running %s make j=odbAccessServer..." % abqExeName)

            # if on CentOS cwitch on latest gcc:
            if 0==os.system(r'(cat /etc/os-release | grep "^NAME=\"CentOS")'
                            r' > /dev/null 2>&1'):
                print ("Running on CentOS, using new GCC from"
                       " 'The Software Collections (SCL)'...")
                rc = 1
                devtoolsetNb = 7
                while rc:
                    rc = os.system("scl enable devtoolset-%d '%s make j=odbAccessServer'"
                          % (devtoolsetNb, abqExeName))
                    if rc:
                        print ("Failed to run <%s make> with devtoolset-%d."
                               % (abqExeName, devtoolsetNb))
                        if devtoolsetNb<15:
                            devtoolsetNb += 1
                            print ("Now trying devtoolset-%d..." % devtoolsetNb)
                        else:
                            print ("ERROR: Failed to compile odbAccessServer.")
                            os.chdir(oldWd)  # change back to old working directory
                            return

            else:
                # no special treatment for others
                os.system("%s make j=odbAccessServer" % abqExeName)

            print ("Renaming the result to %s (cwd: %s)"
                   % (targetPath, os.getcwd()))
            # shutil.move can move the file across file system boundaries,
            # different to os.rename
            shutil.move("odbAccessServer", targetPath)
            # Special for cockatoo-deb-vb: Compile directory is on a filesystem
            # that doesn't permit chmod x. So we do it afterwards in the
            # installed directory.
            os.chmod(targetPath, 64+8+1)  # chmod a+x
            print ("Finished preparing %s." % targetName)

    # change back to old working directory
    os.chdir(oldWd)

# register the puild-procedure with the setup-process:
# append to extraInstallList
extraInstallList.append(_build_odbAccessServer)
