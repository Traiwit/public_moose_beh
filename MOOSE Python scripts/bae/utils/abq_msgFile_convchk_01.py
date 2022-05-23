"""
parsing Abq/Stnd msg-file for convergence checks.
output looks like this::

  STEP INCR ATTP ITER       rmax      cmax      cest     dumax      qbar    qtilde  qbar/qtilde  rmax/qtilde <=Rn <eps <=Rl <=RP  cmax/dumax <=Cn <=Ce  cest/dumax <=Cn  CHK1 CHK2 CHK3 CHK4
     1    1    1    1  +1.00E+00 +1.00E+00 +1.00E+00 +1.00E+00 +1.00E+00 +1.00E+00    +1.00E+00    +1.00E+00 FAIL FAIL FAIL FAIL   +1.00E+00 FAIL FAIL   +1.00E+00 FAIL  FAIL  --  FAIL  --
"""

_version_history_ = """
    201x-xx-xx  AF  initial script development
    2016-11-24  AF  added optional output file object to write to
"""

import re
from bae.log_01 import msg


class AbqStnd_MsgFile_Parser(object):

    def __init__(self):

        self.resetStepVars()
        self.resetIterVars()


    def resetStepVars(self):
        self.kStep = None
        self.kIncr = None
        self.kAttp = None

        self.cp01Rn = 0.005
        self.cp02Cn = 0.010
        #self.cp03q0 = None
        #self.cp04qu = None
        self.cp05Rp = 0.020
        self.cp06eps = 1.00E-05
        self.cp07Ceps = 1.00E-05
        self.cp08Rl = 1.00E-08
        #self.cp09Cf = 1.0
        #self.cp10epsl = 1.00E-05
        #self.cp11epsd = None


    def resetIterVars(self):
        self.kIter = None
        self.rmax = None
        self.cmax = None
        self.cest = None
        self.dumax = None
        self.qmax = None
        self.qtildemax = None
        self.qbar = None
        self.qtilde = None


    def parse(self, msgFilename, output=None):

        reStrFloat = r"([+-.0-9E]*)"
        reStrInt   = r"([0-9]+)"
        reNewStep = re.compile(r"                        S T E P\s+"+reStrInt)
        #reNewStep = re.compile(r"\s*STEP\s+"+reStrInt+r"\s+INCREMENT\s+"+reStrInt)
        reNewIncr = re.compile(r"\s*INCREMENT\s+"+reStrInt+
                               r"\s+STARTS\. ATTEMPT NUMBER\s+"+reStrInt)
        reNewIter = re.compile(r"\s*CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION\s+"+reStrInt)
        reAvgVolFluxes = re.compile(r"\s*AVERAGE VOL\. FLUX\s+"+reStrFloat+
                                    r"\s+TIME AVG\. VOL\. FLUX\s+"+reStrFloat)
        reRmax  = re.compile(r"\s*LARGEST RESIDUAL VOL\. FLUX\s+"+reStrFloat)
        reDumax  = re.compile(r"\s*LARGEST INCREMENT OF P\.PRESS\.\s+"+reStrFloat)
        reCest  = re.compile(r"\s*ESTIMATE OF P\.PRESS\. CORRECTION\s+"+reStrFloat)
        reCmax = re.compile(r"\s*LARGEST CORRECTION TO P\.PRESS\.\s+"+reStrFloat)

        fr = open(msgFilename,"rb")
        self._writeFormattedOutputLine(True,output)
        for line in fr:

            chkNewStep = reNewStep.search(line)
            if chkNewStep:
                newStep = int(chkNewStep.group(1))
                if self.kStep is not None and self.kStep<>newStep:
                    #msg("  => write data for STEP %d INCR %d ATTP %d LAST ITER %d\n"
                    #    %(self.kStep, self.kIncr, self.kAttp, self.kIter) )
                    #    #%(self.kStep, 0, 0, 0) )
                    self._writeFormattedOutputLine(False,output)
                    self.resetStepVars()
                    self.resetIterVars()
                self.kStep = newStep
                msg("- start reading data for STEP %d"
                    %(self.kStep) )

            if self.kIncr is None:
                pass    # check for step data

            chkNewIncr = reNewIncr.search(line)
            if chkNewIncr:
                self.kIncr = int(chkNewIncr.group(1))
                self.kAttp = int(chkNewIncr.group(2))
                #msg("  INFO: a new increment INCR %d ATTP %d"
                #    %(self.kIncr, self.kAttp))

            chkNewIter = reNewIter.search(line)
            if chkNewIter:
                newIter = int(chkNewIter.group(1))
                if self.kIncr is None:
                    self.kIncr = 1
                if self.kAttp is None:
                    self.kAttp = 1
                if self.kIter is not None and self.kIter<>newIter:
                    #msg("  => write data for STEP %d INCR %d ATTP %d ITER %d"
                    #    %(self.kStep, self.kIncr, self.kAttp, self.kIter) )
                    self._writeFormattedOutputLine(False,output)
                    self.resetIterVars()
                #msg("        INFO: a new iteration ITER %d"%newIter)
                self.kIter = newIter

            chkAvgVolFluxes = reAvgVolFluxes.search(line)
            if chkAvgVolFluxes:
                self.qbar = float(chkAvgVolFluxes.group(1))
                self.qtilde = float(chkAvgVolFluxes.group(2))

            chkRmax = reRmax.search(line)
            if chkRmax:
                self.rmax = float(chkRmax.group(1))

            chkCmax = reCmax.search(line)
            if chkCmax:
                self.cmax = float(chkCmax.group(1))

            chkCest = reCest.search(line)
            if chkCest:
                self.cest = float(chkCest.group(1))

            chkDumax = reDumax.search(line)
            if chkDumax:
                self.dumax = float(chkDumax.group(1))


        fr.close()
        #msg("  => write data for  L A S T  STEP %d INCR %d ATTP %d ITER %d\n"
        #    %(self.kStep, self.kIncr, self.kAttp, self.kIter) )
        self._writeFormattedOutputLine(False,output)


    def _writeFormattedOutputLine(self, headerOnly=False, outObj=None):
        """
output format::

  STEP INCR ATTP ITER       rmax      cmax      cest     dumax      qbar    qtilde  qbar/qtilde  rmax/qtilde <=Rn <eps <=Rl <=RP  cmax/dumax <=Cn <=Ce  cest/dumax <=Cn  CHK1 CHK2 CHK3 CHK4
     1    1    1    1  +1.00E+00 +1.00E+00 +1.00E+00 +1.00E+00 +1.00E+00 +1.00E+00    +1.00E+00    +1.00E+00 FAIL FAIL FAIL FAIL   +1.00E+00 FAIL FAIL   +1.00E+00 FAIL  FAIL  --  FAIL  -- 
"""
        def printValue(fmt,val):
            if val is not None:
                return fmt%val
            else:
                return "NONE"

        def printCheck(val):
            if val is None:
                return " -- "
            else:
                if val:
                    return " OK "
                else:
                    return "FAIL"

        output = [
            ("STEP", "%4s"%printValue("%4d",self.kStep)),
            ("INCR", "%4s"%printValue("%4d",self.kIncr)),
            ("ATTP", "%4s"%printValue("%4d",self.kAttp)),
            ("ITER", "%4s"%printValue("%4d",self.kIter)),
            (" ", " "),
            #("     rmax", "%9s"%printValue("%.2e",self.rmax)),
            #("     cmax", "%9s"%printValue("%.2e",self.cmax)),
            #("     cest", "%9s"%printValue("%.2e",self.cest)),
            #("    dumax", "%9s"%printValue("%.2e",self.dumax)),
            ("     qbar", "%9s"%printValue("%.2e",self.qbar)),
            ("   qtilde", "%9s"%printValue("%.2e",self.qtilde)),
            #(" ", " "),
            ]
        zeroflxRatio = None
        zeroflxPrint = "N/A"
        zeroflxCheck = None
        if self.qbar is not None and self.qtilde is not None:
            zeroflxRatio = self.qbar/self.qtilde
            if abs(zeroflxRatio)<self.cp06eps:
                zeroflxPrint = "ZEROVF"
                zeroflxCheck = True
            else:
                zeroflxPrint = "NRML  "
                #zeroflxPrint = "NORMAL"
                zeroflxCheck = False
        output.extend([
            ("qbar/qtilde", "%11s"%printValue("%.2e",zeroflxRatio)),
            ("<eps  ", "%6s"%zeroflxPrint),
            (" ", " "),
            ])
        rmaxqtildeRatio = None
        rmaxqtildeChk1 = None
        rmaxqtildeChk2 = None
        rmaxqtildeChk3 = None
        rmaxqtildeChk4 = None
        if self.rmax is not None and self.qtilde is not None:
            rmaxqtildeRatio = self.rmax/self.qtilde
            rmaxqtildeRatioAbs = abs(rmaxqtildeRatio)
            rmaxqtildeChk1 = (rmaxqtildeRatioAbs<=self.cp01Rn)
            rmaxqtildeChk2 = (rmaxqtildeRatioAbs<=self.cp06eps)
            rmaxqtildeChk3 = (rmaxqtildeRatioAbs<=self.cp08Rl)
            rmaxqtildeChk4 = (rmaxqtildeRatioAbs<=self.cp05Rp)
        output.extend([
            ("     rmax", "%9s"%printValue("%.2e",self.rmax)),
            ("rmax/qtilde", "%11s"%printValue("%.2e",rmaxqtildeRatio)),
            ("<=Rn", printCheck(rmaxqtildeChk1)),
            ("<eps", printCheck(rmaxqtildeChk2)),
            ("<=Rl", printCheck(rmaxqtildeChk3)),
            ("<=Rp", printCheck(rmaxqtildeChk4)),
            (" ", " "),
            ])
        cmaxdumaxRatio = None
        cmaxdumaxChk14 = None
        cmaxdumaxChk2  = None
        if self.cmax is not None and self.dumax is not None:
            cmaxdumaxRatio = self.cmax/self.dumax
            cmaxdumaxRatioAbs = abs(cmaxdumaxRatio)
            cmaxdumaxChk14 = (cmaxdumaxRatioAbs<=self.cp02Cn)
            cmaxdumaxChk2  = (cmaxdumaxRatioAbs<=self.cp07Ceps)
        output.extend([
            ("     cmax", "%9s"%printValue("%.2e",self.cmax)),
            #("     cest", "%9s"%printValue("%.2e",self.cest)),
            ("    dumax", "%9s"%printValue("%.2e",self.dumax)),
            ("cmax/dumax", "%10s"%printValue("%.2e",cmaxdumaxRatio)),
            ("<=Cn", printCheck(cmaxdumaxChk14)),
            ("Ceps", printCheck(cmaxdumaxChk2)),
            (" ", " "),
            ])
        cestdumaxRatio = None
        cestdumaxChk1  = None
        if self.cest is not None and self.dumax is not None:
            cestdumaxRatio = self.cest/self.dumax
            cestdumaxRatioAbs = abs(cestdumaxRatio)
            cestdumaxChk1 = (cestdumaxRatioAbs<=self.cp02Cn)
        output.extend([
            ("     cest", "%9s"%printValue("%.2e",self.cest)),
            ("cest/dumax", "%10s"%printValue("%.2e",cestdumaxRatio)),
            ("<=Cn", printCheck(cestdumaxChk1)),
            ("  ", "  "),
            ])
        chkCondition1 = None
        chkCondition2 = None
        chkCondition3 = None
        chkCondition4 = None
        if not zeroflxCheck:
            if cestdumaxChk1 is not None:
                chkCondition1 = rmaxqtildeChk1 and (cmaxdumaxChk14 or cestdumaxChk1)
            else:
                chkCondition1 = rmaxqtildeChk1 and  cmaxdumaxChk14
        else:
            chkCondition2 = rmaxqtildeChk2 or cmaxdumaxChk2
        chkCondition3 = rmaxqtildeChk3
        if self.kIter>=9:
            chkCondition4 = rmaxqtildeChk4 and cmaxdumaxChk14
        output.extend([
            ("CHK1", printCheck(chkCondition1)),
            ("CHK2", printCheck(chkCondition2)),
            ("CHK3", printCheck(chkCondition3)),
            ("CHK4", printCheck(chkCondition4)),
            ])


        if headerOnly:
            if outObj is None:
                msg(" ".join([ out[0] for out in output ]))
            else:
                outObj.write(" ".join([ out[0] for out in output ])+"\n")
        else:
            if outObj is None:
                msg(" ".join([ out[1] for out in output ]))
            else:
                outObj.write(" ".join([ out[1] for out in output ])+"\n")

if __name__=="__main__":
    cc = AbqStnd_MsgFile_Parser()
    cc.parse("G01_R05H/HCAVE2016_R05HicCPv0_G01H_Q03CAVE00_M06H100_D00H.msg")
