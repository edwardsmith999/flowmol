from mdpostproc import MD_PostProc
from cfdpostproc import CFD_PostProc
from cplpostproc import CPL_PostProc
from postproc import NoResultsInDir

class All_PostProc:
    
    def __init__(self,fdir):

        self.plotlist = {}
        try:
            MD_PP = MD_PostProc(fdir)
            self.plotlist.update(MD_PP.plotlist)
            print(MD_PP)
        except:
            pass
        try:
            CFD_PP = CFD_PostProc(fdir)
            self.plotlist.update(CFD_PP.plotlist)
            print(CFD_PP)
        except:
            pass
        try:
            CPL_PP = CPL_PostProc(fdir)
            self.plotlist.update(CPL_PP.plotlist)
            print(CPL_PP)
        except:
            pass

        if (len(self.plotlist) == 0):
            raise NoResultsInDir

