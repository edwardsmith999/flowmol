from MDRun import MDRun
import multiprocessing
import subprocess as sp
#import time

class MDThread(multiprocessing.Process):

    """ 
        When instatiated, MDThread will create a new thread of 
        sequential MD runs.
        
        Constructor arguments:
        
            semaphore - semaphore object from which MDThread will
                        request a license during execution in run()
            
            info - placeholder

        Example usage from a higher level:

            semaphore = multiprocessing.Semaphore(nprocs)
            threads = []
            for info in range(infolist):
                thread = MDThread(semaphore,info)
                threads.append(thread)
                thread.start()

            for thread in threads:
                thread.join()

    """


    def __init__(self,semaphore,howtosetuprunlist):

        multiprocessing.Process.__init__(self)

        srcdir='./../MD_dCSE/src_code'
        basedir='base'
        rundir='run'+str(howtosetuprunlist)
        executable='./creamsoda'
        inputfile='./rootbeer.in'
        outputfile='./rootbeer.out'
        restartfile=None #'./r0'
        cylinderfile=None
        appendoutput = False
        self.sema = semaphore

        run1 = MDRun(srcdir,basedir,rundir,executable,inputfile,'out1')
        run2 = MDRun(srcdir,rundir,rundir,executable,inputfile,'out2',
                     restartfile='results/final_state')

        self.runlist = [run1,run2]

    	#self.runlist = []
    	#for info in howtosetuprunlist:
    	#	self.runlist.append(extractedinfo)

    def run(self):

        # Get license from semaphore       
        self.sema.acquire()

        for run in self.runlist:

            run.setup_directory(existscheck=False)
            run.change_inputs({})
            run.execute(blocking=True)
            #Run.post_process()

        # Release license from semaphore       
        self.sema.release()

        return
