#! /usr/bin/env python
from MDRun import MDRun
import multiprocessing
import subprocess as sp
import time

class MDThread(multiprocessing.Process):

    """ 
        When instatiated, MDThread will create a new thread of 
        sequential MD runs.
        
        Constructor arguments:
        
            semaphore - semaphore object from which MDThread will
                        request a license during execution in run()
            
            runlist - list of MDRun objects to be executed sequentially

        Example usage from a higher level:

            semaphore = multiprocessing.Semaphore(maxlicenses)

            jobs = []
            for runlist in list_of_runlists:
                job = MDThread(semaphore,runlist)
                jobs.append(job)
                job.start()

            for job in jobs:
                job.join()

    """


    def __init__(self,semaphore,runlist):

        multiprocessing.Process.__init__(self)
        self.sema = semaphore
        self.runlist = runlist

    def run(self):

        # Perform runs per thread sequetially
        for run in self.runlist:

            #Setup run directory
            run.setup_directory(existscheck=False)

            # Check number of processors required for this run
            # and wait until all are avialable using Multiphore  
            nproc = run.get_nprocs()
            self.sema.acquire(nproc)

#             time.sleep(10.0/nproc) #Priority for larger runs
#             acquired = [False for n in range(0,nproc)]

#             while False in acquired:
#                 # Attempt to get all licenses from semaphore      
#                 for n in range(0,len(acquired)):
#                     avail = self.sema.acquire(False)
#                     acquired[n] = avail
#                     #print(nproc,acquired,False in acquired)

#                 #If all licenses not available, release and wait
#                 if False in acquired:
#                     for n in acquired:
#                         if n == True:
#                             self.sema.release()
#                     acquired = [False for n in range(0,nproc)]
#                     time.sleep(5.0)

            run.execute(blocking=True)
            run.finish()
            #run.post_process()

#             for n in range(0,nproc): 
#                 self.sema.release()

            # Release all licenses from MultiPhore
            self.sema.release(nproc)

        return
