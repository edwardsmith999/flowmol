#! /usr/bin/env python
from mdrun import MDRun
import multiprocessing
import subprocess as sp
import time

class Thread(multiprocessing.Process):

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
        self.semaphore = semaphore
        self.runlist = runlist

    def run(self):

        # Perform runs per thread sequetially
        for run in self.runlist:

            #Setup run
            run.setup()

            # Check number of processors required for this run
            # and wait until all are avialable using Multiphore  
            runprocs = run.get_nprocs()
            self.semaphore.acquire(runprocs)

            # Execute and finish the run once license acquired
            run.execute(blocking=True)
            run.finish()

            # Release all licenses from Multiphore
            self.semaphore.release(runprocs)

        return
