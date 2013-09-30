import multiprocessing
import subprocess as sp

from Platform import get_platform
from MDThread import MDThread
from Dummy import DummySemaphore

class MDStudy:

    def __init__(self,threadlist,maxlicenses):

        """
            A single study of multiple MDThreads, each carrying
            the responsibility of executing a series of sequential
            MDRun objects. 
 
            Constructor arguments:
            
                threadlist - a list of runlists. Each runlist contains a
                             series of MDRun objects to be executed 
                             sequentially
                
                maxlicenses - maximum number of licenses the semaphore may
                              release at once. This is likely to be best 
                              set equal to ncpus / nprocsperrun.

        """

        if (get_platform() == 'local'):
            self.semaphore = multiprocessing.Semaphore(maxlicenses)
        else:
            print('Semaphore not available, creating dummy instead.')
            self.semaphore = DummySemaphore()

        jobs = []
        for runlist in threadlist:
            job = MDThread(self.semaphore,runlist)
            jobs.append(job)
            job.start()

        for job in jobs:
            job.join()
