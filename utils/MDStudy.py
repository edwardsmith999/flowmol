import multiprocessing
import subprocess as sp

from Platform import get_platform
from MDThread import MDThread
from Dummy import DummySemaphore
from MultiPhore import MultiPhore

class MDStudy:

    def __init__(self,threadlist,maxproc=None):

        """
            A single study of multiple MDThreads, each carrying
            the responsibility of executing a series of sequential
            MDRun objects. 
 
            Constructor arguments:
            
                threadlist - a list of runlists. Each runlist contains a
                             series of MDRun objects to be executed 
                             sequentially
                
                maxproc - maximum number of licenses the semaphore may
                          release at once. 

        """

        #Get semaphore is possible, otherwise use dummy routine
        if (get_platform() == 'local'):

            # Get number of cpus on computer if not specified
            if maxproc == None:
                print("Maximum number concurrent processors not specified, " +
                      "attempting to use total number of CPUs...")
                try:
                    ncpus = multiprocessing.cpu_count()
                except NotImplementedError:
                    raise

            self.semaphore = MultiPhore(maxproc)

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
