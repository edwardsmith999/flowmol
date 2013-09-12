from MDThread import MDThread
import multiprocessing

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

        self.semaphore = multiprocessing.Semaphore(maxlicenses)

        jobs = []
        for runlist in threadlist:
            job = MDThread(self.semaphore,runlist)
            jobs.append(job)
            job.start()

        for job in jobs:
            job.join()
