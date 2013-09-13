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

        # Get license from semaphore       
        self.sema.acquire()

        # Perform runs sequetially
        for run in self.runlist:

            run.setup_directory(existscheck=False)
            run.execute(blocking=True)
            run.finish()
            #run.post_process()

        # Release license from semaphore       
        self.sema.release()

        return
