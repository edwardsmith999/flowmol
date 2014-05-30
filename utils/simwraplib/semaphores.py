#! /usr/bin/env python2.7
from multiprocessing import BoundedSemaphore
import time 

class DummySemaphore:
    
    def acquire(self, nprocs):
        return 
    
    def release(self, nprocs):
        return

class Multiphore():

    """ 
        Generalisation of the Semaphore to ensure threads calling multiple cpus
        subproccesses only run when they can acquire all required licenses for 
        number of cpus. Corresponding multiple release also provided.

    """

    def __init__(self,maxprocs):

        self.maxprocs = maxprocs
        self.sema = BoundedSemaphore(maxprocs)

    def acquire(self,nproc=1,blocking=True,wait=5.0):

        #Check requested processes not greater than maximum
        if nproc > self.maxprocs:
            raise ValueError("Requested acquire greater than semaphore maximum")

        #Priority for larger runs
        time.sleep(wait/nproc) 

        #Keep looping until all cpus have a license
        acquired = [False for n in range(0,nproc)]
        while False in acquired:

            # Attempt to get all licenses from semaphore      
            for n in range(0,len(acquired)):
                avail = self.sema.acquire(False)
                acquired[n] = avail
                #print(nproc,acquired,False in acquired)

            if blocking == False:
                return acquired

            #If all licenses not available, release all acquired and wait
            if False in acquired:
                for n in acquired:
                    if n == True:
                        self.sema.release()
                acquired = [False for n in range(0,nproc)]
                time.sleep(wait)


    def release(self,nproc=1):

        # Release all licenses from semaphore
        for n in range(0,nproc): 
            self.sema.release()

    def get_value(self):

        return self.sema.get_value()
