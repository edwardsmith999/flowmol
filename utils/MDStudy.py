from MDThread import MDThread
import multiprocessing

nruns = 18
semaphore = multiprocessing.Semaphore(8)

threads = []
for i in range(nruns):
    thread = MDThread(semaphore,i)
    threads.append(thread)
    thread.start()

for thread in threads:
    thread.join()
