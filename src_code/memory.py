import psutil
import time
import sys
import subprocess as sp
import getpass
import itertools
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

command = sys.argv[1]
whoami = getpass.getuser() 
plist = psutil.get_process_list()

#rotatechars = itertools.cycle(['-','\\','|','/'])
#rotatechars = itertools.cycle(['+','x'])
#rotatechars = itertools.cycle(['_','-'])
barlist = []
maximum = 40 
mouth = itertools.cycle(['<']*10+['-']*10)
for i in range(maximum):
    frame =  '|' + ' '*i + r'   _\\\_  ' + ' '*(maximum-i) + '|\n'
    frame += '|' + ' '*i + r' \/    o\ ' + ' '*(maximum-i) + '|\n'
    frame += '|' + ' '*i + r' /\._((_X ' + ' '*(maximum-i) + '|\n'
    frame += '|' + ' '*i + r'     >    ' + ' '*(maximum-i) + '|\n'
    barlist.append(frame.replace('X',mouth.next()))
mouth = itertools.cycle(['>']*10+['-']*10)
for i in range(maximum,0,-1):
    frame =  '|' + ' '*i + r'  _///_   ' + ' '*(maximum-i) + '|\n'
    frame += '|' + ' '*i + r' /o    \/ ' + ' '*(maximum-i) + '|\n'
    frame += '|' + ' '*i + r' X ))_./\ ' + ' '*(maximum-i) + '|\n'
    frame += '|' + ' '*i + r'    <     ' + ' '*(maximum-i) + '|\n'
    barlist.append(frame.replace('X',mouth.next()))
rotatechars = itertools.cycle(barlist)     

watchlist = []
for p in plist:
    if command in p.name() and whoami == p.username(): 
        watchlist.append(p)

if (len(watchlist) == 0):
    quit('Process not found')

print('Found {0:d} processes owned by {1:s} with name containing "{2:s}":\n'.format(len(watchlist),whoami,command))
for p in watchlist:
    pid = p.pid
    name = p.name()
    print('\t{0:d}\t{1:s}'.format(pid,name))

message = ''
while True:

    rss_total = 0
    vss_total = 0
    
    spinner = rotatechars.next()
    prevmessage = message

    try:

        nlines = message.count('\n')
        message = ("\033[F"*nlines)
        
        message += ('\n'+spinner+'\n')
        message += ('\n{0:>10s}:{1:>20s}{2:>20s}\n'.format('PID','RSS','VSS'))
        for p in watchlist:
            rss, vss = p.get_memory_info()
            rss_total += rss
            vss_total += vss
            message += ('\n{0:10d}:{1:20d}{2:20d}'.format(p.pid, rss, vss))
        message += ('\n\n{0:>10s}:{1:20d}{2:20d}\n'.format('TOTAL', rss_total, vss_total))
        message += ('\n'+'-'*51+'\n')
        message += '\n'

        sys.stdout.write(message);
        time.sleep(0.1)

    except psutil.NoSuchProcess:
   
        quit('Process no longer exists.')
