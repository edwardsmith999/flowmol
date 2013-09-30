import subprocess as sp

def get_platform():
    
    # All require different commands to find their name!        
    s = sp.Popen(['hostname'],stdout=sp.PIPE).communicate()[0]
    s += sp.Popen(['domainname'],stdout=sp.PIPE).communicate()[0]
    s += sp.Popen(['dnsdomainname'],stdout=sp.PIPE).communicate()[0]
  
    if ('meflow' in s): 
        platform = 'local'
    elif ('cx1' in s):
        platform = 'cx1'
    elif ('cx2' in s): 
        platform = 'cx2'
    else:
        platform = 'local'

    return platform
