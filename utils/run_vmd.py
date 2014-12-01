import sys
import argparse
import os

import postproclib as ppl
from misclib import Chdir

def run_vmd(parent_parser=argparse.ArgumentParser(add_help=False)):

    def print_fieldlist():
        outstr = 'Type of field to overlay with vmd \n'
        try:
            ppObj = ppl.MD_PostProc('../MD_dCSE/src_code/results/')
            outstr = outstr + str(ppObj)
        except:
            print(' \n')
            pass

        outstr = outstr + '\n N.B. Make sure to include quotes if there is a space in field name \n'
        return outstr

    #Keyword arguments
    parser = argparse.ArgumentParser(description='run_vmd vs. master jay -- Runs VMD with overlayed field',
                                     parents=[parent_parser])

    try:
        argns, unknown = parser.parse_known_args()
        print('Using directory defined as ', argns.fdir)
    except AttributeError:
        parser.add_argument('-d','--fdir',dest='fdir', nargs='?', 
                            help='Directory with vmd file and field files',     
                            default=None)

    parser.add_argument('-f', '--field', dest='field', 
                        help=print_fieldlist(), default=None)
    parser.add_argument('-c', '--comp', dest='comp', 
                        help='Component name', default=None)
    args = vars(parser.parse_args())

    #Static arguments
    if args['fdir'] == None:
        scriptdir = os.path.dirname(os.path.realpath(__file__))
        args['fdir'] = scriptdir + '/../MD_dCSE/src_code/results/'
    if args['field'] == None:
        print("No field type specified -- using default value of no field")
        args['field'] = None
        component = 0
    elif(sys.argv[1] in ['--help', '-help', '-h']):
        ppObj = ppl.MD_PostProc(args['fdir'])
        print("Available field types include")
        print(ppObj)
        sys.exit()
    if args['comp'] == None:
        print("No components direction specified, setting default = 0")
        args['comp'] = 0

    #Plane field case
    if args['field'] == None:
        fobj = ppl.MD_dummyField(args['fdir'])
        fobj.plotfreq = 1
        vmdobj = ppl.VMDFields(fobj,args['fdir'])
        vmdobj.reformat()

        #Open vmd 
        with Chdir(args['fdir']):
            command = "vmd " + "./vmd_out.dcd"
            os.system(command)
    #Overlayed field case
    else:
        try:
            ppObj = ppl.MD_PostProc(args['fdir'])
            fobj = ppObj.plotlist[args['field']]
        except KeyError:
            print("Field not recognised -- available field types include:")
            print(ppObj)
            sys.exit()
        except:
            raise

        vmdobj = ppl.VMDFields(fobj,args['fdir'])
        vmdobj.reformat()
        vmdobj.write_vmd_header()
        vmdobj.write_vmd_intervals()
        vmdobj.write_dx_range(component=args['comp'])
        vmdobj.writecolormap('RdYlBu')

        #Open vmd 
        with Chdir(args['fdir'] + './vmd/'):
            command = "vmd -e " + "./plot_MD_field.vmd"
            os.system(command)

if __name__ == "__main__":
    run_vmd()
