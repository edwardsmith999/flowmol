import sys
import argparse
import os

import postproclib as ppl
from misclib import Chdir

if __name__ == "__main__":

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
    parser = argparse.ArgumentParser(description='run_vmd vs. master jay -- Runs VMD with overlayed field')
    parser.add_argument('d','-fdir', nargs='?', help='Directory with vmd file and field files', default='../MD_dCSE/src_code/results/')
    #parser.add_argument('-h','--help', help='Field type',  default='vbins')
    parser.add_argument('-f','--field', help=print_fieldlist(),  default=None)
    parser.add_argument('-c','--comp', help='Component name', default=None)
    args = vars(parser.parse_args())

    #Static arguments
    ppObj = ppl.MD_PostProc(args['fdir'])
    if args['field'] == None:
        print("No field type specified -- available field types include:")
        print(ppObj)
        print("Using default value of vbins")
        args['field'] = 'vbins'
        component = 0
    elif(sys.argv[1] in ['--help', '-help', '-h']):
        print("Available field types include")
        print(ppObj)
        sys.exit()
    if args['comp'] == None:
        print("No components direction specified, setting default = 0")
        args['comp'] = 0

    try:
        fobj = ppObj.plotlist[args['field']]
    except KeyError:
        print("Field not recognised -- available field types include:")
        print(ppObj)
        sys.exit()
    except:
        raise

    vmdobj = ppl.VMDFields(fobj,args['fdir'])
    vmdobj.write_vmd_header()
    vmdobj.write_vmd_intervals()
    vmdobj.write_dx_range(component=args['comp'])
    vmdobj.writecolormap('RdYlBu')
    with Chdir(args['fdir'] + './vmd/'):
        print(args['fdir'])
        command = "vmd -e " + "./plot_MD_field.vmd"
        os.system(command)
