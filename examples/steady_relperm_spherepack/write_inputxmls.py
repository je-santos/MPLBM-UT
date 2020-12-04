from shutil import copyfile


dom_size = 200
saturations = [.1, .2, .3, .4, .5, .6, .7, .8, .9]

import fileinput
import sys
import os

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)


for sat in saturations:
    
    s_nw = int(sat*100)
    fname = f'input_spherepack_S{s_nw}.xml'
    copyfile('input_spherepack_S25.xml',fname)
    replaceAll(fname,'050',f'{int(dom_size*sat)}')
    replaceAll(fname,'051',f'{int(dom_size*sat)+1}')
    replaceAll(fname,'0.00001',f'{0.00001*10}')
    replaceAll(fname,'tmp_25',f'tmp_{s_nw}')
    replaceAll(fname,'<x> true </x> <y> true </y> <z> true </z>',
                     '<x> true </x> <y> false </y> <z> false </z>')
    try:
        os.mkdir(f'tmp_{s_nw}')
    except:
        pass
    
    bname = f'batchjob_{s_nw}.sh'
    copyfile('batchjob_4TACC.sh', bname)
    replaceAll(bname, '4', '10')
    replaceAll(bname, 'S25', f'S{s_nw}')
    replaceAll(bname, '108', '12')
    #replaceAll(bname, 'pge-fracture', f'S{s_nw}')
    
    
        
    