import os
import glob


def is_int(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False


with open('input_perm_general.xml','r') as f:
    filedata = f.read()


all_dirs = glob.glob('input/*data_4NN*')

for domain in all_dirs:
    domain = domain.replace('input/','')
    domain = domain.replace('_data_4NN_258_128_3.dat','')
    
    if is_int(domain)==True:
        try:
            os.mkdir(f'lbm_results/{domain}')
        except:
            print('Directory already exists')
        
        new_filedata=filedata.replace('X', domain)
        
    
        with open('input_perm.xml','w') as f:
            f.write(new_filedata)
        
        os.system('bash run1-phase.sh')


