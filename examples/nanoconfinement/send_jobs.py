import os
import glob

dir_loc  = 'domains_real'
send_sims = 40 # max num of sims to send


sent_sims = 0 # so far
for item_name in os.listdir(dir_loc):
    if os.path.isdir(f'{dir_loc}/{item_name}')  & \
                                (item_name.find('Psi') == -1):
                                    
        if item_name.find('480')>-1:
            continue
      
        num_sims = len( sims:= glob.glob(
                    f'{dir_loc}/{item_name}/{item_name}_uz/uz_0_*') )
        if num_sims > 1:
            print('Multiple saved states found, send help')
        elif num_sims == 1:
            it = int( sims[0].replace(
                f'{dir_loc}/{item_name}/{item_name}_uz\\uz_0_','') )
            if it == 500000:
                print('checkpointed, might want to send help \n')
                # num_sims = 0
            else:
                print('sim looks succesfull')
            
        elif num_sims == 0:
            # no sims found :(
            p = item_name.split('_')
            if len (os.listdir( glob.glob(
                f'{dir_loc}/{p[0]}_{p[1]}_{p[2]}_*{p[3]}_Psi_mfp' )[0] ))>20:
                #os.system(f'sbatch JOB_{item_name}')
                print(f'sbatch JOB_{item_name}')
                sent_sims += 1
                print(f'I have sent {sent_sims} sims')
                if sent_sims == send_sims:
                    raise Exception('Sim limit reached')
            else:
                print(f'{item_name} domains are not here')
                            
              
    
    