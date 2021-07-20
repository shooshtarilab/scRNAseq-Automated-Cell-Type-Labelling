import pandas as pd 
from os import listdir
import sys


# get a list of files and make an empty dataframe
time_files = listdir(f'./{sys.argv[1]}')
df = pd.DataFrame(columns=['dataset', 'algorithm', 'seconds'])

for file in time_files:
    # skip anything but the txt files
    if file[-4:] != '.txt':
        print(file)
        continue
    else:
        # get the dataset and algorithm names from the filename
        name = file[:-4].split('_')
        alg = name[1]
        data = name[0]

        # open the file and read the time out of it
        f = open(f'./{sys.argv[1]}/{file}', 'r')
        lines = f.readlines()
        for line in lines:
            #only read the line that has the real time
            if line[:4] != 'real':
                continue
            else:
                # strip whitespace, get the minutes and seconds, and save the total time in seconds
                line = line.strip('s\n')
                t = line.split('\t')[1]
                t = t.split('m')
                minutes = float(t[0]) 
                seconds = float(t[1])
                time = minutes * 60 + seconds
                #print(f'{file} = {time}')
        
        # add to the dataframe and close the file
        row = {'dataset': data, 'algorithm': alg, 'seconds': time }
        df = df.append(row, ignore_index=True)
        f.close()

# write to disk        
df.to_csv(f'./{sys.argv[1]}/erik_times.tsv', sep='\t', index=False)
