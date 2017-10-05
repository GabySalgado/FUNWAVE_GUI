import numpy as np
import subprocess
import ipywidgets as widgets
import matplotlib.pyplot as plt
import os

# import pertaining variables:
from pyFiles2.PrincipalTab_2 import  processors_text, time_text #,output_text # output directory & number of processors
from pyFiles2.PrincipalTab_4 import RunFunwave_load                          # funwave progress load bar

originalSize = 0.0

def runFunLoad(logFolder):             # loading bar function
    
    filename = os.path.join(logFolder,'LOG.txt')     
    
    # test
    print("LOG path is: ", filename)
    # test end
    
    total_time = time_text.value      # total time of model (which is given by the user)
    RunFunwave_load.min = 0.0         # min value in progress bar = 0.0
    RunFunwave_load.max = total_time  # max value in progress bar = total time of model

    originalSize = 0.0                # original size of log.txt = 0
    actual_time = 0.0                 # initial actual time of model = 0
    RunFunwave_load.value = actual_time
    
    while actual_time < total_time:
        fileStats=os.stat(filename)          
        fileSizeInBytes=fileStats.st_size     # check size in bytes of log.txt
        
        if fileSizeInBytes > originalSize:
            originalSize = fileSizeInBytes
            lastLines = subprocess.check_output(['tail', '-3', filename])    # take the final line in log.txt           

            for line in lastLines:
                if "TIME/TOTAL" in str(line):    
                    #test
                    print('last line in question: ', line)
                    # test end

                    rowDataList = line.strip().split()                      # convert to list
                    time_total_index = rowDataList.index("TIME/TOTAL:")         # find index of TIME/TOTAL         
                    actual_time = float(rowDataList[time_total_index + 1])      # actual_time's index is +1 after TIME/TOTAL's 
                    RunFunwave_load.value = actual_time                         # value of progress bar is actual_time
                    RunFunwave_load.description = "%4.2f / %4.2f" % (actual_time,total_time)
                
        
def runFUN_function(variable):
    """ this function runs the funwave model"""

    # NOTE: depth.txt, input.txt and mytvd must be in the gui folder!!!
    
    postprocessDir = str(output_text.value)    # output directory (this directory must be inside the gui folder)
    folder_dir = os.path.join(os.path.split(postprocessDir)[0],'data')    # directory where input.txt, depth.txt and mytvd are located (gui folder) 
    os.chdir(folder_dir)                       # move to folder dir (NOT NEEDED)
    
    
    EXEC = 'mytvd'
    BIN = 'FUNWAVE-TVD-master'
    HOME = os.environ['HOME']
    fun_dir = os.path.join(HOME,BIN,EXEC) # path to mytvd 
    inputDir =  os.path.join(os.path.split(postprocessDir)[0],'data','input.txt')  # directory to input.txt 
    
    # run funwave terminal command:
    run_fun = "cd %s && /usr/local/mpich/bin/mpiexec -n %d %s %s > LOG.txt &"%(folder_dir,
                                                                    processors_text.value,
                                                                    fun_dir,
                                                                    inputDir)
 
    
    os.system(run_fun) # run funwave
    runFunLoad(folder_dir)
    

def abortFun_function(variable):     
    """ this function aborts the funwave model"""
    os.system('pkill -9 mytvd-ver3.0')             # in this computer the EXEC is mytvd-ver3.0
    
    actual_time = total_time                       # to stop while loop for progress bar
    RunFunwave_load.description = 'Model Aborted'
