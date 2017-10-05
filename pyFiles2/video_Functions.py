import numpy as np
import ipywidgets as widgets
import matplotlib.pyplot as plt
from IPython.display import display, clear_output
import os

# import pertaining variables: 
from pyFiles2.PRINCIPAL_TAB import principal_tab,title_text     # import tab structure & projects name
from pyFiles2.PrincipalTab_1 import THL as Lt                   # total horizontal length
from pyFiles2.PrincipalTab_2 import time_text, plotInt_text     # total time and plot interval in (s)
from pyFiles2.PrincipalTab_5 import TimeLimit_text,video_load,plot_time    # number of eta files and video loadbar widget


def ffmpeg_run(postprocessDir):
    #os.chdir(postprocessDir) # move to directory where the eta images are located (output directory)
    os.system("ffmpeg -r 12 -y -i surf_%5d.png -s 815x735 surfaceMovie.mp4") # create video with ffmpeg terminal command
    
    
def readETAData(num,numOfETA,postprocessDir,pts):
    
    """Function that reads in the eta_00.. ASCII file and returns
    it in the 2D NumPy array of size [numOfRows,pts]."""

    assert num in range(1,numOfETA+1), "File Index:%d is not in range of station numbers." %(num,)
    if num < 10:
        fileIndex = str(num)
        fileName = postprocessDir+'eta_0000{0:s}'.format(fileIndex)

    elif num <100 and num > 9:
        fileIndex = str(num)
        fileName = postprocessDir+'eta_000{0:s}'.format(fileIndex)

    elif num < 1000 and num>99: 
        fileIndex = str(num)
        fileName = postprocessDir+'eta_00{0:s}'.format(fileIndex)

    else: 
        fileIndex = str(i)
        fileName = postprocessDir+'eta_0{0:s}'.format(fileIndex)

    fin = open(fileName, "r")
    lines = fin.readlines()
    
    numpyFileArray = np.zeros((len(lines),pts))

    for lineIndex,line in enumerate(lines):
        rowDataList = line.strip().split()          # dropping all whitespace (leading & trailing)
        numpyData = np.array(rowDataList).astype('float')
        numpyFileArray[lineIndex,:] = numpyData[:]

    fin.close()
    return numpyFileArray

    
def runVID_function(variable):
    
    """Function that plots the eta files""" 
    
    # move to directory where depth file is located
    pwd = os.getcwd()  # get current path
    s = title_text.value 
    wrd_lst = s.split()
    space = '_'
    folder_name = space.join(wrd_lst) # substitute ' ' space in project title with '_'
    
    postprocessDir=os.path.join(pwd,folder_name+'/output/') # create output folder path
    #os.chdir(postprocessDir+'..')
    
    Depthtext = os.path.join(pwd,folder_name+'/depth.txt')
    depth = (np.loadtxt(Depthtext))*-1    # Change Depth to (-) under MWL and (+) over MWL. 
    pts  = len(depth.T)                   # num of Colums = number of points 
    
    
    data_text = os.path.join(pwd,folder_name+'/data.txt') # create path to open data.txt in project folder
    fin= open(data_text,'r')  # upload data.txt; which has the bathy dx data
    val = fin.read()            
    val=val.split()
    dx = float(val[5])
    
    Lt = np.ceil(dx*(pts+1)) # compute total horizonta length
    
    numOfRows = len(depth)
    numOfETA = int(TimeLimit_text.value/plotInt_text.value) + 1
    
    video_load.max = numOfETA
    plotInt = plotInt_text.value  
    
    
    # data placeholder
    eta = np.zeros((numOfRows,pts,numOfETA))
    freeSurface = np.zeros((numOfRows,numOfETA))
   
    # Plot and save every ETA file:
    for i in range(1, numOfETA+1):
        surf = readETAData(i,numOfETA,postprocessDir,pts)   # Y axis = eta data
        x = np.linspace(0, Lt, pts)                     # X axis = Total Length (m)

        fig  = plt.figure(figsize=(18,4), dpi=200)
        ax = fig.add_subplot(1,1,1)
        plt.plot(x, surf[1,:], 'c-', linewidth = 0.2)
        plt.axis([-1.5,Lt,min(depth[1,:])-.05,5])

        plt.xlabel('X (m)', fontsize = 12, fontweight = 'bold')
        plt.ylabel(r'$\eta$'+' (m)', fontsize = 12, fontweight = 'bold')
        
         # Water Fill:
        plt.fill_between(x, depth[1,:], surf[1,:],
                         where = surf[1,:] > depth[1,:],
                         facecolor = 'cyan', interpolate =True)

        # Bottom Fill:
        plt.fill_between(x, min(depth[1,:])-.05, depth[1,:], 
                                   where= depth[1,:] > (depth[1,:]-.05),facecolor = '0.35',
                                   hatch = 'X')

        # Time Annotations:
        Time = plt.annotate('Time: %4.2f sec'%(i*plotInt-plotInt), xy=(Lt*.85,3),fontsize = 16,
                                     ha='center', va='center')

        if i < 10:
            fileIndex = '0'+str(i)
            fileName = postprocessDir+'surf_000{0:s}'.format(fileIndex)
            plt.savefig(fileName, ext="png", bbox_inches='tight')
            plt.close()
            
            # video load bar progress description:
            video_load.min = 0 
            video_load.max = numOfETA
            video_load.value = i
            video_load.description = '%d/%d'%(i,numOfETA)
            
        elif i < 100:
            fileIndex = str(i)
            fileName = postprocessDir+'surf_000{0:s}'.format(fileIndex)
            plt.savefig(fileName, ext="png", bbox_inches='tight')
            plt.close()
            
            # video load bar progress description:
            video_load.min = 0 
            video_load.max = numOfETA
            video_load.value = i
            video_load.description = '%d/%d'%(i,numOfETA)
            
        elif i < 1000: 
            fileIndex = str(i)
            fileName = postprocessDir+'surf_00{0:s}'.format(fileIndex)
            plt.savefig(fileName, ext="png", bbox_inches='tight')
            plt.close()
            
            # video load bar progress description:
            video_load.min = 0 
            video_load.max = numOfETA
            video_load.value = i
            video_load.description = '%d/%d'%(i,numOfETA)
            
        else: 
            fileIndex = str(i)
            fileName = postprocessDir+'surf_0{0:s}'.format(fileIndex)
            plt.savefig(fileName,ext="png", bbox_inches='tight')
            plt.close()
            
            # video load bar progress description:
            video_load.min = 0 
            video_load.max = numOfETA
            video_load.value = i
            video_load.description = '%d/%d'%(i,numOfETA)       
    #ffmpeg_run(postprocessDir)
        
            


def plotResult_function(variable):
    fig, ax = plt.subplots(figsize=(15,5), dpi=600)
    plt.close(fig)
    files = [(plot_time.value/plotInt_text.value) + 1]
    
    pwd = os.getcwd()  # get current path
    s = title_text.value 
    wrd_lst = s.split()
    space = '_'
    folder_name = space.join(wrd_lst) # substitute ' ' space in project title with '_'
    
    postprocessDir=os.path.join(pwd,folder_name+'/output/') # create output folder path
    #os.chdir(postprocessDir+'..')
    
    Depthtext = os.path.join(pwd,folder_name+'/depth.txt')
    depth = (np.loadtxt(Depthtext))*-1 
    pts  = len(depth.T)     # number of points
    
    data_text = os.path.join(pwd,folder_name+'/data.txt') # create path to open data.txt in project folder
    fin= open(data_text,'r')  # upload data.txt; which has the bathy dx data
    val = fin.read()            
    val=val.split()
    dx = float(val[5])
    
    Lt = np.ceil(dx*(pts+1)) # compute total horizontal length
    
    x = np.linspace(0, Lt, pts)
    X = np.linspace(0,Lt,Lt)
    
    ## Plot Output
    ax.clear()
    for num in range(len(files)):
        fnum = '%.5d' % files[num]


        eta = np.loadtxt(postprocessDir+'eta_'+fnum)

        ax.clear()
        ax.plot(np.asarray(x),depth[0,:],'k',np.asarray(X),eta[0,:],'c',linewidth=1)
        ax.set_xlabel('X (m)',fontsize = 12, fontweight = 'bold')
        ax.set_ylabel(r'$\eta$'+' (m)',fontsize = 12, fontweight = 'bold')
        ax.axis([-1.5,Lt,min(depth[1,:])-.05,5])

        # Water Fill:
        ax.fill_between(x, depth[0,:], eta[0,:],
                         where = eta[0,:] > depth[0,:],
                         facecolor = 'cyan', interpolate =True)

        # Bottom Fill:
        ax.fill_between(x, min(depth[0,:])-.05, depth[0,:], 
                                   where= depth[0,:] > (depth[0,:]-.05),facecolor = '0.35',
                                   hatch = 'X')
        # Time Annotations:
        ax.set_title('Time: %4.2f sec'%(files[num]*plotInt_text.value-plotInt_text.value),fontsize = 16)
    display(fig)

def saveResult_function(variable):
    fig, ax = plt.subplots(figsize=(15,5), dpi=200)
    plt.close(fig)
    files = [(plot_time.value/plotInt_text.value) + 1]
    
    pwd = os.getcwd()  # get current path
    s = title_text.value 
    wrd_lst = s.split()
    space = '_'
    folder_name = space.join(wrd_lst) # substitute ' ' space in project title with '_'
    postprocessDir=os.path.join(pwd,folder_name+'/output/') # create output folder path
    #os.chdir(postprocessDir+'..')

    Depthtext = os.path.join(pwd,folder_name+'/depth.txt')
    depth = (np.loadtxt(Depthtext))*-1 
    pts  = len(depth.T)     # number of points
    
    data_text = os.path.join(pwd,folder_name+'/data.txt') # create path to open data.txt in project folder
    fin= open(data_text,'r')  # upload data.txt; which has the bathy dx data
    val = fin.read()            
    val=val.split()
    dx = float(val[5])
    
    Lt = np.ceil(dx*(pts+1)) # compute total horizonta length
    
    x = np.linspace(0, Lt, pts)
    X = np.linspace(0,Lt,Lt)
    
    ## Plot Output
    ax.clear()
    for num in range(len(files)):
        fnum = '%.5d' % files[num]

        eta = np.loadtxt(postprocessDir+'eta_'+fnum)

        ax.clear()
        ax.plot(np.asarray(x),depth[0,:],'k',np.asarray(X),eta[0,:],'c',linewidth=1)
        ax.set_xlabel('X (m)',fontsize = 12, fontweight = 'bold')
        ax.set_ylabel(r'$\eta$'+' (m)',fontsize = 12, fontweight = 'bold')
        ax.axis([-1.5,Lt,min(depth[1,:])-.05,5])

        # Water Fill:
        ax.fill_between(x, depth[0,:], eta[0,:],
                         where = eta[0,:] > depth[0,:],
                         facecolor = 'cyan', interpolate =True)

        # Bottom Fill:
        ax.fill_between(x, min(depth[0,:])-.05, depth[0,:], 
                                   where= depth[0,:] > (depth[0,:]-.05),facecolor = '0.35',
                                   hatch = 'X')
        # Time Annotations:
        ax.set_title('Time: %4.2f sec'%(files[num]*plotInt_text.value-plotInt_text.value),fontsize = 16)       
    # save figure
    folder_path=os.path.join(pwd,folder_name+'/') # path to save image in project folder
    fileName = folder_path+'eta_1d_Time:%4.2fsec.png'%(files[num]*plotInt_text.value-plotInt_text.value)
    fig.savefig(fileName, dpi=fig.dpi) # save figure
    
def plotSurface_clicked(variable):    # plot results
    clear_output(wait=True)
    display(principal_tab)
    plotResult_function(variable)        

def saveSurface_clicked(variable):    # save plot results
    saveResult_function(variable) 
    
###--------------------------------------------       
# activate functions    
from pyFiles2.PrincipalTab_5 import plot_results_button,save_plot_results_button,Video_button
# ^ import pertinent variables of PrincipalTab_5

plot_results_button.on_click(plotSurface_clicked)             # activate plot result button
save_plot_results_button.on_click(saveSurface_clicked)    # activate save plot button
Video_button.on_click(runVID_function)             # activate run video button 