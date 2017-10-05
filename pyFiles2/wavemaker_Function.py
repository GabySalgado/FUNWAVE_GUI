import ipywidgets as widgets
from IPython.display import display, clear_output
import os
import shutil
import numpy as np

##### This py file contains all the functions that show/hide
##### some widgets on the GUI

## Show GUI function
# this function shows the entire GUI once the "generate project" button is pressed
# this function works with the generate bathy widgets on PRINCIPAL_TAB
from pyFiles2.PRINCIPAL_TAB import GUI_CONT, title_text
def project_clicked(variable):
       
    s = title_text.value 
    pwd = os.getcwd()  # get current path
    wrd_lst = s.split()
    space = '_'
    folder_name = space.join(wrd_lst) # substitute ' ' space in title with '_'
    
    if folder_name == '': # verify that the user gave a project/folder name
        pass
        warning = "Please specify the Project Title."
        raise Exception(warning)
    else:        
        # Step 1 - create project folder        
        project_text = os.path.join(pwd,folder_name +'/') # create project folder path 

        check_dir = os.path.exists(project_text) # check if project folder exist

        if check_dir == True:         # if project folder exist, print warning
            pass
            warning = "Project title '%s' already exists, please use a different title." %(s)
            raise Exception(warning)
            
        else:                         # else, create project folder
            os.mkdir(project_text)
            # Step 2 - Show GUI
            clear_output(wait=True)
            display(GUI_CONT)
            
## generate bathy widgets function
# this function works with the generate bathy widgets on PrincipalTab_1
from pyFiles2.PrincipalTab_1 import Box_upload,Box_SlopeBathy,Box_FlatBathy,MWL,space_box2
def toggle_choose_bathy(change):
    if change['new'] == 'Upload File':
        #show
        Box_upload.layout.display=''
        #hide
        Box_SlopeBathy.layout.display='none'
        Box_FlatBathy.layout.display='none'       
        
    elif change['new'] == 'Slope':
        #show
        Box_SlopeBathy.layout.display=''
        MWL.layout.width = '100%'
        #hide
        Box_upload.layout.display='none'
        Box_FlatBathy.layout.display='none'
        
    elif change['new'] == 'Flat':
        #show
        Box_FlatBathy.layout.display=''
        MWL.description='Depth'
        MWL.layout.width = '50%'
        space_box2.layout.height = '15px'
        #hide
        Box_upload.layout.display='none' 
        Box_SlopeBathy.layout.display='none'
        
    else: 
        Box_upload.layout.display='none'
        Box_SlopeBathy.layout.display='none'
        Box_FlatBathy.layout.display='none'

## generate initial condition widgets function
# this function works with the widgets on PrincipalTab_2b
from pyFiles2.PrincipalTab_2 import init
def toggle_initial(change):
    if change['new']:
        #show
        init.layout.display=''
    else:
        # hide
        init.layout.display='none'

# This function generates new initial conditions files by changing the uploaded
# files to FUNWAVE's format
# this function works with the widgets on PrincipalTab_2b
from pyFiles2.PrincipalTab_2 import iniElev_text,Uvel_text,Vvel_text
def update_initial_conditions(variable):
    results_Z = []
    results_U = []
    results_V = []
    
    pwd = os.getcwd()  # get current path
    s = title_text.value 
    wrd_lst = s.split()
    space = '_'
    folder_name = space.join(wrd_lst) # substitute ' ' space in project title with '_'
    folder_path = file_pathZ = os.path.join(pwd,folder_name+'/') # folder path, where the initial files will be saved
    
    filenameZ = iniElev_text.value
    file_pathZ = os.path.join(pwd,folder_name+'/'+filenameZ) 
    
    filenameU = Uvel_text.value
    file_pathU= os.path.join(pwd,folder_name+'/'+filenameU) 
    
    filenameV = Vvel_text.value
    file_pathV = os.path.join(pwd,folder_name+'/'+filenameU) 
    
    if filenameZ == '':
        pass

    else:
        
        f = open(file_pathZ,'r')
        for line in f:
            all = [float(val) for val in line.split()] 
            results_Z.append(all)
        IniZ = np.asarray(results_Z[0])
        fin = open(folder_path+'Ini_Z.txt','w')
        for i in range(3): 
            # need three rows for Funwave to run a 1D case (really a 2D case)
            for zpoint in IniZ:
                fin.write(str(zpoint)+' ')
            fin.write('\n')
        fin.close()
        
    if filenameU == '':
        pass
        
    else:
        f = open(file_pathU,'r')
        for line in f:
            all = [float(val) for val in line.split()] 
            results_U.append(all)
        IniU = np.asarray(results_U[0])
        fin = open(folder_path+'Ini_U.txt','w')
        for i in range(3): 
            # need three rows for Funwave to run a 1D case (really a 2D case)
            for upoint in IniU:
                fin.write(str(upoint)+' ')
            fin.write('\n')
        fin.close()
        
    if filenameV == '':
        pass
        
    else:
        f = open(file_pathV,'r')
        for line in f:
            all = [float(val) for val in line.split()] 
            results_V.append(all)
        IniV = np.asarray(results_V[0])
        fin = open(folder_path+'Ini_V.txt','w')
        for i in range(3): 
            # need three rows for Funwave to run a 1D case (really a 2D case)
            for vpoint in IniV:
                fin.write(str(vpoint)+' ')
            fin.write('\n')
        fin.close()                

## generate WAVEMAKER variables widgets function
# this function works with the widgets on PrincipalTab_2c
from pyFiles2.PrincipalTab_2 import container_IniRec,container_Gauss,container_LefSol,container_IniSol,container_WkReg
from pyFiles2.PrincipalTab_2 import container_JON1D,container_JON2D,container_WKIRR,container_TMA_1D, GammaTMA
                                
def toggle_waveMaker(change):
    if change['new'] == 'INI_REC':
        #show
        container_IniRec.layout.display=''
        #hide
        container_Gauss.layout.display='none'
        container_LefSol.layout.display='none'
        container_IniSol.layout.display='none'
        container_WkReg.layout.display='none'
        container_JON1D.layout.display='none'
        container_JON2D.layout.display='none'
        container_WKIRR.layout.display='none'
        container_TMA_1D.layout.display='none'
        
    elif change['new'] == 'GAUSIAN':
        #show
        container_Gauss.layout.display=''
        #hide
        container_IniRec.layout.display='none'
        container_LefSol.layout.display='none'
        container_IniSol.layout.display='none'
        container_WkReg.layout.display='none'
        container_JON1D.layout.display='none'
        container_JON2D.layout.display='none'
        container_WKIRR.layout.display='none'
        container_TMA_1D.layout.display='none'
        
    elif change['new'] == 'LEF_SOL':
        #show
        container_LefSol.layout.display=''
        #hide
        container_IniRec.layout.display='none' 
        container_Gauss.layout.display='none'
        container_IniSol.layout.display='none'
        container_WkReg.layout.display='none'
        container_JON1D.layout.display='none'
        container_JON2D.layout.display='none'
        container_WKIRR.layout.display='none'
        container_TMA_1D.layout.display='none'
        
    elif change['new'] == 'INI_SOL':
        #show
        container_IniSol.layout.display=''
        #hide
        container_IniRec.layout.display='none' 
        container_Gauss.layout.display='none'
        container_LefSol.layout.display='none'
        container_WkReg.layout.display='none'
        container_JON1D.layout.display='none'
        container_JON2D.layout.display='none'
        container_WKIRR.layout.display='none'
        container_TMA_1D.layout.display='none'
        
    elif change['new'] == 'WK_REG':
        #show
        container_WkReg.layout.display=''
        #hide
        container_IniRec.layout.display='none' 
        container_Gauss.layout.display='none'
        container_LefSol.layout.display='none'
        container_IniSol.layout.display='none'
        container_JON1D.layout.display='none'
        container_JON2D.layout.display='none'
        container_WKIRR.layout.display='none'
        container_TMA_1D.layout.display='none'
    
    elif change['new'] == 'JON_1D':
        #show
        container_JON1D.layout.display=''
        #hide
        container_IniRec.layout.display='none' 
        container_Gauss.layout.display='none'
        container_LefSol.layout.display='none'
        container_IniSol.layout.display='none'
        container_WkReg.layout.display='none'
        container_JON2D.layout.display='none'
        container_WKIRR.layout.display='none'
        container_TMA_1D.layout.display='none'
    
        # gamma value for Jonwsap = 3.3
        GammaTMA.step = 0.01
        GammaTMA.min = 3.3
        GammaTMA.value = 3.3
    
    elif change['new'] == 'JON_2D':
        #show
        container_JON2D.layout.display=''
        #hide
        container_IniRec.layout.display='none' 
        container_Gauss.layout.display='none'
        container_LefSol.layout.display='none'
        container_IniSol.layout.display='none'
        container_WkReg.layout.display='none'
        container_JON1D.layout.display='none'
        container_WKIRR.layout.display='none'
        container_TMA_1D.layout.display='none'
        
        # gamma value for Jonwsap = 3.3
        GammaTMA.step = 0.01
        GammaTMA.min = 3.3
        GammaTMA.value = 3.3
    
    elif change['new'] == 'WK_IRR':
        #show
        container_WKIRR.layout.display=''
        #hide
        container_IniRec.layout.display='none' 
        container_Gauss.layout.display='none'
        container_LefSol.layout.display='none'
        container_IniSol.layout.display='none'
        container_WkReg.layout.display='none'
        container_JON1D.layout.display='none'
        container_JON2D.layout.display='none'
        container_TMA_1D.layout.display='none'
        
        # gamma value for TMA = 5.0
        GammaTMA.value = 5.0
        GammaTMA.step = 0.01
        GammaTMA.min = 5.0
        
    elif change['new'] == 'TMA_1D':
        #show
        container_TMA_1D.layout.display=''
        #hide
        container_IniRec.layout.display='none' 
        container_Gauss.layout.display='none'
        container_LefSol.layout.display='none'
        container_IniSol.layout.display='none'
        container_WkReg.layout.display='none'
        container_JON1D.layout.display='none'
        container_JON2D.layout.display='none'
        container_WKIRR.layout.display='none'
        
        # gamma value for TMA = 5.0
        GammaTMA.value = 5.0
        GammaTMA.step = 0.01
        GammaTMA.min = 5.0
        
    else: 
        container_Gauss.layout.display='none'
        container_IniRec.layout.display='none'
        container_LefSol.layout.display='none'
        container_IniSol.layout.display='none'
        container_WkReg.layout.display='none'
        container_JON1D.layout.display='none'
        container_JON2D.layout.display='none'
        container_WKIRR.layout.display='none'
        container_TMA_1D.layout.display='none'
        
## generate sponge layer variables widgets function
# this function works with the widgets on PrincipalTab_2d
from pyFiles2.PrincipalTab_2 import container_csp,container_CDsponge,container_R_sponge,container_A_sponge 
def toggle_DifSponge(change):  # show/hide diffusion sponge variables
    if change['new']:
        #show
        container_csp.layout.display=''
        
    else:
        # hide
        container_csp.layout.display='none'
        
           
def toggle_FricSponge(change): # show/hide friction sponge variables
    if change['new']:
        #show
        container_CDsponge.layout.display=''
        
    else:
        # hide
        container_CDsponge.layout.display='none'
        
        
def toggle_DirSponge(change): # show/hide direct sponge variables
    if change['new']:
        #show
        container_R_sponge.layout.display=''
        container_A_sponge.layout.display=''
        
    else:
        # hide
        container_R_sponge.layout.display='none'
        container_A_sponge.layout.display='none'
            
    
## generate wave_Height variables widgets function
# this function works with the widgets on PrincipalTab_3 (Output tab)
from pyFiles2.PrincipalTab_3 import container_col3
def toggle_waveH(change):
    if change['new']:
        #show
        container_col3.layout.display=''
    else:
        # hide
        container_col3.layout.display='none'
        
        
###--------------------------------------------       
# activate functions

## PRINCIPAL_TAB functions
from pyFiles2.PRINCIPAL_TAB import project_button # import pertinent variables of PRINCIPAL_TAB
project_button.on_click(project_clicked) # activate GUI

## PrincipalTab_1 functions
from pyFiles2.PrincipalTab_1 import bathy_list # import pertinent variables of PrincipalTab_1

bathy_list.observe(toggle_choose_bathy, 'value')   # activate bathy type variables
toggle_choose_bathy({'new': bathy_list.value})

## PrincipalTab_2 functions
from pyFiles2.PrincipalTab_2 import show_initial,ini_button,wave_maker,Dir,fric,dif,update_input_button,inputFile_button
# ^ import pertinent variables of PrincipalTab_2

show_initial.observe(toggle_initial, 'value')   # activate initial conditions input variables
toggle_initial({'new': show_initial.value})

ini_button.on_click(update_initial_conditions)  # activate generate inintial condition files button

wave_maker.observe(toggle_waveMaker, 'value')   # activate wavemaker variables
toggle_waveMaker({'new': wave_maker.value})

dif.observe(toggle_DifSponge, 'value')   # activate diffusion sponge variables
toggle_DifSponge({'new': dif.value})

fric.observe(toggle_FricSponge, 'value')   # activate friction sponge variables
toggle_FricSponge({'new': fric.value})

Dir.observe(toggle_DirSponge, 'value')   # activate direct sponge variables
toggle_DirSponge({'new': Dir.value})

## PrincipalTab_3 functions
from pyFiles2.PrincipalTab_3 import WaveHeight # import waveHeight variable from principalTab_3

WaveHeight.observe(toggle_waveH, 'value')   # activate Wave Height Output variables
toggle_waveH({'new': WaveHeight.value})

