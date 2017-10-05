import numpy as np
import matplotlib.pyplot as plt
from pyFiles2.PrincipalTab_1 import bathy_list, vert2_elev,vert2_loc,vert3_elev,THL,NumSeg,MWL,dom,upload_bathy_name
from pyFiles2.PRINCIPAL_TAB import principal_tab,title_text,container_title,space_box2, project_button
from IPython.display import display, clear_output
import os
import shutil

### Plot and Save Bathy Functions  ###
# globals (for slider replotting)
fig, ax = plt.subplots(figsize=(15,5), dpi=600)
plt.close(fig)


def update_bathy(variable):
# This function generates a new depth file by changing the uploaded bathy file to FUNWAVE's format
    results = []
    
    
    pwd = os.getcwd()  # get current path
    
    s = title_text.value 
    wrd_lst = s.split()
    space = '_'
    folder_name = space.join(wrd_lst) # substitute ' ' space in project title with '_'
    
    filename = upload_bathy_name.value   # uploaded file name
    file_path = os.path.join(pwd,folder_name+'/'+filename) # create path to the uploaded depth file in project folder
    f = open(file_path,'r')

    for line in f:
        all = [float(val) for val in line.split()] 
        results.append(all)

    depth = np.asarray(results[0])

    
    depth_text = os.path.join(pwd,folder_name+'/depth.txt') # create path to save depth.txt in project folder
    fin = open(depth_text,'w')
    for i in range(3): 
        # need three rows for Funwave to run a 1D case (really a 2D case)
        for depthPoint in depth:
            depthFunwaveFormat = depthPoint*-1
            fin.write(str(depthFunwaveFormat)+' ')
        fin.write('\n')
    fin.close()

    points = len(depth)
    dx = THL.value/(float(points)+1)
    
    data_text = os.path.join(pwd,folder_name+'/data.txt') # create path to save data.txt in project folder
    DataFile = open(data_text,'w')         # save text file with points and dx for future use
    dataText = """points = %d
    dx = %f"""%(points,dx)
    DataFile.write(dataText)
    DataFile.close()
 

def on_plot_uploaded_file(variable):      # plot uploaded bathy file function
    results = []
    
    pwd = os.getcwd()  # get current path
    s = title_text.value 
    wrd_lst = s.split()
    space = '_'
    folder_name = space.join(wrd_lst) # substitute ' ' space in project title with '_'
    
    filename = upload_bathy_name.value   # uploaded file name
    file_path = os.path.join(pwd,folder_name+'/'+filename) # create path to the uploaded depth file in project folder
    f = open(file_path,'r')
    
    for line in f:
        all = [float(val) for val in line.split()] 
        results.append(all)

    depth = np.asarray(results[0])
    points = len(depth)
    x = np.linspace(0,THL.value,points)
   
    ax.clear()
    ax.plot(x,depth,color='k',linewidth=2.0)
    ax.set_xlabel('Length (m)', fontsize = 15)
    ax.set_ylabel('Height (m)', fontsize = 15)
    ax.axis([0,THL.value,min(depth)-.55,max(depth)+1])
    
    # land fill
    ax.fill_between(x,min(depth)-0.55,depth,where=depth>(min(depth)-0.5),
                facecolor='0.35',hatch='X')
    
    # water fill
    waterLevel =  np.zeros(int(points))
    ax.fill_between(x, depth, waterLevel, 
                            where=waterLevel>depth,facecolor='cyan', interpolate=True)    
    display(fig)

def on_plot_flat(variable): # plot flat bathymetry function
    dx = dom.value
    if dx > 1:
        points = 0
        warning = 'dx must be less than one.'
        raise Exception(warning)
                
    if dx >= 0.01:
        numOfSegments = int(np.ceil(THL.value / dx))
        points = int(numOfSegments - 1)

    x = np.linspace(0,THL.value,points)
    depth = np.zeros(len(x))+ MWL.value

    ax.clear()
    ax.plot(x,depth,color='k',linewidth=2.0)
    ax.set_xlabel('Length (m)', fontsize = 15)
    ax.set_ylabel('Height (m)', fontsize = 15)
    ax.axis([0,THL.value,min(depth)-.55,max(depth)*-1])
    
    # land fill
    ax.fill_between(x,min(depth)-0.55,depth,where=depth>(min(depth)-0.5),
                facecolor='0.35',hatch='X')
    
    # water fill
    waterLevel =  np.zeros(int(points))
    ax.fill_between(x, depth, waterLevel, 
                            where=waterLevel>depth,facecolor='cyan', interpolate=True)    
    display(fig)
    

def compute_flat(variable):
# this function computes the flat bathymetry    
    dx = dom.value
    if dx > 1:
        points = 0
        warning = 'dx must be less than one.'
        raise Exception(warning)
                
    if dx >= 0.01:
        numOfSegments = int(np.ceil(THL.value / dx))
        points = int(numOfSegments - 1)

    x = np.linspace(0,THL.value,points)
    depth = np.zeros(len(x))+ MWL.value
    
    pwd = os.getcwd()  # get current path
    s = title_text.value 
    wrd_lst = s.split()
    space = '_'
    folder_name = space.join(wrd_lst) # substitute ' ' space in project title with '_'
    data_text = os.path.join(pwd,folder_name+'/data.txt') # create path to save data.txt in project folder
    DataFile = open(data_text,'w')         # save text file with points and dx for future use
    dataText = """points = %d
    dx = %f"""%(points,dx)
    DataFile.write(dataText)
    DataFile.close()
    
    return depth
    
    
def compute_high_difference(variable):   
# this function computes the slope bathymetry, it is continued by "setDepth" function
    # compute length difference between vertices
    Vertex_Loc = [0.0,vert2_loc.value,THL.value]
    Vertex_Loc.sort()
    SegLen = ([Vertex_Loc[num]- Vertex_Loc[num-1]
                              for num in range(len(Vertex_Loc)) if num > 0])
    
    # compute high difference between vertices
    VertexElev = [MWL.value,vert2_elev.value,vert3_elev.value]
    
    
    Hi_dif = [0.0]*NumSeg
    Hi_dif = ([VertexElev[num]- VertexElev[num-1]
                              for num in range(len(VertexElev)) if num > 0])
    
    
    dx = dom.value
    if dx > 1:
        points = 0
        warning = 'dx must be less than one.'
        raise Exception(warning)
                
    if dx >= 0.01:
        numOfSegments = int(np.ceil(THL.value / dx))
        points = int(numOfSegments - 1)
    
    pwd = os.getcwd()  # get current path
    s = title_text.value 
    wrd_lst = s.split()
    space = '_'
    folder_name = space.join(wrd_lst) # substitute ' ' space in project title with '_'
    
    data_text = os.path.join(pwd,folder_name+'/data.txt') # create path to save data.txt in project folder
    DataFile = open(data_text,'w')         # save text file with points and dx for future use
    dataText = """points = %d
    dx = %f"""%(points,dx)
    DataFile.write(dataText)
    DataFile.close()
    
    [depth, points, VertexElev, Vertex_Loc] = setDepth(points,dx,VertexElev,Hi_dif,SegLen,Vertex_Loc)
    return [depth, points, VertexElev, Vertex_Loc]


def setDepth(points,dx,VertexElev,Hi_dif,SegLen,Vertex_Loc):    # compute depth function
    # depth NumPy array instantiation
    depth = np.zeros(points)
    
    # lambda (inline) functions for each segment height
    seg = []
    for index in range(NumSeg):
        # Create Line function (y = b + m*x)
        seg.append(lambda x: VertexElev[index]+ Hi_dif[index]/SegLen[index]*x)

    # set segment breakpoints/ends
    segmentEnds = [0.0]*NumSeg   #placeholder
    for outerIndex in range(NumSeg):
        for innerIndex in range(outerIndex+1):
            segmentEnds[outerIndex] += SegLen[innerIndex]
    
    # determine which segment point belongs to and set its depth
    for point in range(points):
        for index,endPoint in enumerate(segmentEnds):
            location = point*dx
            if index is 0:
                    if location <= endPoint:
                        shiftedToOrigin = location
                        depth[point] = seg[index] (shiftedToOrigin)
            else:
                    if location<=endPoint and location>segmentEnds[index-1]:
                        shiftedToOrigin = location-segmentEnds[index-1]
                        depth[point] = seg[index] (shiftedToOrigin)       
    return [depth, points, VertexElev, Vertex_Loc]

    
def on_plot(variable):      # plot slope bathymetry function
    [depth, points, VertexElev, Vertex_Loc] = compute_high_difference(variable)

    # plotting the bathymetry with surface at z = 0
    x = np.linspace(0,THL.value,points)

    # clearing output/plot so we can plot on same graph
    #clear_output(wait=True)
    #fig, ax = plt.subplots(figsize=(15,5), dpi=600)
    ax.clear()

    ax.plot(x,depth,color='k',linewidth=2.0)
    
    ax.axis([0,THL.value,min(depth)-.55,max(depth)+1])

    ax.fill_between(x,min(depth)-0.55,depth,where=depth>(min(depth)-0.5),
                facecolor='0.35',hatch='X')
    ax.set_xlabel('Length (m)', fontsize = 15)
    ax.set_ylabel('Height (m)', fontsize = 15)

    # water fill
    waterLevel =  np.zeros(int(points))
    ax.fill_between(x, depth, waterLevel, 
                            where=waterLevel>depth,facecolor='cyan', interpolate=True)    
    
    ## identify each vertix in the plot using scatter plot:
    #vert1
    ax.scatter([0.0, ],[VertexElev[0], ], 100, color='red')
    ax.annotate(r'$\#1$',
             xy=(0.0,VertexElev[0]), xycoords='data',
             xytext=(+10, +30), textcoords='offset points', fontsize=16,
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
    #vert2
    ax.scatter([Vertex_Loc[1], ],[VertexElev[1], ], 100, color='red')
    ax.annotate(r'$\#2$',
             xy=(Vertex_Loc[1],VertexElev[1]), xycoords='data',
             xytext=(+10, +30), textcoords='offset points', fontsize=16,
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"))
    #vert3             
    ax.scatter([Vertex_Loc[2], ],[VertexElev[2], ], 100, color='red')
    ax.annotate(r'$\#3$',
             xy=(Vertex_Loc[2],VertexElev[2]), xycoords='data',
             xytext=(+10, +30), textcoords='offset points', fontsize=16,
             arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=0.2"))
    
    display(fig)
    
    
def on_save_clicked(variable): # save bathy text file function              
    if bathy_list.value == 'Upload File': # if the chosen bathy type is Upload, go to update_bathy function
        update_bathy(variable)
        
    elif bathy_list.value =='Slope': # if the chosen bathy type is Slope, save depth.txt file
        [depth, points, VertexElev, Vertex_Loc] = compute_high_difference(variable)
        
        pwd = os.getcwd()  # get current path
        s = title_text.value 
        wrd_lst = s.split()
        space = '_'
        folder_name = space.join(wrd_lst) # substitute ' ' space in project title with '_'
        
        depth_text = os.path.join(pwd,folder_name+'/depth.txt') # create path to save depth.txt in project folder
        fin = open(depth_text,'w')
        
        for i in range(3): 
        # need three rows for Funwave to run a 1D case (really a 2D case)
            for depthPoint in depth:
                depthFunwaveFormat = depthPoint*-1
                fin.write(str(depthFunwaveFormat)+' ')
            fin.write('\n')
        fin.close() 
        
    elif bathy_list.value == 'Flat': # if the chosen bathy type is Flat, save depth.txt file
        depth = compute_flat(variable)
        
        pwd = os.getcwd()  # get current path
        s = title_text.value 
        wrd_lst = s.split()
        space = '_'
        folder_name = space.join(wrd_lst) # substitute ' ' space in project title with '_'
        
        depth_text = os.path.join(pwd,folder_name+'/depth.txt') # create path to save depth.txt in project folder
        fin = open(depth_text,'w')
        
        for i in range(3): 
        # need three rows for Funwave to run a 1D case (really a 2D case)
            for depthPoint in depth:
                depthFunwaveFormat = depthPoint*-1
                fin.write(str(depthFunwaveFormat)+' ')
            fin.write('\n')
        fin.close() 
        
    else:
        warning = 'Please choose Bathy type.'
        raise Exception(warning)

        
def on_plot_clicked(variable):    # plot bathymetry fucntion
    clear_output(wait=True)
    display(space_box2,container_title,space_box2) # display project text box and button
    display(principal_tab) # display GUI
    
    if bathy_list.value == 'Slope':  # if the chosen bathy type is Slope, go to on_plot function
        on_plot(variable)
        
    elif bathy_list.value == 'Upload File': # if the chosen bathy type is Upload, go to on_plot_uploaded_file function
        on_plot_uploaded_file(variable)
    
    elif bathy_list.value =='Flat':
        on_plot_flat(variable)
    
    else:
        warning = 'Please choose Bathy type.'
        raise Exception(warning)

    


# -------------------------------------------------------------------------------------
# activate bathymetry functions 
from pyFiles2.PrincipalTab_1 import plot_button, save_button
from pyFiles2.PrincipalTab_1 import MWL, vert2_elev, vert3_elev
from pyFiles2.PrincipalTab_1 import vert2_loc

plot_button.on_click(on_plot_clicked) # user-initialized plotting via button
save_button.on_click(on_save_clicked) # active save depth.txt button


