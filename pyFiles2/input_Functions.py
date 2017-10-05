import numpy as np
import matplotlib.pyplot as plt
import os
import shutil

# import pertinent variables from principal tab
from pyFiles2.PRINCIPAL_TAB import title_text

# import pertinent variables from principal tab 2
from pyFiles2.PrincipalTab_2 import processors_text, time_text,plotInt_text,show_initial, dif, fric, Dir, wave_maker,csp_text,CDsponge_text,LSW_text, RSW_text,BSW_text,TSW_text,R_sponge_text,A_sponge_text,input_verification

# import wavemaker variables from principal tab 2 (continued)
from pyFiles2.PrincipalTab_2 import xc, yc, wid, amp, dep, LagTime, Xwavemaker, xc_wk, yc_wk, tPeriod, ampWK, depWK, thetaWK, TimeRamp, ywidth_wk, deltaWK, FreqPeak, FreqMin, FreqMax, HMO, GammaTMA, ThetaPeak, NFreq, NTheta 

# import pertinent variables from principal tab 3
from pyFiles2.PrincipalTab_3 import DEPTH_OUT,U,V,ETA,Hmax,Hmin,MFmax,Umax,VORmax,Umean,Vmean,ETAmean,MASK,MASK9,SourceX, SourceY,P,Q,Fx,Fy,Gx,Gy,AGE,WaveHeight,steady_time,T_INTV_MEAN

 

def boolean_function(variable):
    
    
# This function changes the checkbox variables True & False to T & F
# since that is the format of FUNWAVE's input file. 
    if show_initial.value == True:
        turn_on_IC = 'T'
    else:
        turn_on_IC = 'F'
        
    if dif.value == True:
        dif_text = 'T'
    else:
        dif_text = 'F'   
    
    if fric.value == True:
        fric_text = 'T'
    else:
        fric_text = 'F'
        
    if Dir.value == True:
        dir_text = 'T'
    else:
        dir_text = 'F'
        
    if DEPTH_OUT.value == True:
        depth_out = 'T'
    else:
        depth_out = 'F'    
        
    if U.value == True:
        u = 'T'
    else:
        u = 'F'
    
    if V.value == True:
        v = 'T'
    else:
        v = 'F' 
        
    if ETA.value == True:
        eta = 'T'
    else:
        eta = 'F' 
    
    if Hmax.value == True:
        hmax = 'T'
    else:
        hmax = 'F'
    
    if Hmin.value == True:
        hmin = 'T'
    else:
        hmin = 'F'
        
    if MFmax.value == True:
        mfmax = 'T'
    else:
        mfmax = 'F'    
        
    if Umax.value == True:
        umax = 'T'
    else:
        umax = 'F'
    
    if VORmax.value == True:
        vormax = 'T'
    else:
        vormax = 'F'
        
    if Umean.value == True:
        umean = 'T'
    else: 
        umean = 'F'
        
    if Vmean.value == True:
        vmean = 'T'
    else:
        vmean = 'F'
        
    if ETAmean.value == True:
        etamean = 'T'
    else:
        etamean = 'F'
        
    if MASK.value == True:
        mask = 'T'
    else:
        mask = 'F'
    
    if MASK9.value == True:
        mask9 = 'T'
    else:
        mask9 = 'F'
    
    if SourceX.value == True:
        sourcex = 'T'
    else:
        sourcex = 'F'
        
    if SourceY.value == True:
        sourcey = 'T'
    else:
        sourcey = 'F'    
    
    if P.value == True:
        p = 'T'
    else:
        p = 'F'
        
    if Q.value == True:
        q = 'T'
    else:
        q = 'F'  
        
    if Fx.value == True:
        fx = 'T'
    else:
        fx = 'F'
        
    if Fy.value == True:
        fy = 'T'
    else:
        fy = 'F'
        
    if Gx.value == True:
        gx = 'T'
    else:
        gx = 'F'
        
    if Gy.value == True:
        gy = 'T'
    else:
        gy = 'F'    
    
    if AGE.value == True:
        age = 'T'
    else:
        age = 'F'
    
    if WaveHeight.value == True:
        waveheight = 'T'
    else:
        waveheight = 'F'
    
    
    generate_input_file(turn_on_IC,dif_text,fric_text,dir_text,depth_out,
                       u,v,eta,hmax,hmin,mfmax,umax,vormax,umean,vmean,etamean,mask,
                       mask9,sourcex,sourcey,p,q,fx,fy,gx,gy,age,waveheight,steady_time,T_INTV_MEAN) 
    
def generate_input_file(turn_on_IC,dif_text,fric_text,dir_text,depth_out,
                       u,v,eta,hmax,hmin,mfmax,umax,vormax,umean,vmean,etamean,mask,
                       mask9,sourcex,sourcey,p,q,fx,fy,gx,gy,age,waveheight,steady_time,T_INTV_MEAN):
# this function generates FUNWAVE's input.txt
        
        # Step 1 - take bathy points and dx data saved in data.txt 
        pwd = os.getcwd()  # get current path
        s = title_text.value 
        wrd_lst = s.split()
        space = '_'
        folder_name = space.join(wrd_lst) # substitute ' ' space in project title with '_'
        
        data_text = os.path.join(pwd,folder_name+'/data.txt') # create path to open data.txt in project folder
        fin= open(data_text,'r')  
        val = fin.read()            
        val=val.split()
        points = val[2]
        dx = float(val[5])
        
        # Step 2 - create output folder and input file
        output_text=os.path.join(pwd,folder_name+'/output/') # create output path 
        input_path = os.path.join(pwd,folder_name+'/input.txt') # create path to save input.txt in project folder
        
        fin = open(input_path,'w') # create input file
        inputFile ="""!INPUT FILE FOR FUNWAVE_TVD
! NOTE: all input parameter are capital sensitive
! --------------------TITLE-------------------------------------
! title only for log file
TITLE = %s
! -------------------HOT START---------------------------------
HOT_START = F
FileNumber_HOTSTART = 0 
! -------------------PARALLEL INFO-----------------------------
! 
!    PX,PY - processor numbers in X and Y
!    NOTE: make sure consistency with mpirun -np n (px*py)
!    
PX = %s
PY = 1
! --------------------DEPTH-------------------------------------
! Depth types, DEPTH_TYPE=DATA: from depth file
!              DEPTH_TYPE=FLAT: idealized flat, need depth_flat
!              DEPTH_TYPE=SLOPE: idealized slope, 
!                                 need slope,SLP starting point, Xslp
!                                 and depth_flat
DEPTH_TYPE = DATA
! Depth file
! depth format NOD: depth at node (M1xN1), ELE: depth at ele (MxN) 
! where (M1,N1)=(M+1,N+1)  
DEPTH_FILE = depth.txt
DepthFormat = ELE
! if depth is flat and slope, specify flat_depth
DEPTH_FLAT = 30.0
if depth is slope, specify slope and starting point
SLP = 0.05
Xslp = 10.0

! -------------------PRINT---------------------------------
! PRINT*,
! result folder
RESULT_FOLDER = %s
! ------------------DIMENSION-----------------------------
! global grid dimension
Mglob = %s
Nglob = 3
! ----------------- TIME----------------------------------
! time: total computational time/ plot time / screen interval 
! all in seconds
TOTAL_TIME = %4.2f
PLOT_INTV = %4.2f
PLOT_INTV_STATION = 0.5
SCREEN_INTV = 10.0
HOTSTART_INTV = 360000000.0
! -----------------GRID----------------------------------
! if use spherical grid, in decimal degrees
StretchGrid = F
Lon_West = 120.0
Lat_South = 0.0
Dphi = 0.0042
Dtheta = 0.0042 
!DX_FILE = ../input/dx_str.txt
!DY_FILE = ../input/dy_str.txt
CORIOLIS_FILE = ../input/cori_str.txt
! cartesian grid sizes
DX = %4.2f
DY = %4.2f
! --------------- INITIAL UVZ ---------------------------
! INI_UVZ - initial UVZ e.g., initial deformation
!         must provide three (3) files 
INI_UVZ = %s
! if true, input eta u and v file names
ETA_FILE = Ini_Z.txt
U_FILE = Ini_U.txt
V_FILE = Ini_V.txt
!MASK_FILE = ../result_eddy_long/mask_0340
HotStartTime = 0.0
OutputStartNumber = 0
! ----------------WAVEMAKER------------------------------
WAVEMAKER = %s
!
AMP = %4.2f
DEP = %4.2f
LAGTIME = %4.2f
XWAVEMAKER = %4.2f
Xc = %4.2f
Yc = %4.2f
WID = %4.2f
Xc_WK = %4.2f
Yc_WK = %4.2f
DEP_WK = %4.2f
Time_ramp = %4.2f
Delta_WK = %4.2f
FreqPeak = %4.2f
FreqMin = %4.2f
FreqMax = %4.2f
Hmo = %4.2f
GammaTMA = %4.2f
ThetaPeak = %4.2f
Nfreq = %4.2f
Ntheta = %4.2f
Tperiod = %4.2f
AMP_WK = %4.2f
Theta_WK = %4.2f
Ywidth_WK = %4.2f
!
! ---------------- PERIODIC BOUNDARY CONDITION ---------
! South-North periodic boundary condition
!
PERIODIC = F
! ---------------- SPONGE LAYER ------------------------
! DHI type sponge layer
! need to specify widths of four boundaries and parameters
! set width=0.0 if no sponge
! R_sponge: decay rate
! A_sponge: maximum decay rate
! e.g., sharp: R=0.85
!       mild:  R=0.90, A=5.0
!       very mild, R=0.95, A=5.0
DIFFUSION_SPONGE = %s
FRICTION_SPONGE = %s
DIRECT_SPONGE = %s
Csp = %4.2f
CDsponge = %4.2f
Sponge_west_width =  %4.2f
Sponge_east_width =  %4.2f
Sponge_south_width = %4.2f
Sponge_north_width = %4.2f
R_sponge = %4.2f
A_sponge = %4.2f
ETA_LIMITER = F
CrestLimit = 6.0
TroughLimit = -6.0
! ----------------OBSTACLES-----------------------------
! obstacle structures using mask_struc file
! mask_struc =0 means structure element
! give a file contains a mask array with Mloc X Nloc
!OBSTACLE_FILE= breakwaters.txt
! ----------------PHYSICS------------------------------
! parameters to control type of equations
! dispersion: all dispersive terms
! gamma1=1.0,gamma2=0.0: NG's equations
! gamma1=1.0,gamma2=1.0: Fully nonlinear equations
DISPERSION = T
Gamma1 = 1.0
Gamma2 = 1.0    
Gamma3 = 1.0    ! If set to zero, full linear case.
Beta_ref=-0.531
SWE_ETA_DEP = 0.80
VISCOSITY_BREAKING = T
!----------------Friction-----------------------------
Friction_Matrix = F
Cd_file= fric.txt
Cd = 0.0
! Cd = 0.01
! ----------------NUMERICS----------------------------
! time scheme: runge_kutta for all types of equations
!              predictor-corrector for NSWE
! space scheme: second-order
!               fourth-order
! construction: HLLC
! cfl condition: CFL
! froude number cap: FroudeCap

Time_Scheme = Runge_Kutta
!Time_Scheme = Predictor_Corrector
! spacial differencing

HIGH_ORDER = THIRD
CONSTRUCTION = HLLC
! CFL
CFL = 0.5
! Froude Number Cap (to avoid jumping drop, set 10)
FroudeCap = 10.0
! --------------WET-DRY-------------------------------
! MinDepth for wetting-drying
MinDepth=0.05
! -----------------
! MinDepthfrc to limit bottom friction
MinDepthFrc = 0.05
!---------------breaking-----------------------------
!  there are two options for breaking algorithm 
!  1: shock-capturing breaking, need SWE_ETA_DEP
!  2: eddy-viscosity breaking, when VISCOSITY_BREAKING = T
!     the shock-capturing breaking is invalid
!     Cbrk1 and Cbrk2 are parameters defined in Kennedy et al 2000
!     suggested in this model Cbrk1=0.65, Cbrk2=0.15
!     WAVEMAKER_Cbrk is to avoid breaking inside wavemaker 
SWE_ETA_DEP = 0.80
VISCOSITY_BREAKING = T
Cbrk1 = 0.65
Cbrk2 = 0.35
WAVEMAKER_Cbrk = 0.65
! ----------------- MIXING ---------------------------
! if use smagorinsky mixing, have to set -DMIXING in Makefile
! and set averaging time interval, T_INTV_mean, default: 20s
STEADY_TIME = %4.2f
T_INTV_mean = %4.2f
C_smg = 0.2
! ----------------- COUPLING -------------------------
! if do coupling, have to set -DCOUPLING in Makefile
COUPLING_FILE = coupling.txt
! -----------------OUTPUT-----------------------------
! stations 
! if NumberStations>0, need input i,j in STATION_FILE
NumberStations = 0
STATIONS_FILE = stations.txt
! output variables, T=.TRUE, F = .FALSE.
DEPTH_OUT = %s
U = %s
V = %s
ETA = %s
Hmax = %s
Hmin = %s
MFmax = %s
Umax = %s
VORmax = %s
Umean = %s
Vmean = %s
ETAmean = %s
MASK = %s
MASK9 = %s
SourceX = %s
SourceY = %s
P = %s
Q = %s
Fx = %s
Fy = %s
Gx = %s
Gy = %s
AGE = %s
WaveHeight = %s""" % (folder_name,int(processors_text.value),
                              str(output_text),points,time_text.value,plotInt_text.value,dx,dx,turn_on_IC,
                             wave_maker.value,amp.value,dep.value,LagTime.value,Xwavemaker.value,xc.value,
                              yc.value,wid.value,xc_wk.value,yc_wk.value,depWK.value,TimeRamp.value,
                              deltaWK.value,FreqPeak.value,FreqMin.value,FreqMax.value,HMO.value,
                              GammaTMA.value,ThetaPeak.value,NFreq.value,NTheta.value,tPeriod.value,
                              ampWK.value,thetaWK.value,ywidth_wk.value,dif_text, fric_text, dir_text,
                             csp_text.value,CDsponge_text.value,LSW_text.value,
                             RSW_text.value,BSW_text.value,TSW_text.value,
                             R_sponge_text.value,float(A_sponge_text.value),steady_time.value,T_INTV_MEAN.value, 
                              depth_out,u,v,eta,hmax,hmin,mfmax,umax,vormax,umean,vmean,
                              etamean,mask,mask9,sourcex,sourcey,p,q,fx,
                              fy,gx,gy,age,waveheight)
        fin.write(inputFile)
        fin.close()
     
def update_input_function(variable):
    pwd = os.getcwd()  # get current path
    s = title_text.value 
    wrd_lst = s.split()
    space = '_'
    folder_name = space.join(wrd_lst) # substiture ' ' space in project title with '_'
    
    output_text=os.path.join(pwd,folder_name+'/output/') # create output path 
    
### this function lets the user update the input.txt before generating the file

# Step 1: change the checkbox variables True & False to T & F
    if show_initial.value == True:
        turn_on_IC = 'T'
    else:
        turn_on_IC = 'F'
        
    if dif.value == True:
        dif_text = 'T'
    else:
        dif_text = 'F'   
    
    if fric.value == True:
        fric_text = 'T'
    else:
        fric_text = 'F'
        
    if Dir.value == True:
        dir_text = 'T'
    else:
        dir_text = 'F'
        
    if DEPTH_OUT.value == True:
        depth_out = 'T'
    else:
        depth_out = 'F'    
        
    if U.value == True:
        u = 'T'
    else:
        u = 'F'
    
    if V.value == True:
        v = 'T'
    else:
        v = 'F' 
        
    if ETA.value == True:
        eta = 'T'
    else:
        eta = 'F' 
    
    if Hmax.value == True:
        hmax = 'T'
    else:
        hmax = 'F'
    
    if Hmin.value == True:
        hmin = 'T'
    else:
        hmin = 'F'
        
    if MFmax.value == True:
        mfmax = 'T'
    else:
        mfmax = 'F'    
        
    if Umax.value == True:
        umax = 'T'
    else:
        umax = 'F'
    
    if VORmax.value == True:
        vormax = 'T'
    else:
        vormax = 'F'
        
    if Umean.value == True:
        umean = 'T'
    else: 
        umean = 'F'
        
    if Vmean.value == True:
        vmean = 'T'
    else:
        vmean = 'F'
        
    if ETAmean.value == True:
        etamean = 'T'
    else:
        etamean = 'F'
        
    if MASK.value == True:
        mask = 'T'
    else:
        mask = 'F'
    
    if MASK9.value == True:
        mask9 = 'T'
    else:
        mask9 = 'F'
    
    if SourceX.value == True:
        sourcex = 'T'
    else:
        sourcex = 'F'
        
    if SourceY.value == True:
        sourcey = 'T'
    else:
        sourcey = 'F'    
    
    if P.value == True:
        p = 'T'
    else:
        p = 'F'
        
    if Q.value == True:
        q = 'T'
    else:
        q = 'F'  
        
    if Fx.value == True:
        fx = 'T'
    else:
        fx = 'F'
        
    if Fy.value == True:
        fy = 'T'
    else:
        fy = 'F'
        
    if Gx.value == True:
        gx = 'T'
    else:
        gx = 'F'
        
    if Gy.value == True:
        gy = 'T'
    else:
        gy = 'F'    
    
    if AGE.value == True:
        age = 'T'
    else:
        age = 'F'
    
    if WaveHeight.value == True:
        waveheight = 'T'
    else:
        waveheight = 'F'

#Step 2: upload data.txt to obtain point and dx values
    data_text = os.path.join(pwd,folder_name+'/data.txt') # create path to open data.txt in project folder
    fin= open(data_text,'r')   
    val = fin.read()            
    val=val.split()
    points = val[2]
    dx = float(val[5])
    
#Step  3: update FUNWAVE's input.txt before generating it      
    inputFile ="""!INPUT FILE FOR FUNWAVE_TVD </br>
! NOTE: all input parameter are capital sensitive </br>
! --------------------TITLE-------------------------------------</br>
! title only for log file</br>
TITLE = %s</br>
! -------------------HOT START---------------------------------</br>
HOT_START = F</br>
FileNumber_HOTSTART = 0 </br>
! -------------------PARALLEL INFO-----------------------------</br>
! </br>
!    PX,PY - processor numbers in X and Y</br>
!    NOTE: make sure consistency with mpirun -np n (px*py)</br>
!    </br>
PX = %s</br>
PY = 1</br>
! --------------------DEPTH-------------------------------------</br>
! Depth types, DEPTH_TYPE=DATA: from depth file</br>
!              DEPTH_TYPE=FLAT: idealized flat, need depth_flat</br>
!              DEPTH_TYPE=SLOPE: idealized slope, </br>
!                                 need slope,SLP starting point, Xslp</br>
!                                 and depth_flat</br>
DEPTH_TYPE = DATA</br>
! Depth file</br>
! depth format NOD: depth at node (M1xN1), ELE: depth at ele (MxN) </br>
! where (M1,N1)=(M+1,N+1) </br> 
DEPTH_FILE = depth.txt</br>
DepthFormat = ELE</br>
</br>
! -------------------PRINT---------------------------------</br>
! PRINT*,</br>
! result folder</br>
RESULT_FOLDER = %s</br>
! ------------------DIMENSION-----------------------------</br>
! global grid dimension</br>
Mglob = %s</br>
Nglob = 3</br>
! ----------------- TIME----------------------------------</br>
! time: total computational time/ plot time / screen interval </br>
! all in seconds</br>
TOTAL_TIME = %4.2f</br>
PLOT_INTV = %4.2f</br>
PLOT_INTV_STATION = 0.5</br>
SCREEN_INTV = 10.0</br>
HOTSTART_INTV = 360000000.0</br>
! -----------------GRID----------------------------------</br>
! if use spherical grid, in decimal degrees</br>
StretchGrid = F</br>
Lon_West = 120.0</br>
Lat_South = 0.0</br>
Dphi = 0.0042</br>
Dtheta = 0.0042 </br>
! cartesian grid sizes</br>
DX = %4.2f</br>
DY = %4.2f</br>
! --------------- INITIAL UVZ ---------------------------</br>
! INI_UVZ - initial UVZ e.g., initial deformation</br>
!         must provide three (3) files </br>
INI_UVZ = %s</br>
! if true, input eta u and v file names</br>
ETA_FILE = Ini_Z.txt</br>
U_FILE = Ini_U.txt</br>
V_FILE = Ini_V.txt</br>
! ----------------WAVEMAKER------------------------------</br>
WAVEMAKER = %s</br>
!</br>
AMP = %4.2f</br>
DEP = %4.2f</br>
LAGTIME = %4.2f</br>
XWAVEMAKER = %4.2f</br>
Xc = %4.2f</br>
Yc = %4.2f</br>
WID = %4.2f</br>
Xc_WK = %4.2f</br>
Yc_WK = %4.2f</br>
DEP_WK = %4.2f</br>
Time_ramp = %4.2f</br>
Delta_WK = %4.2f</br>
FreqPeak = %4.2f</br>
FreqMin = %4.2f</br>
FreqMax = %4.2f</br>
Hmo = %4.2f</br>
GammaTMA = %4.2f</br>
ThetaPeak = %4.2f</br>
Nfreq = %4.2f</br>
Ntheta = %4.2f</br>
Tperiod = %4.2f</br>
AMP_WK = %4.2f</br>
Theta_WK = %4.2f</br>
Ywidth_WK = %4.2f</br>
!</br>
! ---------------- PERIODIC BOUNDARY CONDITION ---------</br>
! South-North periodic boundary condition</br>
!</br>
PERIODIC = F</br>
! ---------------- SPONGE LAYER ------------------------</br>
! DHI type sponge layer</br>
! need to specify widths of four boundaries and parameters</br>
! set width=0.0 if no sponge</br>
! R_sponge: decay rate</br>
! A_sponge: maximum decay rate</br>
! e.g., sharp: R=0.85</br>
!       mild:  R=0.90, A=5.0</br>
!       very mild, R=0.95, A=5.0</br>
DIFFUSION_SPONGE = %s</br>
FRICTION_SPONGE = %s</br>
DIRECT_SPONGE = %s</br>
Csp = %4.2f</br>
CDsponge = %4.2f</br>
Sponge_west_width =  %4.2f</br>
Sponge_east_width =  %4.2f</br>
Sponge_south_width = %4.2f</br>
Sponge_north_width = %4.2f</br>
R_sponge = %4.2f</br>
A_sponge = %4.2f</br>
ETA_LIMITER = F</br>
CrestLimit = 6.0</br>
TroughLimit = -6.0</br>
! ----------------PHYSICS------------------------------</br>
! parameters to control type of equations</br>
! dispersion: all dispersive terms</br>
! gamma1=1.0,gamma2=0.0: NG's equations</br>
! gamma1=1.0,gamma2=1.0: Fully nonlinear equations</br>
DISPERSION = T</br>
Gamma1 = 1.0</br>
Gamma2 = 1.0 </br>   
Gamma3 = 1.0    ! If set to zero, full linear case.</br>
Beta_ref=-0.531</br>
SWE_ETA_DEP = 0.80</br>
VISCOSITY_BREAKING = T</br>
!----------------Friction-----------------------------</br>
Friction_Matrix = F</br>
Cd_file= fric.txt</br>
Cd = 0.0</br>
! Cd = 0.01</br>
! ----------------NUMERICS----------------------------</br>
! time scheme: runge_kutta for all types of equations</br>
!              predictor-corrector for NSWE</br>
! space scheme: second-order</br>
!               fourth-order</br>
! construction: HLLC</br>
! cfl condition: CFL</br>
! froude number cap: FroudeCap</br>
Time_Scheme = Runge_Kutta</br>
!Time_Scheme = Predictor_Corrector</br>
! spacial differencing</br>
HIGH_ORDER = THIRD</br>
CONSTRUCTION = HLLC</br>
! CFL</br>
CFL = 0.5</br>
! Froude Number Cap (to avoid jumping drop, set 10)</br>
FroudeCap = 10.0</br>
! --------------WET-DRY-------------------------------</br>
! MinDepth for wetting-drying</br>
MinDepth=0.05</br>
! -----------------</br>
! MinDepthfrc to limit bottom friction</br>
MinDepthFrc = 0.05</br>
!---------------breaking-----------------------------</br>
!  there are two options for breaking algorithm </br>
!  1: shock-capturing breaking, need SWE_ETA_DEP</br>
!  2: eddy-viscosity breaking, when VISCOSITY_BREAKING = T</br>
!     the shock-capturing breaking is invalid</br>
!     Cbrk1 and Cbrk2 are parameters defined in Kennedy et al 2000</br>
!     suggested in this model Cbrk1=0.65, Cbrk2=0.15</br>
!     WAVEMAKER_Cbrk is to avoid breaking inside wavemaker </br>
SWE_ETA_DEP = 0.80</br>
VISCOSITY_BREAKING = T</br>
Cbrk1 = 0.65</br>
Cbrk2 = 0.35</br>
WAVEMAKER_Cbrk = 0.65</br>
! ----------------- MIXING ---------------------------</br>
! if use smagorinsky mixing, have to set -DMIXING in Makefile</br>
! and set averaging time interval, T_INTV_mean, default: 20s</br>
STEADY_TIME = %4.2f</br>
T_INTV_mean = %4.2f</br>
C_smg = 0.2</br>
! -----------------OUTPUT-----------------------------</br>
! stations </br>
! if NumberStations>0, need input i,j in STATION_FILE</br>
NumberStations = 0</br>
STATIONS_FILE = stations.txt</br>
! output variables, T=.TRUE, F = .FALSE.</br>
DEPTH_OUT = %s</br>
U = %s</br>
V = %s</br>
ETA = %s</br>
Hmax = %s</br>
Hmin = %s</br>
MFmax = %s</br>
Umax = %s</br>
VORmax = %s</br>
Umean = %s</br>
Vmean = %s</br>
ETAmean = %s</br>
MASK = %s</br>
MASK9 = %s</br>
SourceX = %s</br>
SourceY = %s</br>
P = %s</br>
Q = %s</br>
Fx = %s</br>
Fy = %s</br>
Gx = %s</br>
Gy = %s</br>
AGE = %s</br>
WaveHeight = %s""" % (folder_name,int(processors_text.value),
                              str(output_text),points,time_text.value,plotInt_text.value,dx,dx,turn_on_IC,
                              wave_maker.value,amp.value,dep.value,LagTime.value,Xwavemaker.value,xc.value,
                              yc.value,wid.value,xc_wk.value,yc_wk.value,depWK.value,TimeRamp.value,
                              deltaWK.value,FreqPeak.value,FreqMin.value,FreqMax.value,HMO.value,
                              GammaTMA.value,ThetaPeak.value,NFreq.value,NTheta.value,tPeriod.value,
                              ampWK.value,thetaWK.value,ywidth_wk.value,dif_text, fric_text, dir_text,
                             csp_text.value,CDsponge_text.value,LSW_text.value,
                             RSW_text.value,BSW_text.value,TSW_text.value,
                             R_sponge_text.value,float(A_sponge_text.value),steady_time.value,T_INTV_MEAN.value, 
                              depth_out,u,v,eta,hmax,hmin,mfmax,umax,vormax,umean,vmean,
                              etamean,mask,mask9,sourcex,sourcey,p,q,fx,
                              fy,gx,gy,age,waveheight)

    input_verification.value = inputFile

    
    
###--------------------------------------------       
# activate functions

from pyFiles2.PrincipalTab_2 import update_input_button,inputFile_button
# ^ import pertinent variables of PrincipalTab_2

update_input_button.on_click(update_input_function)  # activate update input values button
inputFile_button.on_click(boolean_function)        # activate generate input file button