import ipywidgets as widgets
from traitlets import link

#####################################################
#### Tabs for Principal Tab 2: Funwave Input txt ####
#####################################################


################################
#### Tab #2a: Input Intro ####
################################

space_box1 = widgets.Box(layout=widgets.Layout(height ='25px', width='90%')) # box created to have space among widgets 

# label with input.txt intro
label_input_INTRO = widgets.HTML("""This Notebook creates <b>FUNWAVE's
input.txt</b> according to the user's specifications.<br> It consists of the following six steps:""",
layout = widgets.Layout(height = '50px', width = '90%'))

space_box2 = widgets.Box(layout=widgets.Layout(width='10%')) # box created to have space among widgets

input_step1 = widgets.HTML("""<b>1- Project Intro:</b> In this tab the user submits the  
<b>Number of Processors</b>, the simulation's <b>Total Time</b>, and the <b>Plot Interval.</b>""")

input_step2 = widgets.HTML("""<b>2- Initial Conditions:</b> In this tab the user determines if the project has 
initial conditions and specifies their names. These files must be uploaded in the Jupyter 
Notebook Home folder.""")

input_step3 = widgets.HTML("""<b>3- Wave Maker:</b> The user chooses the wave maker and inputs it's respective 
specifications.""")

input_step4 = widgets.HTML("""<b>4- Sponge Layer:</b> The user chooses the sponge layers and inputs their respective 
specifications.""")

input_step5 = widgets.HTML("""<b>5- Output Options:</b> The user chooses all the variables desired in the
 output text file.""")

input_step6 = widgets.HTML("""<b>6- Generate Input:</b> The user verifies and generates the project's input.txt""")

step_box = widgets.VBox([input_step1,input_step2,input_step3,input_step4,input_step5,input_step6],
                        layout = widgets.Layout(width = '80%'))

space_step_box = widgets.HBox([space_box2,step_box]) # box containing the steps of Input tab


# Number of processors widget container (label & textbox)
label_processors = widgets.Label('Number of Processors:')
processors_text = widgets.BoundedFloatText(max = '32', min = '4',layout = widgets.Layout(width = "50%"))
container_processors = widgets.VBox(children=[label_processors,processors_text],
                                    layout = widgets.Layout(width = "30%"))

# Total project time widget container (label & textbox)
label_time = widgets.Label('Total Time (sec):')
time_text = widgets.BoundedFloatText(max=3600, layout=processors_text.layout)
container_time = widgets.VBox(children=[label_time,time_text],
                              layout = widgets.Layout(width = "30%"))

# Total project time widget container (label & textbox)
label_plotInt = widgets.Label('Plot Interval (sec):')
plotInt_text = widgets.BoundedFloatText(max=3600, layout=processors_text.layout)
container_pltint = widgets.VBox(children=[label_plotInt,plotInt_text],
                                layout = widgets.Layout(width = "30%"))

container_proc_time = widgets.HBox([space_box2,container_processors,container_time,container_pltint],
                                  layout = widgets.Layout(height = '75px'))

# tab 2a container box
page_inputIntro = widgets.VBox(children=[label_input_INTRO,space_step_box,space_box1,space_box1,
                              container_proc_time],layout = widgets.Layout(width = "90%",height = '465px',
                                                              align_items = 'stretch' ))

#####################################
#### Tab #2b: Initial Conditions ####
#####################################

label_iniCond = widgets.HTML("""Click the checkbox if you want to initialize the simulation with
none-zero elevation and velocity values.""")

show_initial = widgets.Checkbox(description='Activate initial conditions') # turn on inicial condition checkbox

label_iniCond1 = widgets.HTML("""The Initial surface (Z) file dictates if there is
        a perturbation on the originally flat water surface. The initial U and V velocity files specifies 
        the velocities in the x and y directions. They must be uploaded in 
        the Jupyter Notebook Home folder before running the simulation. Once you have identified their names
        (e.g. Ini_Z.txt), press "Generate Initial Condition Files" to format the uploaded files to the
        FUNWAVE format.""",
        layout = widgets.Layout(width = '90%',
                    size = '20'))

space_box = widgets.Box(layout=widgets.Layout(height ='20px', width='90%')) 
# this ^ box is created to have space among widgets 

    # initial elevation widget container (label and textbox)
label_iniElev = widgets.Label('Initial Elevation Text File:',
                layout =  widgets.Layout(width = '18%'))
iniElev_text = widgets.Text(layout=widgets.Layout(width = "50%",height = '50px'))
container_iniElev = widgets.HBox(children=[label_iniElev,iniElev_text])

    # initial u vel widget container (label and textbox)
label_Uvel = widgets.Label('Initial U Velocity Text File:', 
             layout = widgets.Layout(width = '18%'))
Uvel_text = widgets.Text(layout = iniElev_text.layout)
container_Uvel = widgets.HBox(children=[label_Uvel,Uvel_text])

    # initial v vel widget container (label and textbox)
label_Vvel = widgets.Label('Initial V Velocity Text File:', 
             layout = widgets.Layout(width = '18%'))
Vvel_text = widgets.Text(layout = iniElev_text.layout)
container_Vvel = widgets.HBox(children=[label_Vvel,Vvel_text])

label_iniCond2=widgets.HTML("""<b>NOTE:</b> The data format must be the same as depth file (1 row).""",
                           layout = widgets.Layout(width = "60%",height = '40px'))

# button to format the files to FUNWAVE format
ini_button = widgets.Button(description = 'Generate Initial Condition Files',     
                            layout = widgets.Layout(width = "30%",height = '40px'))   

# box containing the Note and the button
note_button_container = widgets.HBox([label_iniCond2,ini_button],layout = widgets.Layout(width = "90%"))

# initial condition box that appears if the show_initial checkbox is turned on
# the function that shows/hides this continer is at wavemaker_Function.py
init = widgets.VBox(children=[label_iniCond1,space_box,container_iniElev,
                                container_Uvel,container_Vvel,note_button_container])  

# tab 2b container box
page_iniCond = widgets.VBox([label_iniCond,space_box,show_initial, init],
                            layout = widgets.Layout(width = "90%",height = '465px',align_items = 'stretch' ))
  

############################
#### Tab #2c: Wavemaker ####
############################

# wavemaker dropdown widget
label_wave = widgets.Label("""Specification of Wave maker:""",layout = widgets.Layout(height = '35px'))
wave_options = ('Select Wave Maker','LEF_SOL','INI_SOL','INI_REC','WK_REG',
                'JON_1D','TMA_1D','GAUSIAN')  # add WK_IRR and JON_2D for 2D cases
wave_maker = widgets.Dropdown(options=wave_options)
wave_container = widgets.VBox([label_wave,wave_maker],layout = widgets.Layout(height = '100px'))

label_waveMaker = widgets.HTML("""Click <a href="http://udel.edu/~fyshi/FUNWAVE/definition.html" target="_blank">here</a> for
                               the parameter's definitions.""",
                                layout = widgets.Layout(height = '35px'))

#  wavemaker variables
xc_label = widgets.HTML('Xc',layout = widgets.Layout(width = "20%"))
xc = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_xc = widgets.HBox(children=[xc_label,xc],
                            layout = widgets.Layout(height = '45px'))

yc_label = widgets.Label('Yc',layout = widgets.Layout(width = "20%"))
yc = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_yc = widgets.HBox(children=[yc_label,yc],
                            layout = widgets.Layout(height = '45px'))

wid_label = widgets.Label('WID',layout = widgets.Layout(width = "20%"))
wid = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_wid = widgets.HBox(children=[wid_label,wid],
                             layout = widgets.Layout(height = '45px'))

amp_label = widgets.Label('AMP',layout = widgets.Layout(width = "20%"))
amp = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_amp = widgets.HBox(children=[amp_label,amp],
                             layout = widgets.Layout(height = '45px'))

dep_label = widgets.Label('DEP',layout = widgets.Layout(width = "20%"))
dep = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_dep = widgets.HBox(children=[dep_label,dep],
                             layout = widgets.Layout(height = '45px'))

LagTime_label = widgets.Label('LagTime',layout = widgets.Layout(width = "20%"))
LagTime = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_LagTime = widgets.HBox(children=[LagTime_label,LagTime],
                                 layout = widgets.Layout(height = '45px'))

Xwavemaker_label = widgets.Label('X_wavemaker',layout = widgets.Layout(width = "20%"))
Xwavemaker = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_Xwavemaker = widgets.HBox(children=[Xwavemaker_label,Xwavemaker],
                                 layout = widgets.Layout(height = '45px'))

xc_wk_label = widgets.Label('Xc_WK',layout = widgets.Layout(width = "20%"))
xc_wk = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_xc_wk = widgets.HBox(children=[xc_wk_label,xc_wk],
                               layout = widgets.Layout(height = '45px'))

yc_wk_label = widgets.Label('Yc_WK',layout = widgets.Layout(width = "20%"))
yc_wk = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_yc_wk = widgets.HBox(children=[yc_wk_label,yc_wk],
                               layout = widgets.Layout(height = '45px'))

tPeriod_label = widgets.Label('Tperiod',layout = widgets.Layout(width = "20%"))
tPeriod = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_tPeriod = widgets.HBox(children=[tPeriod_label,tPeriod],
                               layout = widgets.Layout(height = '45px'))

ampWK_label = widgets.Label('AMP_WK',layout = widgets.Layout(width = "20%"))
ampWK = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_ampWK = widgets.HBox(children=[ampWK_label,ampWK],
                             layout = widgets.Layout(height = '45px'))

depWK_label = widgets.Label('DEP_WK',layout = widgets.Layout(width = "20%"))
depWK = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_depWK = widgets.HBox(children=[depWK_label,depWK],
                             layout = widgets.Layout(height = '45px'))

thetaWK_label = widgets.Label('Theta_WK',layout = widgets.Layout(width = "20%"))
thetaWK = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_thetaWK = widgets.HBox(children=[thetaWK_label,thetaWK],
                             layout = widgets.Layout(height = '45px'))

TimeRamp_label = widgets.Label('Time_ramp',layout = widgets.Layout(width = "20%"))
TimeRamp = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_TimeRamp = widgets.HBox(children=[TimeRamp_label,TimeRamp],
                             layout = widgets.Layout(height = '45px'))

ywidth_wk_label = widgets.Label('Ywidth_WK',layout = widgets.Layout(width = "20%"))
ywidth_wk = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_ywidth_wk = widgets.HBox(children=[ywidth_wk_label,ywidth_wk],
                               layout = widgets.Layout(height = '45px'))

deltaWK_label = widgets.Label('DELTA_WK',layout = widgets.Layout(width = "20%"))
deltaWK = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_deltaWK = widgets.HBox(children=[deltaWK_label,deltaWK],
                             layout = widgets.Layout(height = '45px'))

FreqPeak_label = widgets.Label('FreqPeak',layout = widgets.Layout(width = "20%"))
FreqPeak = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_FreqPeak = widgets.HBox(children=[FreqPeak_label,FreqPeak],
                             layout = widgets.Layout(height = '45px'))

FreqMin_label = widgets.Label('FreqMin',layout = widgets.Layout(width = "20%"))
FreqMin = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_FreqMin = widgets.HBox(children=[FreqMin_label,FreqMin],
                             layout = widgets.Layout(height = '45px'))

FreqMax_label = widgets.Label('FreqMax',layout = widgets.Layout(width = "20%"))
FreqMax = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_FreqMax = widgets.HBox(children=[FreqMax_label,FreqMax],
                             layout = widgets.Layout(height = '45px'))

HMO_label = widgets.Label('Hmo',layout = widgets.Layout(width = "20%"))
HMO = widgets.BoundedFloatText(layout = widgets.Layout(width = "30%"))
container_HMO = widgets.HBox(children=[HMO_label,HMO],
                             layout = widgets.Layout(height = '45px'))

GammaTMA_label = widgets.Label('GammaTMA',layout = widgets.Layout(width = "20%"))
GammaTMA = widgets.BoundedFloatText(value = 3.3,layout = widgets.Layout(width = "30%"))
container_GammaTMA = widgets.HBox(children=[GammaTMA_label,GammaTMA],
                                 layout = widgets.Layout(height = '45px'))

ThetaPeak_label = widgets.Label('ThetaPeak',layout = widgets.Layout(width = "20%"))
ThetaPeak = widgets.BoundedFloatText(value = 0.0,layout = widgets.Layout(width = "30%"))
container_ThetaPeak = widgets.HBox(children=[ThetaPeak_label,ThetaPeak],
                                 layout = widgets.Layout(height = '45px'))

NFreq_label = widgets.Label('NFreq',layout = widgets.Layout(width = "20%"))
NFreq = widgets.BoundedFloatText(value = 45,layout = widgets.Layout(width = "30%"))
container_NFreq = widgets.HBox(children=[NFreq_label,NFreq],
                                 layout = widgets.Layout(height = '45px'))

NTheta_label = widgets.Label('NTheta',layout = widgets.Layout(width = "20%"))
NTheta = widgets.BoundedFloatText(value = 24,layout = widgets.Layout(width = "30%"))
container_NTheta = widgets.HBox(children=[NTheta_label,NTheta],
                                 layout = widgets.Layout(height = '45px'))

SigmaTheta_label = widgets.Label('Sigma_Theta',layout = widgets.Layout(width = "20%"))
Sigma_Theta = widgets.BoundedFloatText(value = 24,layout = widgets.Layout(width = "30%"))
container_SigmaTheta = widgets.HBox(children=[SigmaTheta_label,Sigma_Theta],
                                 layout = widgets.Layout(height = '45px'))

### containers for each wavemaker's respective variables:
# the function that shows/hides this continers is at wavemaker_Function.py

container_IniRec = widgets.VBox(children=[container_xc,container_yc,container_wid],
                         layout = widgets.Layout(align_items = 'stretch' ))

container_Gauss = widgets.VBox(children=[container_xc,container_yc,container_wid,container_amp],
                         layout = widgets.Layout(align_items = 'stretch' ))

container_LefSol = widgets.VBox(children=[container_amp,container_dep,container_LagTime],
                         layout = widgets.Layout(align_items = 'stretch' ))

container_IniSol = widgets.VBox(children=[container_amp,container_dep,container_LagTime,container_Xwavemaker],
                         layout = widgets.Layout(align_items = 'stretch' ))

container_WkReg = widgets.VBox(children=[container_xc_wk,container_yc_wk,container_tPeriod,container_ampWK,
                         container_depWK,container_thetaWK,container_TimeRamp],
                         layout = widgets.Layout(align_items = 'stretch' ))

container_JON2D = widgets.VBox(children=[container_xc_wk,container_yc_wk,container_ywidth_wk],
                         layout = widgets.Layout(align_items = 'stretch' ))

container_JON1D_col1 = widgets.VBox(children=[container_xc_wk,container_yc_wk,container_ywidth_wk,
                                             container_depWK,container_TimeRamp,container_deltaWK],
                         layout = widgets.Layout(width = '50%'))
container_JON1D_col2 = widgets.VBox(children=[container_FreqPeak,container_FreqMin,container_FreqMax,container_HMO,
                         container_GammaTMA,container_NFreq],
                         layout = widgets.Layout(width = '50%'))
container_JON1D = widgets.HBox([container_JON1D_col1,container_JON1D_col2],
                              layout = widgets.Layout(width = '90%'))

container_WkIrr_col1 = widgets.VBox(children=[container_xc_wk,container_yc_wk,container_depWK,container_TimeRamp,
                         container_ywidth_wk,container_deltaWK,container_FreqPeak,container_SigmaTheta],
                         layout = widgets.Layout(width = '50%'))
container_WkIrr_col2 = widgets.VBox(children=[container_FreqMin,container_FreqMax,container_HMO,
                         container_GammaTMA,container_ThetaPeak,container_NFreq,container_NTheta],
                         layout = widgets.Layout(width = '50%'))
container_WKIRR = widgets.HBox([container_WkIrr_col1,container_WkIrr_col2],
                              layout = widgets.Layout(width = '90%'))

container_TMA_1D_col1 = widgets.VBox(children=[container_xc_wk,container_yc_wk,container_depWK,
                         container_TimeRamp,container_ywidth_wk,container_deltaWK],
                         layout = widgets.Layout(width = '50%'))
container_TMA_1D_col2 = widgets.VBox(children=[container_FreqPeak,container_FreqMin,
                         container_FreqMax,container_HMO,container_GammaTMA,container_NFreq],
                         layout = widgets.Layout(width = '50%'))
container_TMA_1D = widgets.HBox([container_TMA_1D_col1,container_TMA_1D_col2],
                               layout = widgets.Layout(width = '90%'))

# tab 2c container box
page_waveMaker = widgets.VBox(children=[wave_container,label_waveMaker,container_IniRec,
                             container_Gauss,container_LefSol,container_IniSol,
                             container_WkReg,container_JON1D,container_JON2D,
                             container_WKIRR,container_TMA_1D],
                             layout = widgets.Layout(height = '465px',width = '90%'))

###############################
#### Tab #2d: Sponge layer ####
###############################

# sponge layer tab intro label
sponge_label = widgets.HTML("""FUNWAVE possess a DHI type sponge layer. 
                        The user needs to specify the widths of four boundaries 
                        and parameters.""",layout = widgets.Layout(width = '90%',height = '55px')) 

## column 1 of sponge layer tab:

# diffusion sponge widget container (label and checkbox) 
dif = widgets.Checkbox(description='Diffusion Sponge',value=False,layout=widgets.Layout(height = '45px'))

# friction sponge widget container (label and checkbox) 
fric = widgets.Checkbox(description='Friction Sponge',value=False,layout=dif.layout)

# direct sponge widget container (label and checkbox) 
Dir = widgets.Checkbox(description='Direct Sponge',value=False,layout=dif.layout)

container_SpongeLayer_column1 = widgets.VBox(children=[dif,fric,Dir],
                                layout = widgets.Layout(width = '42%',align_items = 'flex-start'))

## column 2 of sponge layer tab:
# diffusion coefficient widget container (label and textbox)
label_CDsponge = widgets.Label('Friction Coefficient', layout=widgets.Layout(width = "70%"))
CDsponge_text = widgets.BoundedFloatText(layout = widgets.Layout(width = "55%"))
container_CDsponge = widgets.VBox(children=[label_CDsponge,CDsponge_text],
                                 layout=widgets.Layout(height = '90px'))

# friction coefficient widget container (label and textbox)
label_csp = widgets.Label('Diffusion Coefficient', layout = label_CDsponge.layout)
csp_text = widgets.BoundedFloatText(layout = widgets.Layout(width = "55%"))
container_csp = widgets.VBox(children=[label_csp,csp_text],
                            layout=widgets.Layout(height = '90px'))

# decay rate widget container (label and textbox)
label_R_sponge = widgets.Label('Decay Rate', layout = label_CDsponge.layout)
R_sponge_text = widgets.BoundedFloatText(layout=widgets.Layout(width = '55%'),value ='0.85', max = '0.95',
                                    min='0.85',step = '0.01')
container_R_sponge = widgets.VBox(children=[label_R_sponge,R_sponge_text],
                                 layout=widgets.Layout(height = '90px'))

# max decay rate widget container (label and textbox)
label_A_sponge = widgets.Label('Max Decay Rate', layout = label_CDsponge.layout)
A_sponge_text = widgets.BoundedFloatText(layout = widgets.Layout(width = "55%"),value ='5',min='5',step = '0.1')
container_A_sponge = widgets.VBox(children=[label_A_sponge,A_sponge_text],
                                 layout=widgets.Layout(height = '90px'))

container_SpongeLayer_column2 = widgets.VBox(children=[container_csp,container_CDsponge,
                                container_R_sponge,container_A_sponge],
                                layout = widgets.Layout(width = '30%'))

## Column 3 of sponge layer tab:
# sponge layer note label
sponge_note = widgets.HTML("""<b>NOTE: </b>Set width = 0.0 if no sponge.""",
                            layout=sponge_label.layout)

# Left sponge width (LSW) widget container (label and textbox)
label_LSW = widgets.Label('Left Width', layout=widgets.Layout(width = "39%"))
LSW_text = widgets.BoundedFloatText(layout = widgets.Layout(width = "50%"),min ='0.0',step = '0.01')
container_LSW = widgets.HBox(children=[label_LSW,LSW_text],
                             layout=widgets.Layout(width = '100%',height = '65px'))

# Right sponge width (RSW) widget container (label and textbox)
label_RSW = widgets.Label('Right Width', layout = label_LSW.layout)
RSW_text = widgets.BoundedFloatText(layout = LSW_text.layout,min ='0.0',step = '0.01')
container_RSW = widgets.HBox(children=[label_RSW,RSW_text],layout=container_LSW.layout)

# Top sponge width (TSW) widget container (label and textbox)
label_TSW = widgets.Label('Top Width', layout = label_LSW.layout)
TSW_text = widgets.BoundedFloatText(layout = LSW_text.layout,min ='0.0',step = '0.01')
container_TSW = widgets.HBox(children=[label_TSW,TSW_text],layout=container_LSW.layout)

# Bottom sponge width (BSW) widget container (label and textbox)
label_BSW = widgets.Label('Bottom Width', layout = label_LSW.layout)
BSW_text = widgets.BoundedFloatText(layout = LSW_text.layout,min ='0.0',step = '0.01')
container_BSW = widgets.HBox(children=[label_BSW,BSW_text],layout=container_LSW.layout)

container_SpongeLayer_column3 = widgets.VBox(children=[sponge_note,container_LSW,container_RSW,container_TSW,
                                                       container_BSW],layout = widgets.Layout(width = '30%'))

# create box with all columns
container_totalColumns = widgets.HBox(children=[container_SpongeLayer_column1,
                                                container_SpongeLayer_column2,
                                                container_SpongeLayer_column3],
                                      layout = widgets.Layout(height = '550px',width = '90%',
                                                              align_items = 'stretch' ))

# tab 2d container box
page_spongeLayer = widgets.VBox([sponge_label,container_totalColumns],
                                layout = widgets.Layout(width = "90%",height = '465px'))

#######################################################################
#### NOTE: Tab #2e: Output Options is located in PrincipalTab_3.py ####
#######################################################################


#######################################
#### Tab #2f: Input File Generator ####
#######################################

inputFile_label = widgets.HTML("""Press button to review/update <b>input.txt values</b>
 before generating the file.""",layout = widgets.Layout(width = "50%",height = '25px'))

update_input_button = widgets.Button(description = 'Review Input Values',layout = widgets.Layout(width = "50%",
                                  height = '40px'))

inputUpdate_box = widgets.HBox([inputFile_label,update_input_button],
                               layout = widgets.Layout(width = "90%",height = '55px'))

# print input.txt for verification 
input_verification = widgets.HTML() 

inputFile_box = widgets.Box([input_verification],layout = widgets.Layout(height = '295px',width = '90%',
                                                                        border='solid 2px grey'))
# generate input file button
inputFile_label2 = widgets.HTML("""If you are satisfied with your input values, press the Generate Input File button.""",
                                layout = widgets.Layout(width = "50%",height = '25px'))
inputFile_button = widgets.Button(description = 'Generate Input File',layout = widgets.Layout(width = "50%",
                                  height = '40px'))

inputGen_box = widgets.HBox([inputFile_label2,inputFile_button],
                               layout = widgets.Layout(width = "90%",height = '55px'))

page_GenInput = widgets.VBox([inputUpdate_box,space_box1,inputFile_box,space_box1,inputGen_box],
                                 layout = page_inputIntro.layout)
