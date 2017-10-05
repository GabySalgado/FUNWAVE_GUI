import ipywidgets as widgets
from traitlets import link


#######################################
#### Principal Tab 4: Run Funwave  ####
#######################################



#intro label
runFUNWAVE_label = widgets.Label('Press Button to Run FUNWAVE Model.',layout = widgets.Layout(height = '90px'))

# run funwave button
RunFunwave_button = widgets.Button(description="Run FUNWAVE Model",                           
                             layout = widgets.Layout(width = '80%',height = '45px'))
# column 1: run funwave buttom
fun_cont = widgets.HBox(children = [RunFunwave_button],
                       layout = widgets.Layout(width = '50%',height = '100px'))         

# stop fun button
abortFun_button = widgets.Button(description="Abort Model",
                                 layout = widgets.Layout(width = '80%',height = '45px',button_style='danger'))

# column 2: stop funwave buttom
abort_cont = widgets.HBox(children = [abortFun_button],layout = widgets.Layout(width = '50%',height = '100px'))        

# columns container
runfun_container = widgets.HBox(children = [fun_cont,abort_cont],
                                layout = widgets.Layout(width = '100%',height = '145px'))   


# run funwave load progress bar 
Prog_label = widgets.Label('Model Progress:',layout = widgets.Layout(width = '12%'))
RunFunwave_load = widgets.FloatProgress(value=0.0,min=5.0,max=10.0,step=0.1,
                                        layout = widgets.Layout(width = '88%',height = '35px'))
Load_cont = widgets.HBox([Prog_label,RunFunwave_load],layout = widgets.Layout(width = '90%',height = '137px'))

# create run_funwave tab
RUNfunwave_tabs = widgets.VBox(children=[runFUNWAVE_label,runfun_container,Load_cont])
