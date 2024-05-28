# coding:utf-8
#!/usr/bin/env python
'''Runs the experiment.
Author: MdB 04/2017 s
'''

import pygame, time, pickle, csv,sys
import os.path
from psychopy import core, gui, event
from task_sustAttnWM import Task
from datetime import datetime
from shutil import copyfile
import settings_sustAttnWM as settings
import helper_functions_sustAttnWM as hf

if sys.getdefaultencoding() != 'utf-8':
    reload(sys)
    sys.setdefaultencoding('utf-8')

#get basic specifications for this experiment from GUI
gui_specs = hf.get_specs()
print "Subject: ", gui_specs['subj_name']

#set up the experimental display (psychopy window specs)
#expt_design = settings.ExperimentalDesign(expt_display,gui_specs)

#set up the experimental design
if not gui_specs['restarting'] or gui_specs['debug']:
    #create experiment design
    expt_design = settings.ExperimentalDesign(gui_specs)
else:
    filename_expdat = gui_specs['save_dir'] + gui_specs['subj_name'] + '_expdat.p'
    filename_expdesign = gui_specs['save_dir']+ gui_specs['subj_name'] + '_expdesign.p'
    if os.path.isfile(filename_expdat):
        expt_design = pickle.load( open( filename_expdat, "rb" ) )
    elif os.path.isfile(filename_expdesign):
        expt_design = pickle.load( open( filename_expdesign, "rb" ) )
    else:
        expt_design = settings.ExperimentalDesign(gui_specs)
    
    #overwrite the specs to match what you specified in the gui for this time
    expt_design.restarting = gui_specs['restarting']
    expt_design.env = gui_specs['env']
    expt_design.exptDate = gui_specs['expt_date']
    expt_design.debug = gui_specs['debug']

expt_display = settings.ExperimentalDisplay(expt_design,gui_specs)

#design experiment
exp = Task(expt_design,expt_display)

#save files
filename = gui_specs['save_dir'] + gui_specs['subj_name'] + '_' + 'expdat.p'
pickle.dump( exp.dsgn, open( filename, "wb" ) )

#open csv files to save
f = hf.open_csv_data_file(gui_specs,gui_specs['subj_name'] + '_explog')
f_cd = hf.open_csv_data_file(gui_specs,gui_specs['subj_name'] + '_explog_cd')

#calculate the number of first trial
filename = gui_specs['save_dir'] + gui_specs['subj_name'] + '_trialnum.p'
if gui_specs['restarting'] and os.path.isfile(filename):
    trial0 = pickle.load( open( filename, "rb" ) )
    exp.dsgn.counter_trial = trial0

if __name__ == '__main__':
    
    #initialize the clock
    clock = core.Clock()
    clock.reset()
    
    #set mouse to not be visible 
    if not gui_specs['debug']:
        exp.disp.win.setMouseVisible(False)

    #welcome screen
    exp.welcomeToExperiment()
    
    ######################################################################
    ## Encoding/WM Phase
    
    if gui_specs['sustattn']:
        
        #working memory long-term memory instructions
        if exp.dsgn.counter_trial[0]==0:
            exp.instructionsSustAttn(clock)
            exp.instructionsWholeReport(clock)
            exp.instructionsBoth(clock)
        
        #block loop
        block0 = 0 #FIX THIS TO RESTART BLOCKS?
        for iblock in range(block0,exp.dsgn.nblocks):
            
            exp.startOfBlock(iblock,clock,f=f)
            if iblock == 0:
                hf.write_to_csv_data_file(f,0,0,exp.dsgn,headerLine=True)
                
            #set mouse to not be visible 
            if not gui_specs['debug']:
                exp.disp.win.setMouseVisible(False)
            
            #trial loop
            for itrial in range(exp.dsgn.ntrials_perblock):
                exp.encodingArray(iblock,itrial,clock)
                if exp.dsgn.probe_trials[iblock,itrial]==1:
                    exp.retentionInterval(iblock,itrial,clock)
                    exp.memoryProbe(iblock,itrial,clock,f=f)
                    exp.itiWM(iblock,itrial,clock)
                    
                #update counter
                exp.dsgn.counter_trial[iblock] = exp.dsgn.counter_trial[iblock]+1
               
                #save files
                filename = gui_specs['save_dir'] + gui_specs['subj_name'] + '_' + 'expdat.p'
                pickle.dump( exp.dsgn, open( filename, "wb" ) )
                hf.write_to_csv_data_file(f,iblock,itrial,exp.dsgn) 
            exp.endOfBlock(iblock,clock,f=f)
        
        
        #close and copy explog for this section of the experiment
        today = datetime.now()
        filename_old = gui_specs['save_dir'] + gui_specs['subj_name'] + '_explog.csv'
        filename_new = gui_specs['save_dir'] + gui_specs['subj_name'] + '_explog_' + today.strftime('%Y%m%d_%H%M%S') + '.csv'
        copyfile(filename_old,filename_new)
    
    
    ######################################################################
    ## Change Detection Phase
    
    
    # end of experiment screen 
    exp.endOfExperiment(f=f)