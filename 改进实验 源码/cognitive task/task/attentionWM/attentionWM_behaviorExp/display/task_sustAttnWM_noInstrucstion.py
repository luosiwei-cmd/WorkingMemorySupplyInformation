# coding:utf-8
import sys
if 'psychopy' in sys.modules:
    from psychopy import visual, core, tools, event, monitors, logging #, parallel #, gui  
from helper_functions_sustAttnWM import *
import numpy as np
from pylink import *
from scipy.spatial.distance import pdist

import serial
ser = serial.Serial('COM3', 9600)
#from EyeLinkCoreGraphicsPsychoPy import EyeLinkCoreGraphicsPsychoPy

if sys.getdefaultencoding() != 'utf-8':
    reload(sys)
    sys.setdefaultencoding('utf-8')
    
#totolcorrecitem = np.zeros(60)
WRacc = 0

#all the data files 
class Task(object):

    def __init__(self,expt_design,expt_display):

        self.dsgn = expt_design
        self.disp = expt_display             
        
    def welcomeToExperiment(self,f=None):

        self.disp.text_welcome.draw()           #welcome text at top of screen
        self.disp.text_advance.draw()           #advance text at bottom of screen
        self.disp.fix.draw()
        drawpath(self.disp,self.dsgn.nblocks+1)   #draw the path of the experiment
        self.disp.win.flip()                    #flip the window
        
        #if not self.dsgn.debug:
        wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        #空屏
        self.disp.fix.draw()
        self.disp.win.flip()
        core.wait(.5)



    def startOfBlock(self,iblock,clock,f=None):

        #display the start of block screen
        #self.disp.text_block_start.setText('Block ' + str(iblock+1) + ' of ' + str(self.dsgn.nblocks+1))   #set the block text to be the right block number
        self.disp.text_block_start.setText('第 ' + str(iblock+1) + ' 组，共 ' + str(self.dsgn.nblocks) + ' 组')
        self.disp.text_block_start.draw()                       #draw the block text
        self.disp.fix.color = (-1,-1,-1)
        self.disp.fix.draw()
        drawpath(self.disp,self.dsgn.nblocks+1,i_fill=iblock-1)   #display progress
        self.disp.text_advance.pos = (0,-6.5)
        self.disp.text_advance.draw()                           #text advance
        self.disp.win.flip()                                    #flip the window
        
        
        
        #wait for the response to start the block
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.dsgn.time_blockstart[iblock]=clock.getTime()


        self.disp.fix.draw()
        self.disp.win.flip()
        core.wait(self.dsgn.t_poststim)
        ser.write('1')


    def encodingArray(self,iblock,itrial,clock,f=None):

        self.disp.fix.color = (-1,-1,-1)
        self.disp.fix.draw()
        if self.dsgn.freq_trials[iblock,itrial]==1:
            for i in range(self.dsgn.setsize):
                self.disp.circle.fillColor = self.dsgn.trials_colorrgb[iblock,itrial,i]
                self.disp.circle.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                self.disp.circle.draw()
        else:
            for i in range(self.dsgn.setsize):
                self.disp.square.fillColor = self.dsgn.trials_colorrgb[iblock,itrial,i]
                self.disp.square.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                self.disp.square.draw()            
        self.disp.win.flip()

        #save timing information
        self.dsgn.time_encarray_onset[iblock,itrial] = clock.getTime()
        
        response = None
        event.clearEvents(eventType='keyboard')
        while clock.getTime()<(self.dsgn.time_encarray_onset[iblock,itrial]+self.dsgn.t_stim):
            if response is None:
                response = poll_for_response(self.disp.win,self.dsgn.key_freq+self.dsgn.key_infreq+self.dsgn.key_escape)
                if response is not None:
                    #log time resposne was made
                    self.dsgn.time_response[iblock,itrial] = clock.getTime()
                    self.dsgn.rts[iblock,itrial] = self.dsgn.time_response[iblock,itrial]-self.dsgn.time_encarray_onset[iblock,itrial]

                    #log response that was made
                    self.dsgn.actual_response_string[iblock][itrial] = response[0]
                    
                    if response==self.dsgn.key_freq[0]:
                        self.dsgn.actual_response[iblock,itrial]=1
                    elif response==self.dsgn.key_infreq[0]:
                        self.dsgn.actual_response[iblock,itrial]=0

                    #see if the response was correct or not
                    self.dsgn.acc[iblock,itrial]=(self.dsgn.correct_response[iblock,itrial]==self.dsgn.actual_response[iblock,itrial])*1

                    #redraw the display
                    self.disp.fix.color = (1,1,1)
                    self.disp.fix.draw()
                    if self.dsgn.freq_trials[iblock,itrial]==1:
                        for i in range(self.dsgn.setsize):
                            self.disp.circle.fillColor = self.dsgn.trials_colorrgb[iblock,itrial,i]
                            self.disp.circle.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                            self.disp.circle.draw()
                    else:
                        for i in range(self.dsgn.setsize):
                            self.disp.square.fillColor = self.dsgn.trials_colorrgb[iblock,itrial,i]
                            self.disp.square.pos = (self.dsgn.positions_x[i], self.dsgn.positions_y[i])
                            self.disp.square.draw()            
                    self.disp.win.flip()                    


    def retentionInterval(self,iblock,itrial,clock):
        self.disp.fix.color = (-1,-1,-1)
        self.disp.fix.draw()
        
#        if self.disp.env=='eeg':
#            parallel.setData(41)
            
        #ser.write('2')
        
        self.disp.win.flip()

        #save timing information
        self.dsgn.time_encarray_offset[iblock,itrial] = clock.getTime()

        core.wait(self.dsgn.t_poststim)
        
#        if self.disp.env=='eeg':
#            parallel.setData(0)
        

    def memoryProbe(self,iblock,itrial,clock,f=None):

        global WRacc
        
        #draw response screen                   
        self.disp.fix.draw()
        for i in range(self.dsgn.setsize):
            for j in range(9):
                self.disp.wholereportsquare[i,j].draw()
        self.disp.win.flip()
        
        ser.write('2')

        #save timing information
        self.dsgn.time_memprobe_onset[iblock,itrial] = clock.getTime()

        #display mouse at fixation
        self.disp.mouse.setPos(newPos=(0,0))
        self.disp.mouse.setVisible(1)

        count = 0
        #sum = 0
        event.clearEvents(eventType='keyboard')
        #totolcorrecitem = np.zeros(self.dsgn.setsize)
        while np.any(np.isnan(self.dsgn.wholereport_resp[iblock,itrial])):
            response = poll_for_response(self.disp.win,self.dsgn.key_escape)
            for i in range(6):
                if np.isnan(self.dsgn.wholereport_resp[iblock,itrial][i]):
                    for j in range(9):
                        if self.disp.mouse.isPressedIn(self.disp.wholereportsquare[i,j],buttons=[0]):
                            #log the time that they responded to that item
                            self.dsgn.wholereport_rts[iblock,itrial][i] = clock.getTime()-self.dsgn.time_memprobe_onset[iblock,itrial]

                            #which color within that square did they click
                            self.dsgn.wholereport_respcolorind[iblock,itrial][i] = j
                            self.dsgn.wholereport_resprgb[iblock,itrial][i] = self.dsgn.colors_rgb[j]
                            self.dsgn.wholereport_respacc[iblock,itrial][i] = (self.dsgn.wholereport_respcolorind[iblock,itrial][i]==
                                                                            self.dsgn.trials_colorind[iblock,itrial,i])*1
                            if self.dsgn.wholereport_respcolorind[iblock,itrial][i]==self.dsgn.trials_colorind[iblock,itrial,i]:
                                WRacc = WRacc+1
                                #totolcorrecitem[i] = 1
                            #WRacc = int(np.sum(totolcorrecitem)) + WRacc
                            #sum = WRacc + sum
                            #log that they have responded to that item
                            self.dsgn.wholereport_resp[iblock,itrial][i] = 1

                            #log the order that they responded to that item
                            count = count+1
                            self.dsgn.wholereport_resporder[iblock,itrial][i] = count

                            #redraw the screen
                            self.disp.fix.draw()
                            for i in range(self.dsgn.setsize):
                                for j in range(9):
                                    self.disp.wholereportsquare[i,j].draw()
                                    #add black lines around the responded item
                                    if self.dsgn.wholereport_resp[iblock,itrial][i]==1:
                                        self.disp.square_pos[i].draw()
                            self.disp.win.flip()

        #remove memory probe after response
        self.disp.fix.draw()
        self.disp.win.flip()
        
        ser.write('3')

        #save timing information
        self.dsgn.time_memprobe_offset[iblock,itrial] = clock.getTime()
        self.disp.mouse.setVisible(0)

    def itiWM(self,iblock,itrial,clock):
        self.disp.fix.draw()
        self.disp.win.flip()

#        if self.disp.env=='eeg':
#            parallel.setData(71)

        #save timing information
        self.dsgn.time_iti_onset[iblock,itrial] = clock.getTime()

        core.wait(self.dsgn.t_iti)
        
#        if self.disp.env=='eeg':
#            parallel.setData(0)

        #save timing information
        self.dsgn.time_iti_offset[iblock,itrial] = clock.getTime()

        
    def endOfBlock(self,iblock,clock,f=None):

        #display the start of block screen
        #self.disp.text_block_start.setText('完成第一部分')   #set the block text to be the right block number
        #print int(np.sum(totolcorrecitem))
        self.disp.text_block_end.setText('第 ' + str(iblock+1) + ' 组结束，请休息一下')
        #self.disp.text_block_start.draw()
        self.disp.text_block_end.draw()
        #self.disp.text_blockAvgCorrect.setText('你在该部分的正确率为' + str(int(np.round(np.nanmean(self.dsgn.wm_acc[iblock])*100))) + ' %')
        self.disp.text_blockAvgCorrect.setText('您的累计得分为：' + str(WRacc) + ' 分')
        self.disp.text_blockAvgCorrect.draw()
        self.disp.fix.color = (-1,-1,-1)
        self.disp.fix.draw()
        drawpath(self.disp,self.dsgn.nblocks+1,i_fill=iblock-1)   #display progress
        self.disp.text_advance.pos = (0,-6.5)
        self.disp.text_advance.draw()                           #text advance
        self.disp.win.flip()                                    #flip the window
        
        #wait for the response to start the block
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.dsgn.time_blockstart[iblock]=clock.getTime()
        
        self.disp.fix.draw()
        self.disp.win.flip()
        #self.dsgn.time_blockend[iblock] = clock.getTime()
        #core.wait(.5)
        core.wait(self.dsgn.t_poststim)


    def changeDetectionInstructionsIntro(self,f=None):
        
        ##################################################
        ##INSTRUCTIONS SCREEN 0
        #display text
        self.disp.text_CDinstructions0.draw()
        
        #draw_example square
        #self.disp.square.pos = (0,0)
        self.disp.square.pos = (10,10)
        color_ind = np.random.randint(0,self.dsgn.cd_ncolors-3)
        #self.disp.square.fillColor = self.dsgn.cd_colors_rgb[color_ind]
        self.disp.square.fillColor = (0,0,0)
        self.disp.square.draw()
        #self.disp.text_advance.draw()
        self.disp.win.flip()
        
        #wait for response
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.disp.win.flip()
        core.wait(.5)
        
        ##################################################
        ##INSTRUCTIONS SCREEN 1
        #display text
        self.disp.text_CDinstructions1.draw()
        
        #draw example_square
        self.disp.square.draw()
        self.disp.text_advance.draw()
        self.disp.win.flip()

        #wait for response
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.disp.win.flip()
        core.wait(.5)
        
        ##################################################
        ##INSTRUCTIONS SCREEN 2
        #display text
        self.disp.text_CDinstructions2.draw()
        
        #draw example_square
        new_color_ind = np.random.choice(np.setdiff1d(np.arange(self.dsgn.cd_ncolors-3),color_ind),replace=False)
        self.disp.square.fillColor = self.dsgn.cd_colors_rgb[new_color_ind]
        self.disp.square.draw()
        self.disp.text_advance.draw()
        self.disp.win.flip()

        #wait for response
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.disp.win.flip()
        core.wait(.5)
        
        ##################################################
        ##INSTRUCTIONS SCREEN 3
        #display text
        self.disp.text_CDinstructions3.draw()
        
        #draw example_square
        self.disp.square.fillColor = self.dsgn.cd_colors_rgb[color_ind]
        self.disp.square.draw()
        
        #draw other things
        self.disp.text_advance.draw()
        self.disp.win.flip()

        #wait for response
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.disp.win.flip()
        core.wait(.5)

    def changeDetectionInstructionsPracticeSS1(self,f=None):
        #display instructions screen
        self.disp.text_CDinstructions4.draw()
        #drawpath(self.disp,self.dsgn.nblocks)
        self.disp.text_advance.draw()
        self.disp.win.flip()

        #wait for response
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.disp.win.flip()
        core.wait(.5)
        
        acc = 0
        
        while acc == 0:
        
            #encoding array
            self.disp.fix.draw()
            ss = 1
            r =  self.dsgn.cd_rad_min + ( self.dsgn.cd_rad_max - self.dsgn.cd_rad_min)*np.random.random(ss) #radius of eccentricity for square position
            theta = np.random.random(ss)*90+90*np.random.randint(0,4) #[degrees] for square position
            x, y = tools.coordinatetools.pol2cart(theta,r,units='deg')
            self.disp.square.pos = (x[0],y[0])
            color_ind = np.random.randint(0,self.dsgn.cd_ncolors-3)
            self.disp.square.fillColor = self.dsgn.cd_colors_rgb[color_ind]
            self.disp.square.draw()
            self.disp.win.flip()
            core.wait(self.dsgn.cd_t_encarray)
            
            #retention interval
            self.disp.fix.draw()
            self.disp.win.flip()
            core.wait(self.dsgn.cd_t_retinterval)
            
            #memory probe
            self.disp.fix.draw()
            self.disp.square.draw()
            self.disp.text_diff.draw()  #draw response for different trials
            self.disp.text_same.draw()  #draw response for same trials
            self.disp.win.flip()        #flip window
            
            #wait for responses
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_same+self.dsgn.key_diff+self.dsgn.key_escape,f=f)
            
            if [response] == self.dsgn.key_same:
                self.disp.text_CDinstructionsCorrect.draw()
                acc = 1
            else:
                self.disp.text_CDinstructionsIncorrect.draw()
            self.disp.text_advance.draw()
            self.disp.win.flip()
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
            self.disp.win.flip()
            core.wait(.5)

    def changeDetectionInstructionsPracticeSS2Diff(self,f=None):
        
        #display instructions screen
        self.disp.text_CDinstructions5.draw()
        #drawpath(self.disp,self.dsgn.nblocks)
        self.disp.text_advance.draw()
        self.disp.win.flip()

        #wait for response
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
        self.disp.win.flip()
        core.wait(.5)
        
        acc = 0
        
        while acc == 0:
        
            #encoding array
            self.disp.fix.draw()
            ss = 2
            current_min_dist = 0
            while current_min_dist<self.dsgn.cd_stim_mindist:
                    
                #polar angle for square position
                r =  self.dsgn.cd_rad_min + ( self.dsgn.cd_rad_max - self.dsgn.cd_rad_min)*np.random.random(ss) #radius of eccentricity for square position
                theta = np.random.random(ss)*90+90*np.random.choice(np.arange(4),size=ss,replace=False)#[degrees] for square position
               
                #cartesian coordinates corresponding to polar position for square position
                x, y = tools.coordinatetools.pol2cart(theta,r,units='deg')
                
                #check that they all satisfy min_dist
                current_min_dist = np.min(pdist(np.transpose(np.vstack((x,y)))))
                
            color_ind = np.random.choice(np.arange(self.dsgn.cd_ncolors),size=ss,replace=False)
            for i in range(ss):
                self.disp.square.pos = (x[i],y[i])
                self.disp.square.fillColor = self.dsgn.cd_colors_rgb[color_ind[i]]
                self.disp.square.draw()
            self.disp.win.flip()
            core.wait(self.dsgn.cd_t_encarray)
            
            #retention interval
            self.disp.fix.draw()
            self.disp.win.flip()
            core.wait(self.dsgn.cd_t_retinterval)
            
            #memory probe
            self.disp.fix.draw()
            self.disp.square.fillColor = self.dsgn.cd_colors_rgb[color_ind[0]]
            self.disp.square.draw()
            self.disp.text_diff.draw()  #draw response for different trials
            self.disp.text_same.draw()  #draw response for same trials
            self.disp.win.flip()        #flip window
            
            #wait for responses
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_same+self.dsgn.key_diff+self.dsgn.key_escape,f=f)
            
            if [response] == self.dsgn.key_diff:
                self.disp.text_CDinstructionsCorrect.draw()
                acc = 1
            else:
                self.disp.text_CDinstructionsIncorrect.draw()
            self.disp.text_advance.draw()
            self.disp.win.flip()
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
            self.disp.win.flip()
            core.wait(.5)

    def changeDetectionInstructionsPracticeSS6Same(self,f=None):
            
            #display instructions screen
            self.disp.text_CDinstructions6.draw()
            self.disp.text_advance.draw()
            self.disp.win.flip()

            #wait for response
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
            self.disp.win.flip()
            core.wait(.5)
            
            acc = 0
            
            while acc == 0:
            
                #encoding array
                self.disp.fix.draw()
                ss = 6
                current_min_dist = 0
                while current_min_dist<self.dsgn.cd_stim_mindist:
                        
                    #polar angle for square position
                    r =  self.dsgn.cd_rad_min + ( self.dsgn.cd_rad_max - self.dsgn.cd_rad_min)*np.random.random(ss) #radius of eccentricity for square position
                    theta = np.random.random(ss)*90+90*np.random.choice(np.append(np.arange(4),np.arange(4)),size=ss,replace=False)#[degrees] for square position
                   
                    #cartesian coordinates corresponding to polar position for square position
                    x, y = tools.coordinatetools.pol2cart(theta,r,units='deg')
                    
                    #check that they all satisfy min_dist
                    current_min_dist = np.min(pdist(np.transpose(np.vstack((x,y)))))
                    
                color_ind = np.random.choice(np.arange(self.dsgn.cd_ncolors),size=ss,replace=False)
                for i in range(ss):
                    self.disp.square.pos = (x[i],y[i])
                    self.disp.square.fillColor = self.dsgn.cd_colors_rgb[color_ind[i]]
                    self.disp.square.draw()
                self.disp.win.flip()
                core.wait(self.dsgn.cd_t_encarray)
                
                #retention interval
                self.disp.fix.draw()
                self.disp.win.flip()
                core.wait(self.dsgn.cd_t_retinterval)
                
                #memory probe
                self.disp.fix.draw()
                self.disp.square.draw()
                self.disp.text_diff.draw()  #draw response for different trials
                self.disp.text_same.draw()  #draw response for same trials
                self.disp.win.flip()        #flip window
                
                #wait for responses
                response = wait_for_response(self.disp.win,keylist=self.dsgn.key_same+self.dsgn.key_diff+self.dsgn.key_escape,f=f)
                
                if [response] == self.dsgn.key_same:
                    self.disp.text_CDinstructionsCorrect.draw()
                    acc = 1
                else:
                    self.disp.text_CDinstructionsIncorrect.draw()
                self.disp.text_advance.draw()
                self.disp.win.flip()
                response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
                self.disp.win.flip()
                core.wait(.5)
                
            #display instructions screen
            self.disp.text_CDinstructions7.draw()
            self.disp.win.flip()

            #wait for response
            response = wait_for_response(self.disp.win,keylist=self.dsgn.key_advance+self.dsgn.key_escape,f=f)
            self.disp.win.flip()
            core.wait(.5)

    def startOfChangeDetection(self,iblock,clock,f=None):
        
        #display the start of block screen
        self.disp.text_block_start.setText('Start of last block')   #set the block text to be the right block number
        self.disp.text_block_start.draw()                       #draw the block text
        drawpath(self.disp,self.dsgn.nblocks+1,i_fill=iblock-1)   #display progress
        self.disp.text_either_key.draw()                        #text advance
        self.disp.win.flip()                                    #flip the window
        
        #save timing information
        self.dsgn.cd_time_blockstart = clock.getTime()

        #wait for the response to start the block
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_same+self.dsgn.key_diff+self.dsgn.key_escape,f=f)

        core.wait(1)

    def changeDetectionEncodingArray(self,iblock,itrial,clock,f=None):
        
        #draw memory array
        self.disp.fix.draw()    #draw fixation
        for i in range(self.dsgn.cd_setsize):
            self.disp.square.pos = (self.dsgn.cd_encarray_x[itrial][i],self.dsgn.cd_encarray_y[itrial][i])
            self.disp.square.fillColor = self.dsgn.cd_encarray_colorrgb[itrial][i]
            self.disp.square.draw()
        self.disp.win.flip()    #flip the window
        
        #save timing information
        self.dsgn.cd_time_encarray_onset[itrial] = clock.getTime()
        
        #wait for the duration that stimuli should be on the screen
        core.wait(self.dsgn.cd_t_encarray)

    def changeDetectionRetentionInterval(self,iblock,itrial,clock,f=None):
        
        #draw blank retention interval
        self.disp.fix.draw()
        self.disp.win.flip()
        
        #save timing information
        self.dsgn.cd_time_encarray_offset[itrial] = clock.getTime()
        
        #wait for the duration of the retention interval
        core.wait(self.dsgn.cd_t_retinterval)


    def changeDetectionMemoryProbe(self,iblock,itrial,clock,f=None):
        
        #display probe image
        self.disp.fix.draw()        #fixation dot
        self.disp.square.pos = (self.dsgn.cd_memprobe_x[itrial],self.dsgn.cd_memprobe_y[itrial])
        self.disp.square.fillColor = self.dsgn.cd_memprobe_probecolorrgb[itrial]
        self.disp.square.draw()
        self.disp.text_diff.draw()  #draw response for different trials
        self.disp.text_same.draw()  #draw response for same trials
        self.disp.win.flip()        #flip window
        
        #save timing information
        self.dsgn.cd_time_memprobe_onset[itrial] = clock.getTime()

        #wait for responses
        response = wait_for_response(self.disp.win,keylist=self.dsgn.key_same+self.dsgn.key_diff+self.dsgn.key_escape,f=f)
        
        #collect responses and RT
        self.dsgn.cd_time_memprobe_offset[itrial] = clock.getTime()
        self.dsgn.cd_rt[itrial]=self.dsgn.cd_time_memprobe_offset[itrial]-self.dsgn.cd_time_memprobe_onset[itrial]
        self.dsgn.cd_actual_response_string[itrial] = response
        if [self.dsgn.cd_actual_response_string[itrial]] == self.dsgn.key_same:
            self.dsgn.cd_actual_response[itrial] = 1
        else:
            self.dsgn.cd_actual_response[itrial] = 0
        
        self.dsgn.cd_acc[itrial] = (self.dsgn.cd_actual_response[itrial]==self.dsgn.cd_correct_response[itrial])*1 

        #remove memory probe after response
        self.disp.fix.draw()
        self.disp.win.flip()
        
        #save timing information
        self.dsgn.cd_time_memprobe_offset[itrial] = clock.getTime()

    def changeDetectionITI(self,iblock,itrial,clock,f=None):
        
        #draw blank iti
        self.disp.fix.draw()
        self.disp.win.flip()
        
        #log time
        self.dsgn.cd_time_iti_onset[itrial] = clock.getTime()
        
        #wait for the duration of the iti
        core.wait(self.dsgn.cd_t_iti[itrial])
        
        #log time
        self.dsgn.cd_time_iti_offset[itrial] = clock.getTime()

    def endOfExperiment(self,f=None):
        self.disp.text_end_expt.draw()
        self.disp.win.flip()
        keys = event.waitKeys(keyList=['space','escape']) 
        