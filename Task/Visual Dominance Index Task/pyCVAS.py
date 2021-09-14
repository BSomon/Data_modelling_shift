#!/usr/bin/env python
# -*- coding: utf-8 -*-
from psychopy import visual, core, event
from pygame import mixer
#parallel
import numpy as np  # whole numpy lib is available, prepend 'np.'
from random import random, shuffle
import os  # handy system and path functions
import sys  # to get file system encoding
import csv
from datetime import datetime
from easygui import msgbox, multenterbox
# import pdb




### GUI ### 

fieldNames = ["Subject code","number"]
fieldValues = list(multenterbox(msg='Fill in values for the fields.', title='Enter', fields=(fieldNames)))
msgbox(msg=(fieldValues), title = "Is this ok?")
name, num = fieldValues[0], int(fieldValues[1])    #To be used later for target modality randomization
lineDict = {}       #initiate the csv file
lineDict['ID'] = name

#### VARIABLES ####

condRepetition = 9   # will be multipicated by 16 (4 AV + 4 VA + 4 AS + 4 VS) for the total number of trials 
rep_duration=0.8
circleDuration = 0.2
#port = parallel.ParallelPort(address=0x378)        #for port parallel triggers


### VISUAL STUFF ###

win = visual.Window(fullscr=True, size=[1920, 1080], color = [-1,-1,-1])
screenWidth = win.size[0]
instructionWidth = 0.1*screenWidth

crossFixation = visual.TextStim(win=win,text="+", font='times', pos = (0,0),color = [0.5,0.5,0.5], height=0.3)

taskInstruction = visual.TextStim(win=win,text="Welcome.\n Visual targets (white circles) \
on the left or the right of the screen and sound targets (bips) in the left or the right ear will be presented. \n\
You have to give their presentation side with 'a' (for left) or 'p' (for right) keys. \n Press space bar to continu." ,pos = (0,0), height=None, color = [0.5,0.5,0.5], wrapWidth=1.8)

auditoryInstruction = visual.TextStim(win=win,text="Sound target: press 'a' for left sounds and 'p' for right sounds while looking at the fixation cross. \n Press space bar to continu.", pos = (0,0), color = [0.5,0.5,0.5], height=None)

visualInstruction = visual.TextStim(win=win,text="Visual target: press 'a' for left circle and 'p' for right circle while looking at the fixation cross. \n Press space bar to continu." , pos = (0,0),height=None)
TrainVisuInstruction = visual.TextStim(win=win,text="Training: Visual target: press 'a' for left circle and 'p' for right circle while looking at the fixation cross. \n Press space bar to continu." , pos = (0,0),height=None)
TrainAuditInstruction = visual.TextStim(win=win,text="Training: Auditory target: press 'a' for left sound and 'p' for right sound while looking at the fixation cross. \n Press space bar to continu." , pos = (0,0),height=None)

# STIM #
VD = visual.Polygon(win=win, edges=90, size=(0.15, 0.25), ori=0, pos=(0.7, 0), lineWidth=1, lineColor=[1,1,1], fillColor=[1,1,1])
VG = visual.Polygon(win=win, edges=90, size=(0.15, 0.25), ori=0, pos=(-0.7, 0), lineWidth=1, lineColor=[1,1,1], fillColor=[1,1,1])


#### SOUND STUFF ###
mixer.pre_init(44100, -16, 2, 1024)  #permet d'eviter le lag du son
mixer.init()
AD = mixer.Sound("750D.wav")
AG = mixer.Sound("750G.wav")


### experimental design ### Latin square
targetSequence = [['AV','VA','AS','VS'],['VA','AS','VS', 'AV'],['AS','VS', 'AV','VA'],['VS', 'AV','VA','AS']][num%4]

#targetSequence = ['A','V'] if num%2 else ['V','A']
conditions = ['CG','CD','IG','ID']*condRepetition

blocSequence = []
for thisTask in targetSequence:
    shuffle(conditions)
    blocSequence.append([thisTask+cond for cond in conditions])
print(blocSequence)
fieldsname = ['ID','condition', 'marker', 'response', 'RT', 'accuracy']



with open('Perf_'+name+'_'+datetime.today().strftime('%Y%m%d')+'_CVAS.csv', 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldsname, delimiter=';', lineterminator = '\n')
    writer.writeheader()
    
    jitter = round(random(),2)  
    timer = core.Clock()
    
    #### START TASK ####
    win.flip()
    taskInstruction.draw()
    win.flip()
    event.waitKeys(maxWait=25, keyList='space')

            ## Training ###
    # win.flip()
    # event.clearEvents()
    # win.flip()
    # TrainVisuInstruction.draw()
    # win.flip()
    # event.waitKeys(maxWait=10, keyList='space')
    # crossFixation.draw()
    # event.waitKeys(maxWait=20, keyList='space')
    # VD.draw()
    # event.waitKeys(maxWait=1.5, keyList='p')
    # VG.draw()
    # event.waitKeys(maxWait=1.5, keyList='a')
    # VD.draw()
    # event.waitKeys(maxWait=1.5, keyList='p')
    # VG.draw()
    # event.waitKeys(maxWait=1.5, keyList='a')
    # win.flip()


    
    for thisBloc in blocSequence:    ####C'est là!!!!!
        
        if thisBloc[0][0:2] in ['AV', 'AS'] : auditoryInstruction.draw()
        if thisBloc[0][0:2] in ['VA', 'VS'] : visualInstruction.draw()
        win.flip()
        event.waitKeys(maxWait=10, keyList='space')
        #core.wait(1)
        crossFixation.draw()
        win.flip()
        core.wait(2)
        
        for thisTrial in thisBloc:
            lineDict['condition'] = thisTrial
            lineDict['marker'] = 'nan'
            jitter = round(random(),2)
            #port.setData(0) 
            
         
            
            ### Rules for audiovisual stimulus presentation ###            
            if thisTrial[0:2] == 'AV' and thisTrial[2] == 'I':  #NB: l'entree est incluse, la sortie est exclusive (donc [0:2] = pos1 et pos 2 de la liste)
                sound = AG if thisTrial[3] == 'G' else AD
                circle = VD if thisTrial[2] == 'G' else VG
            elif thisTrial[0:2] == 'AV' and thisTrial[2] == 'C':
                sound = AG if thisTrial[3] == 'G' else AD
                circle = VG if thisTrial[3] == 'G' else VD
            elif thisTrial[0:2] == 'VA' and thisTrial[2] == 'I':
                circle = VG if thisTrial[3] == 'G' else VD
                sound = AD if thisTrial[3] == 'G' else AG
            elif thisTrial[0:2] == 'VA' and thisTrial[2] == 'C':    
                circle = VG if thisTrial[3] == 'G' else VD
                sound = AG if thisTrial[3] == 'G' else AD 
            
            elif thisTrial[0:2] == 'AS':
                sound = AG if thisTrial[3] == 'G' else AD    
                
            elif thisTrial[0:2] == 'VS':   
                circle = VG if thisTrial[3] == 'G' else VD
                
            soundPlayed = False
            responseGiven = False
            
            ### next steps to prepare the csv file ###
            lineDict['response'] = 'nan'
            lineDict['RT'] = 'nan'
            lineDict['accuracy'] = 'miss'
            
            ### trial timing ###
            event.clearEvents()
            timer.reset()
            while 0 < timer.getTime() < circleDuration + rep_duration:
                
                crossFixation.draw()
                
                if 0 < timer.getTime() < circleDuration:
                   if not thisTrial[0:2] == 'AS': circle.draw() # Circle
                    # pdb.set_trace()
                if 0 < timer.getTime() and not soundPlayed:
                    if not thisTrial[0:2] == 'VS': sound.play()
                    soundPlayed = True
                
                ### Get answers ###
                if not responseGiven:
                    key = event.getKeys(timeStamped=timer)
                    
                    # port.setData(8)
                    if key:
                        if str(key[0][0]).lower() == 'a':
                            #pdb.set_trace()
                            lineDict['response'] = 'G'
                            lineDict['RT'] = round(key[0][1],3)
                            responseGiven = True
                        if str(key[0][0]).lower() == 'p':
                            lineDict['response'] = 'D'
                            lineDict['RT'] = round(key[0][1],3)
                            responseGiven = True
                            
                        if lineDict['response'] == thisTrial[3]:
                            lineDict['accuracy'] = 'hit'
                        elif lineDict['response'] != thisTrial[3]:
                            lineDict['accuracy'] = 'false'  
                        elif lineDict['response'] == 'nan':
                            lineDict['accuracy'] = 'miss'
                        
                win.flip()
            
            core.wait(jitter)
                
            writer.writerow(lineDict)
            
               
    win.close()
    core.quit()
    # make sure everything is closed down

