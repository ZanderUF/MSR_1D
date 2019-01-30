#!/usr/bin/python
import os,sys,random,subprocess,fileinput

####################################################################################
# USAGE:
# 
####################################################################################


def CreateDirectory(DirectoryID,BaseFileName, PreviousParm, ParmToVary,ExePath) :

    #---Make the directory we are going to put everything in 
    os.system('mkdir ' + DirectoryID )
    
    #---Copy supporting files into new director
    BetaFileName = 'beta_flow.txt'
    Dif3dFileName = 'dif3d_values.txt'
    PlottingDirectory = 'plotting'
    RunScriptName = 'run_msr.sh'
    #---Cp 
    os.system('cp '+ BaseFileName      + ' ' + DirectoryID )
    os.system('cp '+ BetaFileName      + ' ' + DirectoryID ) 
    os.system('cp '+ Dif3dFileName     + ' ' + DirectoryID ) 
    os.system('cp -R '+ PlottingDirectory + ' ' + DirectoryID ) 
    os.system('cp '+ RunScriptName     + ' ' + DirectoryID ) 
    
    #---Enter new directory
    CurrentDirectory = os.getcwd() 
    os.chdir(CurrentDirectory + '/' + DirectoryID)
    
    #---Open 'new' input file to modify
    
    # Read in the file
    with open(BaseFileName, 'r') as file :
      filedata = file.read()
    
    # Replace the target string
    filedata = filedata.replace(PreviousParm, ParmToVary)
    
    # Write the file out again
    with open(BaseFileName, 'w') as file:
      file.write(filedata)

    #---Run program 
    p = subprocess.Popen( ExePath +'/msr1d' )
    p.wait()

    #---Return to the starting location 
    os.chdir(CurrentDirectory) 
    
#----Setup the problem

#---Starting problem from the 'base' input
BaseFileName = 'input_t'
#---What we are naming each of the directories
DirectoryBaseName = 'MassFlow'
#---Path to MSR1D executable
ExePath = '/Users/zanderm/Documents/Dissertation/MSR_1D/PrecursorSolve/src/'

NumberCases = 2
ParameterInverval = 500

#---Select parameter to vary
StartingParmID = 'mflow='
#---This starting parm should match original input_t starting value
StartingParmValue = 10000 
StartingParm = StartingParmID + str(StartingParmValue)

for i in range(NumberCases) :
    #---New value
    Parm = StartingParmValue - (i+1)*ParameterInverval
    ParmToVary = StartingParmID + str(Parm)
    DirectoryID = str(i) + '_' + DirectoryBaseName + '_' + str(Parm) 

    CreateDirectory(DirectoryID,BaseFileName,StartingParm,ParmToVary,ExePath)

