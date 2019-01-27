#!/usr/bin/python
import os,sys,random,subprocess,fileinput
import numpy as np

####################################################################################
# USAGE: python BasicRegressionTests.py 
# 
####################################################################################

#  Check if a command exists
def CmdExists(cmd):
    return subprocess.call("type " + cmd, shell=True, 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

#   This goes into a new problem directory and runs + compares to REFERENCE output
def RunProblem(DirectoryID,ReferenceFileName) :

    #---Get current directory
    CurrentDirectory = os.getcwd() 
    ProblemDirectory = CurrentDirectory + '/' + DirectoryID 
    print 'Current working directory:', str(ProblemDirectory) 
    #---Enter new problem directory
    os.chdir(ProblemDirectory)

    #---Get number of lines for REFERENCE file
    ReferenceFile       = np.loadtxt(ReferenceFileName,skiprows=1) 
    PowerReference      = ReferenceFile[:,1]
    TotalReferenceLines = len(PowerReference)

    print '# of lines in Reference file: ', str(TotalReferenceLines)    

    #---Run program 
    p = subprocess.Popen( ExePath +'/msr1d' )
    p.wait()
    
    #---Newly calculated file information
    CalculationFile       = np.loadtxt('power_amp_soln.txt',skiprows=1)
    PowerCalculation      = CalculationFile[:,1]
    TotalCalculationLines = len(PowerCalculation)
   
    print '# of lines in calculated file: ', str(TotalCalculationLines)
   
    #---Check the number of lines between files
    if TotalCalculationLines - TotalReferenceLines > 0:
        print 'Potentailly a problem - mismatched number of output lines '
    else:
        print 'Same number of lines in both Reference and calculated '

    #---Compare values between reference and test case
    Tolerance = 1E-3
    j = 0
    for j in range(0,TotalCalculationLines):
        RefValue  = PowerReference[j]
        CalcValue = PowerCalculation[j]
        Difference = abs(RefValue - CalcValue)
        q=0
        #---Count # of times the power differs greater than the tolerance
        if Difference > Tolerance:
            q=1+1
        
        j = j+1
    
    print 'Reference and Calculated differed by more than tolerance: ', Tolerance,' this many times: ',q

    print '****************************************************'
    print '*** Move on to next test case ***'

    #---Return to the starting location 
    os.chdir(CurrentDirectory) 


#****************************************************************

# Begin main part of program

#---Select the directories to be read in
ListDirs="list_of_dirs_for_comp.txt"
DirectoryNames = [line.rstrip('\n') for line in open(ListDirs)]
SummaryFileName='SummaryRegressionTests.txt'

#---Reference output to compare to 
ReferenceFileName = 'REFERENCE_SOLUTION.TXT'

#---Location of MSR1D executable
ExePath = '/Users/zanderm/Documents/Dissertation/MSR_1D/PrecursorSolve/src' 
DoesExecExist = CmdExists(ExePath +'/msr1d')

#---Test if exectutable is actually there
if DoesExecExist == True:
    print 'Using executable in: ', str(ExePath)
else:
    print 'MSR1D executable is not present in directory specified'
    sys.exit()

print '*** Entering regression sequence ***'

i=0
#---Go into each directory, run problem and compare to REFERENCE solution
while i < len(DirectoryNames):
    #---Run problem
    RunProblem(DirectoryNames[i],ReferenceFileName)
    i = i+1
