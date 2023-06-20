#!/usr/bin/env python3 
import matplotlib.pyplot as plt
from matplotlib.figure import figaspect
import os 
import pandas as pd
import mplcursors
import pyperclip
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "Helvetica"
# })

# Script created by Rachel Sun 
# Last modified: 2021-08-10
# This script plots and extracts bandgap data from the dispersion relation of a Bloch Wave Analysis simulation.

######## DEFINE FUNCTIONS ########
def findBandGaps(freqs): 
    # function that finds bandgaps from a list of frequencies in a dispersion relation 
    # remove duplicate frequencies
    freqs = list(set(freqs))
    # sort frequencies from small to large
    freqs.sort() 

    # find bandgaps based on a set minimum width 
    BGW = 0.163 # bandgap width in MHz
    bandgaps = [] # list to store bandgap information
    for freqInd, freq in enumerate(freqs[0:-1]):
        if abs(freqs[freqInd] - freqs[freqInd+1]) > BGW:
            bandgapWidth = abs(freqs[freqInd] - freqs[freqInd+1])
            bandgapCenter = (freqs[freqInd] + freqs[freqInd+1])/2 
            bandgaps.append([freqs[freqInd+1], freqs[freqInd], bandgapWidth, bandgapCenter]) # append bandgap information to list
    return bandgaps

def extractlist(myList): 
    # converts a list of lists into a list of values
    newlist = []
    for item in myList: 
        value = item[0]*1e-6
        newlist.append(value)
    return newlist


######## MAIN SCRIPT ########
directory = 'TBD' # set directory to folder containing dispersion relation data 
fileInd = 0 # index for each file
for filename in os.listdir(directory): # loop through each file in directory
    fileID = filename[0:-5] # extract file ID from filename
    if filename.endswith(".xlsx"): # only process excel files
        file = os.path.join(directory, filename) # create path to file
        wavenum = pd.read_excel(file, sheet_name='Sheet1', usecols=[0]) # extract wavenumbers
        realfreqs = pd.read_excel(file, sheet_name='Sheet1', usecols=[2]) # extract real-valued frequencies
        usablefreqs = extractlist(realfreqs.values.tolist()) # convert frequencies to list

        # find bandgaps
        extractedBGs = findBandGaps(usablefreqs)

        # plot dispersion relation
        # plt.figure(fileInd)
        plt.rcParams["figure.figsize"] = [5.50, 7.50] # set figure size
        plt.rcParams["figure.autolayout"] = True # fit plot to window
        fig, ax = plt.subplots(figsize=figaspect(4/3))
        plt.scatter(wavenum, realfreqs*1e-6, s=10) 

        # plot bandgaps
        print('Bandgaps: ' + fileID) 
        BGInfo = [] # list to store bandgap width and center
        for BG in extractedBGs: 
            print(BG)
            plt.axhspan(BG[0], BG[1], alpha=0.5, color='red') # plot bandgap as shaded rectangle
            BGInfo.append(BG[2]) # append bandgap width to list
            BGInfo.append(BG[3]) # append bandgap center to list
        print(BGInfo)

        # plot formatting
        csfont = {'fontname':'Helvetica', 'fontsize': 20} # set font
        plt.title(fileID.replace('p','.'), **csfont) # set title
        plt.yticks(**csfont) # set y-axis ticks
        plt.xlabel('Reduced Wavevector', **csfont) # set x-axis label
        plt.ylabel('Frequency (MHz)', **csfont) # set y-axis label
        fig.canvas.draw() # draw figure
        labels = ['Gamma', 'X', 'M', 'Gamma','R','M','X','R','Gamma'] # name x-axis tick labels
        ax.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8]) # set x-axis ticks
        ax.set_xticklabels(labels, **csfont) # set x-axis tick labels


        # uncomment for interactive plot values upon hover 
        cursor = mplcursors.cursor(hover=True)
        plt.show(block=True)
        pyperclip.copy(str(BGInfo)) # copy bandgap width and center to clipboard
        ########

        # uncomment for fast processing 
        # plt.show(block=False)
        # # save plot for viewing later 
        # plt.savefig(fileID + '.png')
        ########
    fileInd += 1 
