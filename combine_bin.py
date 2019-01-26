#!/usr/bin/python
import os
import math
import matplotlib.pyplot as plt
import numpy
from numpy import sin, pi, arange
import astropy.io
from astropy.io import fits
from PIL import Image

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script contrasts photon index and luminosity for various observations. #
# In order to reduce error, the data is binned together                       #
#                                                                             #
# Usage: The script can automatically handle any number of galaxies of which  #
# the fits files should be avaliable in one folder (named "Final-fits" here). #
# Bin sizes and error should be set as needed. All the changes that the user  #
# should make are identifyable with the tag "NOTE". Lastly, change the first  #
# line of this script to the appropriate path to your python.
#                                                                             #
# Author: April Walker                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#This function does the binning and graph making
def BinData():
    gamma_avg = 0
    fits_file = []

    #NOTE: set to the directory containing all data you wish to plot
    fits_dir = os.listdir(os.environ["RESEARCH_DATA_PATH"] + "/Final-fits/")

    #specify the path each data file
    for i in range(len(fits_dir)):
        fits_file.append(os.environ["RESEARCH_DATA_PATH"] + "/Final-fits/" + fits_dir[i])

    fig = plt.figure(figsize=(20, 10))
    ax1 = fig.add_subplot(111)

    x = []
    x_error = []
    x_error_mean = []
    y = []

    for iterator, value in enumerate(fits_file):
        hdu_list = fits.open(value)

        #our data is now put into columns
        data = hdu_list[1].data
        cols = hdu_list[1].columns
        hdu_list.close()

        #take the values you need. If you need a new field, use print(cols) and use
        gamma = data.field('gamma')
        gamma_error = data.field('gamma_err')
        lh8 = data.field('lh8')

        #gamma of value 1.7 is an error as is NaN, so we only append indexes wihtout
        #those values
        gamma_indices = numpy.where(gamma != 1.7)

        for i, value in enumerate(gamma_indices[0]):
            if (gamma[value] == gamma[value] and lh8[value] == lh8[value] and gamma_error[value] > 0):
                x.append(gamma[value])
                x_error.append(gamma_error[value])
                y.append(lh8[value])

    #this guy holds our data set sorted by luminosity
    data_set = []
    for i in range(len(x)):
        temp = [y[i],x[i],x_error[i]]
        data_set.append(temp)

    #sort it
    data_set.sort(key=lambda x: x[0])

    Ay = []
    Ax = []
    Ax_error = []

    for i in range(len(x)):
        Ay.append(data_set[i][0])
        Ax.append(data_set[i][1])
        Ax_error.append(data_set[i][2])

    gamma_avg = sum(Ax)/float(len(Ax))
    # print(gamma_avg)


    bin_minimum = Ay[0]

    #NOTE: set these as needed to define bin limits
    bin_upper_limit = 50
    bin_lower_limit = 15

    y_binned = []
    x_binned = []
    x_error_binned = []

    y_unbinned = []
    x_unbinned = []
    x_error_unbinned = []

    counter = []

    # THE BELOW CODE AUTOBINS NORMALLY #
    j = 0
    for i in range(len(Ay)):
        if(j == len(x_error_binned) - 1):
            #NOTE: set the first portion of this if statement as needed (upon developement set to 0.1)
            if((numpy.sqrt(x_error_binned[j])/(counter[j]) < 0.1) and (counter[j] >= bin_lower_limit)):
                j += 1
            elif(counter[j] >= bin_upper_limit):
                j += 1
            else:
                counter[j] += 1
                y_binned[j] += Ay[i]
                x_binned[j] += Ax[i]
                x_error_binned[j] += Ax_error[i]**2
        else:
            y_binned.append(0)
            x_binned.append(0)
            x_error_binned.append(0)
            counter.append(0)

            counter[j] += 1
            y_binned[j] += Ay[i]
            x_binned[j] += Ax[i]
            x_error_binned[j] += Ax_error[i]**2

    #calculates the mean error as sqrt(sum(errors^2))/sqrt(n)
    for j in range(len(y_binned)):
        if value == value:
            y_binned[j] = y_binned[j]/counter[j]
            x_binned[j] = x_binned[j]/counter[j]
            x_error_binned[j] = numpy.sqrt(x_error_binned[j])/(counter[j])

    no = '_nolegend_'

    #all following code is for plotting to our figure
    ax1.scatter(x, y, color = '#aabaff', zorder = 2,  s=5)
    ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#231193', fmt = 'o', marker = "s", zorder = 4)
    ax1.errorbar(x_unbinned, y_unbinned, xerr = x_error_unbinned, color = '#3219cd', fmt = 'o', zorder = 3, label = no, ms=3)

    #test script
    # ax1.scatter(x, y, color = '#e6e6e6', zorder = 1, label = no, s=5)
    # ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#e6e6e6', fmt = '.', marker = "s", zorder = 1, label = no)
    # ax1.errorbar(x_unbinned, y_unbinned, xerr = x_error_unbinned, color = '#e6e6e6', fmt = '.', zorder = 1, label = no)

    plt.yscale('log');
    ax1.set_xlabel("Photon Index")
    ax1.set_ylabel("Luminosity")
    ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) from Observations in 13 Galaxies')
    legend_labels = ["Single Data Point", "Binned Data"]

    ax1.legend(legend_labels)

    plt.grid(True)
    plt.draw()
    ax1.apply_aspect()

    fig.savefig('spectral_analysis.eps', dpi=fig.dpi)
    Image.open('spectral_analysis.eps').save('spectral_analysis.png','png')
    plt.show()

    return gamma_avg

def main():
    gamma_avg = BinData()
    return 0
main()
