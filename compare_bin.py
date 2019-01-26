#!/usr/bin/python
import os
import re
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
# Usage: The script was strictly developed for the 13 galaxies avaliable at   #
# the time of development, and thus is more of a script to be used for       #
# reference.
#                                                                             #
# Author: April Walker                                                        #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def main():
    fits_file = []

    #set to the directory containing all data you wish to plot
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

        bin_minimum = Ay[0]
        #Set this in case you're binning too many values

        if iterator == 4:
            bin_upper_limit = 16
            bin_lower_limit = 8

        if iterator == 5 or 6 or 7 or 10:
            bin_upper_limit = 14
            bin_lower_limit = 6

        if iterator == 9:
            bin_lower_limit = 4
            bin_upper_limit = 10

        if iterator == 8:
            bin_upper_limit = 25
            bin_lower_limit = 10

        else:
            bin_upper_limit = 15
            bin_lower_limit = 10

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
            # if(Ax_error[i] < 0.17):
            #     y_unbinned.append(Ay[i])
            #     x_unbinned.append(Ax[i])
            #     x_error_unbinned.append(Ax_error[i])
            if(j == len(x_error_binned) - 1):
                if((numpy.sqrt(x_error_binned[j])/(counter[j]) < 0.13) and (counter[j] >= bin_lower_limit)):
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

        if iterator == 0 or iterator == 5 or iterator == 10:
            if iterator == 0:
                ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#231193', fmt = 'o', marker = "s", zorder = 4)
            elif iterator == 5:
                ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#231193', fmt = 'o', marker = "o", zorder = 4)
            elif iterator == 8:
                ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#231193', fmt = 'o', marker = "^", zorder = 4)
            ax1.scatter(x, y, color = '#e6e6e6', zorder = 1, label = no, s=5)
            ax1.errorbar(x_unbinned, y_unbinned, xerr = x_error_unbinned, color = '#3219cd', fmt = 'o', zorder = 3, label = no, ms=3)

            # The following lines can be uncommented (and the above commented) to develop graphs for individual galaxies
            # ax1.scatter(x, y, color = '#aabaff', zorder = 2, label = no, s=5)
            # ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) of NGC4382')
            # legend_labels = ["NGC4382"]

            # ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) of M63')
            # legend_labels = ["M63"]

            # ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) of NGC3184')
            # legend_labels = ["NGC3184"]

            # The following ignores this galaxy
            # ax1.scatter(x, y, color = '#e6e6e6', zorder = 1, label = no, s=5)
            # ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#e6e6e6', fmt = '.', marker = "s", zorder = 1, label = no)
            # ax1.errorbar(x_unbinned, y_unbinned, xerr = x_error_unbinned, color = '#e6e6e6', fmt = '.', zorder = 1, label = no)

            dummy_variable = 1
        if iterator == 1 or iterator == 6 or iterator == 11:
            if iterator == 1:
                ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#096612', fmt = 'o', marker = "s", zorder = 4)
            elif iterator == 6:
                ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#096612', fmt = 'o', marker = "o", zorder = 4)
            elif iterator == 11:
                ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#096612', fmt = 'o', marker = "^", zorder = 4)
            ax1.scatter(x, y, color = '#e6e6e6', zorder = 1, label = no, s=5)
            ax1.errorbar(x_unbinned, y_unbinned, xerr = x_error_unbinned, color = '#0c8718', fmt = 'o', zorder = 3, label = no, ms=3)

            # The following lines can be uncommented (and the above commented) to develop graphs for individual galaxies
            # ax1.scatter(x, y, color = '#90e097', zorder = 2, label = no, s=5)
            # ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) of NGC0628')
            # legend_labels.append(["NGC0628"])

            # ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) of M94')
            # legend_labels = ["M94"]

            # ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) of NGC3198')
            # legend_labels = ["NGC3198"]

            # The following ignores this galaxy
            # ax1.scatter(x, y, color = '#e6e6e6', zorder = 1, label = no, s=5)
            # ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#e6e6e6', fmt = '.', marker = "s", zorder = 1, label = no)
            # ax1.errorbar(x_unbinned, y_unbinned, xerr = x_error_unbinned, color = '#e6e6e6', fmt = '.', zorder = 1, label = no)

            dummy_variable = 1
        if iterator == 2 or iterator ==  7 or iterator == 12:
            if iterator == 2:
                ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#6a0a7a', fmt = 'o', marker = "s", zorder = 4)
            elif iterator == 7:
                ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#6a0a7a', fmt = 'o', marker = "o", zorder = 4)
            elif iterator == 12:
                ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#6a0a7a', fmt = 'o', marker = "^", zorder = 4)
            ax1.scatter(x, y, color = '#e6e6e6', zorder = 1, label = no, s=5)
            ax1.errorbar(x_unbinned, y_unbinned, xerr = x_error_unbinned, color = '#880d9d', fmt = 'o', zorder = 3, label = no, ms=3)

            # The following lines can be uncommented (and the above commented) to develop graphs for individual galaxies
            # ax1.scatter(x, y, color = '#edddff', zorder = 2, label = no, s=5)
            # ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) of NGC2403')
            # legend_labels = ["NGC2403"]

            # ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) of M95')
            # legend_labels = ["M95"]

            # ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) of NGC4559')
            # legend_labels = ["NGC4559"]

            # The following lines ignore this galaxy
            # ax1.scatter(x, y, color = '#e6e6e6', zorder = 1, label = no, s=5)
            # ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#e6e6e6', fmt = '.', marker = "s", zorder = 1, label = no)
            # ax1.errorbar(x_unbinned, y_unbinned, xerr = x_error_unbinned, color = '#e6e6e6', fmt = '.', zorder = 1, label = no)

            dummy_variable = 1
        if iterator == 3 or iterator == 8:
            if iterator == 3:
                ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#b50636', fmt = 'o', marker = "s", zorder = 4)
            elif iterator == 8:
                ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#b50636', fmt = 'o', marker = "o", zorder = 4)
            ax1.scatter(x, y, color = '#e6e6e6', zorder = 1, label = no, s=5)
            ax1.errorbar(x_unbinned, y_unbinned, xerr = x_error_unbinned, color = '#ec0040', fmt = 'o', zorder = 3, label = no, ms=3)

            # The following lines can be uncommented (and the above commented) to develop graphs for individual galaxies
            # ax1.scatter(x, y, color = '#f2bcc5', zorder = 2, label = no, s=5)
            # ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) of NGC6946')
            # legend_labels = ["NGC6946"]

            # ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) of M100')
            # legend_labels = ["M100"]

            # The following lines ignore this galaxy
            # ax1.scatter(x, y, color = '#e6e6e6', zorder = 1, label = no, s=5)
            # ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#e6e6e6', fmt = '.', marker = "s", zorder = 1, label = no)
            # ax1.errorbar(x_unbinned, y_unbinned, xerr = x_error_unbinned, color = '#e6e6e6', fmt = '.', zorder = 1, label = no)

            dummy_variable = 1
        if iterator == 4 or iterator == 9:
            if iterator == 4:
                ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#bb9407', fmt = 'o', marker = "s", zorder = 4)
            elif iterator == 9:
                ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#bb9407', fmt = 'o', marker = "*", zorder = 4)
            ax1.scatter(x, y, color = '#e6e6e6', zorder = 1, label = no, s=5)
            ax1.errorbar(x_unbinned, y_unbinned, xerr = x_error_unbinned, color = '#f1bc00', fmt = 'o', zorder = 3, label = no, ms=3)

            # The following lines can be uncommented (and the above commented) to develop graphs for individual galaxies
            # ax1.scatter(x, y, color = '#f1d9ac', zorder = 2, label = no, s=5)
            # legend_labels = ["NGC7793"]
            # ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) of NGC7793')

            # legend_labels = ["NGC2841"]
            # ax1.set_title(r'Luminosity vs Photon Index ($\Gamma$) of NGC2841')

            # The following lines ignore this galaxy
            # ax1.scatter(x, y, color = '#e6e6e6', zorder = 1, label = no, s=5)
            # ax1.errorbar(x_binned, y_binned, xerr = x_error_binned, color = '#e6e6e6', fmt = '.', marker = "s", zorder = 1, label = no)
            # ax1.errorbar(x_unbinned, y_unbinned, xerr = x_error_unbinned, color = '#e6e6e6', fmt = '.', zorder = 1, label = no)

            dummy_variable = 1


    plt.yscale('log');
    ax1.set_xlabel("Photon Index")
    ax1.set_ylabel("Luminosity")
    ax1.set_title('Comparitive Chart')

    legend_labels = []
    for i in range(len(fits_dir)):
        legend_labels.append("".join(fits_dir[i].rsplit("-final-sample.fits")))

    ax1.legend(legend_labels)

    plt.grid(True)
    plt.draw()
    ax1.apply_aspect()
    fig.savefig('comparitive_auto.eps', dpi=fig.dpi)
    Image.open('comparitive_auto.eps').save('comparitive_auto.png','png')

    # The following lines can be uncommented (and the above commented) to develop graphs for individual galaxies

    # if (legend_labels == ["NGC4382", "NGC0628", "NGC2403", "NGC6946", "NGC7793", "M63", "M94", "M95", "M100", "NGC2841", "NGC3184", "NGC3198", "NGC4559"]):
    #     fig.savefig('comparitive_auto_update.eps', dpi=fig.dpi)
    #     Image.open('comparitive_auto_update.eps').save('comparitive_auto_update.png','png')
    #
    # elif (legend_labels == ["NGC4382"]):
    #     fig.savefig('NGC4382.eps', dpi=fig.dpi)
    #     Image.open('NGC4382.eps').save('NGC4382.png','png')
    #
    # elif (legend_labels == ["NGC0628"]):
    #     fig.savefig('NGC0628.eps', dpi=fig.dpi)
    #     Image.open('NGC0628.eps').save('NGC0628.png','png')
    #
    # elif (legend_labels == ["NGC2403"]):
    #     fig.savefig('NGC2403.eps', dpi=fig.dpi)
    #     Image.open('NGC2403.eps').save('NGC2403.png','png')
    #
    # elif (legend_labels == ["NGC6946"]):
    #     fig.savefig('NGC6946.eps', dpi=fig.dpi)
    #     Image.open('NGC6946.eps').save('NGC6946.png','png')
    #
    # elif (legend_labels == ["NGC7793"]):
    #     fig.savefig('NGC7793.eps', dpi=fig.dpi)
    #     Image.open('NGC7793.eps').save('NGC7793.png','png')
    #
    # if (legend_labels == ["M63"]):
    #     fig.savefig('M63.eps', dpi=fig.dpi)
    #     Image.open('M63.eps').save('M63.png','png')
    #
    # elif (legend_labels == ["M94"]):
    #     fig.savefig('M94.eps', dpi=fig.dpi)
    #     Image.open('M94.eps').save('M94.png','png')
    #
    # elif (legend_labels == ["M95"]):
    #     fig.savefig('M95.eps', dpi=fig.dpi)
    #     Image.open('M95.eps').save('M95.png','png')
    #
    # elif (legend_labels == ["M100"]):
    #     fig.savefig('M100.eps', dpi=fig.dpi)
    #     Image.open('M100.eps').save('M100.png','png')
    #
    # elif (legend_labels == ["NGC2841"]):
    #     fig.savefig('NGC2841.eps', dpi=fig.dpi)
    #     Image.open('NGC2841.eps').save('NGC2841.png','png')
    #
    # elif (legend_labels == ["NGC3184"]):
    #     fig.savefig('NGC3184.eps', dpi=fig.dpi)
    #     Image.open('NGC3184.eps').save('NGC3184.png','png')
    #
    # elif (legend_labels == ["NGC3198"]):
    #     fig.savefig('NGC3198.eps', dpi=fig.dpi)
    #     Image.open('NGC3198.eps').save('NGC3198.png','png')
    #
    # elif (legend_labels == ["NGC4559"]):
    #     fig.savefig('NGC4559.eps', dpi=fig.dpi)
    #     Image.open('NGC4559.eps').save('NGC4559.png','png')

    plt.show()
    return 0
main()
