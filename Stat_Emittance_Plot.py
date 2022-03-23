import os

import numpy as np
import pandas as pd

# Note to the user, this program -MAY- throw an error when saving the files for each BPM on windows, there exist ...
# solutions online regarding chaning folder types, however i reccomend not bothing with this, simply delete any ...
# files beggining with .png in the folder specific for the BPMS, this then should have no issues and removes ...
# the need to worry about setting folder permissions to numeric (not simple). On linux this shpuld not be an ...
# issue - I point blank refuse to use a Mac to test due to chubby fingers not liking the command button. Regards \S

dir = "D:\Downloads\BPMS"  # File directory containing the BPM files
dats = ["x", "y"]# Flick between x or y at convinience - does need formatting as a text string.
shift_data = True # Centers the ellipse back on 0,0 for the purposes of plotting, depending on how G4BL is used ...
#Sometimes the data can shift away from 0,0 despite it being centered on 0,0 from the beamline perspective - be careful.
show_plots = True  # Boolean for displaying plots or not.
calc_emit = True

def emit(data,flag):
    if flag == "x":

        sxd = data["x"]
        sxprimed = data["xprime"]
        sxxprimed = data["xxprime"]
        x2 = np.power(sxd, 2)
        xprime2 = np.power(sxprimed, 2)
        avgx2 = np.average(x2)
        avgxprime2 = np.average(xprime2)

        avgxxprime = np.average(sxxprimed)
        avgxxprime2 = np.power(avgxxprime, 2)
        emitx = np.sqrt((avgx2 * avgxprime2) - avgxxprime2)
        return (emitx)
    elif flag== "y":

        syd = data["y"]
        syprimed = data["yprime"]
        syyprimed = data["yyprime"]
        y2 = np.power(syd, 2)
        yprime2 = np.power(syprimed, 2)
        avgy2 = np.average(y2)
        avgyprime2 = np.average(yprime2)
        # here we average then square
        avgyyprime = np.average(syyprimed)
        avgyyprime2 = np.power(avgyyprime, 2)
        emity = np.sqrt((avgy2 * avgyprime2) - avgyyprime2)
        return emity

for filename in os.listdir(dir):
    if filename.endswith(".txt"):
        File = os.path.join(dir, filename)  # Setting the file to this iteration
        filename = filename.strip(
            ".txt")  # Removing the .txt for later formatting - if your file isnt a .txt change to the extension here.
        pd.set_option('display.max_columns', None)  # set display columns pandas
        cols = ["x", "y", "z", "Px", "Py", "Pz", "t", "PDGid", "EventID", "TrackID", "ParentID", "Weight"]
        data = pd.read_csv(File, sep=" ", skiprows=3, header=None, names=cols)  # read data with above vols
        if shift_data == True: #Shifting the data to 0,0 if needed
            avgx = data["x"].mean()  # X can sometimes be shifted, global coords - just moving it back to roughly 0 for the purposes of plotting. Draw attention to this shift if required
            data["x"] = data["x"] - avgx
            avgy = data["y"].mean()
            data["y"] = data["y"] - avgy

        xprime = data["Px"] / data["Pz"]  # Calculating x' using for small angles (x1-x2/L) ~ (Px/Pz) assuming paraxial approximation
        data["xprime"] = xprime
        yprime = data["Py"] / data["Pz"]
        data["yprime"] = yprime
        data["xxprime"] = data["x"] * data["xprime"]
        data["yyprime"] = data["y"] * data["yprime"]

        import matplotlib.pyplot as plt
        for xory in dats:
            from matplotlib.patches import Ellipse
            x = data[xory]  # setting our "x" data to "y"
            y = data[xory + "prime"]  # and y data to y'
            xmax = max(abs(x)) * 1.012
            ymax = max(abs(y)) * 1.012
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            plt.xlim(-xmax, xmax)
            plt.ylim(-ymax, ymax)

            ax.spines['left'].set_position(
                'zero')  # ax setting to center the axis, comment any ax.spines to return to edge axis.
            ax.spines['bottom'].set_position('zero')

            # Eliminate upper and right axes
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')

            # Show ticks in the left and lower axes only
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            ax.set_xlabel(xory, fontsize=15)  # Setting the labels for x and y
            ax.set_ylabel(xory + "'", rotation=0, fontsize=15)
            ax.xaxis.set_label_coords(0.98, 0.45)  # the X label position
            ax.yaxis.set_label_coords(0.55,
                                      0.95)  # And the y label position (0,0 bottom left, 0.5,0.5 center, 1,1 top right)
            def eigso(cov):  # Function to return eigenvalues of covariant matrix
                vs, vc = np.linalg.eigh(cov)
                ord = vs.argsort()[::-1]
                return vs[ord], vc[:, ord]
            plt.scatter(x, y, marker="+")  # scatter before ellipse, makes it clearer
            stdevs = [1, 2, 3]  # Standard deviations
            devs = ["68%", "95%", "99.7%"]  # Each stdev amount
            colours = ["red", "black", "orange"]  # colours of ellipses
            for i in range(len(stdevs)):  # calculating and plotting each ellipse.
                ax = plt.subplot(111)
                cov = np.cov(x, y)
                vs, vc = eigso(cov)
                t = np.degrees(np.arctan2(*vc[:, 0][::-1]))
                w, h = 2 * stdevs[i] * np.sqrt(vs)
                if calc_emit == True:
                    e = w * h
                    str = devs[i] +(" ($\epsilon_%s = %1.2e$)" % (xory,e))
                else:
                    str = devs[i]
                ell = Ellipse(xy=(np.mean(x), np.mean(y)),
                              width=w, height=h,
                              angle=t, color=colours[i], label=str)
                ell.set_facecolor('none')
                ax.add_artist(ell)
            if calc_emit == True:
                erms = emit(data,xory)
                plt.legend(title='$\epsilon_{%s rms} = %.2e$' % (xory, erms))
            else:
                plt.legend()
            plt.title(filename)
            plt.savefig(dir + "\\" + filename + "_emittance_" + xory + ".png")
            if show_plots == True:
                plt.show()
