import matplotlib
#This allows us to create png files without a display object (for remote severs without displays)
matplotlib.use('Agg') 
import sys, getopt
from pylab import *
import matplotlib.pyplot  as pyplot
import os.path
import shlex


def main(argv):

    y_data = ""
    x_data = ""
    x_limits = [0,0]
    y_limits = [0,0]
    error_data = ""
    x_label = "x data"
    y_label = "y data"
    legend = ""
    save_plot_file = ""
    title = ""
    error_data = ""

    #print(argv)

    try:

        opts, args = getopt.getopt(argv, "h", ["xdata=", "ydata=", "errordata=", "xlabel=","ylabel=", "legend=", "saveplot=", "title=", "ylimits=", "xlimits="])
        

    except getopt.GetoptError:

        print argv[0] + ' --xdata="1 2 4 ..." --ydata="3 6 7 ..."'
        sys.exit(2)

    for opt, arg in opts:

        

        if opt == '-h':

            print 'test.py -i <inputfile> -o <outputfile>'
            sys.exit()
            
        
        elif opt in "--xdata":

            x_data = arg

        elif opt in "--ydata":

            y_data = arg

        elif opt in "--errordata":

            error_data = arg

        elif opt in "--ylabel":

            y_label = arg

        elif opt in "--xlabel":

            x_label = arg

        elif opt in "--legend":

            legend = arg

        elif opt in "--saveplot":

            save_plot_file = arg

        elif opt in "--title":

            title  = arg
        
        elif opt in "--ylimits":

            y_limits  =  [ float(x) for x in arg.split() ]
            
        elif opt in "--xlimits":

            x_limits  =  [ float(x) for x in arg.split() ]

    if len(y_data.split()) != len(x_data.split()) :

        print '--x_data and --y_data don\'t have the same number of parameters ' + str(len(x_data)) + " " + str(len(y_data))
        sys.exit();

    if( error_data == "" ):    

        print("plotting")
        plot(x_data,y_data,x_label,y_label, legend, title, save_plot_file, x_limits, y_limits)
        
    else:
        
        error_plot(x_data,y_data,error_data,x_label,y_label, legend, title, save_plot_file, x_limits, y_limits)

def error_plot(x_data,y_data,error_data,x_label,y_label, legend, title, save_plot_file, xlimits, ylimits):

    colors = ['blue', 'red', 'green', 'purple', 'orange', 'black', '#ffff00', '#999999', '#ff00ff', '#00ffff', '#555555', '#ff9900']

    standard_font_size = 32
    standard_line_width = 10
    standard_figure_size_width = 25
    standard_figure_size_height = 15
    smaller_text = 24


    if legend != "" :

        standard_font_size = 70
        standard_line_width = 15
        standard_figure_size_width = 50
        standard_figure_size_height = 30
        smaller_text = 55

    font = {'size'   : standard_font_size}
    matplotlib.rc('font', **font)

    fig = pyplot.figure(figsize=(standard_figure_size_width, standard_figure_size_height))
    ax1 = fig.add_subplot(1,1,1)
    ax1.grid(lw=standard_line_width)
    plt.xlabel(x_label,fontsize=standard_font_size)
    plt.ylabel(y_label,fontsize=standard_font_size)
    legend_name = ""
    legend_names = []

    if(legend != ""):

        legend_names = legend.split()
        legend_name = legend[0]

    x_data_sets = x_data.split("#")
    y_data_sets = y_data.split("#")
    line_data_sets = error_data.split("#")

    for i_index in range( len(y_data_sets)):

        y_coordinates = [ float(x) for x in y_data_sets[i_index].split() ]
        error_data_sets = [float(x) for x in line_data_sets[i_index].split() ]
        

        if( i_index > len( x_data_sets) ) :

            x_coordinates = [ float(x) for x in x_data_sets[0].split() ]

        else:

            x_coordinates = [ float(x) for x in x_data_sets[i_index].split() ]

        if len(legend_names) > i_index:

            legend_name = legend_names[i_index]

        print error_data_sets
        print y_coordinates
        print x_coordinates
        ax1.errorbar( x_coordinates, y_coordinates, yerr=error_data_sets, fmt="-o", color=colors[ i_index % len(colors)  ], lw=standard_line_width, label=legend_name)

    if legend != "" :

        ax1.legend( loc=1,ncol=2,fontsize=smaller_text)

    if title != "" :

        plt.suptitle(title, fontsize=standard_font_size)
    
    if xlimits[0] != xlimits[1] :
    
        ax1.set_xlim(xlimits)
    
    if ylimits[0] != ylimits[1] :

        ax1.set_ylim(ylimits)

    if save_plot_file == "":

        plt.show()

    else:

        plt.savefig(save_plot_file)

def plot(x_data, y_data,x_label, y_label, legend, title, save_plot_file, xlimits, ylimits):

    colors = ['blue', 'red', 'green', 'purple', 'orange', 'black', '#ffff00', '#999999', '#ff00ff', '#00ffff', '#555555', '#ff9900']

    standard_font_size = 32
    standard_line_width = 10
    standard_figure_size_width = 25
    standard_figure_size_height = 15
    smaller_text = 24


    if legend != "" :

        standard_font_size = 70
        standard_line_width = 15
        standard_figure_size_width = 50
        standard_figure_size_height = 30
        smaller_text = 55

    font = {'size'   : standard_font_size}
    matplotlib.rc('font', **font)

    fig = pyplot.figure(figsize=(standard_figure_size_width, standard_figure_size_height))
    ax1 = fig.add_subplot(1,1,1)
    ax1.grid(lw=standard_line_width)
    plt.xlabel(x_label,fontsize=standard_font_size)
    plt.ylabel(y_label,fontsize=standard_font_size)
    legend_name = ""
    legend_names = []

    if(legend != ""):

        legend_names = legend.split()
        legend_name = legend[0]

    x_data_sets = x_data.split("#")
    y_data_sets = y_data.split("#")

    for i_index in range( len(y_data_sets)):

        y_coordinates = y_data_sets[i_index].split()

        if( i_index > len( x_data_sets) ) :

            x_coordinates = x_data_sets[0].split()

        else:

            x_coordinates = x_data_sets[i_index].split()

        if len(legend_names) > i_index:

            legend_name = legend_names[i_index]


        ax1.plot( x_coordinates, y_coordinates, color=colors[ i_index % len(colors)  ], lw=standard_line_width, label=legend_name)

    if legend != "" :

        ax1.legend( loc=1,ncol=2,fontsize=smaller_text)

    if title != "" :

        plt.suptitle(title, fontsize=standard_font_size)
        
    if xlimits[0] != xlimits[1] :
    
        ax1.set_xlim(xlimits)
    
    if ylimits[0] != ylimits[1] :

        ax1.set_ylim(ylimits)

    if save_plot_file == "":

        plt.show()

    else:

        plt.savefig(save_plot_file)

def runFileArguments(file_name):
    
    
    f = open(file_name,"r") 
    
    for line in f:

        if line.strip() != "" :
                        
            arguments = shlex.split(line)
            main(arguments[2:])
            pyplot.close("all")
        
 

if __name__ == "__main__":


    input_file = ""

    print(sys.argv[2])

    if sys.argv[1] == "inputfile" :
    
        
        print("Input file")
        runFileArguments(sys.argv[2])

    
    else:    
        
        print("No input file");
        main(sys.argv[1:])
        
        