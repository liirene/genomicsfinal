def pathway_parse(path_to_script, args):
    import subprocess  # for running r through later
    command = 'C:/Program Files/R/R-3.3.2/bin/Rscript.exe'
    path_to_script = 'pathwayscript.R'
    args = ['clean_data.csv']
    cmd = [command, path_to_script] + args

    x = subprocess.check_output(cmd, universal_newlines=True)

###### PLOT MAKING
 # Read plot name from the pathway?


# Do this for each group of pathways that you're interested in (1 figure for
# each group of pathways)
#pathwayname = 'pathhere'
#title = str(pathwayname)

# get each protein name into the tuple form


#array = np.ndarray()
#f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3,
                                                                      #squeeze=False,
                                                                      #sharey='row')

#plt.suptitle(title, fontsize=12)

# include however many you want w/in the tuple

# Include subplot (and making into multiple plots) here

# Add pathway title name at top of the entire figure

#path = 'path'
# name this as the name of the pathway itself
#plt.savefig('path' + title + '.png', bbox_inches='tight')
# Make multiple graphs at once?
# Make 8 at once to be default

#plt.close('all')  # closes all open figures
