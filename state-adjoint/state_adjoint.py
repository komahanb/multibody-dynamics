# Import numpy
import numpy as np

# Import matplotlib
import matplotlib.pylab as plt

# Configure 
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'

# Optionally set font to Computer Modern to avoid common missing font errors
params = {
  'axes.labelsize': 24,
  'legend.fontsize': 14,
  'xtick.labelsize': 24,
  'ytick.labelsize': 24,
  'text.usetex': True}
plt.rcParams.update(params)

# Latex math
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath}']
plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = 'courier'
plt.rcParams['font.size'] = 18
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['lines.color'] = 'r'

# Make sure everything is within the frame
plt.rcParams.update({'figure.autolayout': True})

# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)

def readInputFile(filename):
    """
    Method to read the input file and return the data as string.
    """
    try:
        inpFile = open(filename, "r")
        content = list(inpFile.readlines())
        inpFile.close()
    except:
        raise IOError
    return content

def parse_data(fil):
    data = readInputFile(fil)
    t = []
    state1 = []
    state2 = []
    state3 = []
    for line in data:
        cols = line.split()
        t.append(float(cols[0]))
        state1.append(float(cols[1]))
        state2.append(float(cols[2]))
        state3.append(float(cols[2]))
    return t, state1, state2, state3

def make_plot_abm(time_steps, error, name):
    plt.figure()
    # order, error, time_steps
    #plt.axis([0, 12, 1.0e-16, 1.0e0])
    plt.loglog(time_steps, error[0,0,:], '--v', label='absolute error -- order 1', ms=12.0, color=tableau20[0])
    plt.loglog(time_steps, error[0,1,:], '-v' , label='relative error -- order 1', ms=12.0, color=tableau20[0])
    
    plt.loglog(time_steps, error[1,0,:], '--v', label='absolute error -- order 2', ms=12.0, color=tableau20[2])
    plt.loglog(time_steps, error[1,1,:], '-v' , label='relative error -- order 2', ms=12.0, color=tableau20[2])

    plt.loglog(time_steps, error[2,0,:], '--v', label='absolute error -- order 3', ms=12.0, color=tableau20[4])
    plt.loglog(time_steps, error[2,1,:], '-v' , label='relative error -- order 3', ms=12.0, color=tableau20[4])
        
    plt.xlabel('Number of Time Steps')
    plt.ylabel('Error')
    
    plt.legend(loc=2)
    plt.savefig(name)

def make_plot_bdf(time_steps, error, name):
    plt.figure()
    # order, error, time_steps
    #plt.axis([0, 12, 1.0e-16, 1.0e0])
    plt.loglog(time_steps, error[0,0,:], '--D', label='absolute error -- order 1', ms=12.0, color=tableau20[0])
    plt.loglog(time_steps, error[0,1,:], '-D' , label='relative error -- order 1', ms=12.0, color=tableau20[0])
    
    plt.loglog(time_steps, error[1,0,:], '--D', label='absolute error -- order 2', ms=12.0, color=tableau20[2])
    plt.loglog(time_steps, error[1,1,:], '-D' , label='relative error -- order 2', ms=12.0, color=tableau20[2])

    plt.loglog(time_steps, error[2,0,:], '--D', label='absolute error -- order 3', ms=12.0, color=tableau20[4])
    plt.loglog(time_steps, error[2,1,:], '-D' , label='relative error -- order 3', ms=12.0, color=tableau20[4])
    
    plt.xlabel('Number of Time Steps')
    plt.ylabel('Error')
    plt.legend(loc=2)
    plt.savefig(name)


def make_plot_nbg(time_steps, error, name):
    plt.figure()
    # order, error, time_steps
    #plt.axis([0, 12, 1.0e-16, 1.0e0])
    plt.loglog(time_steps, error[0,0,:], '--*', label='absolute error -- order 2', ms=12.0, color=tableau20[2])
    plt.loglog(time_steps, error[0,1,:], '-*' , label='relative error -- order 2', ms=12.0, color=tableau20[2])
    
    plt.loglog(time_steps, error[1,0,:], '--*', label='absolute error -- order 2', ms=12.0, color=tableau20[2])
    plt.loglog(time_steps, error[1,1,:], '-*' , label='relative error -- order 2', ms=12.0, color=tableau20[2])

    plt.loglog(time_steps, error[2,0,:], '--*', label='absolute error -- order 2', ms=12.0, color=tableau20[2])
    plt.loglog(time_steps, error[2,1,:], '-*' , label='relative error -- order 2', ms=12.0, color=tableau20[2])

    plt.loglog(time_steps, error[3,0,:], '--*', label='absolute error -- order 3', ms=12.0, color=tableau20[4])
    plt.loglog(time_steps, error[3,1,:], '-*' , label='relative error -- order 3', ms=12.0, color=tableau20[4])
    
    plt.xlabel('Number of Time Steps')
    plt.ylabel('Error')
    plt.legend(loc=2)
    plt.savefig(name)

def make_plot_dirk(time_steps, error, name):
    plt.figure()
    # order, error, time_steps
    #plt.axis([0, 12, 1.0e-16, 1.0e0])
    plt.loglog(time_steps, error[0,0,:], '--o', label='absolute error -- order 2', ms=12.0, color=tableau20[2])
    plt.loglog(time_steps, error[0,1,:], '-o' , label='relative error -- order 2', ms=12.0, color=tableau20[2])
    
    plt.loglog(time_steps, error[1,0,:], '--o', label='absolute error -- order 3', ms=12.0, color=tableau20[4])
    plt.loglog(time_steps, error[1,1,:], '-o' , label='relative error -- order 3', ms=12.0, color=tableau20[4])

    plt.loglog(time_steps, error[2,0,:], '--o', label='absolute error -- order 4', ms=12.0, color=tableau20[6])
    plt.loglog(time_steps, error[2,1,:], '-o' , label='relative error -- order 4', ms=12.0, color=tableau20[6])
    
    plt.xlabel('Number of Time Steps')
    plt.ylabel('Error')
    plt.legend(loc=2)
    plt.savefig(name)

def plotall_order2(time_steps, error, name):
    plt.figure()
    # method, error, time_steps        
    #plt.axis([0, 12, 1.0e-16, 1.0e0])
    plt.loglog(time_steps, error[0,0,:], '--D', label='absolute error -- BDF', ms=12.0, color=tableau20[2])
    plt.loglog(time_steps, error[0,1,:], '-D' , label='relative error -- BDF', ms=12.0, color=tableau20[2])
    
    plt.loglog(time_steps, error[1,0,:], '--v', label='absolute error -- ABM', ms=12.0, color=tableau20[2])
    plt.loglog(time_steps, error[1,1,:], '-v' , label='relative error -- ABM', ms=12.0, color=tableau20[2])

    plt.loglog(time_steps, error[2,0,:], '--*', label='absolute error -- NBG', ms=12.0, color=tableau20[2])
    plt.loglog(time_steps, error[2,1,:], '-*' , label='relative error -- NBG', ms=12.0, color=tableau20[2])

    plt.loglog(time_steps, error[3,0,:], '--o', label='absolute error -- DIRK', ms=12.0, color=tableau20[2])
    plt.loglog(time_steps, error[3,1,:], '-o' , label='relative error -- DIRK', ms=12.0, color=tableau20[2])
    
    plt.xlabel('Number of Time Steps')
    plt.ylabel('Error')
    
    plt.legend(loc=2)
    plt.savefig(name)

def plotall_order3(time_steps, error, name):
    plt.figure()
    # method, error, time_steps        
    #plt.axis([0, 12, 1.0e-16, 1.0e0])
    plt.loglog(time_steps, error[0,0,:], '--D', label='absolute error -- BDF', ms=12.0, color=tableau20[4])
    plt.loglog(time_steps, error[0,1,:], '-D' , label='relative error -- BDF', ms=12.0, color=tableau20[4])
    
    plt.loglog(time_steps, error[1,0,:], '--v', label='absolute error -- ABM', ms=12.0, color=tableau20[4])
    plt.loglog(time_steps, error[1,1,:], '-v' , label='relative error -- ABM', ms=12.0, color=tableau20[4])

    plt.loglog(time_steps, error[2,0,:], '--*', label='absolute error -- NBG', ms=12.0, color=tableau20[4])
    plt.loglog(time_steps, error[2,1,:], '-*' , label='relative error -- NBG', ms=12.0, color=tableau20[4])

    plt.loglog(time_steps, error[3,0,:], '--o', label='absolute error -- DIRK', ms=12.0, color=tableau20[4])
    plt.loglog(time_steps, error[3,1,:], '-o' , label='relative error -- DIRK', ms=12.0, color=tableau20[4])
    
    plt.xlabel('Number of Time Steps')
    plt.ylabel('Error')
    
    plt.legend(loc=2)
    plt.savefig(name)

[t, q, qdot, qddot] = parse_data("nbg_adjoint.dat")
num_steps = len(t)
data = np.asarray([t,q,qdot,qddot])
print data
