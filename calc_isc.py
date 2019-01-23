import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, pearsonr, zscore


# Create simple simulated data with high intersubject correlation
def simulated_timeseries(n_subjects, n_TRs, n_voxels=1, noise=1):
    signal = np.random.randn(n_TRs, n_voxels // 100)
    data = [zscore(np.repeat(signal, 100, axis=1) +
                 np.random.randn(n_TRs, n_voxels) * noise,
                 axis=0)
          for subject in np.arange(n_subjects)]
    
    return data

def plot_isc(c,fout):
    """Given a 2D matrix c creates a 2D colormap from -1 to 1"""

    plt.figure()
    plt.imshow(c, vmin = -1, vmax = 1, cmap = 'PRGn')
    plt.xlabel('# Voxels')
    plt.ylabel('# Voxels')
    plt.title('Leave-one-out ISC')
    plt.savefig(fout,dpi = 300,fmt=fig_fmt)
    plt.close('all')
    
    return
    
    
def leave_one_out_isc(data,n_subjects,froot):
    
    data_average = np.average(data,axis = 0)
        
    for participant in range(2):
    
        # I exclude the participant from the total average
        corr_average = data_average - data[participant] / (1.0 * n_subjects)
        
        # Calculate the  correlation matrix
        c = np.corrcoef(data[participant],corr_average, rowvar = False)[n_voxels:,:n_voxels] # to select corr between part and av
        # And save it
        fout = froot + 'c_{}.txt'.format(participant)
        np.savetxt(fout,c)

        # Perform the Fisher z-transformation
        z = np.arctanh(c)
        # And save it
        fout = froot + 'z_{}.txt'.format(participant)
        np.savetxt(fout,z)

        # And save a plot as well
        fout = froot + 'c_{}.{}'.format(participant,fig_fmt)        
        plot_isc(c,fout)             


    return

##### MAIN CODE ####
if __name__ == '__main__':
    

    fig_fmt = 'png'    
    froot = 'M:/data/'
    ### IMPORT INPUT DATA ###
    
    
    # Set parameters for toy time series data
    n_subjects = 20
    n_TRs = 300
    n_voxels = 1000
    
    
    # Go fo simulated data
    data = np.array(simulated_timeseries(n_subjects, n_TRs, n_voxels=n_voxels))
    
    # Compute leave-one-out isc
    leave_one_out_isc(data,n_subjects,froot)