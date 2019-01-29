import numpy as np
import matplotlib.pyplot as plt
import configparser, os, re
import h5py
from scipy.stats import zscore
import nibabel as nib
#%%
# Create simple simulated data with high intersubject correlation
def simulated_timeseries(n_subjects, n_TRs, n_voxels=1, noise=1):
    signal = np.random.randn(n_TRs, n_voxels // 100)
    data = [zscore(np.repeat(signal, 100, axis=1) +
                 np.random.randn(n_TRs, n_voxels) * noise,
                 axis=0)
          for subject in np.arange(n_subjects)]
    
    return data

def plot_isc(c,participant,foutroot):
    """Given a 2D matrix c creates a 2D colormap from -1 to 1"""

    plt.figure(figsize = (5,4))
    plt.imshow(c, vmin = -1, vmax = 1, cmap = 'PRGn')
    plt.xlabel('# Voxels')
    plt.ylabel('# Voxels')
    plt.colorbar()
    plt.title('Leave-one-out ISC // Participant {}'.format(participant))
    fout = foutroot + 'c_{}.{}'.format(participant,fig_fmt)        
    plt.savefig(fout,dpi = 300,fmt=fig_fmt)
    plt.close('all')
    
    return
    
    
def leave_one_out_isc(data,n_subjects,froot,n_voxels):
    
    data_average = np.average(data,axis = 0)
    foutroot = froot + '/ISC/'
    if not os.path.exists(foutroot):
        os.makedirs(foutroot)
        
    for participant in range(2):
    
        # I exclude the participant from the total average
        corr_average = data_average - data[participant] / (1.0 * n_subjects)
        
        # Calculate the  correlation matrix
        c = np.corrcoef(data[participant],corr_average, rowvar = False)[n_voxels:,:n_voxels] # to select corr between part and av
        # And save it
        fout = foutroot + 'c_{}.txt'.format(participant)
        np.savetxt(fout,c)

        # Perform the Fisher z-transformation
        z = np.arctanh(c)
        # And save it
        fout = foutroot + 'z_{}.txt'.format(participant)
        np.savetxt(fout,z)

        # And save a plot as well
        plot_isc(c,participant,foutroot)             


    return

def fill_n(string,n):
    
    beg, temp = string.split('{')
    to_rep, final = temp.split('}')
    num_digits = len(to_rep)
    
    return beg + '{num:0{ndigits}}'.format(num=n,ndigits=num_digits) + final

def import_input_file(fin,sub_sample = False):
    
    input_ext = fin.split('.')[-1]
    
    if input_ext == 'mat':
        
        h5file = h5py.File(fin)
        data_pointer = h5file['data']

        data = data_pointer[:,:] # TODO might require a lot of memory
        
    elif input_ext == 'nii':
        
        img = nib.load(fin)            
        data = img.get_fdata()
        
        vx,vy,vz, n_TRs = data.shape
        n_voxels = np.prod([vx,vy,vz])
        
        data = data.reshape((n_voxels,n_TRs)).T # TODO double check that this transpose doesn't hurt
        
        # If you want to play with less data
        if sub_sample == True:
            
            n_sampled_voxels = 100
            r_indices = np.arange(n_sampled_voxels)
#                r_indices = np.random.permutation(n_voxels)[:n_sampled_voxels]
            v_indices = np.sort(r_indices)
        
            data = data[:,v_indices]
    
    return data
    
#%%
########## MAIN CODE ##########
if __name__ == '__main__':
    


    ### IMPORT INPUT DATA ###
    config = configparser.ConfigParser()
    config.read('settings.ini')
    fig_fmt = config['OUTPUT']['Figure format']
    fig_dpi = config['OUTPUT']['Figure dpi']
    froot = config['INPUT']['Folder with input data']
    N_subjects = int(config['INPUT']['N participants'])
    fld_example = config['INPUT']['Folder hierarchy']
    input_fname = config['INPUT']['File name']
    input_ext = input_fname.split('.')[-1]
    
    # Set parameters for toy time series data
    #%%
    subjects = np.arange(N_subjects) + 1
    
    file_list = []
    n_true = 0
    
    for n in subjects:
    
        fin = froot + fill_n(fld_example,n) + fill_n(input_fname,n)
        
        if os.path.isfile(fin):
            file_list.append(fin)
            n_true += 1
        else:
            print('Warning!!! No file found in %s' % fin)
            
    print('Ready to load data for {} participants'.format(n_true))


#%%
    file_list = [file_list[0]]
    # for the first settings

    for k,fin in enumerate(file_list):

        # If it's the first I need to initialize the average
        if k == 0:
            average = import_input_file(fin)
        else:
            average += import_input_file(fin)

    #%%
    
    # Go fo simulated data
    simdata = np.array(simulated_timeseries(2, n_TRs=300, n_voxels=1000))
    #%%    
    # Compute leave-one-out isc
    leave_one_out_isc(data,n_subjects,froot,n_sampled_voxels)
