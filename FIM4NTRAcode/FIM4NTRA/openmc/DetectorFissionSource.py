from .common_imports import np, matlab, plt, sns, csr_matrix, eigs


class DetectorFissionSource: 


    def __init__(self, path, detector_name):
        self.path = path
        self.detector_name = detector_name
        self.matlab_det_extraction() # Internal method



    # Data extraction method
    def matlab_det_extraction(self, normalize: bool = True): 
        eng = matlab.engine.start_matlab()
        

        eng.run(self.path,nargout=0)
        
        
        self.det_t = np.array(eng.workspace[f'DETxy_{self.detector_name}'])[:,10]
        

        self.dic_mesh_dimensions = {
                         'nx' : int(len(eng.workspace["DETxy_fission_sourceX"])),
                         'ny' : int(len(eng.workspace["DETxy_fission_sourceY"])),    
                         'xmin': eng.workspace['DETxy_fission_sourceX'][0][0],
                         'xmax': eng.workspace['DETxy_fission_sourceX'][-1][-2],
                         'ymin': eng.workspace['DETxy_fission_sourceY'][0][0],
                         'ymax': eng.workspace['DETxy_fission_sourceY'][-1][-2]      
        }


        if normalize:
            self.det_t = self.det_t/np.sum(self.det_t) 
        

        self.det_t = self.det_t.reshape(self.dic_mesh_dimensions['nx'],self.dic_mesh_dimensions['ny'])


        self.det_t_err = np.array(eng.workspace[f'DETxy_{self.detector_name}'])[:,11]
    
    

    # Return fmxt tallies
    def get_det_tallies(self):
        return self.det_t
    

    # Return fmtx statistical errors
    def get_det_errors(self):
        return self.det_t_err    
    


    # Plot Fission source
    def plot_fission_source(self, cmap: str = 'inferno'):
        plt.figure()
        

        plt.imshow(self.det_t, extent = [self.dic_mesh_dimensions['xmin'], self.dic_mesh_dimensions['xmax'], self.dic_mesh_dimensions['ymin'], self.dic_mesh_dimensions['ymax']], cmap = cmap)


        cbar = plt.colorbar()


        plt.title("Detector Fission Source")


    
