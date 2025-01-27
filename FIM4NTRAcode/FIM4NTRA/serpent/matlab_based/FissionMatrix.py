from .common_imports import np, matlab, plt, sns, csr_matrix, eigs


class FissionMatrix: 


    def __init__(self,path):
        self.path = path
        self.matlab_fmxt_extraction() # Internal method


    # Data extraction method
    def matlab_fmxt_extraction(self): 
        eng = matlab.engine.start_matlab()
        

        eng.run(self.path,nargout=0)
        

        self.fmxt_dim = {
                         'nx' : int(eng.workspace["fmtx_nx"]),
                         'ny' : int(eng.workspace["fmtx_ny"]),
                         'nz' : int(eng.workspace["fmtx_nz"]),
                           }
        

        self.fmxt_mesh_limits = {'xmin': eng.workspace["fmtx_xmin"], 'xmax': eng.workspace["fmtx_xmax"],
                                 'ymin': eng.workspace["fmtx_ymin"], 'ymax': eng.workspace["fmtx_ymax"],
                                 'zmin': eng.workspace["fmtx_zmin"], 'zmax': eng.workspace["fmtx_zmax"],
                                 }
        

        self.fmtx_t = np.array(eng.workspace['fmtx_t']) 
        
        
        self.fmtx_t_err = np.array(eng.workspace['fmtx_t_err'])


    # Return fmxt dimensions
    def get_fmtx_dim(self):
        return self.fmxt_dim
        

    # Return 
    def get_fmtx_dim(self):
        return self.fmxt_dim        


    # Return fmxt tallies
    def get_fmtx_tallies(self):
        return self.fmtx_t
    

    # Return fmtx statistical errors
    def get_fmtx_errors(self):
        return self.fmtx_t_err    
    

    # Heatmap plot for matrix visulisation  
    @staticmethod
    def matrix_plot(matrix: np.array = None, cmap: str = "inferno"):
        sns.heatmap(matrix, cmap=cmap)
        
        plt.show()
        
    
    # Eigenvalue resolution
    def solve_eigenproblem(self, normalize: bool = True):
        self.k_eigv, self.fission_source = eigs(self.fmtx_t,k=1)


        if normalize:
            self.fission_source = self.fission_source.real/np.sum(self.fission_source.real)

        
        self.fission_source = self.fission_source.reshape(self.fmxt_dim['nx'],self.fmxt_dim['ny'])

        return self.k_eigv , self.fission_source
    

    # Plot Fission source
    @staticmethod
    def plot_fission_source(self, set_threshold: float = None, cmap: str = 'inferno'):
        plt.figure()

        
        xmin, xmax, ymin, ymax = self.fmxt_mesh_limits['xmin'], self.fmxt_mesh_limits['xmax'], self.fmxt_mesh_limits['ymin'], self.fmxt_mesh_limits['ymax']
        

        if set_threshold is not None: 
            self.fission_source = np.where(self.fission_source > set_threshold, self.fission_source, 0.0)
        

        plt.imshow(self.fission_source, extent = [xmin, xmax, ymin, ymax], cmap = cmap)

        cbar = plt.colorbar()

        plt.title("Fission Source by solving eigenvalue problem of the FMTX")


