from .common_imports import np, re, plt, sns, csr_matrix, eigs


class FissionMatrix: 


    def __init__(self,path):
        self.path = path
        self.read_fission_matrix_data() # Internal method


    def read_fission_matrix_data(self):

        
        with open(self.path, 'r') as file:
            lines = file.readlines()

        # Variabili per le dimensioni della matrice
        nx = ny = nz = 0
        fmtx_t = None
        matrix_size = 0

        # Pattern per leggere le dimensioni
        dim_pattern = re.compile(r"fmtx_\s*(\w+)\s*=\s*([\d\.\-E+]+)\s*;")
        value_pattern = re.compile(r"fmtx_t\s*\(\s*(\d+),\s*(\d+)\s*\)\s*=\s*([\d\.\-E+]+)\s*;\s*fmtx_t_err\s*\(\s*\d+,\s*\d+\s*\)\s*=\s*([\d\.\-E+]+)\s*;")

        for line in lines:
            # Legge le dimensioni
            dim_match = dim_pattern.match(line)
            if dim_match:
                key, value = dim_match.groups()
                if key == "nx":
                    nx = int(float(value))
                elif key == "ny":
                    ny = int(float(value))
                elif key == "nz":
                    nz = int(float(value))  # nz non è usato direttamente in 2D

                
            # Crea la matrice quando nx e ny sono noti
            if nx and ny and nz and fmtx_t is None:
                matrix_size = nx * ny * nz
                print(matrix_size)
                fmtx_t = np.zeros((matrix_size, matrix_size))

            # Legge i valori della matrice
            value_match = value_pattern.match(line)
            if value_match and fmtx_t is not None:
                i, j, value, error = value_match.groups()
                i, j = int(i) - 1, int(j) - 1  # Converti gli indici in base 0
                fmtx_t[i, j] = float(value)
        
        self.fmtx_t = np.array(fmtx_t)
        self.nx = nx
        self.ny = ny
        self.nz = nz


        self.solve_eigenproblem()
        
  


    # Return fmxt dimensions
    def get_fmtx_dim(self):
        
        self.fmxt_dim = {'nx':self.nx,'ny':self.ny,'nz':self.nz}

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

        self.get_fmtx_dim()
        self.fission_source = self.fission_source.reshape(self.fmxt_dim['nx'],self.fmxt_dim['ny'],self.fmxt_dim['nz'])

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


    @staticmethod
    def generate_ifc_file(T_distributions: list, material_name: str, reference_density: float,  nx_xmin_xmax: list, ny_ymin_ymax: list, nz_zmin_zmax: list, density_function = None, directory: str = ""):
       
        """
        Generate an IFC (Input File for Calculations) for a material, specifying the temperature and density distribution.

        Parameters:
        - material_name (str): Name of the material.
        - reference_density (float): Reference density for the material.
        - nx_xmin_xmax (list): [nx, xmin, xmax] for the x-dimension.
        - ny_ymin_ymax (list): [ny, ymin, ymax] for the y-dimension.
        - nz_zmin_zmax (list): [nz, zmin, zmax] for the z-dimension.
        - density_function (Callable, optional): User-defined function to compute density distribution.
        """

        flatten_T_dist = T_distributions.flatten()

        # Using user-defined relation (function) to generate density distribution.
        if density_function is not None: 
            rho_distribution = density_function(flatten_T_dist)
        else:
            rho_distribution = np.ones_like(flatten_T_dist)*reference_density
        
        
        # Valori fissi del modello TREAT
        XMIN, XMAX = nx_xmin_xmax[1], nx_xmin_xmax[2] 
        YMIN, YMAX = ny_ymin_ymax[1], ny_ymin_ymax[2] 
        ZMIN, ZMAX = nz_zmin_zmax[1], nz_zmin_zmax[2]

        # Corregge il primo valore di Tdist per il bug TMS di Serpent
        # self.T_distribution[0] += 0.001

        # Apre il file in modalità scrittura
        with open(f"{directory}{material_name}.ifc", 'w', newline='\n') as file:
            # Prima linea: TYPE MAT OUT
            file.write(f"2 {material_name} 0\n")

            # Terza linea: MESHTYPE
            file.write("1\n")

            # Quarta linea: mesh cartesiana con coordinate e dimensioni
            file.write(f"{nx_xmin_xmax[0]} {XMIN:.5f} {XMAX:.5f} {ny_ymin_ymax[0]} {YMIN:.5f} {YMAX:.5f} {nz_zmin_zmax[0]} {ZMIN:.5f} {ZMAX:.5f}\n")

            # Scrive i valori di densità e temperatura in ordine x, y, z
            i = 0
            for iz in range(nz_zmin_zmax[0]):
                for iy in range(ny_ymin_ymax[0]):
                    for ix in range(nx_xmin_xmax[0]):
                        file.write(f"{rho_distribution[i]:.6f} {flatten_T_dist[i]:.6f} ")
                        i += 1
                    file.write("\n") # Aggiunge una nuova riga al termine di ogni y-loop


        print(f"File '{material_name}.ifc' successfully created!")


