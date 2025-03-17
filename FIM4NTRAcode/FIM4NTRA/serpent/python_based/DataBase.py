from .common_imports import np, plt, sns, csr_matrix, re, eigs, os, contextmanager
from typing import Callable


class DataBase:


    def __init__(self, serpent_files_folder: str):
        
        self.serpent_files_folder = serpent_files_folder

        
    def build_database_non_uniform_T(self, list_T_file: list, matrix_xy_dimension: list, find_fisison_source: bool = False):
        
        """
        Build the database using the provided temperature files and dimensions.

        Parameters:
        - list_T_file (np.array): A list containing  are temperatures distribution, and values are file paths.
        - matrix_xy_dimension (list): List with three elements [nx, ny, nz] representing matrix dimensions.
        - find_fisison_source (bool): If True, calculate the fission source for each file.
        """


        self.number_of_data = len(list_T_file)

        self.mesh_dim_dic = {'nx': matrix_xy_dimension[0], 'ny': matrix_xy_dimension[1], 'nz': matrix_xy_dimension[2]}
        
        self.number_of_cells = int(self.mesh_dim_dic['nx']*self.mesh_dim_dic['ny']*self.mesh_dim_dic['nz']) 


        self.database_teperature_profiles = []
        

        for idataset,dataset in enumerate(list_T_file):

            self.database_teperature_profiles.append(dataset[0]) 

        self.database_teperature_profiles = np.array([self.database_teperature_profiles])
        self.database_teperature_profiles = self.database_teperature_profiles[0]

        # Iteration over data files
        self.data_base = np.zeros((self.number_of_data, self.number_of_cells, self.number_of_cells))
        self.k_eig_data_base = np.zeros((self.number_of_data))
        self.fission_source_data_base = np.zeros((self.number_of_data, self.number_of_cells,1))
        self.data_base_error = np.zeros((self.number_of_data, self.number_of_cells, self.number_of_cells))


        for idataset,dataset in enumerate(list_T_file):


            if find_fisison_source:
                results = self.read_fission_matrix_data(dataset[1]) # First input is just the fmtx
                self.data_base[idataset][:,:] = results[0]
                self.k_eig_data_base[idataset] = results[1]
                self.fission_source_data_base[idataset][:,:] = results[2]
                self.data_base_error[idataset][:,:] = results[3]

            else: 
                self.data_base[idataset][:,:] =  self.read_fission_matrix_data(dataset[1])[0]
                self.data_base_error[idataset][:,:] = self.read_fission_matrix_data(dataset[1])[3]

    # Build the database using isothermal cases as dataset: 
    def build_database(self, dic_T_file: dict, matrix_xy_dimension: list, find_fisison_source: bool = False):
        
        """
        Build the database using the provided temperature files and dimensions.

        Parameters:
        - dic_T_file (dict): A dictionary where keys are temperatures, and values are file paths.
        - matrix_xy_dimension (list): List with two elements [nx, ny] representing matrix dimensions.
        - find_fisison_source (bool): If True, calculate the fission source for each file.
        """


        self.number_of_data = len(list(dic_T_file.keys()))

        self.mesh_dim_dic = {'nx': matrix_xy_dimension[0], 'ny': matrix_xy_dimension[1], 'nz': matrix_xy_dimension[2]}
        
        self.number_of_cells = int(self.mesh_dim_dic['nx']*self.mesh_dim_dic['ny']*self.mesh_dim_dic['nz']) 


        self.data_base_teperatures = np.array([Temp for Temp in list(dic_T_file.keys())])
        
        # Iteration over data files
        self.data_base = np.zeros((self.number_of_data, self.number_of_cells, self.number_of_cells))
        self.k_eig_data_base = np.zeros((self.number_of_data))
        self.fission_source_data_base = np.zeros((self.number_of_data, self.number_of_cells,1))
        self.data_base_error = np.zeros((self.number_of_data, self.number_of_cells, self.number_of_cells))

        for iT, T in enumerate(dic_T_file):

            if find_fisison_source:
                results = self.read_fission_matrix_data(dic_T_file[T]) # First input is just the fmtx
                self.data_base[iT][:,:] = results[0]
                self.k_eig_data_base[iT] = results[1]
                self.fission_source_data_base[iT][:,:] = results[2]
                self.data_base_error[iT][:,:] = results[3]
            else: 
                self.data_base[iT][:,:] =  self.read_fission_matrix_data(dic_T_file[T])[0]
                self.data_base_error[iT][:,:] = self.read_fission_matrix_data(dic_T_file[T])[3]

    
    # def read_fission_matrix_data(self, path):
        
    #     self.eng.run(path, nargout=0) 

    #     fmtx_t = np.array(self.eng.workspace['fmtx_t'])

    #     k_eigv, fission_source = self.solve_fission_matrix(fmtx_t)

    #     return fmtx_t, k_eigv, fission_source
    



    def read_fission_matrix_data(self, path):

        
        with open(path, 'r') as file:
            lines = file.readlines()

        # Variabili per le dimensioni della matrice
        nx = ny = nz = None
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
                    self.nx = int(float(value))
                elif key == "ny":
                    self.ny = int(float(value))
                elif key == "nz":
                    self.nz = int(float(value))  # nz non è usato direttamente in 2D

        # Crea la matrice quando nx e ny sono noti
        if self.nx and self.ny and self.nz and fmtx_t is None:
            matrix_size = self.nx * self.ny * self.nz
            fmtx_t = np.zeros((matrix_size, matrix_size))
            fmtx_t_err = np.zeros((matrix_size, matrix_size))

        for line in lines:
            # Legge i valori della matrice
            value_match = value_pattern.match(line)
            if value_match and fmtx_t is not None:
                i, j, value, error = value_match.groups()
                i, j = int(i) - 1, int(j) - 1  # Converti gli indici in base 0
                fmtx_t[i, j] = float(value)
                fmtx_t_err[i, j] = float(error)
        

        
        k_eigv, fission_source = self.solve_fission_matrix(fmtx_t)

        return fmtx_t, k_eigv, fission_source, fmtx_t_err        

        

    
    def linear_interpolation(self, T_distribution):
        
        
        self.T_distribution = T_distribution

        self.T_distribution[0] = self.T_distribution[0] - 10E-5 
        self.T_distribution[-1] = self.T_distribution[-1] + 10E-5 

        fmtx_interpolated = np.zeros_like(self.data_base[0]) #fmtx_interpolated = np.zeros((self.number_of_cells, self.number_of_cells))

        for i_cell, T_cell in enumerate(T_distribution):
            

            if T_cell > self.data_base_teperatures[-1]:
                
                print(f"""Temperature of the cell number {i_cell}, is {T_cell}. It is above data base range.
                        The cell temperature is setted to the maximum temperature of the database.""")

                T_cell = self.data_base_teperatures[-1]

            elif T_cell < self.data_base_teperatures[0]:
                
                print(f"""Temperature of the cell number {i_cell}, is {T_cell}. It is below data base range.
                        The cell temperature is setted to the minimum temperature of the database.""")

                T_cell = self.data_base_teperatures[0]


            else:

                index_database_temperature_higher = (self.data_base_teperatures > T_cell).argmax() #next((T_DB for T_DB,  in self.data_base_teperatures if T_DB >= T_cell), None) 
                index_database_temperature_lower = index_database_temperature_higher - 1
                lower_DB_temperature = self.data_base_teperatures[index_database_temperature_lower] 
                higher_DB_temperature = self.data_base_teperatures[index_database_temperature_higher] 

                w = (T_cell - lower_DB_temperature)/(higher_DB_temperature - lower_DB_temperature)

                fmtx_higher = self.data_base[index_database_temperature_higher][i_cell, :]
                fmtx_lower  = self.data_base[index_database_temperature_lower][i_cell, :] 

                fmtx_interpolated[i_cell][:] = fmtx_lower*(1-w) + fmtx_higher*(w)

        self.fmtx_interpolated = fmtx_interpolated
        self.k_fmtx_interpolated, self.fission_source_fmtx_interpolated = self.solve_fission_matrix(fmtx_interpolated) 



    def linear_interpolation_non_uniform_T(self, T_distribution):

        """
        Perform linear interpolation of the fission matrix for a given temperature distribution.

        Parameters:
        - T_distribution (list): List of temperatures for each cell in the system.
        """


        
        self.T_distribution = T_distribution

        # self.T_distribution[0] = self.T_distribution[0] + 10E-5
        # self.T_distribution[-1] = self.T_distribution[-1] - 10E-5 

        fmtx_interpolated = np.zeros_like(self.data_base[0]) #fmtx_interpolated = np.zeros((self.number_of_cells, self.number_of_cells))

        fmtx_uncertainty_interpolated = np.zeros_like(self.data_base[0]) #fmtx_interpolated = np.zeros((self.number_of_cells, self.number_of_cells))


        fmtx_cell_counter = 0
        for i_cell in range(T_distribution.shape[0]):
            
            for j_cell in range(T_distribution.shape[1]):

                for k_cell in range(T_distribution.shape[2]):

                    if T_distribution[i_cell,j_cell,k_cell] > np.max(self.database_teperature_profiles[:,i_cell,j_cell,k_cell]):

                        print(f"""Temperature of the cell number {i_cell}, is { T_distribution[i_cell,j_cell,k_cell]}. It is above data base range.
                                The cell temperature is setted to the maximum temperature of the database.""")

                        T_distribution[i_cell,j_cell,k_cell] = np.max(self.database_teperature_profiles[:,i_cell,j_cell,k_cell])

                    elif T_distribution[i_cell,j_cell,k_cell] < np.min(self.database_teperature_profiles[:,i_cell,j_cell,k_cell]):

                        print(f"""Temperature of the cell number {i_cell}, is { T_distribution[i_cell,j_cell,k_cell]}. It is below data base range.
                                The cell temperature is setted to the minimum temperature of the database.""")

                        T_distribution[i_cell,j_cell,k_cell] = np.min(self.database_teperature_profiles[:,i_cell,j_cell,k_cell])


                    else:

                        index_database_temperature_higher = (self.database_teperature_profiles[:,i_cell,j_cell,k_cell] > T_distribution[i_cell,j_cell,k_cell]).argmax() #next((T_DB for T_DB,  in self.data_base_teperatures if T_DB >= T_distribution[i_cell,j_cell]), None) 
                        index_database_temperature_lower = index_database_temperature_higher - 1
                        lower_DB_temperature = self.database_teperature_profiles[:,i_cell,j_cell,k_cell][index_database_temperature_lower] 
                        higher_DB_temperature = self.database_teperature_profiles[:,i_cell,j_cell,k_cell][index_database_temperature_higher] 

                        w = (T_distribution[i_cell,j_cell,k_cell] - lower_DB_temperature)/(higher_DB_temperature - lower_DB_temperature)

                        fmtx_higher = self.data_base[index_database_temperature_higher][fmtx_cell_counter, :] 
                        fmtx_higher_uncertainty = self.data_base_error[index_database_temperature_lower][fmtx_cell_counter, :]
                        fmtx_lower  = self.data_base[index_database_temperature_lower][fmtx_cell_counter, :]  
                        fmtx_lower_uncertainty = self.data_base_error[index_database_temperature_lower][fmtx_cell_counter, :]

                        fmtx_interpolated[fmtx_cell_counter][:] = fmtx_lower*(1-w) + fmtx_higher*(w)

                        fmtx_uncertainty_interpolated[fmtx_cell_counter][:]  = np.sqrt((fmtx_lower_uncertainty*(1-w) + fmtx_higher_uncertainty*(w))**2)

                    fmtx_cell_counter += 1

            

        
        self.fmtx_interpolated = fmtx_interpolated
        self.k_fmtx_interpolated, self.fission_source_fmtx_interpolated = self.solve_fission_matrix(fmtx_interpolated) 
        
        self.fmtx_uncertainty_interpolated = fmtx_uncertainty_interpolated





    def linear_interpolation_non_uniform_T_uncertainity(self, T_distribution):

        """
        Perform linear interpolation of the fission matrix for a given temperature distribution.

        Parameters:
        - T_distribution (list): List of temperatures for each cell in the system.
        """


        
        self.T_distribution = T_distribution

        # self.T_distribution[0] = self.T_distribution[0] + 10E-5
        # self.T_distribution[-1] = self.T_distribution[-1] - 10E-5 

        fmtx_interpolated = np.zeros_like(self.data_base[0]) #fmtx_interpolated = np.zeros((self.number_of_cells, self.number_of_cells))


        fmtx_cell_counter = 0
        for i_cell in range(T_distribution.shape[0]):
            
            for j_cell in range(T_distribution.shape[1]):

                for k_cell in range(T_distribution.shape[2]):

                    if T_distribution[i_cell,j_cell,k_cell] > np.max(self.database_teperature_profiles[:,i_cell,j_cell,k_cell]):

                        print(f"""Temperature of the cell number {i_cell}, is { T_distribution[i_cell,j_cell,k_cell]}. It is above data base range.
                                The cell temperature is setted to the maximum temperature of the database.""")

                        T_distribution[i_cell,j_cell,k_cell] = np.max(self.database_teperature_profiles[:,i_cell,j_cell,k_cell])

                    elif T_distribution[i_cell,j_cell,k_cell] < np.min(self.database_teperature_profiles[:,i_cell,j_cell,k_cell]):

                        print(f"""Temperature of the cell number {i_cell}, is { T_distribution[i_cell,j_cell,k_cell]}. It is below data base range.
                                The cell temperature is setted to the minimum temperature of the database.""")

                        T_distribution[i_cell,j_cell,k_cell] = np.min(self.database_teperature_profiles[:,i_cell,j_cell,k_cell])


                    else:

                        index_database_temperature_higher = (self.database_teperature_profiles[:,i_cell,j_cell,k_cell] > T_distribution[i_cell,j_cell,k_cell]).argmax() #next((T_DB for T_DB,  in self.data_base_teperatures if T_DB >= T_distribution[i_cell,j_cell]), None) 
                        index_database_temperature_lower = index_database_temperature_higher - 1
                        lower_DB_temperature = self.database_teperature_profiles[:,i_cell,j_cell,k_cell][index_database_temperature_lower] 
                        higher_DB_temperature = self.database_teperature_profiles[:,i_cell,j_cell,k_cell][index_database_temperature_higher] 

                        w = (T_distribution[i_cell,j_cell,k_cell] - lower_DB_temperature)/(higher_DB_temperature - lower_DB_temperature)

                        fmtx_higher = self.data_base[index_database_temperature_higher][fmtx_cell_counter, :] # da modificare qua
                        fmtx_lower  = self.data_base[index_database_temperature_lower][fmtx_cell_counter, :]  # da modificare qua

                        fmtx_interpolated[fmtx_cell_counter][:] = fmtx_lower*(1-w) + fmtx_higher*(w)

                    fmtx_cell_counter += 1

            

        self.fmtx_interpolated = fmtx_interpolated
        self.k_fmtx_interpolated, self.fission_source_fmtx_interpolated = self.solve_fission_matrix(fmtx_interpolated) 











    # Method for ifc file generation

    def generate_ifc_file(self, material_name: str, reference_density: float,  nx_xmin_xmax: list, ny_ymin_ymax: list, nz_zmin_zmax: list, density_function = None, directory: str = ""):
       
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

        flatten_T_dist = self.T_distribution.flatten()

        # Using user-defined relation (function) to generate density distribution.
        if density_function is not None: 
            self.rho_distribution = density_function(flatten_T_dist)
        else:
            self.rho_distribution = np.ones_like(flatten_T_dist)*reference_density
        
        
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
                        file.write(f"{self.rho_distribution[i]:.6f} {flatten_T_dist[i]:.6f} ")
                        i += 1
                    file.write("\n") # Aggiunge una nuova riga al termine di ogni y-loop


        print(f"File '{material_name}.ifc' successfully created!")





    @staticmethod
    def solve_fission_matrix(fmtx_to_solve):

        k_eigv, fission_source = eigs(fmtx_to_solve,k=1)

        k_eigv = k_eigv.real

        fission_source = fission_source.real

        fission_source = fission_source /np.sum(fission_source)

        return k_eigv, fission_source
        

    @staticmethod
    def calculate_fission_source_uncertainty(keff,fs_flatten,fm,fm_uncertainty,nxnynz_list):

        fs_uncertainty = np.zeros(nxnynz_list[0]*nxnynz_list[1]*nxnynz_list[2])
        # fs_uncertainty.reshape(nxnynz_list[0],nxnynz_list[1],nxnynz_list[2])
        fm_uncertainty_extended = fm_uncertainty*fm
        
        for i in range(fm.shape[0]):  
            
            error_diff = 1
            sigma_fs_guess = 1.0
            while_iteriation = 0

            while error_diff > 10E-6 and while_iteriation < 10:

                if while_iteriation == 0:
                    sigma_fs_old = sigma_fs_guess
                
                
                sigma_fs = 1/(keff-fm[i,i]) * ( np.sqrt( np.sum( (fs_flatten/keff*fm_uncertainty_extended[i,:])**2 + (fm[i,:]*sigma_fs_old)**2 ) ) )
                

                error_diff = (sigma_fs - sigma_fs_old)/sigma_fs_old
                
                while_iteriation += 1

                sigma_fs_old = sigma_fs
                    
            
            fs_uncertainty[i] =  sigma_fs

        return fs_uncertainty




@contextmanager
def change_directory(new_directory):
    # Save the current directory
    current_directory = os.getcwd()
    try:
        # Change to the new directory
        os.chdir(new_directory)
        print(f"Entered directory: {os.getcwd()}")
        yield  # Pass control to the `with` block
    finally:
        # Return to the original directory
        os.chdir(current_directory)
        print(f"Returned to directory: {os.getcwd()}")