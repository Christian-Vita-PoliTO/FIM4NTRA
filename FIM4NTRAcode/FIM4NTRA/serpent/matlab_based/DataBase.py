from .common_imports import np, matlab, plt, sns, csr_matrix, eigs, os, contextmanager
from typing import Callable


class DataBase:

    """
    This class is responsible for managing and building a database for fission matrix calculations
    using input files and temperature-dependent data.
    """


    def __init__(self, serpent_files_folder: str):
      
        """
        Initialize the database object with the folder containing Serpent output files.

        Parameters:
        - serpent_files_folder (str): Path to the folder containing Serpent files.
        """        
        
        self.serpent_files_folder = serpent_files_folder

        
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

        self.mesh_dim_dic = {'nx': matrix_xy_dimension[0], 'ny': matrix_xy_dimension[1]}
        
        self.number_of_cells = int(self.mesh_dim_dic['nx']*self.mesh_dim_dic['ny']) 

        self.eng = matlab.engine.start_matlab()

        self.data_base_teperatures = np.array([Temp for Temp in list(dic_T_file.keys())])
        
        # Iteration over data files
        self.data_base = np.zeros((self.number_of_data, self.number_of_cells, self.number_of_cells))
        self.k_eig_data_base = np.zeros((self.number_of_data))
        self.fission_source_data_base = np.zeros((self.number_of_data, self.number_of_cells,1))

        for iT, T in enumerate(dic_T_file):

            if find_fisison_source:
                results = self.read_fission_matrix_data(dic_T_file[T]) # First input is just the fmtx
                self.data_base[iT][:,:] = results[0]
                self.k_eig_data_base[iT] = results[1]
                self.fission_source_data_base[iT][:,:] = results[2]
            else: 
                self.data_base[iT][:,:] =  self.read_fission_matrix_data(dic_T_file[T])[0]




    def build_database_non_uniform_T(self, list_T_file: list, matrix_xy_dimension: list, find_fisison_source: bool = False):
        
        """
        Build the database using the provided temperature files and dimensions.

        Parameters:
        - list_T_file (np.array): A list containing  are temperatures distribution, and values are file paths.
        - matrix_xy_dimension (list): List with two elements [nx, ny] representing matrix dimensions.
        - find_fisison_source (bool): If True, calculate the fission source for each file.
        """


        self.number_of_data = len(list_T_file)

        self.mesh_dim_dic = {'nx': matrix_xy_dimension[0], 'ny': matrix_xy_dimension[1]}
        
        self.number_of_cells = int(self.mesh_dim_dic['nx']*self.mesh_dim_dic['ny']) 

        self.eng = matlab.engine.start_matlab()


        self.database_teperature_profiles = []
        

        for idataset,dataset in enumerate(list_T_file):

            self.database_teperature_profiles.append(dataset[0]) 

        self.database_teperature_profiles = np.array([self.database_teperature_profiles])
        self.database_teperature_profiles = self.database_teperature_profiles[0]

        # Iteration over data files
        self.data_base = np.zeros((self.number_of_data, self.number_of_cells, self.number_of_cells))
        self.k_eig_data_base = np.zeros((self.number_of_data))
        self.fission_source_data_base = np.zeros((self.number_of_data, self.number_of_cells,1))


        for idataset,dataset in enumerate(list_T_file):


            if find_fisison_source:
                results = self.read_fission_matrix_data(dataset[1]) # First input is just the fmtx
                self.data_base[idataset][:,:] = results[0]
                self.k_eig_data_base[idataset] = results[1]
                self.fission_source_data_base[idataset][:,:] = results[2]
            else: 
                self.data_base[idataset][:,:] =  self.read_fission_matrix_data(dataset[1])[0]



    def linear_interpolation_non_uniform_T(self, T_distribution):

        """
        Perform linear interpolation of the fission matrix for a given temperature distribution.

        Parameters:
        - T_distribution (list): List of temperatures for each cell in the system.
        """


        
        self.T_distribution = T_distribution

        self.T_distribution[0] = self.T_distribution[0] - 10E-5 
        self.T_distribution[-1] = self.T_distribution[-1] + 10E-5 

        fmtx_interpolated = np.zeros_like(self.data_base[0]) #fmtx_interpolated = np.zeros((self.number_of_cells, self.number_of_cells))

        fmtx_cell_counter = 0
        for i_cell in range(T_distribution.shape[0]):
            
            for j_cell in range(T_distribution.shape[1]):

                if T_distribution[i_cell,j_cell] > np.max(self.database_teperature_profiles[:,i_cell,j_cell]):

                    print(f"""Temperature of the cell number {i_cell}, is { T_distribution[i_cell,j_cell]}. It is above data base range.
                            The cell temperature is setted to the maximum temperature of the database.""")

                    T_distribution[i_cell,j_cell] = np.max(self.database_teperature_profiles[:,i_cell,j_cell])

                elif T_distribution[i_cell,j_cell] < np.min(self.database_teperature_profiles[:,i_cell,j_cell]):

                    print(f"""Temperature of the cell number {i_cell}, is { T_distribution[i_cell,j_cell]}. It is below data base range.
                            The cell temperature is setted to the minimum temperature of the database.""")

                    T_distribution[i_cell,j_cell] = np.min(self.database_teperature_profiles[:,i_cell,j_cell])


                else:

                    index_database_temperature_higher = (self.database_teperature_profiles[:,i_cell,j_cell] > T_distribution[i_cell,j_cell]).argmax() #next((T_DB for T_DB,  in self.data_base_teperatures if T_DB >= T_distribution[i_cell,j_cell]), None) 
                    index_database_temperature_lower = index_database_temperature_higher - 1
                    lower_DB_temperature = self.database_teperature_profiles[:,i_cell,j_cell][index_database_temperature_lower] 
                    higher_DB_temperature = self.database_teperature_profiles[:,i_cell,j_cell][index_database_temperature_higher] 

                    w = (T_distribution[i_cell,j_cell] - lower_DB_temperature)/(higher_DB_temperature - lower_DB_temperature)

                    fmtx_higher = self.data_base[index_database_temperature_higher][fmtx_cell_counter, :] # da modificare qua
                    fmtx_lower  = self.data_base[index_database_temperature_lower][fmtx_cell_counter, :]  # da modificare qua

                    fmtx_interpolated[fmtx_cell_counter][:] = fmtx_lower*(1-w) + fmtx_higher*(w)
                
                fmtx_cell_counter += 1

            

        self.fmtx_interpolated = fmtx_interpolated
        self.k_fmtx_interpolated, self.fission_source_fmtx_interpolated = self.solve_fission_matrix(fmtx_interpolated) 



    
    def read_fission_matrix_data(self, path):

        """
        Read the fission matrix data from the specified file path using MATLAB.

        Parameters:
        - path (str): File path of the Serpent output.

        Returns:
        - tuple: Fission matrix, eigenvalue, and fission source.
        """


        self.eng.run(path, nargout=0) 

        fmtx_t = np.array(self.eng.workspace['fmtx_t'])

        k_eigv, fission_source = self.solve_fission_matrix(fmtx_t)

        return fmtx_t, k_eigv, fission_source
    

    
    def linear_interpolation(self, T_distribution):

        """
        Perform linear interpolation of the fission matrix for a given temperature distribution.

        Parameters:
        - T_distribution (list): List of temperatures for each cell in the system.
        """


        
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





    # Method for ifc file generation

    def generate_ifc_file(self, material_name: str, reference_density: float,  nx_xmin_xmax: list, ny_ymin_ymax: list, nz_zmin_zmax: list, density_function = None ):
       
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


        # Using user-defined relation (function) to generate density distribution.
        if density_function is not None: 
            self.rho_distribution = density_function(self.T_distribution)
        else:
            self.rho_distribution = np.ones_like(self.T_distribution)*reference_density
        
        
        # Valori fissi del modello TREAT
        XMIN, XMAX = nx_xmin_xmax[1], nx_xmin_xmax[2] 
        YMIN, YMAX = ny_ymin_ymax[1], ny_ymin_ymax[2] 
        ZMIN, ZMAX = nz_zmin_zmax[1], nz_zmin_zmax[2]

        # Corregge il primo valore di Tdist per il bug TMS di Serpent
        self.T_distribution[0] += 0.001

        # Apre il file in modalità scrittura
        with open(f"{material_name}.ifc", 'w', newline='\n') as file:
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
                        file.write(f"{self.rho_distribution[i]:.6f} {self.T_distribution[i]:.6f} ")
                        i += 1
                    file.write("\n")  # Aggiunge una nuova riga al termine di ogni y-loop

        print(f"File '{material_name}.ifc' successfully created!")





    @staticmethod
    def solve_fission_matrix(fmtx_to_solve):

        """
        Solve the fission matrix to compute the dominant eigenvalue and fission source.

        Parameters:
        - fmtx_to_solve (numpy.ndarray): Fission matrix to solve.

        Returns:
        - tuple: Dominant eigenvalue and normalized fission source.
        """

        

        k_eigv, fission_source = eigs(fmtx_to_solve,k=1)

        k_eigv = k_eigv.real

        fission_source = fission_source.real

        fission_source = fission_source/np.sum(fission_source)

        return k_eigv, fission_source
        






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