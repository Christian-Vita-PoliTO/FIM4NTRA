from .common_imports import np, plt, sns, csr_matrix, eigs


class DetectorFissionSource: 


    def __init__(self, path, detector_name):
        self.path = path
        self.detector_name = detector_name
        self.matlab_det_extraction() # Internal method




    def matlab_det_extraction(self, normalize: bool = True): 
        """
            Reads a file, searches for the variable DET<detector_name>, extracts the fission source values,
            and determines the mesh dimensions based on the maximum values in columns 8, 9, and 10.

            :param file_path: Path to the file to read.
            :param detector_name: Base name of the detector (e.g., 'xy_fission_source').
            :return: Tuple (fission_source_values, nz, nx, ny), where:
                     - fission_source_values: list of fission source values (float).
                     - nz, nx, ny: deduced mesh dimensions.
        """        

        detector_variable = f"DET{self.detector_name}"  # Variabile da cercare
        reading_data = False
        fission_source_values = []
        fission_source_errors_values = []
        z_values = []  # Per dedurre nz
        y_values = []  # Per dedurre ny
        x_values = []  # Per dedurre nx
        
        with open(self.path, 'r') as file:
            for line in file:
                # Controlla se inizia la sezione della variabile detector
                if line.strip().startswith(detector_variable):
                    reading_data = True
                    continue  # Salta questa riga e inizia a leggere i dati successivi
                
                # Se siamo nella sezione dati, leggi i valori
                if reading_data:
                    if line.strip() == "];":  # Fine della variabile detector
                        break
                    
                    # Salta righe vuote o malformattate
                    if not line.strip() or line.strip().startswith("#"):
                        continue
                    
                    columns = line.split()
                    try:
                        # Colonne: 8 (z), 9 (y), 10 (x), 11 (fission source)
                        z = int(columns[7])
                        y = int(columns[8])
                        x = int(columns[9])
                        fission_source = float(columns[10])
                        fission_source_errors = float(columns[11])
                        
                        z_values.append(z)
                        y_values.append(y)
                        x_values.append(x)
                        fission_source_values.append(fission_source)
                        fission_source_errors_values.append(fission_source_errors)
                    except (IndexError, ValueError):
                        continue  # Salta righe mal formattate
                    
        # Deduzione delle dimensioni della mesh
        self.nz = max(z_values) if z_values else 0
        self.ny = max(y_values) if y_values else 0
        self.nx = max(x_values) if x_values else 0
        

        if self.nz == 0 or self.ny == 0 or self.nx == 0:
            raise ValueError("Non è stato possibile dedurre le dimensioni della mesh.")
        
        if not fission_source_values:
            raise ValueError(f"Nessun dato trovato per il detector {detector_variable}.")
        
        


        self.det_fission_source = np.array(fission_source_values) 
        self.det_fission_source_error = np.array(fission_source_errors_values)





    # Return fmxt tallies
    def get_det_tallies(self):
        return self.det_fission_source
    

    # Return fmtx statistical errors
    def get_det_errors(self):
        return self.det_fission_source_error    
    


    # Plot Fission source
    def plot_fission_source(self, cmap: str = 'inferno'):
        plt.figure()
        

        plt.matshow(self.det_fission_source.reshape(self.nx,self.ny,self.nz), cmap = cmap)


        cbar = plt.colorbar()


        plt.title("Detector Fission Source")


    
