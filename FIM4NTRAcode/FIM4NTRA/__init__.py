# FIM4NTRA/__init__.py

# Inizializza le classi a seconda della configurazione
def configure(code="serpent", use_matlab=False):

    """
    Configure the module to use specific implementations of the classes FissionMatrix, DataBase,
    and DetectorFissionSource based on the selected computational code and backend.

    Parameters:
    - code (str, optional): The computational code to use. Options are:
        * "serpent" (default)
        * "openmc"
    - use_matlab (bool, optional): If True, use MATLAB-based implementations (only for "serpent" code).

    Raises:
    - ValueError: If an unsupported computational code is provided.
    """
    
        
    global FissionMatrix, DataBase, DetectorFissionSource

    if code.lower() == "serpent":
        if use_matlab:
            from .serpent.matlab_based.FissionMatrix import FissionMatrix
            from .serpent.matlab_based.DataBase import DataBase
            from .serpent.matlab_based.DetectorFissionSource import DetectorFissionSource
        else:
            from .serpent.python_based.FissionMatrix import FissionMatrix
            from .serpent.python_based.DataBase import DataBase
            from .serpent.python_based.DetectorFissionSource import DetectorFissionSource

    elif code.lower() == "openmc":
        from .openmc.FissionMatrix import FissionMatrix
        from .openmc.DataBase import DataBase
        from .openmc.DetectorFissionSource import DetectorFissionSource

    else:
        raise ValueError(f"Unsupported code: {code}")