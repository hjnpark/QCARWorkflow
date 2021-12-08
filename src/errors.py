class Error(Exception):
    pass

class NoDatasetError(Error):
    "Dataset can't be found."
    pass

class InvalidCommandError(Error):
    "Provided command is invalid."
    pass 

class QCFractalError(Error):
    "QCFractal error"
    pass

class QCEngineError(Error):
    "QCEngine error"
    pass
