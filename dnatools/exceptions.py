class configFileError (Exception):
    def __init__(self, file_name, message = None):
        if message == None:
            message = "Error: Configuration file {} is not found. Check file name and path.".format(file_name)
        super(configFileError,self).__init__(message)

class dnaFileError (Exception):
    def __init__(self, file_name, message = None):
        if message == None:
            message = "Error: DNA file {} is not found. Check file name and path.".format(file_name)
        super(dnaFileError,self).__init__(message)

class configSectionError (Exception):
    def __init__(self, section_name, message = None):
        if message == None:
            message = "Error: Configuration file section {} not found. " \
                      " This section is mandatory in config file".format(section_name)
        super(configSectionError,self).__init__(message)


class configOptionError (Exception):
    def __init__(self, section_name, option_name,message = None):
        if message == None:
            message = "Error: Configuration file option {1} not found. \n" \
                      " Option {1} in section {0} is mandatory in config file"\
                        .format(section_name, option_name)
        super(configOptionError,self).__init__(message)