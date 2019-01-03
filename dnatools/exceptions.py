class configFileError (Exception):
    def __init__(self, file_name, message = None):
        if message == None:
            message = "Error: Configuration file {0} is not found.\n" \
                      "This program requires a valid configuration file in order to run.\n" \
                      "Check the file name and path of file ({0})."\
                .format(file_name)
        super(configFileError,self).__init__(message)

class dnaFileError (Exception):
    def __init__(self, file_name, message = None):
        if message == None:
            message = "Error: DNA file {0} does not exist.\n" \
                      "This program requires a DNA file in order to run.\n" \
                      "Check the file name and path of file ({0})."\
                        .format(file_name)
        super(dnaFileError,self).__init__(message)

class configSectionError (Exception):
    def __init__(self, section_name, message = None):
        if message == None:
            message = "Error: Configuration file section {0} not found.\n" \
                      "This program requires a section in configuration file labelled {0}\n"\
                        "Please check and correct the configuration file "\
                .format(section_name)
        super(configSectionError,self).__init__(message)


class configOptionError (Exception):
    def __init__(self, section_name, option_name,message = None):
        if message == None:
            message = "Error: Option \"{1}\" not found in configuration file.\n" \
                      "Option {1} in section {0} is required in configuration file.\n" \
                      "Please check and correct the configuration file " \
                .format(section_name, option_name)
        super(configOptionError,self).__init__(message)