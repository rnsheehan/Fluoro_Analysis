# Import libraries
# You should try an import the bare minimum of modules
import sys # access system routines
import os

# add path to our file

import Fluoro_Data_Mangling

# The aim of this project is to perform some basic plotting and filtering of measured fluoresence data
# R. Sheehan 27 - 3 - 2019
 
def main():
    pass

if __name__ == '__main__':
    main()

    pwd = os.getcwd() # get current working directory

    print(pwd)

    Fluoro_Data_Mangling.Analyse_All_Fl_Data(False, True)

    #Fluoro_Data_Mangling.Analyse_Fl_Data()
