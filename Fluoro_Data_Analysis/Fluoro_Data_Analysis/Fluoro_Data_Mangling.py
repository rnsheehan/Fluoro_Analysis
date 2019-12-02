
import os
import glob
import re
import sys # access system routines, including writing console output to file

import math
import scipy
import numpy as np
import matplotlib.pyplot as plt

import pandas
import pprint

import Common
import Plotting

# Take a look at the fluoresence data file format
# First task: Import the data, make some basic plots of the measured spectra: Done
# Second Task: Investigate some basic filtering, try and re-produce the mask filter that MC has
# Can now locate and distinguish scattering and fluoresence peaks
# Need to remove the scattering peaks from the data
# Add a function for removing the scattering peaks
# function will eliminate scattering peaks
# then interpolate missing data based on fluoresence peak
# Make 2D plot of the filtered data
# R. Sheehan 27 - 3 - 2019

MOD_NAME_STR = "Fluoro_Data_Mangling" # use this in exception handling messages

def Analyse_Fl_Data(perform_analysis = True):
    # run the various functions needed to analyse the FL data
    # R. Sheehan 28 - 3 - 2019

    FUNC_NAME = ".Analyse_Fl_Data()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        #DATA_HOME = "C:/Users/Robert/Research/CAPPA/Data/EEM_fluorescence/";
        #DATA_HOME = "C:/Users/Robert/Research/CAPPA/Data/EEM_fluorescence/3_Collection/";
        DATA_HOME = "C:/Users/Robert/Research/CAPPA/Data/EEM_fluorescence/3_Collection/2_Test_high_and_medium_gain/";
        
        if os.path.isdir(DATA_HOME):
            os.chdir(DATA_HOME)

            #MCfile = 'CIT_Water_02Mar2015.txt' # this is the file to be read
            #MCfile = 'UCD_old_samples.txt' # this is the file to be read
            
            # New data from PL 12 - 6 - 2019
            #MCfile = 'Result(16.46.22_05.06.2019).txt' # this is the file to be read
            #MCfile = 'Result(18.28.27_05.06.2019)_2.txt' # this is the file to be read
            #MCfile = 'Result(18.28.27_05.06.2019)_3.txt' # this is the file to be read
            #MCfile = 'Result(18.28.27_05.06.2019)_4.txt' # this is the file to be read
            #MCfile = 'Result(17.03.12_10.06.2019)_5.txt' # this is the file to be read

            # New Data from 30 - 7 - 2019
            #MCfile = 'UCD_3Collection_2nd_set_EXC_200-450.txt' # this is the file to be read
            #MCfile = 'UCD_3Collection_3nd_set_EXC_200-450.txt' # this is the file to be read
            #MCfile = 'UCD_3Collection_EXC_200-450.txt' # this is the file to be read
            MCfile = 'UCD_3Collection_EXC_200-450_medium_gain.txt' # this is the file to be read

            the_data = Read_Fl_File(MCfile) # read the data into memory

            fl_preamble, fl_locs = Parse_Fl_File(the_data, False) # Parse the file for information

            if perform_analysis:
                #sample_name = 'Tap_water'; #sample_name = 'Ground_H20'; #sample_name = 'Inside_Puddle'; #sample_name = 'Puddle_Mud'
                #sample_name = 'RK1'; sample_name = 'OF'; #sample_name = 'US'; #sample_name = 'TS1'; 
                sample_name = 'TS2'; 
                sample_Data_set = Extract_Fl_Data(sample_name, the_data, fl_preamble, fl_locs)
            
                #di_sample_name = 'DIW'
                #di_sample_name = 'DI_water'; 
                di_sample_name = 'DI_Water'; 
                di_Data_set = Extract_Fl_Data(di_sample_name, the_data, fl_preamble, fl_locs)

                # Compute the difference between the measured sample and DI
                sample_Data_set[2] = Compute_Spectral_Difference(sample_Data_set, di_Data_set)
            
                fig_name = ''
                #Plot_All_Spectra(sample_Data_set[0], sample_Data_set[1], sample_Data_set[2], fig_name, True)

                #choice = 5
                #FL_peak_data, RES_peak_data = Locate_Peaks(sample_Data_set[0], sample_Data_set[1], sample_Data_set[2], choice, True)
                #Remove_Scatt_Peaks_Spctrm(sample_Data_set[0], sample_Data_set[1], sample_Data_set[2], choice, RES_peak_data, True)
                #Extract_Peak_Data(sample_Data_set[0], sample_Data_set[1], sample_Data_set[2])

                Remove_Scatt_Peaks(sample_Data_set[0], sample_Data_set[1], sample_Data_set[2])

                #Plot_All_Spectra(sample_Data_set[0], sample_Data_set[1], sample_Data_set[2], fig_name, True)

                # Write the Wavelength and Filtered Data to a File
                #Common.write_data(sample_name+'_Scan_WL.txt', sample_Data_set[0])
                #Common.write_data(sample_name+'_Excit_WL.txt', sample_Data_set[1])
                #Common.write_matrix(sample_name+'_Filt_Data.txt', sample_Data_set[2])

                Plot_Contours(sample_Data_set[0], sample_Data_set[1], sample_Data_set[2], sample_name, True)

                # Compute the integral of the data set
                Integrate_Fl_Data(sample_Data_set, [300.0, 500.0], [225.0, 275.0], loud = True)

                del di_Data_set; del sample_Data_set;
                
            del the_data; 
        else:
            raise EnvironmentError
    except EnvironmentError:
        print(ERR_STATEMENT)
        print("Directory:",DATA_HOME," not found")
    except Exception as e:
        print(ERR_STATEMENT); 
        print(e); 

def Analyse_All_Fl_Data(make_plots, compute_integrals):

    FUNC_NAME = ".Analyse_All_Fl_Data()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        #DATA_HOME = "C:/Users/Robert/Research/CAPPA/Data/EEM_fluorescence/";
        #DATA_HOME = "C:/Users/Robert/Research/CAPPA/Data/EEM_fluorescence/3_Collection/";
        #DATA_HOME = "C:/Users/Robert/Research/CAPPA/Data/EEM_fluorescence/3_Collection/2_Test_high_and_medium_gain/";
        DATA_HOME = "C:/Users/Robert/Research/CAPPA/Data/EEM_fluorescence/4_Collection/";

        if os.path.isdir(DATA_HOME):
            os.chdir(DATA_HOME)

            #MCfile = 'UCD_old_samples.txt' # this is the file to be read

            # New data from PL 12 - 6 - 2019
            #MCfile = 'Result(16.46.22_05.06.2019).txt' # this is the file to be read
            #MCfile = 'Result(18.28.27_05.06.2019)_2.txt' # this is the file to be read
            #MCfile = 'Result(18.28.27_05.06.2019)_3.txt' # this is the file to be read
            #MCfile = 'Result(18.28.27_05.06.2019)_4.txt' # this is the file to be read
            
            #MCfile = 'Result(17.03.12_10.06.2019)_5.txt' # this is the file to be read

            # New Data 30 - 7 - 2019
            #MCfile = 'UCD_3Collection_2nd_set_EXC_200-450.txt' # this is the file to be read
            #MCfile = 'UCD_3Collection_3nd_set_EXC_200-450.txt' # this is the file to be read
            #MCfile = 'UCD_3Collection_EXC_200-450.txt' # this is the file to be read
            #MCfile = 'UCD_3Collection_EXC_200-450_medium_gain.txt' # this is the file to be read
            
            # New Data 11 - 11 - 2019
            MCfile = 'UCD_4_Collection.txt' # this is the file to be read

            the_data = Read_Fl_File(MCfile) # read the data into memory

            fl_preamble, fl_locs = Parse_Fl_File(the_data, False) # Parse the file for information

            # Extract the DIW data set to be extracted from each data set
            di_sample_name = 'DI_water'; 
            #di_sample_name = 'DI_Water'; 
            di_Data_set = Extract_Fl_Data(di_sample_name, the_data, fl_preamble, fl_locs)

            # lists of integration ranges in the form [lambda_em, lambda_ex]
            peak_T = [[250,350],[200, 250]]
            peak_M = [[350,400],[250, 300]]
            peak_A1 = [[400,460],[220, 250]]
            peak_A2 = [[330,370],[235, 260]]
            peak_C1 = [[400,500],[300, 340]]
            peak_C2 = [[400,500],[300, 370]]

            ranges = [peak_T, peak_M, peak_A1, peak_A2, peak_C1, peak_C2]
            range_labels = ['T', 'M', 'A[1]', 'A[2]', 'C[1]', 'C[2]']

            if compute_integrals:
                integrationpath = "%(v1)s_EEM_Integration.txt"%{"v1":MCfile.replace(".txt","")}
                integrationfile = open(integrationpath, "w")
                integrationfile.write("Integration Results %(v1)s\n"%{"v1":MCfile})
                integrationfile.write("Sample name, Peak Label, EM start / nm, EM stop / nm, EX start / nm, EX stop / nm, Integral Value, Ratio of Total Integral\n")

            for x in fl_locs:
                if x != di_sample_name:
                    print('Analysing Dataset:', x, 'at location:', fl_locs[x])

                    # Extract the raw data from the file
                    sample_Data_set = Extract_Fl_Data( x, the_data, fl_preamble, fl_locs)

                    # Compute the difference between the measured sample and DI
                    sample_Data_set[2] = Compute_Spectral_Difference(sample_Data_set, di_Data_set)

                    Remove_Scatt_Peaks(sample_Data_set[0], sample_Data_set[1], sample_Data_set[2], False)

                    if make_plots:
                        # write the filtered data to a file and make a plot of the data
                        Common.write_data(x+'_Scan_WL.txt', sample_Data_set[0])
                        Common.write_data(x+'_Excit_WL.txt', sample_Data_set[1])
                        Common.write_matrix(x+'_Filt_Data.txt', sample_Data_set[2])

                        Plot_Contours(sample_Data_set[0], sample_Data_set[1], sample_Data_set[2], x, False)

                    if compute_integrals:
                        # integrate the filtered data sets over the various regions

                        res_string = "%(v1)s, %(v2)s, %(v3)d, %(v4)d, %(v5)d, %(v6)d, %(v7)0.5f, %(v8)0.5f\n"

                        # first compute the total integral
                        em_range = [sample_Data_set[0][0], sample_Data_set[0][-1]]
                        ex_range = [sample_Data_set[1][0], sample_Data_set[1][-1]]
                        
                        result_total = Integrate_Fl_Data(sample_Data_set, em_range, ex_range)

                        out_string = res_string%{"v1":x, "v2":"Total", "v3":em_range[0], "v4":em_range[1], "v5":ex_range[0], "v6":ex_range[1], "v7":result_total, "v8":result_total/result_total}

                        integrationfile.write(out_string)

                        # then compute the integral over the domain subsets
                        for i in range(0, len(ranges), 1): 
                            em_range = ranges[i][0]
                            ex_range = ranges[i][1]
                            result = Integrate_Fl_Data(sample_Data_set, em_range, ex_range)

                            out_string = res_string%{"v1":x, "v2":range_labels[i], "v3":em_range[0], "v4":em_range[1], "v5":ex_range[0], "v6":ex_range[1], "v7":result, "v8":(result / result_total)}

                            integrationfile.write(out_string)

                    del sample_Data_set;
                    
            if compute_integrals:
                integrationfile.close()

            del di_Data_set; del the_data; del fl_preamble; del fl_locs; 
        else:
            raise EnvironmentError
    except EnvironmentError:
        print(ERR_STATEMENT)
        print("Directory:",DATA_HOME," not found")
    except Exception as e:
        print(ERR_STATEMENT); 
        print(e);

def Integrate_Fl_Data(data_set, scan_range, ex_range, loud = False):
    
    # Compute the integral of a data set that is in memory
    # data_set comprises [scan_wl, ex_wl, fl_values]
    # scan_range and ex_range specify the limits of the integration region
    # they should be inside scan_wl and ex_wl
    # R. Sheehan 7 - 10 - 2019

    FUNC_NAME = ".Integrate_Fl_Data()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        c1 = True if data_set is not None else False
        c2 = True if scan_range is not None else False
        c3 = True if ex_range is not None else False
        c4 = True if scan_range[0] >= data_set[0][0] else False
        c5 = True if scan_range[1] <= data_set[0][-1] else False
        c6 = True if ex_range[0] >= data_set[1][0] else False
        c7 = True if ex_range[1] <= data_set[1][-1] else False
        c8 = True if scan_range[1] > scan_range[0] else False
        c9 = True if ex_range[1] > ex_range[0] else False
        
        c10 = c1 and c2 and c3 and c4 and c5 and c6 and c7 and c8 and c9

        if c10:

            xs = np.where(data_set[0] == scan_range[0])[0][0]
            xe = np.where(data_set[0] == scan_range[1])[0][0]
            ys = np.where(data_set[1] == ex_range[0])[0][0]
            ye = np.where(data_set[1] == ex_range[1])[0][0]

            integral = 0.0
            term = 0.0
            yindx = ys
            while yindx < ye:
                dy = data_set[1][yindx+1] - data_set[1][yindx]
                xindx = xs
                while xindx < xe:
                    dx = data_set[0][xindx+1] - data_set[0][xindx]
                    term = data_set[2][xindx][yindx] * dx * dy
                    integral = integral + term
                    xindx = xindx + 1
                yindx = yindx + 1

            if loud:
                print("The integral of the input data set is",integral)

            return integral
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT)
        print(e)

def Read_Fl_File(fl_filename):
    # Read an entire Fl data file into memory
    # use pandas.read_csv as this is now standard
    # R. Sheehan 27 - 3 - 2019

    FUNC_NAME = ".Read_Fl_File()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        if glob.glob(fl_filename):

            data = pandas.read_csv(fl_filename) # Use pandas to read the entire file
            
            print("\nFile:",fl_filename,"located")
            print("Data from file is saved\n")

            # print the file preamble
            #print(MCdata[:38]) # display first 38 columns
            #print(MCdata.loc[3][0]); print(MCdata.loc[4][0]); print(MCdata.loc[5][0]); print(MCdata.loc[6][0]);
            #print(MCdata.loc[10][0]); print(MCdata.loc[11][0]); print(MCdata.loc[12][0]);
            #print(MCdata.loc[13][0]); print(MCdata.loc[14][0]); print(MCdata.loc[15][0]);
            
            return data
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT); 
        print(e);
    
def Parse_Fl_File(fl_data, loud = False):    
    # extract the fl data file pre-amble, save it as a dictionary
    # extract the names and locations of all the samples within the file
    # return [pre_amble_dict, sample_dict]
    # R. Sheehan 27 - 3 - 2019

    FUNC_NAME = ".Parse_Fl_File()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        if fl_data is not None:
            # Extract the scan parameters from the pre-amble
            # NoSamples r = 3, NoScanDims r = 4, NoScanPoints(1) r = 5, NoScanPoints(2) r = 6
            # Scan Type1 r = 10, Scan Start1 r = 11, Scan End1 r= 12 
            # Scan Type2 r = 13, Scan Start2 r = 14, Scan End2 r= 15
            # Return a dictionary based on the data pre-amble
            # Return a dictionary based on the sample location within the file
                
            param_labels = []; param_values = []; 
            pre_amble_dict = {}
            rows = [3, 4, 5, 6, 10, 11, 12, 13, 14, 15] # relevant data is stored in these rows of fl_data
            for r in rows:
                split_str = fl_data.loc[r][0].split('\t')
                param_labels.append(split_str[0]); 
                param_values.append(int(split_str[1]) if split_str[1].isdigit() else split_str[1]); 

            pre_amble_dict = dict( zip(param_labels, param_values) )

            if loud:
                print("\nData Preamble:")
                pprint.pprint(pre_amble_dict) # print the dictionary contents alphabetically
                print("\n")

            # Find the locations and names of the individual data sets
            # Store the data as a dictionary
            param_labels = []; param_values = []; 
            sample_dict = {}
            
            count = 0; 
            for r in range(0, len(fl_data), 1):
                if "SampleID" in fl_data.loc[r][0]:
                    if loud: print("r:",r,fl_data.loc[r][0])
                    param_values.append(r); 
                    param_labels.append(fl_data.loc[r][0].split('\t')[1])
                    count = count + 1
            
            print("\nStated No. Samples:",pre_amble_dict['NoSamples'])
            print("Found Samples:",count)
            print("Found Samples == Stated No. Samples?:", count == pre_amble_dict['NoSamples'])
            
            sample_dict = dict( zip(param_labels, param_values) )
            
            print("\nSamples:")
            pprint.pprint(sample_dict)
            print("\n")

            return [pre_amble_dict, sample_dict]
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT); 
        print(e);

def Extract_Fl_Data(smpl_name, fl_data, fl_pre_amble, fl_smpl_locs, loud = False):
    # extract the measured scan data from a file
    # get the scan1 wavelengths, array of length NoScanPoints(1), this is the scan WL data
    # get the scan2 wavelengths, array of length NoScanPoints(2), this is the excitation WL data
    # get the excitation data, matrix of size NoScanPoints(1)*NoScanPoints(2), each column of measData corresponds to a different excitation WL 
    # measured data is starts row fl_smpl_locs[smpl_name] + 6
    # measured data stops on row fl_smpl_locs[smpl_name] + 6 + fl_pre_amble['NoScanPoints(1)'] - 1
    # R. Sheehan 28 - 3 - 2019

    FUNC_NAME = ".Extract_Fl_Data()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        c1 = fl_data is not None
        c2 = fl_pre_amble is not None
        c3 = fl_smpl_locs is not None
        c4 = Common.dict_contains_key(fl_smpl_locs, smpl_name)

        if c1 and c2 and c3 and c4:
            # Extract the Scan2 WL
            rWL2 = fl_smpl_locs[smpl_name] + 5 # this is the row that the scan2 WL are stored on
            nWL2 = fl_pre_amble['NoScanPoints(2)'] # this is the number of scan2 WL
            scan2WL = np.zeros([ nWL2 ]) # declarate memory space for numpy array
            WL2list = fl_data.loc[ rWL2 ][0].split('\t') # save the scan2 WL as a list

            if loud:
                print("\nScan2 WL:")
                print("r:", rWL2, fl_data.loc[ rWL2 ][0])
                print(nWL2)
                print(len(WL2list))

            count = 0
            for c in range(0, len(WL2list), 1):
                if 'Int' in WL2list[c]:
                    scan2WL[count] = float( WL2list[c].replace('Int[','').replace(']','') )
                    #print(count, scan2WL[count], float( WL2list[c].replace('Int[','').replace(']','') ) )
                    count = count + 1           

            # Where is the measured data located? 
            # Extract the Scan1 WL and the measured data, more efficient to extract both data sets simultaneously
            nWL1 = fl_pre_amble['NoScanPoints(1)'] # number of WL for scan1 direction
            meas_r_start = fl_smpl_locs[smpl_name] + 6 # row number where data starts
            meas_r_stop = meas_r_start + nWL1 - 1 # row number where data stops

            scan1WL = np.zeros([ nWL1 ]) # declarate memory space for numpy array

            measData = np.zeros([nWL1, nWL2]) # declarate memory space for the numpy array

            if loud:
                print("\nData First:")
                print("r:", meas_r_start, fl_data.loc[ meas_r_start ][0].split('\t'))
                print("\nData Last:")
                print("r:", meas_r_stop, fl_data.loc[ meas_r_stop ][0].split('\t'))

            # loop over the rows of the measured data
            r_count = 0
            for r in range(meas_r_start, meas_r_stop+1, 1):
                fl_row = fl_data.loc[ r ][0].split('\t')
                c_count = 0
                # store the data from the columns of the file
                for c in range(0, len(fl_row)-1, 1):
                    if c == 0:
                        # save the scan1 WL value at each step
                        scan1WL[r_count] = float( fl_row[c] )
                    else:
                        # save the measured data for each of the scan2 WL 
                        measData[r_count][c_count] = float(fl_row[c])
                        c_count = c_count + 1        
                r_count = r_count + 1

            #print(len(scan1WL)); print(len(scan2WL)); 
            #print(len(measData)); print(len(measData[0]));
            #print(measData[0]); print(measData[-1]);

            print(smpl_name,"data extracted")

            return [scan1WL, scan2WL, measData]
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT); 
        print(e);

def Plot_All_Spectra(scanWL, excitationWL, measData, fig_name, loud = False):
    # Make a 2D plot of all the measured spectra
    # all spectra are measured across the WL set scanWL
    # each column of measData corresponds to a spectrum measured at a value excitationWL
    # size(measData) = len(scanWL)*len(excitationWL)
    # This function should be generic enough to plot measData even after its contents have been manipulated
    # R. Sheehan 28 - 3 - 2019

    FUNC_NAME = ".Plot_All_Spectra()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        if loud:
            c1 = scanWL is not None
            c2 = excitationWL is not None
            c3 = measData is not None
            c4 = len(measData) == len(scanWL) if True else False
            c5 = len(measData[0]) == len(excitationWL) if True else False
            c10 = c1 and c2 and c3 and c4 and c5

            if c10:
                plt_data = []; labels = []; marks = []

                for c in range(0, len(excitationWL), 2):
                    plt_data.append([ scanWL, measData[:,c] ] ); 
                    marks.append(Plotting.labs_lins[ c%len(Plotting.labs_lins) ]); 
                    labels.append('$\lambda_{Ex}$ = %(v1)0.1f nm'%{"v1":excitationWL[c]} ); 

                args = Plotting.plot_arg_multiple()

                args.loud = loud
                args.crv_lab_list = labels
                args.mrk_list = marks

                Plotting.plot_multiple_curves(plt_data, args)

                del plt_data; del args; 
            else:
                raise Exception
    except Exception as e:
        print(ERR_STATEMENT); 
        print(e);

def Plot_Spectrum(scanWL, excitationWL, measData, spctrm_choice, derivative_order, loud = False):
    # Make a plot of an individual spectrum with Savitky-Golay approximation
    # scanWL is the array containing all scan wavelengths
    # excitationWL is the array containing all excitation wavelengths
    # measData is the 2D array containing the measured spectra for each excitation wavelength
    # each column of measData corresponds to a measured spectrum at an excitation wavelength
    # spctrm_choice selects the spectrum to be approximated by the SG filter
    # derivative_order specifies what order of derivative to compute
    # R. Sheehan 29 - 3 - 2019

    FUNC_NAME = ".Plot_Spectrum()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        c1 = scanWL is not None
        c2 = excitationWL is not None
        c3 = measData is not None
        c4 = len(measData) == len(scanWL) if True else False
        c5 = len(measData[0]) == len(excitationWL) if True else False
        c6 = True if spctrm_choice > -1 and spctrm_choice < len(excitationWL) else False
        c10 = c1 and c2 and c3 and c4 and c5 and c6
        
        if c10:
            sgappr = Spectrum_Approximation(measData, spctrm_choice, derivative_order)

            # make a plot of the spectrum with its SG approximation
            args = Plotting.plot_arg_multiple()

            args.loud = loud
            args.curve_label = '$\lambda_{Ex}$ = %(v1)0.1f nm'%{"v1":excitationWL[spctrm_choice]}
            args.marker = 'k:'
            args.crv_lab_list = ['$\lambda_{Ex}$ = %(v1)0.1f nm'%{"v1":excitationWL[spctrm_choice]}, 'Savitzky-Golay']
            args.mrk_list = ['r+', 'b-']
            args.x_label = 'Scan Wavelength (nm)'
            args.y_label = 'Intensity (a. u.)'
            args.fig_name = 'Smpl_Spctrm'

            Plotting.plot_multiple_curves([[scanWL, measData[:,spctrm_choice]],[scanWL, sgappr]], args)

            del args; del sgappr; 
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT); 
        print(e);

def Plot_Contours(scanWL, excitationWL, measData, fig_name, loud = False):
    # make a contour plot of the data
    # R. Sheehan 5 - 4 - 2019

    FUNC_NAME = ".Plot_Contours()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME
    
    try:
        c1 = scanWL is not None
        c2 = excitationWL is not None
        c3 = measData is not None
        c4 = len(measData) == len(scanWL) if True else False
        c5 = len(measData[0]) == len(excitationWL) if True else False
        c10 = c1 and c2 and c3 and c4 and c5

        if c10:
            from mpl_toolkits.mplot3d import Axes3D
            from matplotlib import cm
            from matplotlib.ticker import LinearLocator, FormatStrFormatter

            # Compute the cartesian product of the two position sets
            XX, YY = np.meshgrid(scanWL, excitationWL)

            # generate the levels for the contour plot
            nlevels = 50
            vmax = measData.max()
            vmin = 10
            delta = (vmax - vmin) / (nlevels - 1)
            levels = np.arange(vmin, vmax, delta)
            #indx = 0
            ## exclud the zero contour
            #for i in range(0, len(levels), 1):
            #    if abs(levels[i]) < 1e-4 or abs(levels[i]) == 0.0:
            #        indx = i
            #if indx > 0:
            #    levels = np.delete(levels,indx)

            # make the contour plot
            #fig = plt.figure()
            #CS = plt.contour(XX, YY, np.transpose(measData), levels, origin = 'upper', cmap=cm.Reds, linewidths = 2.0)
            #CS = plt.contour(XX, YY, measData)
            #plt.clabel(CS, fontsize=9, inline=1) # add labels to the contours
            #fig.colorbar(CS, shrink=0.5, aspect=5)

            # another method
            fig1, ax = plt.subplots(constrained_layout=True)
            CS = ax.contourf(XX, YY, np.transpose(measData), nlevels, cmap=plt.cm.rainbow, origin='lower')

            # color maps
            # https://matplotlib.org/tutorials/colors/colormaps.html
            
            ax.set_title(fig_name)
            ax.set_xlabel('Scan $\lambda$ (nm)')
            ax.set_ylabel('Excitation $\lambda$ (nm)')        
            
            cbar = fig1.colorbar(CS)
            cbar.ax.set_ylabel('Intensity (a.u.)')
            #cbar.add_lines(CS)

            plt.savefig(fig_name)
            if loud: plt.show()

        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT); 
        print(e);

def Plot_Surface(scanWL, excitationWL, measData, fig_name, loud = False):
    # Make a contour plot of the measData

    FUNC_NAME = ".Plot_Surface()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        c1 = scanWL is not None
        c2 = excitationWL is not None
        c3 = measData is not None
        c4 = len(measData) == len(scanWL) if True else False
        c5 = len(measData[0]) == len(excitationWL) if True else False
        c10 = c1 and c2 and c3 and c4 and c5

        if c10:
            from mpl_toolkits.mplot3d import Axes3D
            from matplotlib import cm
            from matplotlib.ticker import LinearLocator, FormatStrFormatter

            # Compute the cartesian product of the two position sets
            XX, YY = np.meshgrid(scanWL, excitationWL)

            fig = plt.figure()
            ax = fig.gca(projection='3d')
            surf = ax.plot_surface(XX, YY, np.transpose(measData), cmap=cm.coolwarm, linewidth=0, antialiased=False)

            ax.set_xlabel('Excitation $\lambda$ (nm)')
            ax.set_ylabel('Scan $\lambda$ (nm)')
            ax.set_zlabel('Intensity (a.u.)')

            #ax.set_ylim(ydata[0], ydata[len(ydata)-1])
            #ax.set_xlim(xdata[0], xdata[len(xdata)-1])
            #ax.set_zlim(zrange[0], zrange[1])

            # Customize the z axis.
            #ax.set_zlim(-1.01, 1.01)
            ax.zaxis.set_major_locator(LinearLocator(10))
            ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

            # Add a color bar which maps values to colors.
            fig.colorbar(surf, shrink=0.5, aspect=5)

            plt.show()
        else:
            raise Exception
    except Exception:
        print(ERR_STATEMENT)

def Compute_Spectral_Difference(Data_Set_1, Data_Set_2):
    # Compute the difference between two sets of spectral data
    # R. Sheehan 28 - 3 - 2019

    FUNC_NAME = ".Compute_Spectral_Difference()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        # Check that the arrays are the same size and have been measured over the same wavelength ranges
        c1 = len(Data_Set_1[0]) == len(Data_Set_2[0]) if True else False # scanWL have same size
        c2 = len(Data_Set_1[1]) == len(Data_Set_2[1]) if True else False # exWL have same size
        c3 = len(Data_Set_1[2]) == len(Data_Set_2[2]) if True else False # measured Data have same size
        c4 = len(Data_Set_1[2][0]) == len(Data_Set_2[2][0]) if True else False # measured Data have same size
        c5 = Data_Set_1[0][0] == Data_Set_2[0][0] if True else False # scanWL start pos are the same
        c6 = Data_Set_1[0][-1] == Data_Set_2[0][-1] if True else False # scanWL end pos are the same
        c7 = Data_Set_1[1][0] == Data_Set_2[1][0] if True else False # scanWL end pos are the same
        c8 = Data_Set_1[1][-1] == Data_Set_2[1][-1] if True else False # scanWL end pos are the same
        
        c10 = c1 and c2 and c3 and c4 and c5 and c6 and c7 and c8

        if c10:
            # subtract spectrum 2 from spectrum 1
            spctr_Diff = np.zeros([len(Data_Set_1[0]), len(Data_Set_1[1])])
            for i in range(0, len(Data_Set_1[0]), 1):
                for j in range(0, len(Data_Set_1[1]), 1):
                    diff = Data_Set_1[2][i][j] - Data_Set_2[2][i][j]
                    spctr_Diff[i][j] = diff if diff > 0.0 else 0.0 # subtract spectrum 2 from spectrum 1
            return spctr_Diff
        else:
            raise Exception
    except Exception:
        print(ERR_STATEMENT)

def Spectrum_Approximation(measData, spctrm_choice, derivative_order, loud = False):
    # compute a Savitzky-Golay approximation to a measured-spectrum
    # scanWL is the array containing all scan wavelengths
    # excitationWL is the array containing all excitation wavelengths
    # measData is the 2D array containing the measured spectra for each excitation wavelength
    # each column of measData corresponds to a measured spectrum at an excitation wavelength
    # spctrm_choice selects the spectrum to be approximated by the SG filter
    # derivative_order specifies what order of derivative to compute
    # derivative_order = 0 => spectrum approximation
    # derivative_order = 1 => spectrum first derivative approximation
    # derivative_order = 2 => spectrum second derivative approximation
    # R. Sheehan 29 - 3 - 2019

    FUNC_NAME = ".Spectrum_Approximation()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        c3 = measData is not None
        c6 = True if spctrm_choice > -1 and spctrm_choice < len(measData[0]) else False
        c7 = True if derivative_order > -1 and derivative_order < 3 else False
        c10 = c3 and c6 and c7

        if c10:
            from scipy.signal import savgol_filter

            # try a Savitzky-Golay filter
            # window_length is the length of the filter window, i.e. no. data used to construct filter
            # polyorder is the order of the polynomial used to do the filtering
            # deriv is order of derivative to be computed
            # this actually works quite well, it gives a really good approximation to all spectra that I've seen
            sgappr = savgol_filter(measData[:,spctrm_choice], window_length = 5, polyorder = 2, deriv = derivative_order)

            return sgappr
        else:
            raise Exception
    except Exception:
        print(ERR_STATEMENT)

def Locate_Peaks(scanWL, excitationWL, measData, spctrm_choice, loud = False):
    # use in-built Python code to find all peaks in the spectra
    # filter out the scattering peaks based on their FWHM
    # use https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
    # scanWL is the array containing all scan wavelengths
    # excitationWL is the array containing all excitation wavelengths
    # measData is the 2D array containing the measured spectra for each excitation wavelength
    # each column of measData corresponds to a measured spectrum at an excitation wavelength
    # spctrm_choice selects the spectrum to be approximated by the SG filter

    FUNC_NAME = ".Locate_Peaks()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        c1 = scanWL is not None
        c2 = excitationWL is not None
        c3 = measData is not None
        c4 = len(measData) == len(scanWL) if True else False
        c5 = len(measData[0]) == len(excitationWL) if True else False
        c6 = True if spctrm_choice > -1 and spctrm_choice < len(excitationWL) else False
        c10 = c1 and c2 and c3 and c4 and c5 and c6

        # location of scattering resonance seems strongly correlated with value of excitation wavelength
        # check if this is really the case

        if c10:
            from scipy.signal import find_peaks, peak_prominences, peak_widths
            #from scipy.signal import argrelextrema # can be used to locate maxima / minima
            
            # find the peaks in the signal
            peaks, heights = find_peaks(measData[:,spctrm_choice], height = 15)
            prominences = peak_prominences(measData[:,spctrm_choice], peaks)[0] # compute peak prominences
            widths = peak_widths(measData[:,spctrm_choice], peaks, rel_height=0.5)[0] # compute peak FWHM

            fl_dict = {}
            res_dict = {}
            keys = ['indices', 'peak wl', 'peak height', 'peak width']
            fl_indx = []; res_indx = []; 
            fl_WL = []; res_WL = []; 
            fl_height = []; res_height = [];
            fl_width = []; res_width = [];
            for i in range(0, len(peaks), 1):
                if prominences[i] > 12:
                    if widths[i] > 50:
                        if loud: print("Fluoresence Peak:",scanWL[ peaks[i] ], heights['peak_heights'][i], prominences[i], widths[i])
                        fl_indx.append(peaks[i]); fl_WL.append(scanWL[ peaks[i] ]); 
                        fl_height.append(heights['peak_heights'][i]); fl_width.append(widths[i]); 
                    else:
                        if loud: print("Scattering Peak:",scanWL[ peaks[i] ], heights['peak_heights'][i], prominences[i], widths[i])
                        res_indx.append(peaks[i]); res_WL.append(scanWL[ peaks[i] ]); 
                        res_height.append(heights['peak_heights'][i]); res_width.append(widths[i]);

            # make dictionary for fl peak data
            fl_dict['indices'] = fl_indx; fl_dict['peak wl'] = fl_WL; 
            fl_dict['peak height'] = fl_height; fl_dict['peak width'] = fl_width; 

            # make dictionary for resonance peak data
            res_dict['indices'] = res_indx; res_dict['peak wl'] = res_WL; 
            res_dict['peak height'] = res_height; res_dict['peak width'] = res_width;

            del fl_indx; del fl_WL; del fl_height; del fl_width; 
            del res_indx; del res_WL; del res_height; del res_width; 
            
            if loud: Plot_Spectrum(scanWL, excitationWL, measData, spctrm_choice, 0, loud)

            return [fl_dict, res_dict]
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT); 
        print(e);

def Extract_Peak_Data(scanWL, excitationWL, measData, loud = False):
    # Extract all resonance peak data for the data set
    # scan each column of the measData set
    # locate all flouresence and scattering peaks in each data
    # R. Sheehan 29 - 3 - 2018

    FUNC_NAME = ".Extract_Peak_Data()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        c1 = scanWL is not None
        c2 = excitationWL is not None
        c3 = measData is not None
        c4 = len(measData) == len(scanWL) if True else False
        c5 = len(measData[0]) == len(excitationWL) if True else False
        
        c10 = c1 and c2 and c3 and c4 and c5

        if c10:

            peak_WL = np.zeros([len(excitationWL)]); # WL at which Fluoresence peak occurs
            peak_H = np.zeros([len(excitationWL)]); # Height of FL peak
            peak_W = np.zeros([len(excitationWL)]); # Width of FL peak
            
            non_zero_count = 0; 
            for c in range(0, len(excitationWL), 1):
                FL_pk, Res_pk = Locate_Peaks(scanWL, excitationWL, measData, c)
                #print(FL_pk)
                if len(FL_pk['peak wl']) > 0:
                    peak_WL[c] = FL_pk['peak wl'][0]; 
                    peak_H[c] = FL_pk['peak height'][0]; 
                    peak_W[c] = FL_pk['peak width'][0]; 
                    non_zero_count = non_zero_count + 1
                del FL_pk; del Res_pk; 

            if loud:
                # make a plot of the data
                args = Plotting.plot_arg_single(); 

                args.loud = loud

                args.x_label = "Excitation Wavelength (nm)"

                args.y_label = "Fluoresence Peak (nm)"
                Plotting.plot_single_curve(excitationWL[0:non_zero_count], peak_WL[0:non_zero_count], args);

                args.y_label = "Fluoresence Height (a. u.)"
                Plotting.plot_single_curve(excitationWL[0:non_zero_count], peak_H[0:non_zero_count], args);

                args.y_label = "Fluoresence Width (nm)"
                Plotting.plot_single_curve(excitationWL[0:non_zero_count], peak_W[0:non_zero_count], args); 

                # Would this data would make more sense on a plot with two scales?
                # Not sure the data needs to be plotted like this
                # should add a version of the two scale plot to Plotting module though. 
                #fig, ax1 = plt.subplots()
                #color = 'tab:red'
                #ax1.set_xlabel('Excitation Wavelength (nm)')
                #ax1.set_ylabel('Fluoresence Peak (nm)', color = color)
                #ax1.plot(excitationWL[0:non_zero_count], peak_WL[0:non_zero_count], color = color)
                ##ax1.set_ylabel('Fluoresence Width (nm)', color = color)
                ##ax1.plot(excitationWL[0:non_zero_count], peak_W[0:non_zero_count], color = color)
                #ax1.tick_params(axis='y', labelcolor = color)

                #ax2 = ax1.twinx()
                #color = 'tab:blue'
                #ax2.set_ylabel('Fluoresence Height (a. u.)', color = color)
                #ax2.plot(excitationWL[0:non_zero_count], peak_H[0:non_zero_count], color = color)
                #ax2.tick_params(axis='y', labelcolor = color)

                #fig.tight_layout()
                #plt.show()
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT); 
        print(e);

def Remove_Scatt_Peaks_Spctrm(scanWL, excitationWL, measData, spctrm_choice, pk_dict, loud = False):
    # remove the scattering peaks from an individual spectrum
    # pk_dict contains all info pertaining to spectrum peak in a data set
    # it is a dictionary of the form 
    # pk_dict['indices'] = res_indx; pk_dict['peak wl'] = res_WL; 
    # pk_dict['peak height'] = res_height; pk_dict['peak width'] = res_width;
    # R. Sheehan 4 - 4 - 2019

    # It's possible to refine this process further
    # You could search

    FUNC_NAME = ".Remove_Scatt_Peaks_Spctrm()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        c1 = scanWL is not None
        c2 = excitationWL is not None
        c3 = measData is not None
        c4 = len(measData) == len(scanWL) if True else False
        c5 = len(measData[0]) == len(excitationWL) if True else False
        c6 = True if spctrm_choice > -1 and spctrm_choice < len(excitationWL) else False
        
        c10 = c1 and c2 and c3 and c4 and c5 and c6

        if c10:
            c7 = True if pk_dict is not None else False
            c8 = True if len(pk_dict['indices']) > 0 else False
            
            
            if c7 and c8:
                # there is peak data to be removed
                # if there is no peak data then do nothing
                starts = []; ends = []; 
                for v in range(0, len( pk_dict['indices'] ), 1):
                    # only remove a resonance peak if it far from fluoresence
                    # i.e. only remove the first order resonance peaks
                    c9 = False if pk_dict['peak wl'][v] - excitationWL[spctrm_choice] > 20 else True
                    if c9:
                        pk_strt = -1 + pk_dict['peak wl'][v]-0.5*pk_dict['peak width'][v]
                        pk_end = 1 + pk_dict['peak wl'][v]+0.5*pk_dict['peak width'][v]

                        # find indices of values closest to Peak start and Peak end
                        pk_strt_indx = np.abs(scanWL - pk_strt).argmin()
                        pk_end_indx = np.abs(scanWL - pk_end).argmin()

                        starts.append(pk_strt_indx); ends.append(pk_end_indx)

                        if loud:
                            print('Peak location:',pk_dict['indices'][v])
                            print('Peak width:',pk_dict['peak width'][v])                
                            print('Peak start:',pk_strt,', closest:',scanWL[pk_strt_indx])
                            print('Peak wavelength:',pk_dict['peak wl'][v])                
                            print('Peak end:',pk_end,', closest:',scanWL[pk_end_indx])
                            print('')

                        for i in range(pk_strt_indx, pk_end_indx, 1):
                            if loud: print(scanWL[i],',',measData[:,spctrm_choice][i])
                            measData[:,spctrm_choice][i] = 0.0

                # re-assign the value of the peak to be zero
                #for j in range(0, len(starts), 1):
                #    for i in range(starts[j], ends[j], 1):
                #        #if loud: print(scanWL[i],',',measData[:,spctrm_choice][i])
                #        measData[:,spctrm_choice][i] = 0.0

                # make a plot for comparison
                if loud: Plot_Spectrum(scanWL, excitationWL, measData, spctrm_choice, 0, loud)
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT); 
        print(e);

def Remove_Scatt_Peaks_Second(scanWL, excitationWL, measData, spctrm_choice, fl_pk_dict, res_pk_dict, loud = False):
    # remove the second order scattering peaks from an individual spectrum
    # x_pk_dict contains all info pertaining to spectrum peak in a data set
    # it is a dictionary of the form 
    # x_pk_dict['indices'] = res_indx; x_pk_dict['peak wl'] = res_WL; 
    # x_pk_dict['peak height'] = res_height; x_pk_dict['peak width'] = res_width;
    # the idea here is to remove the second order scattering effects from the data
    # if res_pk_dict['peak_height'] > fl_pk_dict['peak_height']
    # the res_pk_dict data will be filtered
    # R. Sheehan 30 - 7 - 2019

    FUNC_NAME = ".Remove_Scatt_Peaks_Second()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        c1 = scanWL is not None
        c2 = excitationWL is not None
        c3 = measData is not None
        c4 = len(measData) == len(scanWL) if True else False
        c5 = len(measData[0]) == len(excitationWL) if True else False
        c6 = True if spctrm_choice > -1 and spctrm_choice < len(excitationWL) else False
        
        c10 = c1 and c2 and c3 and c4 and c5 and c6

        if c10:
            c5 = True if fl_pk_dict is not None else False
            c6 = True if len(fl_pk_dict['indices']) > 0 else False
            c7 = True if res_pk_dict is not None else False
            c8 = True if len(res_pk_dict['indices']) > 0 else False 

            c10 = c5 and c6 and c7 and c8
            
            if c10:
                # there is peak data to be removed
                # if there is no peak data then do nothing
                starts = []; ends = []; 
                for v in range(0, len( res_pk_dict['indices'] ), 1):
                    # only remove a resonance peak if it far from fluoresence
                    # i.e. only remove the first order resonance peaks
                    c8 = True if res_pk_dict['peak wl'][v] - excitationWL[spctrm_choice] > 20 else False
                    c9 = True if res_pk_dict['peak height'][v] - fl_pk_dict['peak height'][0] > 20 else False
                    if c8 and c9:
                        pk_strt = -1 + res_pk_dict['peak wl'][v]-0.5*res_pk_dict['peak width'][v]
                        pk_end = 1 + res_pk_dict['peak wl'][v]+0.5*res_pk_dict['peak width'][v]

                        # find indices of values closest to Peak start and Peak end
                        pk_strt_indx = np.abs(scanWL - pk_strt).argmin()
                        pk_end_indx = np.abs(scanWL - pk_end).argmin()

                        starts.append(pk_strt_indx); ends.append(pk_end_indx)

                        if loud:
                            print('Peak location:',res_pk_dict['indices'][v])
                            print('Peak width:',res_pk_dict['peak width'][v])                
                            print('Peak start:',pk_strt,', closest:',scanWL[pk_strt_indx])
                            print('Peak wavelength:',res_pk_dict['peak wl'][v])                
                            print('Peak end:',pk_end,', closest:',scanWL[pk_end_indx])
                            print('')

                        for i in range(pk_strt_indx, pk_end_indx, 1):
                            if loud: print(scanWL[i],',',measData[:,spctrm_choice][i])
                            measData[:,spctrm_choice][i] = 0.0

                # re-assign the value of the peak to be zero
                #for j in range(0, len(starts), 1):
                #    for i in range(starts[j], ends[j], 1):
                #        #if loud: print(scanWL[i],',',measData[:,spctrm_choice][i])
                #        measData[:,spctrm_choice][i] = 0.0

                # make a plot for comparison
                if loud: Plot_Spectrum(scanWL, excitationWL, measData, spctrm_choice, 0, loud)
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT); 
        print(e);

def Remove_Scatt_Peaks(scanWL, excitationWL, measData, loud = False):
    # locate the scattering peaks in each data set
    # remove the scattering peaks from each column of the data set
    # R. Sheehan 4 - 4 -2019

    FUNC_NAME = ".Remove_Scatt_Peaks()" # use this in exception handling messages
    ERR_STATEMENT = "Error: " + MOD_NAME_STR + FUNC_NAME

    try:
        c1 = scanWL is not None
        c2 = excitationWL is not None
        c3 = measData is not None
        c4 = len(measData) == len(scanWL) if True else False
        c5 = len(measData[0]) == len(excitationWL) if True else False
        
        c10 = c1 and c2 and c3 and c4 and c5

        if c10:
            
            non_zero_count = 0; 
            for c in range(0, len(excitationWL), 1):

                # only look for resonances in region where \lambda_{ex} > \lambda_{scan}
                if excitationWL[c] > scanWL[c]:
                    # locate the fluoresence and scattering peaks in a single spectrum
                    FL_pk, Res_pk = Locate_Peaks(scanWL, excitationWL, measData, c, loud)

                    if Res_pk['peak wl'][0] > excitationWL[c]:
                        print("Excitation WL:",excitationWL[c])
                        print("Fluor WL:",FL_pk['peak wl'],"Resonance H:",FL_pk['peak height'])
                        print("Resonance WL:",Res_pk['peak wl'],"Resonance H:",Res_pk['peak height'])
                        #for i in range(0, len(Res_pk['peak wl']), 1):
                        #    print("Resonance WL - Excitation WL:",Res_pk['peak wl'][i] - excitationWL[c])
                        print("")

                    #if len(Res_pk['peak wl']) > 1:
                    #    print("Excitation WL:",excitationWL[c],"Resonance WL:",Res_pk['peak wl'])
                
                    # remove the scattering peaks from the spectral data
                    Remove_Scatt_Peaks_Spctrm(scanWL, excitationWL, measData, c, Res_pk, loud)

                    # remove the second order scattering peaks
                    Remove_Scatt_Peaks_Second(scanWL, excitationWL, measData, c, FL_pk, Res_pk, loud)

                    del FL_pk; del Res_pk;
        else:
            raise Exception
    except Exception as e:
        print(ERR_STATEMENT); 
        print(e);