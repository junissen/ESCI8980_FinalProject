#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Input functions for Age Calculation, in order to calculate U-Th data using MC-ICP-MS on SEM measurements. 

Created by Julia Nissen, 2017, for use in Edwards' Isotope Lab.

"""

import isofilter
import numpy as np
import pandas as pd

    
class Ucalculation():
    """
    Class Ucalculation functions as the U sheet in the age calculation spreadsheet. Ucalculation gives outputs for 
    both the Thcalculation function and the Agecalculation function. 
    
    U_normalized_forTh output is a list of the following values: 
        [0]: 236/233 measured ratio
        [1]: 236/233 measured ratio error
        [2]: 235/233 normalized ratio
        [3]: 235/233 normalized ratio error
        [4]: 236/233 corrected ratio
        [5]: 236/233 corrected ratio error
    
    U_normalized_forAge output is a list of the following values: 
        [0]: 235/233 normalized ratio
        [1]: 235/233 normalized ratio error
        [2]: 235/234 normalized and corrected ratio
        [3]: 235/234 normalized and corrected ratio error
        [4]: Unfiltered 233 counts
        [5]: Filtered 234/235 counts
        [6]: Unfiltered 233 mean
    
    """
    def __init__ (self, spike_input, AS_input, filename_input, inquiry_input):
        
        spike = str(spike_input)
    
        spike_six_three_dictionary = {"DIII-B":1.008398,"DIII-A": 1.008398,"1I":1.010128,"1H":1.010128}
        
        
        #derives 236/233 value of spike from preset dictionary
        if spike in spike_six_three_dictionary:
            self.spike = float(spike_six_three_dictionary[spike])
        else: 
            print 'ERROR: You did not enter a valid spike option'
    
        
        #allows you the ability to print as you go
        inquiry = str(inquiry_input)
            
        if inquiry.lower() == "y":
            self.inquiry = True
        else:
            self.inquiry = False
            print "Program Ucalculation ran without printing"
        
        #AS is the abundant sensitivity 237/238, measured through the AS method on the ICP-MS    
        self.AS = float(AS_input)
        
        #uses the filename given for your U run
        filename = str(filename_input)
                            
        #236/233 filtered measured mean and 2s error
        working = isofilter.IsoFilter(filename,"G", 44)
        a = working.getMean()
        b = working.getStanddev()
        c = working.getCounts()
        self.six_three_mean_meas = working.Filtered_mean(a,b,c)
        self.six_three_err_meas = working.Filtered_err(a,b,c)
        
        #235/233 filtered measured mean and 2s error
        working_b = isofilter.IsoFilter(filename, "H", 44)
        a = working_b.getMean()
        b = working_b.getStanddev()
        c = working_b.getCounts()
        self.five_three_mean_meas = working_b.Filtered_mean(a,b,c)
        self.five_three_err_meas = working_b.Filtered_err(a,b,c) 
    
        #234/235 filtered measured mean and 2s error
        working_c = isofilter.IsoFilter(filename,"I", 44)
        a = working_c.getMean()
        b = working_c.getStanddev()
        c = working_c.getCounts()
        self.four_five_mean_meas = working_c.Filtered_mean(a,b,c)
        self.four_five_err_meas = working_c.Filtered_err(a,b,c) 
        self.four_five_counts = working_c.Filtered_counts(a,b,c)
        
        #233 unfiltered mean and counts
        working_d = isofilter.IsoFilter(filename, "C", 44)
        self.three_mean_meas = working_d.getMean()
        self.three_counts = working_d.getCounts()
        
        #constants to be used throughout the class
        self.wt_235 = 235.043924
        self.wt_233 = 233.039629
        self.wt_236 = 236.045563
        self.wt_234 = 234.040947
        self.eight_five_rat = 137.83
        self.AS_six_eight = self.AS/5
        self.AS_four_eight = self.AS/20
        self.eight_five_rat_err_rel = 0.0003
        
        
    def U_normalization_forTh(self):
        """
        Function outputs the measured 236/233 ratio and error, the 236/233 ratio and error corrected for the 238 tail, 
        and the 235/233 normalized ratio and error using the 235/233 corrected ratio and further correcting
        235/233 for mass fractionation in the ICP-MS. These values are used later in the Th_normalization function.
        """
        
        #corrects 236/233 ratio for 238 tail 
        self.six_three_corr = self.six_three_mean_meas * ( 1 - (self.AS_six_eight * self.five_three_mean_meas 
                                                                * self.eight_five_rat/self.spike) )
        
        #provides the ratio that will be used to correct for mass fractionation
        rat = float(np.log(self.wt_235/self.wt_233)/np.log(self.wt_236/self.wt_233))
        
        #corrects for mass fractionation in the ICP-MS
        self.five_three_norm = self.five_three_mean_meas * (self.spike/self.six_three_corr)**rat
        
        #provides relative error constants to be used in this function
        AS_six_eight_err_rel = 0.3
        five_three_err_rel = self.five_three_err_meas/self.five_three_mean_meas
        six_three_err_rel = self.six_three_err_meas/self.six_three_mean_meas
       
        #calculculates the 236/233 corrected error
        self.six_three_corr_err = self.six_three_corr * np.sqrt( six_three_err_rel**2 + 
                                                                ( (self.AS_six_eight * self.five_three_mean_meas 
                                                                   * self.eight_five_rat)/self.spike  
                                                                   * np.sqrt( AS_six_eight_err_rel**2 + five_three_err_rel ** 2 + self.eight_five_rat_err_rel**2 ) 
                                                                   / (1 - (self.AS_six_eight * self.five_three_mean_meas * self.eight_five_rat)
                                                                   / self.spike) ) ** 2 ) 
        #calculates the 236/233 relative corrected error
        self.six_three_corr_err_rel = self.six_three_corr_err/self.six_three_corr
        
        #calculates the 235/233 normalized error
        self.five_three_norm_err = self.five_three_norm * np.sqrt( five_three_err_rel**2 
                                                                  + (2 * (self.six_three_corr_err_rel/3))**2  ) 
       
        #if you have chosen to print as you go, this will print when the function is finished
        while self.inquiry:
            print "RESULTS FOR TH FILTERING FROM U:"
            print "236/233 measured ratio: " + str(self.six_three_mean_meas) + " ± " + str(self.six_three_err_meas)
            print "235/233 normalized ratio: " + str(self.five_three_norm) + " ± " + str(self.five_three_norm_err)
            print "236/233 corrected ratio: " + str(self.six_three_corr) + " ± " + str(self.six_three_corr_err)
            break
        
        #a list of your outputs is created and returned, to be used in the Th functions
        lstU_Th = [self.six_three_mean_meas, self.six_three_err_meas, self.five_three_norm, 
                 self.five_three_norm_err, self.six_three_corr,self.six_three_corr_err ]
        
        return lstU_Th
    
    def U_normalization_forAge(self):
        """
        Function outputs the 235/233 normalized ratio and error, the 235/234 normalized and corrected ratio and error, 
        the unfiltered number of cycles for 233 and mean value and the filtered number of cycles 234/235. These values will be used
        later in the Age Calculation.
        """
        #calculates constants that will be used to calculate normalized 234/235
        four_five_err_rel = self.four_five_err_meas / self.four_five_mean_meas
        
        rat = float(np.log(self.wt_234/self.wt_235)/np.log(self.wt_236/self.wt_233))
        
        #normalizes the 234/235 ratio by correcting for mass fractionation and calculates the resulting error
        self.four_five_norm = self.four_five_mean_meas * (self.spike/self.six_three_corr)**rat
        
        self.four_five_norm_err = self.four_five_norm * np.sqrt( four_five_err_rel**2 + 
                                                                (self.six_three_corr_err_rel/3)**2 )
        
        #calculates constants that will be used to calculated corrected 234/235
        AS_four_eight_err_rel = 0.3
        four_five_norm_err_rel = self.four_five_norm_err/self.four_five_norm
        
        #corrects the normalized 234/235 ratio for 238 tail and calculated the resulting error
        self.four_five_normcorr = self.four_five_norm * (1 - ( self.eight_five_rat 
                                                              * self.AS_four_eight/ self.four_five_norm ))

        self.four_five_normcorr_err = self.four_five_normcorr * np.sqrt( four_five_norm_err_rel**2 + 
                                                                        ( (self.eight_five_rat * self.AS_four_eight / self.four_five_norm) *
                                                                         np.sqrt( self.eight_five_rat_err_rel**2 + AS_four_eight_err_rel**2 + four_five_norm_err_rel**2 )
                                                                         / (1 - ( self.eight_five_rat * self.AS_four_eight/ self.four_five_norm)) ) **2 ) 
        
        self.three_mean_meas = int(self.three_mean_meas)
        
        #if you have chosen to print as you go, this will print when the function is finished
        while self.inquiry:
            print "RESULTS FOR AGE CALC FROM U:"
            print "235/233 normalized ratio : " + str(self.five_three_norm) + " ± " + str(self.five_three_norm_err)
            print "234/235 normalized and corrected ratio: " + str(self.four_five_normcorr) + " ± " + str(self.four_five_normcorr_err)
            print "Unfiltered cycles of 233: " + str(self.three_counts)
            print "Filtered cycles of 234/235: " + str(self.four_five_counts)
            print "Unfiltered mean 233U cps : " + str(self.three_mean_meas)
            break
        
        #a list of your outputs is created and returned, to be used in the Age functions
        lstU_Age = [self.five_three_norm, self.five_three_norm_err, self.four_five_normcorr, 
                    self.four_five_normcorr_err, self.three_counts, self.four_five_counts, self.three_mean_meas] 
       
        return lstU_Age
    
class Thcalculation():
    """
    Class Thcalculation functions as the Th sheet in the age calculation spreadsheet. Thcalculation gives outputs for 
    the Agecalculation function, and needs to be provided inputs from the Ucalculation class U_normalization_forTh function.
    
    Th_normalization_forAge output is a list of the following values: 
        [0]: 230/229 corrected and normalized ratio
        [1]: 230/229 corrected and normalized ratio error
        [2]: 232/229 corrected and normalized ratio
        [3]: 232/229 corrected and normalized ratio error
        [4]: Unfiltered 229 mean
        [5]: Unfiltered 229 counts
        
    """
    
    def __init__ (self, spike_input, AS_input, filename_input, inquiry_input, lstU_Th):
        
        spike = str(spike_input)
    
        spike_six_three_dictionary = {"DIII-B":1.008398,"DIII-A": 1.008398,"1I":1.010128,"1H":1.010128}
        
        
        #derives 236/233 value of spike from preset dictionary
        if spike in spike_six_three_dictionary:
            self.spike = float(spike_six_three_dictionary[spike])
        else: 
            print 'ERROR: You did not enter a valid spike option'
    
        
        #allows you the ability to print as you go
        inquiry = str(inquiry_input)
            
        if inquiry.lower() == "y":
            self.inquiry = True
        else:
            self.inquiry = False
            print "Program Thcalculation ran without printing"
            
        #AS is the abundant sensitivity 237/238, measured through the AS method on the ICP-MS    
        self.AS = float(AS_input)
        
        #uses the filename given for your Th run
        filename = str(filename_input)
        
        #Compiles the values of the lstU_Th provided by your U_normalization_forTh function
        self.six_three_mean_meas = lstU_Th[0]
        self.six_three_err_meas = lstU_Th[1]
        self.five_three_norm = lstU_Th[2]
        self.five_three_norm_err = lstU_Th[3]
        self.six_three_corr = lstU_Th[4]
        self.six_three_corr_err = lstU_Th[5]
        
        #Note: Hai's macro only filters 230/229 column
        
        #230/232 filtered measured mean and 2s error
        working = isofilter.IsoFilter(filename,"G", 28)
        self.zero_two_mean_meas = working.getMean()/1.02
        self.zero_two_counts = working.getCounts()
        self.zero_two_standdev_meas = working.getStanddev()
        self.zero_two_rel_err_meas = (2 * self.zero_two_standdev_meas/(self.zero_two_counts**0.5))/self.zero_two_mean_meas
        self.zero_two_rel_err = max(self.zero_two_rel_err_meas, 0.02)
        self.zero_two_err_meas = self.zero_two_mean_meas * self.zero_two_rel_err
        
        #230/229 filtered measured mean and 2s error
        working_b = isofilter.IsoFilter(filename, "E", 28)
        a = working_b.getMean()
        b = working_b.getStanddev()
        c = working_b.getCounts()
        self.zero_nine_mean_meas = working_b.Filtered_mean(a,b,c)
        self.zero_nine_err_meas = working_b.Filtered_err(a,b,c)
        
        #232/229 filtered measured mean and 2s error
        working_c = isofilter.IsoFilter(filename, "F", 28)
        self.nine_two_mean_meas = working_c.getMean()
        self.two_nine_mean_meas = 1 / (self.nine_two_mean_meas/1.02)
        self.two_nine_counts = working.getCounts()
        self.nine_two_standdev_meas = working_c.getStanddev()
        self.nine_two_rel_err_meas = (2 * self.nine_two_standdev_meas/(self.two_nine_counts**0.5))/self.nine_two_mean_meas
        self.two_nine_rel_err = max(self.nine_two_rel_err_meas, 0.02)
        self.two_nine_err_meas = self.two_nine_mean_meas * self.two_nine_rel_err
        
        #229 unfiltered mean and counts
        working_d = isofilter.IsoFilter(filename, "C", 28)
        self.nine_mean_meas = working_d.getMean()
        self.nine_counts = working_d.getCounts()
        
        #constants to be used throughout the class
        self.wt_233 = 233.039629
        self.wt_236 = 236.045563
        self.wt_229 = 229.031756
        self.wt_230 = 230.033128
        self.wt_232 = 232.038051
        self.AS_zero_nine = self.AS
        self.AS_zero_two = self.AS_zero_nine / 5
        self.AS_two_nine = self.AS_zero_two / 3
        self.eight_five_rat = 137.83
        self.eight_five_rat_err_rel = 0.0003
        
    def Th_normalization_forAge(self):
        
        #corrects the 230/229 and 232/229 ratios for both the 232 and 229 tails
        self.zero_nine_corr = self.zero_nine_mean_meas * (1 - self.AS_zero_two/self.zero_two_mean_meas) * (1 - self.AS_zero_nine)
        
        self.two_nine_corr = self.two_nine_mean_meas * (1 / (1 - (self.AS_two_nine * self.two_nine_mean_meas)))
        
        #constants needed for error calculations
        self.zero_nine_rel_err = self.zero_nine_err_meas/self.zero_nine_mean_meas
        self.AS_zero_two_rel_err = 0.3
        self.AS_zero_nine_rel_err = 0.3
        self.AS_two_nine_rel_err = 0.3
        
        #errors for corrected 230/229 and 232/229 ratios
        self.zero_nine_corr_err = self.zero_nine_corr * ( self.zero_nine_rel_err**2 + 
                                                         ( (self.AS_zero_two/self.zero_two_mean_meas) *  
                                                                (self.AS_zero_two_rel_err**2 + self.zero_two_rel_err**2)**0.5
                                                                 / (1 - self.AS_zero_two/self.zero_two_mean_meas))**2 
                                                         + ( self.AS_zero_nine * self.AS_zero_nine_rel_err/(1 - self.AS_zero_nine) )**2) ** 0.5
       
        self.two_nine_corr_err = self.two_nine_corr * ( self.two_nine_rel_err**2 + 
                                                      ( ((self.AS_two_nine_rel_err**2 + self.two_nine_rel_err**2)**0.5)  
                                                       * (self.AS_two_nine * self.two_nine_mean_meas)/
                                                         (1 - (self.AS_two_nine * self.two_nine_mean_meas))**2) ** 2
                                                       ) ** 0.5
        
        #constant needed for normalization
        
        rat_1 = np.log(self.wt_230/self.wt_229) / np.log(self.wt_236/self.wt_233)
        
        rat_2 = np.log(self.wt_232/self.wt_229) / np.log(self.wt_236/self.wt_233)
        

        #normalizes for corrected 230/229 and 232/229 ratios for mass fractionation
        
        self.zero_nine_corrnorm = self.zero_nine_corr * ((self.spike / self.six_three_corr)**rat_1)
        
        self.two_nine_corrnorm = self.two_nine_corr * ((self.spike / self.six_three_corr)**rat_2)
        
        #constants needed for error calculations
        self.zero_nine_corr_rel_err = self.zero_nine_corr_err / self.zero_nine_corr
        self.six_three_corr_rel_err = self.six_three_corr_err / self.six_three_corr
        self.two_nine_corr_rel_err = self.two_nine_corr_err / self.two_nine_corr
        
        #errors for normalized 230/229 and 232/229 ratios 
        self.zero_nine_corrnorm_err = self.zero_nine_corrnorm * ( self.zero_nine_corr_rel_err**2
                                                                 + ( (self.six_three_corr_rel_err/3)**2 )
                                                                 )**0.5
        
        self.two_nine_corrnorm_err = self.two_nine_corrnorm * ( self.two_nine_corr_rel_err**2
                                                               +  self.six_three_corr_rel_err**2  
                                                               )**0.5
        
        self.nine_mean_meas = int(self.nine_mean_meas)
        
        #if you have chosen to print as you go, this will print when the function is finished
        while self.inquiry:
            print "RESULTS FOR AGE CALC FROM TH:"
            print "230/229 corrected and normalized ratio : " + str(self.zero_nine_corrnorm)  + " ± " + str(self.zero_nine_corrnorm_err) 
            print "232/229 corrected and normalized ratio: " + str(self.two_nine_corrnorm) + " ± " + str(self.two_nine_corrnorm_err)
            print "Unfiltered mean 229Th cps: " + str(self.nine_mean_meas)    
            print "Unfiltered cycles of 229: " + str(self.nine_counts)
            break
        
        #a list of your outputs is created and returned, to be used in the Age functions
        lstTh_age = [self.zero_nine_corrnorm, self.zero_nine_corrnorm_err, self.two_nine_corrnorm,
                     self.two_nine_corrnorm_err, self.nine_mean_meas, self.nine_counts]
        
        return lstTh_age
    
    
class background_values():
    
    def __init__(self, U_file, Th_file, inquiry_input):
        
        inquiry = str(inquiry_input)
        
        #allows you the ability to print as you go
        if inquiry.lower() == "y":
            self.inquiry = True
        else:
            self.inquiry = False
            print "Program background_values ran without printing"
        
        #uses the filename given for your U wash
        self.filename_U = str(U_file)
        
        #uses the filename give for your Th wash
        self.filename_Th = str(Th_file)
        
    def U_wash(self):
        """
        U_wash provides a list the following outputs for the Age Calculation: 
            [0]: 233 unfiltered wash in cps
            [1]: 234 unfiltered wash in cps
            [2]: 235 unfiltered wash in cps
            
        """
        #233 wash value
        working_a = isofilter.IsoFilter(self.filename_U,"C", 44)
        self.three_wash = working_a.getMean()
        
        #234 wash value
        working_b = isofilter.IsoFilter(self.filename_U,"D", 44)
        self.four_wash = working_b.getMean()
        
        #235 wash value
        working_c = isofilter.IsoFilter(self.filename_U,"E", 44)
        self.five_wash = working_c.getMean()
        
        while self.inquiry: 
            print "BACKGROUND VALUES FOR AGE CALC:"
            print "233 wash value: " + str(self.three_wash) + ' cps'
            print "234 wash value: " + str(self.four_wash) + ' cps'
            print "235 wash value: " + str(self.five_wash) + ' cps'
            break
        
        lstU_wash = [self.three_wash, self.four_wash, self.five_wash]
        
        return lstU_wash
    
    def Th_wash(self):
        """
        Th_wash provides the following outputs for the Age Calculation: 
            230 unfiltered wash in cpm 
                
        """
        #230 wash value
        working_a = isofilter.IsoFilter(self.filename_Th, "D", 28)
        self.zero_wash = working_a.getMean()
        
        #calculate "darknoise" value for Age calculation
        self.darknoise = self.zero_wash * 60
        
        while self.inquiry:
            print "230 wash value/darknoise: " + str(self.darknoise) + ' cpm'
            break
            
        return self.darknoise

class chemblank_values():
    """
    Calculates chem blank values for use in Age Calculation spreadsheet
    
    lstChemblank provides a list of the following outputs: 
        [0]: 238 chemblank value in pmol
        [1]: 238 chemblank error in pmol
        [2]: 232 chemblank value in pmol
        [3]: 232 chemblank error in pmol
        [4]: 230 chemblank value in fmol
        [5]: 230 chemblank error in fmol
        
    Provides the option of printing excel file with chem blank information for personal use.
    """
    
    def __init__(self, spike_input, chem_spike_wt, U_wash, Th_wash, U_chemblank, Th_chemblank, inquiry_input,
                 chemblank_filename):
        
        spike = str(spike_input)
        
        #derives spike value based off dictionary entries
        spike_six_three_dictionary = {"DIII-B":1.008398,"DIII-A": 1.008398,"1I":1.010128,"1H":1.010128}
        spike_six_three_err_dictionary = {"DIII-B": 0.00015, "DIII-A": 0.00015, "1I": 0.00015, "1H": 0.00015}
        spike_three_dictionary = {"DIII-B": 0.78938, "DIII-A": 0.78933, "1I": 0.61351, "1H": 0.78997}
        spike_nine_dictionary = {"DIII-B": 0.21734, "DIII-A": 0.21705, "1I": 0.177187, "1H": 0.22815}
        spike_zero_nine_dictionary = {"DIII-B": 0.0000625, "DIII-A": 0.0000625, "1I": 0.0000402, "1H": 0.0000402}
        spike_zero_nine_err_dictionary = {"DIII-B": 0.000003, "DIII-A": 0.000003, "1I": 0.0000011, "1H": 0.0000011}
        spike_nine_two_dictionary = {"DIII-B": 0.00, "DIII-A": 0.00, "1I": 0.00, "1H": 0.00}
        spike_nine_two_err_dictionary = {"DIII-B": 0.00, "DIII-A": 0.00, "1I": 0.00, "1H": 0.00}
        spike_four_three_dictionary = {"DIII-B": 0.003195, "DIII-A": 0.003195, "1I":0.003180, "1H": 0.003180}
        spike_four_three_err_dictionary= {"DIII-B": 0.000003, "DIII-A": 0.000003, "1I": 0.000003, "1H": 0.000003}
        spike_five_three_dictionary = {"DIII-B": 0.10532, "DIII-A": 0.10532, "1I": 0.10521, "1H":0.10521}
        spike_five_three_err_dictionary = {"DIII-B": 0.00003, "DIII-A": 0.00003, "1I": 0.00003, "1H": 0.00003}
        spike_eight_three_dictionary = {"DIII-B": 0.01680, "DIII-A": 0.01680, "1I": 0.01700, "1H":0.01700 }
        spike_eight_three_err_dictionary = {"DIII-B": 0.00001, "DIII-A": 0.00001,"1I": 0.00001, "1H": 0.00001}
        
        if spike in spike_six_three_dictionary:
            self.spike_six_three = float(spike_six_three_dictionary[spike]) #spike ratio
        else: 
            print 'ERROR: You did not enter a valid spike option'
        
        if spike in spike_six_three_err_dictionary: 
            self.spike_six_three_err = float(spike_six_three_err_dictionary[spike]) #error of spike ratio
            
        if spike in spike_three_dictionary:
            self.spike_three = float(spike_three_dictionary[spike]) #in pmol/g
        else:pass

        if spike in spike_nine_dictionary:
            self.spike_nine = float(spike_nine_dictionary[spike]) #in pmol/g
        else: pass
    
        if spike in spike_zero_nine_dictionary:
            self.spike_zero_nine = float(spike_zero_nine_dictionary[spike]) #spike ratio
        else: pass
    
        if spike in spike_zero_nine_err_dictionary:
            self.spike_zero_nine_err = float(spike_zero_nine_err_dictionary[spike]) #error of spike ratio
        else: pass
    
        if spike in spike_nine_two_dictionary: 
            self.spike_nine_two = float(spike_nine_two_dictionary[spike]) #spike ratio
        else: pass
    
        if spike in spike_nine_two_err_dictionary:
            self.spike_nine_two_err = float(spike_nine_two_err_dictionary[spike]) #error of spike ratio
        else: pass
    
        if spike in spike_four_three_dictionary:
            self.spike_four_three = float(spike_four_three_dictionary[spike]) #spike ratio
        else: pass
    
        if spike in spike_four_three_err_dictionary:
            self.spike_four_three_err = float(spike_four_three_err_dictionary[spike]) #error of spike ratio
        else: pass
        
        if spike in spike_five_three_dictionary:
            self.spike_five_three = float(spike_five_three_dictionary[spike]) #spike ratio
        else: pass
    
        if spike in spike_five_three_err_dictionary:
            self.spike_five_three_err = float(spike_five_three_err_dictionary[spike]) #error of spike ratio
        else: pass
    
        if spike in spike_eight_three_dictionary:
            self.spike_eight_three = float(spike_eight_three_dictionary[spike]) #spike ratio
        else: pass
    
        if spike in spike_eight_three_err_dictionary:
            self.spike_eight_three_err = float(spike_eight_three_err_dictionary[spike]) #error of spike ratio
        else: pass
    
        inquiry = str(inquiry_input)
         
        #allows you to print as you go
        if inquiry.lower() == "y":
            self.inquiry = True
        else:
            self.inquiry = False
            print "Program chemblank_values ran without printing"
        
        #spike weight used for chem blank
        self.spike_wt = float(chem_spike_wt)
        
        #spike in chem blank (pmol)
        self.spike_three_used = self.spike_three * self.spike_wt
        self.spike_nine_used = self.spike_nine * self.spike_wt
        
        #allows you to write new chem blank file
        while not chemblank_filename:
            self.chemblank_file = False
            break
        else: 
            self.chemblank_file = True
            self.chemblank_filename = str(chemblank_filename)
            
        #file inputs
        self.U_wash = U_wash
        self.Th_wash = Th_wash
        self.U_file = U_chemblank
        self.Th_file = Th_chemblank
        
        #constants to be used throughout the class
        self.wt_229 = 229.031756
        self.wt_230 = 230.033128
        self.wt_232 = 232.038051
        self.wt_233 = 233.039629
        self.wt_234 = 234.040947
        self.wt_235 = 235.043924
        self.wt_236 = 236.045563
        self.wt_238 = 238.050785
    
    def blank_calculate(self):
        """
        Calculates wash and chem blank values for all isotopes. Returns a list of three values:
            1. mean
            2. counts
            3. relative error
        Input needed: file name, column letter, isotope analyzed
        """
        
        """
        Th wash and chem blank values
        """
        #wash 229 Th
        working_a = isofilter.chem_blank(self.Th_wash, "C", "229")
        nine_wash = working_a.calc()
        
        #chem blank 229 Th
        working_b = isofilter.chem_blank(self.Th_file, "C", "229")
        nine = working_b.calc()
        
        #wash 230 Th
        working_c = isofilter.chem_blank(self.Th_wash, "D", "230")
        zero_wash = working_c.calc()
        
        #chem blank 230 Th
        working_d = isofilter.chem_blank(self.Th_file, "D", "230")
        zero = working_d.calc()
        
        #wash 232 Th
        working_e = isofilter.chem_blank(self.Th_wash, "E", "232")
        two_wash = working_e.calc()
        
        #chem blank 232Th
        working_f = isofilter.chem_blank(self.Th_file, "E", "232")
        two = working_f.calc()
        
        """
        U wash and chem blank values
        """
        
        #wash 233U
        working_g = isofilter.chem_blank(self.U_wash, "D", "233")
        three_wash = working_g.calc()
        
        #chem blank 233U
        working_h = isofilter.chem_blank(self.U_file, "D", "233")
        three = working_h.calc()
        
        #wash 234U
        working_i = isofilter.chem_blank(self.U_wash, "E", "234")
        four_wash = working_i.calc()
        
        #chem blank 234U
        working_j = isofilter.chem_blank(self.U_file, "E", "234")
        four = working_j.calc()
        
        #wash 235U
        working_k = isofilter.chem_blank(self.U_wash, "F", "235")
        five_wash = working_k.calc()
        
        #chem blank 235U
        working_l = isofilter.chem_blank(self.U_file, "F", "235")
        five = working_l.calc()
        
        #wash 236U
        working_m = isofilter.chem_blank(self.U_wash, "G", "236")
        six_wash = working_m.calc()
        
        #chem blank 236U
        working_n = isofilter.chem_blank(self.U_file, "G", "236")
        six = working_n.calc()
        
        #wash 238U
        working_o = isofilter.chem_blank(self.U_wash, "H", "238")
        eight_wash = working_o.calc()
        
        #chem blank 238U
        working_p = isofilter.chem_blank(self.U_file, "H", "238")
        eight = working_p.calc()
        
        """
        Calculates signal isotopic ratio and 2s error
        
        Note: [0]: mean, [1]: counts, [2] = 2s rel error
        
        """
        #230/229
        zero_nine = (zero[0] - zero_wash[0]) / (nine[0] - nine_wash[0])
        zero_nine_err = np.sqrt( ((zero[0]*zero[2])**2/(nine[0] - nine_wash[0])**2) + 
                                 ((zero_wash[0]*zero_wash[2])**2/(nine_wash[0]-nine[0])**2) + 
                                 ((nine[0]*nine[2])**2 * ((zero_wash[0]-zero[0])/((nine[0]-nine_wash[0])**2))**2) + 
                                 ((nine_wash[0]*nine_wash[2])**2 * ((zero[0]-zero_wash[0])/((nine[0]-nine_wash[0])**2))**2)
                                 )/zero_nine
        
        #229/232 
        nine_two = (nine[0] - nine_wash[0]) / (two[0] - two_wash[0])
        nine_two_err = np.sqrt( ((nine[0]*nine[2])**2/(two[0] - two_wash[0])**2) + 
                                 ((nine_wash[0]*nine_wash[2])**2/(two_wash[0]-two[0])**2) + 
                                 ((two[0]*two[2])**2 * ((nine_wash[0]-nine[0])/((two[0]-two_wash[0])**2))**2) + 
                                 ((two_wash[0]*two_wash[2])**2 * ((nine[0]-nine_wash[0])/((two[0]-two_wash[0])**2))**2)
                                 )/nine_two
       
        #234/233
        four_three = (four[0] - four_wash[0])/(three[0] - three_wash[0])
        four_three_err = np.sqrt( ((four[0]*four[2])**2/(three[0] - three_wash[0])**2) + 
                                 ((four_wash[0]*four_wash[2])**2/(three_wash[0]-three[0])**2) + 
                                 ((three[0]*three[2])**2 * ((four_wash[0]-four[0])/((three[0]-three_wash[0])**2))**2) + 
                                 ((three_wash[0]*three_wash[2])**2 * ((four[0]-four_wash[0])/((three[0]-three_wash[0])**2))**2)
                                 )/four_three
        
        #235/233
        five_three = (five[0] - five_wash[0])/(three[0] - three_wash[0])
        five_three_err = np.sqrt( ((five[0]*five[2])**2/(three[0] - three_wash[0])**2) + 
                                 ((five_wash[0]*five_wash[2])**2/(three_wash[0]-three[0])**2) + 
                                 ((three[0]*three[2])**2 * ((five_wash[0]-five[0])/((three[0]-three_wash[0])**2))**2) + 
                                 ((three_wash[0]*three_wash[2])**2 * ((five[0]-five_wash[0])/((three[0]-three_wash[0])**2))**2)
                                 )/five_three
        
        #236/233 
        six_three = (six[0] - six_wash[0])/(three[0] - three_wash[0])
        
        #238/233
        eight_three = (eight[0] - eight_wash[0])/(three[0] - three_wash[0])
        eight_three_err = np.sqrt( ((eight[0]*eight[2])**2/(three[0] - three_wash[0])**2) + 
                                 ((eight_wash[0]*eight_wash[2])**2/(three_wash[0]-three[0])**2) + 
                                 ((three[0]*three[2])**2 * ((eight_wash[0]-eight[0])/((three[0]-three_wash[0])**2))**2) + 
                                 ((three_wash[0]*three_wash[2])**2 * ((eight[0]-eight_wash[0])/((three[0]-three_wash[0])**2))**2)
                                 )/eight_three
        
        """
        Corrects signal isotopic ratios for fractionation
    
        """
        #230/229 fract. corrected
        zero_nine_corr = zero_nine * (self.spike_six_three/six_three)**(np.log(self.wt_230/self.wt_229)/np.log(self.wt_236/self.wt_233))
        
        #229/232 fract. corrected
        nine_two_corr = nine_two * (self.spike_six_three/six_three)**(np.log(self.wt_229/self.wt_232)/np.log(self.wt_236/self.wt_233))
        
        #234/233 fract. corrected
        
        four_three_corr = four_three * (self.spike_six_three/six_three)**(np.log(self.wt_234/self.wt_233)/np.log(self.wt_236/self.wt_233))
        
        #235/233 fract. corrected
        
        five_three_corr = five_three * (self.spike_six_three/six_three)**(np.log(self.wt_235/self.wt_233)/np.log(self.wt_236/self.wt_233))
        
        #238/233 fract. corrected
        
        eight_three_corr = eight_three * (self.spike_six_three/six_three)**(np.log(self.wt_238/self.wt_233)/np.log(self.wt_236/self.wt_233))
        
        
        """
        2s relative spike errors
        """
        
        zero_nine_spike_err = self.spike_zero_nine_err/self.spike_zero_nine
        
        nine_two_spike_err = 0 #may need to change for different spikes
        
        four_three_spike_err = self.spike_four_three_err/self.spike_four_three
        
        five_three_spike_err = self.spike_five_three_err/self.spike_five_three
        
        eight_three_spike_err = self.spike_eight_three_err/self.spike_eight_three
        
        """
        Isotopic ratio and 2s error corrected for fractionation and spike signal
        """
        
        #230/229 fract. corrected and spike corrected
        
        zero_nine_corr_spike = zero_nine_corr - self.spike_zero_nine
        
        zero_nine_corr_spike_err = np.sqrt((zero_nine * zero_nine_err)**2 + (self.spike_zero_nine * zero_nine_spike_err)**2)/ abs(zero_nine_corr_spike)
        
        #229/232 fract. corrected and spike corrected
        
        nine_two_corr_spike = nine_two_corr - self.spike_nine_two
        
        nine_two_corr_spike_err = np.sqrt((nine_two * nine_two_err)**2 + (self.spike_nine_two * nine_two_spike_err)**2)/ abs(nine_two_corr_spike)
        
        #234/233 fract. corrected and spike corrected
        
        four_three_corr_spike = four_three_corr - self.spike_four_three
        
        four_three_corr_spike_err = np.sqrt((four_three * four_three_err)**2 + (self.spike_four_three * four_three_spike_err)**2)/ abs(four_three_corr_spike)
        
        #235/233 fract. corrected and spike corrected
        
        five_three_corr_spike = five_three_corr - self.spike_five_three
        
        five_three_corr_spike_err = np.sqrt((five_three * five_three_err)**2 + (self.spike_five_three * five_three_spike_err)**2)/ abs(five_three_corr_spike)
        
        #238/233 fact. corrected and spike corrected 
        
        eight_three_corr_spike = eight_three_corr - self.spike_eight_three
        
        eight_three_corr_spike_err = np.sqrt((eight_three * eight_three_err)**2 + (self.spike_eight_three * eight_three_spike_err)**2)/ abs(eight_three_corr_spike)
        
        
        """
        Option of printing chem blank values for export into Excel
        """
        
        while self.chemblank_file:
            
            sample_name = str(raw_input("What is the name of your chem blank?: " ))
            sample_date = str(raw_input("On what date did you run your chem blank?: "))
            
            self.zero_chemblank = ((self.spike_nine_used * zero_nine_corr_spike)/(10**12))* self.wt_230 * (10**18) #in ag
            self.zero_chemblank_err = abs(self.zero_chemblank * zero_nine_corr_spike_err)
            
            self.two_chemblank = ((self.spike_nine_used / nine_two_corr_spike)/(10**12))* self.wt_232 * (10**15) #in fg
            self.two_chemblank_err = abs(self.two_chemblank * nine_two_corr_spike_err)
            
            self.four_chemblank = ((self.spike_three_used * four_three_corr_spike)/(10**12))* self.wt_234 * (10**18) #in ag
            self.four_chemblank_err = abs(self.four_chemblank * four_three_corr_spike_err)
            
            self.five_chemblank = ((self.spike_three_used * five_three_corr_spike)/(10**12))* self.wt_235 * (10**15) #in fg
            self.five_chemblank_err = abs(self.five_chemblank * five_three_corr_spike_err)
            
            self.eight_chemblank = ((self.spike_three_used * eight_three_corr_spike)/(10**12))* self.wt_238 * (10**15) #in fg
            self.eight_chemblank_err = abs(self.eight_chemblank * eight_three_corr_spike_err)
            
            data = {'1_info': pd.Series([sample_name, sample_date], index = ['1_filename', '2_date']),
                    '230Th': pd.Series([self.zero_chemblank, self.zero_chemblank_err, 'ag'], index = ['3_chemistry blank', '4_2s abs. err', '5_units']),
                    '232Th': pd.Series([self.two_chemblank, self.two_chemblank_err, 'fg'], index = ['3_chemistry blank', '4_2s abs. err', '5_units']),
                    '234U': pd.Series([self.four_chemblank, self.four_chemblank_err, 'ag'], index = ['3_chemistry blank', '4_2s abs. err', '5_units']),
                    '235U': pd.Series([self.five_chemblank, self.five_chemblank_err, 'fg'], index = ['3_chemistry blank', '4_2s abs. err', '5_units']),
                    '238U': pd.Series([self.eight_chemblank, self.eight_chemblank_err, 'fg'], index = ['3_chemistry blank', '4_2s abs. err', '5_units'])}   
            
            df = pd.DataFrame(data)
            
            writer = pd.ExcelWriter(self.chemblank_filename, engine = 'openpyxl')
            
            df.to_excel(writer)
            
            writer.save()
            
            print "Chemblank data file saved in current folder under filename " + str(self.chemblank_filename)
            
            break
            
        """
        Final calculation for chem blank values for input into Age Calculation
        """
        
        #238 chem blank value and error in pmol
        
        self.chem_blank_eight = self.spike_three_used * eight_three_corr_spike 
        
        self.chem_blank_eight_err = abs(self.chem_blank_eight * eight_three_corr_spike_err) 
        
        #232 chem blank value and error in pmol
        
        self.chem_blank_two = self.spike_nine_used * nine_two_corr_spike 
        
        self.chem_blank_two_err = abs(self.chem_blank_two * nine_two_corr_spike_err)
        
        #230 chem blank value and error in fmol
        
        self.chem_blank_zero = (self.spike_nine_used * zero_nine_corr_spike) * 1000
        
        self.chem_blank_zero_err = abs(self.chem_blank_zero * zero_nine_corr_spike_err)
        
        while self.inquiry: 
            print "CHEMBLANK VALUES FOR AGE CALC: " 
            print "238U chemblank: " + str(self.chem_blank_eight) + " ± " + str(self.chem_blank_eight_err) + " pmol"
            print "232Th chemblank: " + str(self.chem_blank_two) + " ± " + str(self.chem_blank_two_err) + " pmol"
            print "230Th chemblank: " + str(self.chem_blank_zero) +  " ± " + str(self.chem_blank_zero_err) + " fmol"
            break
            
        lstChemBlank = [self.chem_blank_eight,self.chem_blank_eight_err, self.chem_blank_two, self.chem_blank_two_err,
                        self.chem_blank_zero, self.chem_blank_zero_err]
        
        return lstChemBlank
    

            
            
            
     
        
        
        
    
        
        
        
        
        
        
        
        
        
        
        
    
        
    
        
        
        
        
            
        

        
        
        
        
        
        

    
                

       