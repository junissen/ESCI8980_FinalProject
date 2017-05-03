#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Age Calculation function to analyze U-Th data using an MC-ICP-MS. 

Created by Julia Nissen, 2017, for use in Edwards' Isotope Lab.
"""

import isocalc
import numpy as np
from scipy.optimize import fsolve
import openpyxl


def Age_Calculation(U_file, Th_file, U_wash, Th_wash, 
                    U_chemblank, Th_chemblank, U_chemblank_wash, Th_chemblank_wash, 
                    AS, sample_wt, spike_wt, chem_spike_wt, 
                    spike_input, inquiry_input, year_input):
    """
    Input variables for Age Calculation
    """
    #spike values for spike used
    spike = str(spike_input)
    
    #derives spike value based off dictionary entries
    spike_six_three_dictionary = {"DIII-B":1.008398,"DIII-A": 1.008398,"1I":1.010128,"1H":1.010128}
    spike_six_three_err_dictionary = {"DIII-B": 0.00015, "DIII-A": 0.00015, "1I": 0.00015, "1H": 0.00015}
    spike_three_dictionary = {"DIII-B": 0.78938, "DIII-A": 0.78933, "1I": 0.61351, "1H": 0.78997}
    spike_three_err_dictionary = {"DIII-B": 0.00002, "DIII-A": 0.00002, "1I": 0.00002, "1H": 0.00002}
    spike_nine_dictionary = {"DIII-B": 0.21734, "DIII-A": 0.21705, "1I": 0.177187, "1H": 0.22815}
    spike_nine_err_dictionary = {"DIII-B": 0.00001, "DIII-A": 0.00002, "1I": 0.00001, "1H": 0.00001}
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
        spike_six_three = float(spike_six_three_dictionary[spike]) #spike ratio
    else: 
        print 'ERROR: You did not enter a valid spike option'
    
    if spike in spike_six_three_err_dictionary: 
        spike_six_three_err = float(spike_six_three_err_dictionary[spike]) #error of spike ratio
        
    if spike in spike_three_dictionary:
        spike_three = float(spike_three_dictionary[spike]) #in pmol/g
    else:pass

    if spike in spike_three_err_dictionary:
        spike_three_err = float(spike_three_err_dictionary[spike]) #in pmol/g
    else:pass

    if spike in spike_nine_dictionary:
        spike_nine = float(spike_nine_dictionary[spike]) #in pmol/g
    else: pass

    if spike in spike_nine_err_dictionary: 
        spike_nine_err = float(spike_nine_err_dictionary[spike]) #in pmol/g
    else: pass

    if spike in spike_zero_nine_dictionary:
        spike_zero_nine = float(spike_zero_nine_dictionary[spike]) #spike ratio
    else: pass

    if spike in spike_zero_nine_err_dictionary:
        spike_zero_nine_err = float(spike_zero_nine_err_dictionary[spike]) #error of spike ratio
    else: pass

    if spike in spike_nine_two_dictionary: 
        spike_nine_two = float(spike_nine_two_dictionary[spike]) #spike ratio
    else: pass

    if spike in spike_nine_two_err_dictionary:
        spike_nine_two_err = float(spike_nine_two_err_dictionary[spike]) #error of spike ratio
    else: pass

    if spike in spike_four_three_dictionary:
        spike_four_three = float(spike_four_three_dictionary[spike]) #spike ratio
    else: pass

    if spike in spike_four_three_err_dictionary:
        spike_four_three_err = float(spike_four_three_err_dictionary[spike]) #error of spike ratio
    else: pass
        
    if spike in spike_five_three_dictionary:
        spike_five_three = float(spike_five_three_dictionary[spike]) #spike ratio
    else: pass
    
    if spike in spike_five_three_err_dictionary:
        spike_five_three_err = float(spike_five_three_err_dictionary[spike]) #error of spike ratio
    else: pass
    
    if spike in spike_eight_three_dictionary:
        spike_eight_three = float(spike_eight_three_dictionary[spike]) #spike ratio
    else: pass
    
    if spike in spike_eight_three_err_dictionary:
        spike_eight_three_err = float(spike_eight_three_err_dictionary[spike]) #error of spike ratio
    else: pass
    
    #sample information
    sample_name = str(raw_input("What is the sample ID?: "))
    
    #year run
    year = float(year_input)

    #option for exporting chem blank values

    chemblank_inquiry = raw_input("Would you like to export your chem blank values? [y/n]: ")
    
    if str(chemblank_inquiry.lower()) == 'y':
        question2 = raw_input("What file would you like to write to (please include type, i.e. '.xlsx')?: ")
        chemblank_filename = str(question2)
    else: chemblank_filename = False    
    
    #option for printing as you go
    if str(inquiry_input.lower()) == 'y': 
        inquiry = True
    else: inquiry = False
    
    #constants needed in calculations
    wt_229 = 229.031756
    wt_230 = 230.033128
    wt_232 = 232.038051
    wt_233 = 233.039629
    wt_234 = 234.040947
    wt_235 = 235.043924
    wt_236 = 236.045563
    wt_238 = 238.050785
    five_counttime = 0.131
    four_counttime = 1.049
    three_counttime = 0.393
    two_nine_counttime = 1.049
    eight_five_rat = 137.82 #why not 137.83? 
    eight_filament_blank = 0.0001
    eight_filament_blank_err = 0.1
    sample_wt_err = 0.000005
    spike_wt_err = 0.000005
    two_nine_spike = 0.00065
    two_nine_spike_err = 0.00005
    AS_1amu = 1.00E-10
    AS_1amu_err = 0.25 * AS_1amu
    AS_2amu = AS_1amu/2.5
    AS_2amu_err = 0.25 * AS_2amu
    lambda_238 = 0.000000000155125
    lambda_234 = 0.0000028263*0.9985
    lambda_230 = 0.0000091577*1.0014
    threefive_four = 1E-11
    fourfour_four = 1E-11
    
    """
    Input functions for U, Th, wash, and chem blank values for use in Age Calculation
    """
    
    wb_U = isocalc.Ucalculation(spike_input, AS, U_file, inquiry_input)
    
    lstU_Th = wb_U.U_normalization_forTh() #provides a list for use in Th normalization

    lstU_Age = wb_U.U_normalization_forAge() #provides a list for use in Age Calculation
    """
        lstU_Age output is a list of the following values: 
            [0]: 235/233 normalized ratio
            [1]: 235/233 normalized ratio error
            [2]: 235/234 normalized and corrected ratio
            [3]: 235/234 normalized and corrected ratio error
            [4]: Unfiltered 233 counts
            [5]: Filtered 234/235 counts
            [6]: Unfiltered 233 mean
    """

    wb_Th = isocalc.Thcalculation(spike_input, AS, Th_file, inquiry_input, lstU_Th)
    
    lstTh_Age = wb_Th.Th_normalization_forAge() #provides a list to use for Age Calculation
    """
        lstTh_Age provides a list of the following outputs for the Age Calculation: 
            [0]: 230/229 corrected and normalized ratio
            [1]: 230/229 corrected and normalized ratio error
            [2]: 232/229 corrected and normalized ratio
            [3]: 232/229 corrected and normalized ratio error
            [4]: Unfiltered 229 mean
            [5]: Unfiltered 229 counts
    """

    wb_wash = isocalc.background_values(U_wash, Th_wash, inquiry_input)

    lstU_wash = wb_wash.U_wash() #provides a list of 233, 234, 235 wash values for use in Age Calculation
    """
        lstU_wash provides a list of the following outputs for the Age Calculation: 
            [0]: 233 unfiltered wash in cps
            [1]: 234 unfiltered wash in cps
            [2]: 235 unfiltered wash in cps
            
    """
    
    Th_wash = wb_wash.Th_wash() #provides the 230 darknoise cpm for use in Age Calculation
    
    wb_chemblank = isocalc.chemblank_values("1H", chem_spike_wt,
                                            U_chemblank_wash, Th_chemblank_wash, 
                                            U_chemblank, Th_chemblank, inquiry_input, chemblank_filename)
    lst_chemblank = wb_chemblank.blank_calculate() #calculates chem blanks for use in Age Calculation
    """
        lst_chemblank provides a list of the following outputs for the Age Calculation: 
            [0]: 238 chemblank value in pmol
            [1]: 238 chemblank error in pmol
            [2]: 232 chemblank value in pmol
            [3]: 232 chemblank error in pmol
            [4]: 230 chemblank value in fmol
            [5]: 230 chemblank error in fmol
    """
    
    """
    Age Calculation equations
    """
    
    #238 ppb
    
    five_three_max_err = ( (lstU_Age[6] * lstU_Age[0]) - lstU_wash[2] ) / (lstU_Age[6] - lstU_wash[0])
    
    eight_nmol = (((five_three_max_err -  spike_five_three) * spike_wt * spike_three * eight_five_rat)/1000) / sample_wt  
    
    chemblank_corr_238 = ((eight_nmol * sample_wt) - (lst_chemblank[0]/1000)) / sample_wt
    
    filament_blank_corr_238 = chemblank_corr_238 * (1 - (eight_filament_blank/ (lstU_Age[6] * five_three_max_err
                                                    * eight_five_rat)))
    
    eight_ppb = filament_blank_corr_238 * wt_238
    
    #238 ppb error
    
    rel_err_1 = (lstU_Age[1]/lstU_Age[0]) 
    
    three_counting_err = 2 / (lstU_Age[6] * lstU_Age[4] * three_counttime)**0.5
    
    five_counting_err = 2 / (lstU_Age[6] * lstU_Age[0] * five_counttime * lstU_Age[4])**0.5
    
    rel_err_2 = np.sqrt( (five_counting_err**2) + (three_counting_err**2) + (three_counting_err**2)*(8.0/9.0) )
    
    rel_err_five_three = max(rel_err_1, rel_err_2)
    
    abs_err_five_three = rel_err_five_three * five_three_max_err

    eight_nmol_err = eight_nmol * np.sqrt( ((np.sqrt((abs_err_five_three**2) + (0.0000527**2)))/(five_three_max_err - spike_five_three))**2 +
                                          (spike_wt_err/spike_wt)**2 + 
                                          (spike_three_err/spike_three)**2 + 
                                          (sample_wt_err/sample_wt)**2 )
    
    eight_nmol_err_rel = eight_nmol_err/eight_nmol
    
    chemblank_corr_238_err = np.sqrt( (eight_nmol_err**2) + ((lst_chemblank[1]/1000)**2) )
    
    chemblank_corr_238_err_rel = chemblank_corr_238_err / chemblank_corr_238
    
    filament_blank_corr_238_err_rel = np.sqrt( (chemblank_corr_238_err_rel**2) + 
                                               ( ((eight_filament_blank/(lstU_Age[6]*lstU_Age[0]*eight_five_rat)) *
                                                 np.sqrt((eight_filament_blank_err/eight_filament_blank)**2 + 
                                                          ((lstU_Age[6]*0.05)/lstU_Age[6])**2 + 
                                                          ((lstU_Age[1]/lstU_Age[0])**2)))
                                                / (1 - (eight_filament_blank/(lstU_Age[6]*lstU_Age[0]*eight_five_rat))))**2)
    
    filament_blank_corr_238_err = filament_blank_corr_238_err_rel * filament_blank_corr_238
    
    eight_ppb_err = filament_blank_corr_238_err * wt_238
    
    #232 ppt
    
    two_nine_max_err = lstTh_Age[2]
    
    two_nine_spike_corr = two_nine_max_err - two_nine_spike
    
    two_nine_chemblank_corr = two_nine_spike_corr - ( lst_chemblank[2]/(spike_wt * spike_nine)  )
    
    two_pmol = two_nine_chemblank_corr * spike_wt * spike_nine/sample_wt
                                               
    two_ppt = two_pmol * wt_232
    
    #232 ppt error
    
    abs_err_two_nine = lstTh_Age[3]
    
    two_nine_spike_corr_err = np.sqrt( (abs_err_two_nine**2) + (two_nine_spike_err **2) )
    
    two_nine_chemblank_corr_err = np.sqrt( (lst_chemblank[2]/spike_wt*spike_nine) * 
                                            np.sqrt( (lst_chemblank[3]/lst_chemblank[2])**2 + 
                                                  (spike_wt_err/spike_wt)**2 + 
                                                  (spike_nine_err/spike_nine)**2)**2 +
                                            two_nine_spike_corr_err**2)
    
    two_pmol_err = two_pmol * np.sqrt( (two_nine_chemblank_corr_err/two_nine_chemblank_corr)**2 + 
                                       (spike_wt_err/spike_wt)**2 + 
                                       (spike_nine_err/spike_nine)**2 + 
                                       (sample_wt_err/sample_wt)**2 )
    
    two_pmol_err_rel = two_pmol_err / two_pmol
    
    two_ppt_err = two_ppt * two_pmol_err_rel
    
    #230 pmol/g
    
    zero_nine_max_err = lstTh_Age[0]
    
    zero_nine_spike_corr = zero_nine_max_err - spike_zero_nine
    
    zero_nine_AS_corr = zero_nine_spike_corr - AS_1amu - (AS_2amu * lstTh_Age[2])
    
    zero_nine_darknoise_corr = zero_nine_AS_corr * (1 - ((Th_wash/60)/(lstTh_Age[4]*zero_nine_AS_corr)) )
    
    zero_nine_chemblank_corr = zero_nine_darknoise_corr - ( lst_chemblank[4]/(spike_wt * spike_nine * 1000) )
    
    zero_pmol = (zero_nine_chemblank_corr * spike_wt * spike_nine) / sample_wt
    
    #230 pmol/g error
    
    zero_nine_counting_err = lstTh_Age[0] * 2 * np.sqrt( (1 / ((lstTh_Age[4]/lstTh_Age[0])*lstTh_Age[5]*two_nine_counttime)) + 
                                                         (1 / (lstTh_Age[4]*lstTh_Age[5]*two_nine_counttime)  ) )
    
    abs_err_zero_nine = max((zero_nine_max_err*0.00001), zero_nine_counting_err, lstTh_Age[1]  )
    
    zero_nine_spike_corr_err = np.sqrt( (abs_err_zero_nine**2) + (0.000003**2) )
    
    zero_nine_AS_corr_err = np.sqrt( (zero_nine_spike_corr_err**2) + (AS_1amu_err**2) + 
                                    ( AS_2amu * lstTh_Age[2] * np.sqrt( (AS_2amu_err/AS_2amu)**2 + 
                                     (lstTh_Age[3]/lstTh_Age[2])**2 ) )**2 )
    
    zero_nine_darknoise_corr_err = zero_nine_darknoise_corr * np.sqrt( (zero_nine_AS_corr_err/zero_nine_AS_corr)**2 + 
                                                                  (((Th_wash/60)/(lstTh_Age[4]*zero_nine_AS_corr)) * 
                                                                  np.sqrt((0.2**2) + (10/lstTh_Age[4])**2 + (zero_nine_AS_corr_err/zero_nine_AS_corr)**2 
                                                                          / (1 - ((Th_wash/60)/lstTh_Age[4]*zero_nine_AS_corr) ))
                                                                          )**2)
    
    zero_nine_chemblank_corr_err = np.sqrt( zero_nine_darknoise_corr_err**2 +
                                           ( (lst_chemblank[4]/(spike_wt * spike_nine * 1000)) * 
                                            np.sqrt( (lst_chemblank[5]/lst_chemblank[4])**2 + 
                                                     (spike_wt_err/ spike_wt)**2 + 
                                                     (spike_nine_err/ spike_nine)**2 ))**2)
    
    zero_pmol_err = zero_pmol * np.sqrt((zero_nine_chemblank_corr_err/zero_nine_chemblank_corr)**2 + 
                                        (spike_wt_err/spike_wt)**2 + 
                                        (spike_nine_err/spike_nine)**2 + 
                                        (sample_wt_err/sample_wt)**2)
    
    zero_pmol_err_rel = zero_pmol_err / zero_pmol
    
    #230/232 atomic ratio
        
    zero_two_atomic = zero_pmol / two_pmol
    
    zero_two_atomic_final = zero_two_atomic * 10**6
    
    #230/232 atomic ratio error 
        
    zero_two_atomic_err_rel = np.sqrt( two_pmol_err_rel**2 + zero_pmol_err_rel**2 )
    
    zero_two_atomic_err = zero_two_atomic_err_rel * zero_two_atomic
    
    zero_two_atomic_err_final = zero_two_atomic_err * 10**6
    
    #d234U measured
        
    zero_nine_measuredU = lstU_Age[2] * (1 - lstU_wash[1]/(lstU_Age[6] * lstU_Age[2] * lstU_Age[0]))
    
    four_five_wt_avg = zero_nine_measuredU

    four_three_max_err = four_five_wt_avg * lstU_Age[0]
    
    four_five_tail_corr = four_five_wt_avg * (1 - ((4.0/9.0 * threefive_four) + (5.0/9.0 * fourfour_four)))
    
    four_five_spike_corr_234 = four_five_tail_corr * (1 - (spike_four_three/four_three_max_err))
    
    four_five_spike_corr_235 = four_five_spike_corr_234 * (1 / (1- (spike_five_three/five_three_max_err)))

    four_eight_ppm = (four_five_spike_corr_235 * 10**6) / eight_five_rat
    
    d234U_m = (( four_eight_ppm / ((lambda_238/lambda_234) * 10**6)) - 1) * 1000
    
    #d234U measured error
    
    zero_nine_measuredU_err_rel = lstU_Age[3] / zero_nine_measuredU
    
    rel_err_1 = np.sqrt(zero_nine_measuredU_err_rel**2 + (lstU_Age[1]/lstU_Age[0])**2)
    
    four_counting_err = 2 / (lstU_Age[6] * four_three_max_err * four_counttime * lstU_Age[4])**0.5
    
    rel_err_2 = np.sqrt(four_counting_err**2 + 2*three_counting_err**2 + (2.0/9.0)*three_counting_err**2)
    
    rel_err_four_three = max(rel_err_1, rel_err_2)
    
    four_five_wt_avg_err_rel = max(zero_nine_measuredU_err_rel**2, 
                                   np.sqrt(four_counting_err**2 + five_counting_err**2 + (2.0/9.0 * three_counting_err**2) ))
    
    four_five_tail_corr_err_rel = np.sqrt((four_five_wt_avg_err_rel**2) + 
                                          (np.sqrt((4.0/9.0 * threefive_four)**2 + (5.0/9.0 * fourfour_four)**2)/
                                           (1 - (4.0/9.0 * threefive_four + 5/9 * fourfour_four)) )**2)
    
    four_five_spike_corr_234_err_rel = np.sqrt((four_five_tail_corr_err_rel**2) + 
                                               ((spike_four_three/four_three_max_err) * np.sqrt(0.002**2 + rel_err_four_three**2) /
                                                (1 - spike_four_three/four_three_max_err))**2)
    
    four_five_spike_corr_235_err_rel = np.sqrt((four_five_spike_corr_234_err_rel**2) + 
                                               ((spike_five_three/five_three_max_err) * np.sqrt(0.0005**2 + (rel_err_five_three/1000)**2) /
                                                (1 - spike_five_three/five_three_max_err))**2)
    
    four_five_spike_corr_235_err = four_five_spike_corr_235_err_rel * four_five_spike_corr_235
    
    four_eight_ppm_err = (four_five_spike_corr_235_err * 10**6)/ eight_five_rat
    
    d234U_m_err = (four_eight_ppm_err / ((lambda_238/lambda_234) * 10**6)) * 1000
    
    #230Th/238U activity ratio
    
    zero_eight_atomic = (zero_pmol/(eight_ppb/wt_238))/1000
    
    zero_eight_activity = zero_eight_atomic * (lambda_230/lambda_238)
    
    #230Th/238U activity ratio error 
    
    zero_eight_atomic_err_rel = np.sqrt(zero_pmol_err_rel**2 + eight_nmol_err_rel **2 )
    
    zero_eight_activity_err = zero_eight_atomic_err_rel * zero_eight_activity
    
    #Uncorrected age calculation and error
    
    age_func = lambda t : zero_eight_activity - (1 - np.exp(-lambda_230*t) + (d234U_m/1000) * 
                                             (lambda_230/(lambda_230-lambda_234)) * 
                                             (1 - np.exp((lambda_234 - lambda_230)*t)))
    
    t_initial_guess = 0
    uncorrected_t = fsolve(age_func, t_initial_guess) #returns the value for t at which the solution is 0. This is true of all fsolve functions following this. 
    
    age_func_ThUmax = lambda t : (zero_eight_activity+zero_eight_activity_err) - (1 - np.exp(-lambda_230*t) + (d234U_m/1000) * 
                                             (lambda_230/(lambda_230-lambda_234)) * 
                                             (1 - np.exp((lambda_234 - lambda_230)*t)))
    
    uncorrected_ThUmax = fsolve(age_func_ThUmax, t_initial_guess)
    
    age_func_ThUmin = lambda t : (zero_eight_activity-zero_eight_activity_err) - (1 - np.exp(-lambda_230*t) + (d234U_m/1000) * 
                                             (lambda_230/(lambda_230-lambda_234)) * 
                                             (1 - np.exp((lambda_234 - lambda_230)*t)))
    
    uncorrected_ThUmin = fsolve(age_func_ThUmin, t_initial_guess)
    
    age_func_d234Umax = lambda t : zero_eight_activity - (1 - np.exp(-lambda_230*t) + ((d234U_m + d234U_m_err)/1000) * 
                                             (lambda_230/(lambda_230-lambda_234)) * 
                                             (1 - np.exp((lambda_234 - lambda_230)*t)))
    
    uncorrected_d234Umax = fsolve(age_func_d234Umax, t_initial_guess)
    
    age_func_d234Umin = lambda t : zero_eight_activity - (1 - np.exp(-lambda_230*t) + ((d234U_m - d234U_m_err)/1000) * 
                                             (lambda_230/(lambda_230-lambda_234)) * 
                                             (1 - np.exp((lambda_234 - lambda_230)*t)))
    
    uncorrected_d234Umin = fsolve(age_func_d234Umin, t_initial_guess)

    uncorrected_t_maxerr = np.sqrt((uncorrected_ThUmax - uncorrected_t)**2 + (uncorrected_d234Umax - uncorrected_t)**2)
    
    uncorrected_t_minerr = np.sqrt((uncorrected_ThUmin - uncorrected_t)**2 + (uncorrected_d234Umin - uncorrected_t)**2)
    
    uncorrected_t_err = (uncorrected_t_maxerr + uncorrected_t_minerr)/2
    
    #Corrected age calculation and error
    
    zero_two_initial = 0.0000044
    zero_two_initial_err = zero_two_initial/2
    
    age_func_corrected_t = lambda t : (((zero_pmol - zero_two_initial*np.exp(-lambda_230*t)*two_pmol) * lambda_230/(filament_blank_corr_238 * 1000 * lambda_238)) - 
                              (1 - np.exp(-lambda_230 * t) + (d234U_m/1000 * (lambda_230/(lambda_230-lambda_234)) * 
                              (1 - np.exp((lambda_234-lambda_230)*t)))))
    
    t_initial_guess = 0
    corrected_t = fsolve(age_func_corrected_t, t_initial_guess)
    
    zero_two_initial_now = zero_two_initial * np.exp(-lambda_230 * corrected_t)
    
    zero_two_initial_now_err = zero_two_initial_now * (zero_two_initial_err / zero_two_initial)
    
    corrected_zero_eight_activity = (zero_pmol - zero_two_initial_now*two_pmol) * lambda_230/(filament_blank_corr_238 * 1000 * lambda_238)
    
    corrected_zero_eight_activity_err = corrected_zero_eight_activity * np.sqrt( 
                                                                        (np.sqrt(((zero_two_initial_now * two_pmol) * np.sqrt((zero_two_initial_now_err/zero_two_initial_now)**2 
                                                                                + (two_pmol_err/two_pmol))**2)**2 + zero_pmol_err **2) / 
                                                                                (zero_pmol - zero_two_initial_now*two_pmol))**2 +
                                                                                (filament_blank_corr_238_err/filament_blank_corr_238)**2)

    
    age_func_ThUmax = lambda t : (corrected_zero_eight_activity+corrected_zero_eight_activity_err) - (1 - np.exp(-lambda_230*t) + (d234U_m/1000) * 
                                             (lambda_230/(lambda_230-lambda_234)) * 
                                             (1 - np.exp((lambda_234 - lambda_230)*t)))
    
    corrected_ThUmax = fsolve(age_func_ThUmax, t_initial_guess)
    
    age_func_ThUmin = lambda t : (corrected_zero_eight_activity-corrected_zero_eight_activity_err) - (1 - np.exp(-lambda_230*t) + (d234U_m/1000) * 
                                             (lambda_230/(lambda_230-lambda_234)) * 
                                             (1 - np.exp((lambda_234 - lambda_230)*t)))
    
    corrected_ThUmin = fsolve(age_func_ThUmin, t_initial_guess)
    
    age_func_d234Umax = lambda t : corrected_zero_eight_activity - (1 - np.exp(-lambda_230*t) + ((d234U_m + d234U_m_err)/1000) * 
                                             (lambda_230/(lambda_230-lambda_234)) * 
                                             (1 - np.exp((lambda_234 - lambda_230)*t)))
    
    corrected_d234Umax = fsolve(age_func_d234Umax, t_initial_guess)
    
    age_func_d234Umin = lambda t : corrected_zero_eight_activity - (1 - np.exp(-lambda_230*t) + ((d234U_m - d234U_m_err)/1000) * 
                                             (lambda_230/(lambda_230-lambda_234)) * 
                                             (1 - np.exp((lambda_234 - lambda_230)*t)))
    
    corrected_d234Umin = fsolve(age_func_d234Umin, t_initial_guess)
    
    age_func_low = lambda t: ((zero_pmol - ((zero_two_initial_now + zero_two_initial_now_err) * np.exp(-lambda_230 * t)) *two_pmol) 
                            * lambda_230/(filament_blank_corr_238 * 1000 * lambda_238)) - (1 - np.exp(-lambda_230*t) + ((d234U_m)/1000) * 
                                         (lambda_230/(lambda_230-lambda_234)) * 
                                         (1 - np.exp((lambda_234 - lambda_230)*t)))
    
    age_func_high = lambda t: ((zero_pmol - ((zero_two_initial_now - zero_two_initial_now_err) * np.exp(-lambda_230 * t)) *two_pmol) 
                            * lambda_230/(filament_blank_corr_238 * 1000 * lambda_238)) - (1 - np.exp(-lambda_230*t) + ((d234U_m)/1000) * 
                                         (lambda_230/(lambda_230-lambda_234)) * 
                                         (1 - np.exp((lambda_234 - lambda_230)*t)))
    
    corrected_age_low = fsolve(age_func_low, t_initial_guess)
    
    corrected_age_high = fsolve(age_func_high, t_initial_guess)
    
    corrected_t_maxerr = np.sqrt((corrected_ThUmax - corrected_t)**2 + (corrected_d234Umax - corrected_t)**2 + (corrected_age_high - corrected_t)**2 )
    
    corrected_t_minerr = np.sqrt((corrected_ThUmin - corrected_t)**2 + (corrected_d234Umin - corrected_t)**2 + (corrected_age_low - corrected_t)**2 )
    
    corrected_t_err = (corrected_t_maxerr + corrected_t_minerr)/2
    
    #Corrected initial d234U and error
    
    d234U_i = d234U_m * np.exp(lambda_234 * corrected_t)
    
    d234U_i_maxerr = np.sqrt( (d234U_m_err * np.exp(lambda_234 * corrected_t))**2 + 
                             (d234U_m * np.exp((lambda_234 * (corrected_t + corrected_t_maxerr)) - d234U_i))**2)
    
    d234U_i_minerr = np.sqrt( (d234U_m_err * np.exp(lambda_234 * corrected_t))**2 + 
                             (d234U_m * np.exp((lambda_234 * (corrected_t - corrected_t_minerr)) - d234U_i))**2)
    
    d234U_i_err = (d234U_i_maxerr + d234U_i_minerr)/2
    
    #Corrected age BP
    
    corrected_t_BP = corrected_t - year
    
    corrected_t_BP_err = corrected_t_err
    
    while inquiry: 
            print "AGE CALCULATION VALUES:"
            print "238 ppb: " + str(eight_ppb) + " ± " + str(eight_ppb_err)
            print "232 ppt: " + str(two_ppt) + " ± " + str(two_ppt_err)
            print "230/232 atomic (10*6) ratio: " + str(zero_two_atomic_final) + " ± " + str(zero_two_atomic_err_final)
            print "d234U measured: " + str(d234U_m) + " ± " + str(d234U_m_err)
            print "230/238 activity ratio: " + str(zero_eight_activity) + " ± " + str(zero_eight_activity_err)
            print "230Th Age uncorrected: %f" % uncorrected_t + " ± %f" % uncorrected_t_err + " yrs"
            print "230Th Age corrected: %f" % corrected_t + " ± %f" % corrected_t_err + " yrs"
            print "d234U initial corrected: %f" % d234U_i + " ± %f" % d234U_i_err
            print "230Th Age corrected: %f" % corrected_t_BP + " ± %f" %corrected_t_BP_err + " yrs BP"
            print "Age Calculation has finished"
            break
    while not inquiry:
            print "Program Age Calculation ran without printing"
            break
    
    
    workbook = str(raw_input("Into what file should the age data be written? Please include file format (i.e. xlsx): " ))
    row = str(raw_input("Into what row should the age data be written?: "))
    
    age_file = openpyxl.load_workbook(workbook)
    
    sheet = age_file.get_sheet_by_name('Sheet1')
    
    sheet['B' + row] = sample_name
    sheet['C' + row] = "{0:.1f}".format(eight_ppb)
    sheet['D' + row] = "± " + "{0:.1f}".format(eight_ppb_err)
    sheet['E' + row] = "{0:.0f}".format(two_ppt)
    sheet['F' + row] = "± " + "{0:.0f}".format(two_ppt_err)
    sheet['G' + row] = "{0:.1f}".format(zero_two_atomic_final)
    sheet['H' + row] = "± " + "{0:.1f}".format(zero_two_atomic_final)
    sheet['I' + row] = "{0:.1f}".format(d234U_m)
    sheet['J' + row] = "± " + "{0:.1f}".format(d234U_m_err)
    sheet['K' + row] = "{0:.5f}".format(zero_eight_activity)
    sheet['L' + row] = "± " + "{0:.5f}".format(zero_eight_activity_err)
    sheet['M' + row] = "%.0f" % uncorrected_t
    sheet['N' + row] = "± %.0f" % uncorrected_t_err
    sheet['O' + row] = "%.0f" % corrected_t
    sheet['P' + row] = "± %.0f" % corrected_t_err
    sheet['Q' + row] = "%.1f" % d234U_i
    sheet['R' + row] = "± %.1f" % d234U_i_err
    sheet['S' + row] = "%.0f" % corrected_t_BP
    sheet['T' + row] = "± %.0f" % corrected_t_BP_err
    
    age_file.save(workbook)
    
    print "Spreadsheet has been updated. Age Calculation is finished"
    
    
   
    
    
    
    
    
        
        