#C file created using an r4ss function
#C file write time: 2026-11-13  17:31:20
#
0 # 0 means do not read wtatage.ss; 1 means read and usewtatage.ss and also read and use growth parameters
1 #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern
4 # recr_dist_method for parameters
1 # not yet implemented; Future usage:Spawner-Recruitment; 1=global; 2=by area
1 # number of recruitment settlement assignments 
0 # unused option
# for each settlement assignment:
#_GPattern	month	area	age
1	1	1	0	#_recr_dist_pattern1
#
#_Cond 0 # N_movement_definitions goes here if N_areas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
0 #_Nblock_Patterns
#_Cond 0 #_blocks_per_pattern
# begin and end years of blocks
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#
# AUTOGEN
0 0 0 0 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement
#
3 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=Maunder_M;_6=Age-range_Lorenzen
#_ #_Age_natmort_by sex x growthpattern
#_Age_0	Age_1	Age_2	Age_3	Age_4	Age_5	Age_6	Age_7	Age_8	Age_9	Age_10	Age_11	Age_12	Age_13	Age_14	Age_15
0.42	0.37	0.32	0.27	0.22	0.22	0.22	0.22	0.22	0.22	0.22	0.22	0.22	0.22	0.22	0.22	#_natM1
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr;5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0.5 #_Age(post-settlement)_for_L1;linear growth below this
15 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0 #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
2 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn
   50	     200	   80.1	   97.7	 99	0	 -4	0	0	0	0	0	0	0	#_L_at_Amin_Fem_GP_1 
  100	     300	  226.3	  226.3	 99	0	 -2	0	0	0	0	0	0	0	#_L_at_Amax_Fem_GP_1 
 0.05	     0.3	  0.246	   0.25	 99	0	 -4	0	0	0	0	0	0	0	#_VonBert_K_Fem_GP_1 
 0.01	     0.5	    0.1	    0.1	 99	0	 -3	0	0	0	0	0	0	0	#_CV_young_Fem_GP_1  
 0.01	     0.5	    0.1	   0.15	 99	0	 -3	0	0	0	0	0	0	0	#_CV_old_Fem_GP_1    
    0	       3	1.3e-05	1.3e-05	 99	0	 -3	0	0	0	0	0	0	0	#_Wtlen_1_Fem_GP_1   
    0	       4	   3.07	   3.07	 99	0	 -3	0	0	0	0	0	0	0	#_Wtlen_2_Fem_GP_1   
    1	     200	 143.68	 143.68	 99	0	 -3	0	0	0	0	0	0	0	#_Mat50%_Fem_GP_1    
   -3	       3	-0.1034	-0.1034	 99	0	 -3	0	0	0	0	0	0	0	#_Mat_slope_Fem_GP_1 
    0	       3	      1	      1	 99	0	 -3	0	0	0	0	0	0	0	#_Eggs_alpha_Fem_GP_1
    0	       3	      0	      0	 99	0	 -3	0	0	0	0	0	0	0	#_Eggs_beta_Fem_GP_1 
  0.1	      10	      1	      1	  1	0	 -1	0	0	0	0	0	0	0	#_CohortGrowDev      
1e-06	0.999999	    0.5	    0.5	0.5	0	-99	0	0	0	0	0	0	0	#_FracFemale_GP_1    
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; 2=Ricker (2 parms); 3=std_B-H(2); 4=SCAA(2);5=Hockey(3); 6=B-H_flattop(2); 7=Survival(3);8=Shepard(3);9=Ricker_Power(3);10=B-H_a,b(4)
1 # 0/1 to use steepness in initial equ recruitment calculation
1 # future feature: 0/1 to make realized sigmaR a function of SR curvature
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn # parm_name
  3	50	6.57807	9.3	99	0	 1	0	0	0	0	0	0	0	#_SR_LN(R0)  
0.2	 1	    0.9	0.9	99	0	-4	0	0	0	0	0	0	0	#_SR_BH_steep
  0	 2	   0.42	0.6	99	0	-3	0	0	0	0	0	0	0	#_SR_sigmaR  
 -5	 5	      0	  0	99	0	-1	0	0	0	0	0	0	0	#_SR_regime  
  0	 0	      0	  0	99	0	-1	0	0	0	0	0	0	0	#_SR_autocorr
#_no timevary SR parameters
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1985 # first year of main recr_devs; early devs can preceed this era
2021 # last year of main recr_devs; forecast devs start in following year
3 #_recdev phase
1 # (0/1) to read 13 advanced options
1965 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
4 #_recdev_early_phase
0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
1 #_lambda for Fcast_recr_like occurring before endyr+1
1967.78 #_last_yr_nobias_adj_in_MPD; begin of ramp
1995.25 #_first_yr_fullbias_adj_in_MPD; begin of plateau
2020.9 #_last_yr_fullbias_adj_in_MPD
2021.73 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
0.9366 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
0 #_period of cycles in recruitment (N parms read below)
-5 #min rec_dev
5 #max rec_dev
0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
#Fishing Mortality info
0.5 # F ballpark
-1960 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
4 # max F or harvest rate, depends on F_Method
5 # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
0	3	0.342968	1	99	0	2	#_InitF_seas_1_flt_6F6_JPN_WCNPO_OSDWLL_early_Area1
#
#_Q_setup for fleets with cpue or survey data
#_fleet	link	link_info	extra_se	biasadj	float  #  fleetname
   20	1	0	0	0	1	#_S1_JPN_WCNPO_OSDWLL_early_Area1         
   21	1	0	0	0	1	#_S2_JPN_WCNPO_OSDWCOLL_late_Area1        
   22	1	0	0	0	1	#_S3_JPN_WCNPO_OSDWLL_early_Area2         
   23	1	0	0	0	1	#_S4_JPN_WCNPO_OSDWLL_late_Area2          
   24	1	0	0	0	1	#_S5_TWN_WCNPO_DWLL_late                  
   25	1	0	0	0	1	#_S6_US_WCNPO_LL_deep                     
   26	1	0	0	0	1	#_S7_US_WCNPO_LL_shallow_early            
   27	1	0	0	0	1	#_S8_US_WCNPO_LL_shallow_late             
   30	1	0	0	0	1	#_S2_JPN_WCNPO_OSDWCOLL_late_Area1_1994_97
-9999	0	0	0	0	0	#_terminator                              
#_Q_parms(if_any);Qunits_are_ln(q)
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
-15	0	-7.33837	-1	99	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S1_JPN_WCNPO_OSDWLL_early_Area1(20)         
-15	0	-7.38069	-1	99	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S2_JPN_WCNPO_OSDWCOLL_late_Area1(21)        
-15	0	-7.30484	-1	99	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S3_JPN_WCNPO_OSDWLL_early_Area2(22)         
-15	0	-7.58938	-1	99	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S4_JPN_WCNPO_OSDWLL_late_Area2(23)          
-15	0	-5.89548	-1	99	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S5_TWN_WCNPO_DWLL_late(24)                  
-15	0	-8.49161	-1	99	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S6_US_WCNPO_LL_deep(25)                     
-15	0	 -4.2936	-1	99	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S7_US_WCNPO_LL_shallow_early(26)            
-15	0	-3.95014	-1	99	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S8_US_WCNPO_LL_shallow_late(27)             
-15	0	-7.17105	-1	99	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S2_JPN_WCNPO_OSDWCOLL_late_Area1_1994_97(30)
#_no timevary Q parameters
#
#_size_selex_patterns
#_Pattern	Discard	Male	Special
24	0	0	 0	#_1 F1_JPN_WCNPO_OSDWCOLL_late_Area1         
 1	0	0	 0	#_2 F2_TWN_WCNPO_DWLL_late                   
24	0	0	 0	#_3 F3_US_WCNPO_LL_shallow_late              
 1	0	0	 0	#_4 F4_IATTC                                 
15	0	0	 4	#_5 F5_JPN_EPO_OSDWLL                        
24	0	0	 0	#_6 F6_JPN_WCNPO_OSDWLL_early_Area1          
24	0	0	 0	#_7 F7_JPN_WCNPO_CODF                        
24	0	0	 0	#_8 F8_US_WCNPO_LL_shallow_early             
15	0	0	 8	#_9 F9_US_WCNPO_LL_deep                      
15	0	0	 6	#_10 F10_JPN_WCNPO_OSDF                      
15	0	0	 6	#_11 F11_JPN_WCNPO_Other_early               
15	0	0	 1	#_12 F12_JPN_WCNPO_Other_late                
15	0	0	 2	#_13 F13_TWN_WCNPO_DWLL_early                
15	0	0	 1	#_14 F14_TWN_WCNPO_Other                     
15	0	0	 2	#_15 F15_US_WCNPO_GN                         
15	0	0	 2	#_16 F16_US_WCNPO_Other                      
15	0	0	 3	#_17 F17_JPN_WCNPO_OSDWLL_early_Area2        
15	0	0	 3	#_18 F18_JPN_WCNPO_OSDWLL_late_Area2         
15	0	0	 2	#_19 F19_WCPFC                               
15	0	0	 6	#_20 S1_JPN_WCNPO_OSDWLL_early_Area1         
15	0	0	 1	#_21 S2_JPN_WCNPO_OSDWCOLL_late_Area1        
15	0	0	 3	#_22 S3_JPN_WCNPO_OSDWLL_early_Area2         
15	0	0	 3	#_23 S4_JPN_WCNPO_OSDWLL_late_Area2          
15	0	0	 2	#_24 S5_TWN_WCNPO_DWLL_late                  
 0	0	0	 0	#_25 S6_US_WCNPO_LL_deep                     
15	0	0	 8	#_26 S7_US_WCNPO_LL_shallow_early            
15	0	0	 3	#_27 S8_US_WCNPO_LL_shallow_late             
24	0	0	 0	#_28 F1_JPN_WCNPO_OSDWCOLL_late_Area1_1994_97
15	0	0	28	#_29 F12_JPN_WCNPO_Other_late_1994_97        
15	0	0	28	#_30 S2_JPN_WCNPO_OSDWCOLL_late_Area1_1994_97
 1	0	0	 0	#_31 F2_TWN_WCNPO_DWLL_late_2004_09          
15	0	0	31	#_32 F15_US_WCNPO_GN_2004_09                 
15	0	0	31	#_33 F16_US_WCNPO_Other_2004_09              
15	0	0	31	#_34 F19_WCPFC_2004_09                       
#
#_age_selex_patterns
#_Pattern	Discard	Male	Special
0	0	0	0	#_1 F1_JPN_WCNPO_OSDWCOLL_late_Area1         
0	0	0	0	#_2 F2_TWN_WCNPO_DWLL_late                   
0	0	0	0	#_3 F3_US_WCNPO_LL_shallow_late              
0	0	0	0	#_4 F4_IATTC                                 
0	0	0	0	#_5 F5_JPN_EPO_OSDWLL                        
0	0	0	0	#_6 F6_JPN_WCNPO_OSDWLL_early_Area1          
0	0	0	0	#_7 F7_JPN_WCNPO_CODF                        
0	0	0	0	#_8 F8_US_WCNPO_LL_shallow_early             
0	0	0	0	#_9 F9_US_WCNPO_LL_deep                      
0	0	0	0	#_10 F10_JPN_WCNPO_OSDF                      
0	0	0	0	#_11 F11_JPN_WCNPO_Other_early               
0	0	0	0	#_12 F12_JPN_WCNPO_Other_late                
0	0	0	0	#_13 F13_TWN_WCNPO_DWLL_early                
0	0	0	0	#_14 F14_TWN_WCNPO_Other                     
0	0	0	0	#_15 F15_US_WCNPO_GN                         
0	0	0	0	#_16 F16_US_WCNPO_Other                      
0	0	0	0	#_17 F17_JPN_WCNPO_OSDWLL_early_Area2        
0	0	0	0	#_18 F18_JPN_WCNPO_OSDWLL_late_Area2         
0	0	0	0	#_19 F19_WCPFC                               
0	0	0	0	#_20 S1_JPN_WCNPO_OSDWLL_early_Area1         
0	0	0	0	#_21 S2_JPN_WCNPO_OSDWCOLL_late_Area1        
0	0	0	0	#_22 S3_JPN_WCNPO_OSDWLL_early_Area2         
0	0	0	0	#_23 S4_JPN_WCNPO_OSDWLL_late_Area2          
0	0	0	0	#_24 S5_TWN_WCNPO_DWLL_late                  
0	0	0	0	#_25 S6_US_WCNPO_LL_deep                     
0	0	0	0	#_26 S7_US_WCNPO_LL_shallow_early            
0	0	0	0	#_27 S8_US_WCNPO_LL_shallow_late             
0	0	0	0	#_28 F1_JPN_WCNPO_OSDWCOLL_late_Area1_1994_97
0	0	0	0	#_29 F12_JPN_WCNPO_Other_late_1994_97        
0	0	0	0	#_30 S2_JPN_WCNPO_OSDWCOLL_late_Area1_1994_97
0	0	0	0	#_31 F2_TWN_WCNPO_DWLL_late_2004_09          
0	0	0	0	#_32 F15_US_WCNPO_GN_2004_09                 
0	0	0	0	#_33 F16_US_WCNPO_Other_2004_09              
0	0	0	0	#_34 F19_WCPFC_2004_09                       
#
#_SizeSelex
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
  53	257	 156.799	  90	0.05	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F1_JPN_WCNPO_OSDWCOLL_late_Area1(1)         
 -10	  3	-9.68831	 0.1	0.05	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_2_F1_JPN_WCNPO_OSDWCOLL_late_Area1(1)         
  -4	 12	 7.44528	   6	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_3_F1_JPN_WCNPO_OSDWCOLL_late_Area1(1)         
  -2	 20	  8.2649	 0.1	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_4_F1_JPN_WCNPO_OSDWCOLL_late_Area1(1)         
-999	999	    -999	-999	  99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_5_F1_JPN_WCNPO_OSDWCOLL_late_Area1(1)         
-999	999	    -999	-999	  99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_6_F1_JPN_WCNPO_OSDWCOLL_late_Area1(1)         
  10	150	 145.714	  70	0.01	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_1_F2_TWN_WCNPO_DWLL_late(2)                   
0.01	150	 28.2468	  15	0.01	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_2_F2_TWN_WCNPO_DWLL_late(2)                   
  53	250	 167.869	  90	0.05	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F3_US_WCNPO_LL_shallow_late(3)              
 -10	  3	 -9.5107	 0.1	0.05	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_2_F3_US_WCNPO_LL_shallow_late(3)              
  -4	 12	 7.87504	   6	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_3_F3_US_WCNPO_LL_shallow_late(3)              
  -6	 10	 8.17408	 0.1	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_4_F3_US_WCNPO_LL_shallow_late(3)              
-999	999	    -999	-999	  99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_5_F3_US_WCNPO_LL_shallow_late(3)              
-999	999	    -999	-999	  99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_6_F3_US_WCNPO_LL_shallow_late(3)              
  10	150	  139.19	  70	0.01	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_1_F4_IATTC(4)                                 
0.01	150	 78.7609	  15	0.01	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_2_F4_IATTC(4)                                 
  53	257	  144.38	  90	0.05	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F6_JPN_WCNPO_OSDWLL_early_Area1(6)          
 -10	  3	-8.46536	 0.1	0.05	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_2_F6_JPN_WCNPO_OSDWLL_early_Area1(6)          
  -4	 12	 7.37286	   6	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_3_F6_JPN_WCNPO_OSDWLL_early_Area1(6)          
  -2	 10	 9.49211	 0.1	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_4_F6_JPN_WCNPO_OSDWLL_early_Area1(6)          
-999	999	    -999	-999	  99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_5_F6_JPN_WCNPO_OSDWLL_early_Area1(6)          
-999	999	    -999	-999	  99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_6_F6_JPN_WCNPO_OSDWLL_early_Area1(6)          
  50	257	 232.753	  90	0.05	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F7_JPN_WCNPO_CODF(7)                        
 -10	  3	 2.07615	 0.1	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_2_F7_JPN_WCNPO_CODF(7)                        
  -4	 12	 8.36293	   6	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_3_F7_JPN_WCNPO_CODF(7)                        
  -2	 15	 4.48318	 0.1	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_4_F7_JPN_WCNPO_CODF(7)                        
-999	999	    -999	-999	  99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_5_F7_JPN_WCNPO_CODF(7)                        
-999	999	    -999	-999	  99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_6_F7_JPN_WCNPO_CODF(7)                        
  53	250	 157.182	  90	0.05	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F8_US_WCNPO_LL_shallow_early(8)             
 -10	  3	-7.03737	 0.1	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_2_F8_US_WCNPO_LL_shallow_early(8)             
  -4	 12	 7.54159	   6	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_3_F8_US_WCNPO_LL_shallow_early(8)             
  -6	 10	 8.09351	 0.1	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_4_F8_US_WCNPO_LL_shallow_early(8)             
-999	999	    -999	-999	  99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_5_F8_US_WCNPO_LL_shallow_early(8)             
-999	999	    -999	-999	  99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_6_F8_US_WCNPO_LL_shallow_early(8)             
  53	257	 141.693	  90	0.05	0	 2	0	0	0	0	0	0	0	#_SizeSel_P_1_F1_JPN_WCNPO_OSDWCOLL_late_Area1_1994_97(28)
 -10	  3	-9.33759	 0.1	0.05	0	-4	0	0	0	0	0	0	0	#_SizeSel_P_2_F1_JPN_WCNPO_OSDWCOLL_late_Area1_1994_97(28)
  -4	 12	 6.67587	   6	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_3_F1_JPN_WCNPO_OSDWCOLL_late_Area1_1994_97(28)
  -2	 20	 16.8203	 0.1	0.05	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_4_F1_JPN_WCNPO_OSDWCOLL_late_Area1_1994_97(28)
-999	999	    -999	-999	  99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_5_F1_JPN_WCNPO_OSDWCOLL_late_Area1_1994_97(28)
-999	999	    -999	-999	  99	0	-2	0	0	0	0	0	0	0	#_SizeSel_P_6_F1_JPN_WCNPO_OSDWCOLL_late_Area1_1994_97(28)
  10	150	 142.243	  70	0.01	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_1_F2_TWN_WCNPO_DWLL_late_2004_09(31)          
0.01	150	 27.9581	  15	0.01	0	 3	0	0	0	0	0	0	0	#_SizeSel_P_2_F2_TWN_WCNPO_DWLL_late_2004_09(31)          
#_AgeSelex
#_No age_selex_parm
#_no timevary selex parameters
#
0 #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
# Tag loss and Tag reporting parameters go next
0 # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# Input variance adjustments factors: 
#_factor	fleet	value
    1	24	0.11	#_Variance_adjustment_list1
-9999	 0	   0	#_terminator               
#
3 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 16 changes to default Lambdas (default value is 1.0)
#_like_comp	fleet	phase	value	sizefreq_method
    1	24	1	  0	0	#_Surv_S5_TWN_WCNPO_DWLL_late_Phz1                                      
    1	25	1	  1	1	#_Surv_S6_US_WCNPO_LL_deep_Phz1                                         
    1	27	1	  0	0	#_Surv_S8_US_WCNPO_LL_shallow_late_Phz1                                 
    4	 1	1	0.5	1	#_length_F1_JPN_WCNPO_OSDWCOLL_late_Area1_sizefreq_method_1_Phz1        
    4	 2	1	0.5	1	#_length_F2_TWN_WCNPO_DWLL_late_sizefreq_method_1_Phz1                  
    4	 3	1	0.5	1	#_length_F3_US_WCNPO_LL_shallow_late_sizefreq_method_1_Phz1             
    4	 4	1	0.5	1	#_length_F4_IATTC_sizefreq_method_1_Phz1                                
    4	 5	1	  0	0	#_length_F5_JPN_EPO_OSDWLL_sizefreq_method_0_Phz1                       
    4	 6	1	0.5	1	#_length_F6_JPN_WCNPO_OSDWLL_early_Area1_sizefreq_method_1_Phz1         
    4	 7	1	0.5	1	#_length_F7_JPN_WCNPO_CODF_sizefreq_method_1_Phz1                       
    4	 8	1	0.5	1	#_length_F8_US_WCNPO_LL_shallow_early_sizefreq_method_1_Phz1            
    4	 9	1	0.5	1	#_length_F9_US_WCNPO_LL_deep_sizefreq_method_1_Phz1                     
    4	28	1	0.5	1	#_length_F1_JPN_WCNPO_OSDWCOLL_late_Area1_1994_97_sizefreq_method_1_Phz1
    4	31	1	0.5	1	#_length_F2_TWN_WCNPO_DWLL_late_2004_09_sizefreq_method_1_Phz1          
    9	 6	1	  0	0	#_init_equ_catch_F6_JPN_WCNPO_OSDWLL_early_Area1_Phz1                   
   10	 1	1	  1	1	#_recrdev_Phz1                                                          
-9999	 0	0	  0	0	#_terminator                                                            
#
0 # 0/1 read specs for more stddev reporting
#
999
