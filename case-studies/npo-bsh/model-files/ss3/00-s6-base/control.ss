#V3.30
#C file created using the SS_writectl function in the R package r4ss
#C file write time: 2022-04-01 19:48:10
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
1	7	1	0	#_recr_dist_pattern1
#
#_Cond 0 # N_movement_definitions goes here if N_areas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
2 #_Nblock_Patterns
1 2 #_blocks_per_pattern
#_begin and end years of blocks
1970 1970
1973 1981 1982 1992
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#
# AUTOGEN
1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement
#
3 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
#_ #_Age_natmort_by sex x growthpattern
#_Age_0	Age_1	Age_2	Age_3	Age_4	Age_5	Age_6	Age_7	Age_8	Age_9	Age_10	Age_11	Age_12	Age_13	Age_14	Age_15	Age_16	Age_17	Age_18	Age_19	Age_20	Age_21	Age_22	Age_23	Age_24
0.787	0.489	0.371	0.307	0.267	0.240	0.221	0.207	0.196	0.187	0.180	0.175	0.170	0.167	0.164	0.161	0.159	0.157	0.156	0.154	0.153	0.152	0.151	0.151	0.150	#_natM1
0.726	0.491	0.382	0.320	0.279	0.251	0.230	0.214	0.202	0.192	0.184	0.178	0.172	0.168	0.164	0.160	0.158	0.155	0.153	0.151	0.150	0.148	0.147	0.146	0.145	#_natM2
2 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr;5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
1 #_Age(post-settlement)_for_L1;linear growth below this
20 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0 #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
4 #_First_Mature_Age
2 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn
 1e+01	120.000000	 6.430e+01	 6.50000e+01	10.0	0	 -4	0	0	0	0	0.5	0	0	#_L_at_Amin_Fem_GP_1       
 4e+01	410.000000	 2.452e+02	 4.00000e+02	10.0	0	 -2	0	0	0	0	0.5	0	0	#_L_at_Amax_Fem_GP_1       
 1e-01	  0.250000	 1.460e-01	 1.50000e-01	 0.8	0	 -4	0	0	0	0	0.5	0	0	#_VonBert_K_Fem_GP_1       
-1e+01	 10.000000	 1.000e+00	 1.00000e+00	 0.8	0	 -4	0	0	0	0	0.5	0	0	#_Richards_Fem_GP_1        
 1e-02	  1.000000	 2.500e-01	 8.34877e-02	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_CV_young_Fem_GP_1        
-3e+00	  3.000000	 1.000e-01	 0.00000e+00	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_CV_old_Fem_GP_1          
-3e+00	  3.000000	 5.388e-06	 5.38800e-06	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_Wtlen_1_Fem_GP_1         
-3e+00	  3.500000	 3.102e+00	 3.10200e+00	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_Wtlen_2_Fem_GP_1         
-3e+00	300.000000	 1.566e+02	 5.50000e+01	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_Mat50%_Fem_GP_1          
-3e+00	  3.000000	-1.600e-01	-1.60000e-01	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_Mat_slope_Fem_GP_1       
-3e+00	 50.000000	 4.500e+01	 4.50000e+01	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_Eggs/kg_inter_Fem_GP_1   
-3e+00	  3.000000	 0.000e+00	 0.00000e+00	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_Eggs/kg_slope_wt_Fem_GP_1
 1e+01	120.000000	 6.850e+01	 0.00000e+00	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_L_at_Amin_Mal_GP_1       
 4e+01	410.000000	 2.615e+02	 0.00000e+00	 0.8	0	 -2	0	0	0	0	0.5	0	0	#_L_at_Amax_Mal_GP_1       
 1e-01	  0.250000	 1.170e-01	 0.00000e+00	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_VonBert_K_Mal_GP_1       
-1e+01	 10.000000	 1.000e+00	 0.00000e+00	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_Richards_Mal_GP_1        
 1e-02	  1.000000	 2.500e-01	 0.00000e+00	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_CV_young_Mal_GP_1        
-3e+00	  3.000000	 1.000e-01	 0.00000e+00	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_CV_old_Mal_GP_1          
-3e+00	  3.000000	 3.293e-06	 3.29300e-06	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_Wtlen_1_Mal_GP_1         
-3e+00	  3.500000	 3.225e+00	 3.22500e+00	 0.8	0	 -3	0	0	0	0	0.5	0	0	#_Wtlen_2_Mal_GP_1         
 1e-01	 10.000000	 1.000e+00	 1.00000e+00	 1.0	0	 -1	0	0	0	0	0.5	0	0	#_CohortGrowDev            
 1e-06	  0.999999	 5.000e-01	 5.00000e-01	 0.5	0	-99	0	0	0	0	0.0	0	0	#_FracFemale_GP_1          
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; 2=Ricker; 3=std_B-H; 4=SCAA;5=Hockey; 6=B-H_flattop; 7=survival_3Parm;8=Shepard_3Parm
1 # 0/1 to use steepness in initial equ recruitment calculation
1 # future feature: 0/1 to make realized sigmaR a function of SR curvature
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn # parm_name
 3.0	20	9.66115	9.000	10.00	0	 1	0	0	0	0	0	0	0	#_SR_LN(R0)  
 0.2	 1	0.61300	0.613	 0.05	0	-4	0	0	0	0	0	0	0	#_SR_BH_steep
 0.0	 2	0.40000	0.600	 0.80	0	-3	0	0	0	0	0	0	0	#_SR_sigmaR  
-5.0	 5	0.00000	0.000	 1.00	0	-1	0	0	0	0	0	1	1	#_SR_regime  
 0.0	 0	0.00000	0.000	99.00	0	-1	0	0	0	0	0	0	0	#_SR_autocorr
# timevary SR parameters
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
-5	5	0.00492594	0	2.5	0	4	#_SR_regime_BLK1add_1970
2 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
1971 # first year of main recr_devs; early devs can preceed this era
2020 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase
1 # (0/1) to read 13 advanced options
-10 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
2 #_recdev_early_phase
-1 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
1 #_lambda for Fcast_recr_like occurring before endyr+1
1969.3 #_last_yr_nobias_adj_in_MPD; begin of ramp
2014 #_first_yr_fullbias_adj_in_MPD; begin of plateau
2019.5 #_last_yr_fullbias_adj_in_MPD
2020 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
0.5 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
0 #_period of cycles in recruitment (N parms read below)
-10 #min rec_dev
10 #max rec_dev
0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
#Fishing Mortality info
0.2 # F ballpark
2013 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
5 # max F or harvest rate, depends on F_Method
4 # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
0.001	5	0.249539	0.01	99	0	1	#_InitF_seas_1_flt_4F4_JPN_KK_SH
#
#_Q_setup for fleets with cpue or survey data
#_fleet	link	link_info	extra_se	biasadj	float  #  fleetname
   21	1	0	0	0	1	#_S1_HW_DP         
   23	1	0	0	0	1	#_S3_TAIW_LG       
   25	1	0	0	0	1	#_S5_JPN_EARLY     
   26	1	0	0	0	1	#_S6_JPN_LATE      
   27	1	0	0	0	1	#_S7_JPN_RTV       
   29	1	0	0	0	1	#_S9_SPC_OBS_TROPIC
   30	1	0	0	0	1	#_S10_MEX          
   31	1	0	0	0	1	#_S11_DFA_LATE     
-9999	0	0	0	0	0	#_terminator       
#_Q_parms(if_any);Qunits_are_ln(q)
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
-25	25	-8.15168	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S1_HW_DP(21)         
-25	25	-7.51272	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S3_TAIW_LG(23)       
-25	25	-8.05065	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S5_JPN_EARLY(25)     
-25	25	-8.15599	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S6_JPN_LATE(26)      
-25	25	-7.60901	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S7_JPN_RTV(27)       
-25	25	-8.29024	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S9_SPC_OBS_TROPIC(29)
-25	25	-8.21581	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S10_MEX(30)          
-25	25	-7.58984	0	1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_S11_DFA_LATE(31)     
#_no timevary Q parameters
#
#_size_selex_patterns
#_Pattern	Discard	Male	Special
24	0	4	 0	#_1 F1_MEX              
 5	0	0	 1	#_2 F2_CAN              
24	0	4	 0	#_3 F3_CHINA            
24	0	4	 0	#_4 F4_JPN_KK_SH        
24	0	4	 0	#_5 F5_JPN_KK_DP        
 5	0	0	 4	#_6 F6_JPN_ENY_SHL      
24	0	4	 0	#_7 F7_JPN_ENY_DP       
24	0	4	 0	#_8 F8_JPN_LG_MESH_EARLY
24	0	4	 0	#_9 F9_JPN_LG_MESH_LATE 
 5	0	0	 9	#_10 F10_JPN_CST_Oth    
24	0	0	 0	#_11 F11_JPN_SM_MESH    
 5	0	0	 1	#_12 F12_IATTC          
 5	0	0	 3	#_13 F13_KOREA          
24	0	4	 0	#_14 F14_NON_ISC        
24	0	4	 0	#_15 F15_USA_GIILL      
 5	0	0	15	#_16 F16_USA_SPORT      
24	0	4	 0	#_17 F17_USA_Lonline_DP 
27	0	2	 5	#_18 F18_USA_Lonline_SH 
24	0	4	 0	#_19 F19_TAIW_LG        
24	0	4	 0	#_20 F20_TAIW_SM        
 5	0	0	17	#_21 S1_HW_DP           
 5	0	0	18	#_22 S2_HW_SH           
 5	0	0	19	#_23 S3_TAIW_LG         
 5	0	0	20	#_24 S4_TAIW_SM         
 5	0	0	 4	#_25 S5_JPN_EARLY       
 5	0	0	 4	#_26 S6_JPN_LATE        
 5	0	0	 7	#_27 S7_JPN_RTV         
 5	0	0	14	#_28 S8_SPC_OBS         
 5	0	0	14	#_29 S9_SPC_OBS_TROPIC  
 5	0	0	 1	#_30 S10_MEX            
 5	0	0	 7	#_31 S11_DFA_LATE       
#
#_age_selex_patterns
#_Pattern	Discard	Male	Special
0	0	0	0	#_1 F1_MEX              
0	0	0	0	#_2 F2_CAN              
0	0	0	0	#_3 F3_CHINA            
0	0	0	0	#_4 F4_JPN_KK_SH        
0	0	0	0	#_5 F5_JPN_KK_DP        
0	0	0	0	#_6 F6_JPN_ENY_SHL      
0	0	0	0	#_7 F7_JPN_ENY_DP       
0	0	0	0	#_8 F8_JPN_LG_MESH_EARLY
0	0	0	0	#_9 F9_JPN_LG_MESH_LATE 
0	0	0	0	#_10 F10_JPN_CST_Oth    
0	0	0	0	#_11 F11_JPN_SM_MESH    
0	0	0	0	#_12 F12_IATTC          
0	0	0	0	#_13 F13_KOREA          
0	0	0	0	#_14 F14_NON_ISC        
0	0	0	0	#_15 F15_USA_GIILL      
0	0	0	0	#_16 F16_USA_SPORT      
0	0	0	0	#_17 F17_USA_Lonline_DP 
0	0	0	0	#_18 F18_USA_Lonline_SH 
0	0	0	0	#_19 F19_TAIW_LG        
0	0	0	0	#_20 F20_TAIW_SM        
0	0	0	0	#_21 S1_HW_DP           
0	0	0	0	#_22 S2_HW_SH           
0	0	0	0	#_23 S3_TAIW_LG         
0	0	0	0	#_24 S4_TAIW_SM         
0	0	0	0	#_25 S5_JPN_EARLY       
0	0	0	0	#_26 S6_JPN_LATE        
0	0	0	0	#_27 S7_JPN_RTV         
0	0	0	0	#_28 S8_SPC_OBS         
0	0	0	0	#_29 S9_SPC_OBS_TROPIC  
0	0	0	0	#_30 S10_MEX            
0	0	0	0	#_31 S11_DFA_LATE       
#
#_SizeSelex
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
   35.000	  250.000	 1.57527e+02	 50	 0.000	0	  2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F1_MEX(1)                       
  -15.000	   15.000	-6.00000e+00	  0	 0.000	0	 -4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F1_MEX(1)                       
  -15.000	   15.000	 7.36579e+00	  0	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_3_F1_MEX(1)                       
  -15.000	   15.000	 6.80182e+00	  0	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_4_F1_MEX(1)                       
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_5_F1_MEX(1)                       
 -999.000	 -999.000	-9.99000e+02	  0	 5.000	0	 -2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_6_F1_MEX(1)                       
  -20.000	  200.000	-6.68529e+00	125	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_1_F1_MEX(1)                 
  -15.000	   15.000	 1.23966e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_2_F1_MEX(1)                 
  -15.000	   15.000	 2.20324e-02	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_3_F1_MEX(1)                 
  -15.000	   15.000	 0.00000e+00	  4	50.000	0	 -4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_4_F1_MEX(1)                 
    0.000	    5.000	 3.36280e-01	  4	50.000	0	  5	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_5_F1_MEX(1)                 
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F2_CAN(2)                       
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F2_CAN(2)                       
   35.000	  250.000	 1.71126e+02	 50	 0.000	0	  2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F3_CHINA(3)                     
  -15.000	   15.000	-6.00000e+00	  0	 0.000	0	 -4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F3_CHINA(3)                     
  -15.000	   15.000	 7.28940e+00	  0	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_3_F3_CHINA(3)                     
  -15.000	   15.000	 7.81161e+00	  0	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_4_F3_CHINA(3)                     
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_5_F3_CHINA(3)                     
 -999.000	 -999.000	-9.99000e+02	  0	 5.000	0	 -2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_6_F3_CHINA(3)                     
  -20.000	  200.000	-9.03685e+00	125	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_1_F3_CHINA(3)               
  -15.000	   15.000	 1.94023e-02	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_2_F3_CHINA(3)               
  -15.000	   15.000	-8.67034e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_3_F3_CHINA(3)               
  -15.000	   15.000	 0.00000e+00	  4	50.000	0	 -4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_4_F3_CHINA(3)               
    0.000	    5.000	 5.81209e-01	  4	50.000	0	  5	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_5_F3_CHINA(3)               
   35.000	  250.000	 1.51059e+02	 50	 0.000	0	  2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F4_JPN_KK_SH(4)                 
  -15.000	   15.000	-6.00000e+00	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F4_JPN_KK_SH(4)                 
  -15.000	   15.000	 6.86811e+00	  0	 0.000	0	  3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_3_F4_JPN_KK_SH(4)                 
  -15.000	   15.000	 7.84481e+00	  0	 0.000	0	  3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_4_F4_JPN_KK_SH(4)                 
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_5_F4_JPN_KK_SH(4)                 
 -999.000	 -999.000	-9.99000e+02	  0	 5.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_6_F4_JPN_KK_SH(4)                 
  -20.000	  200.000	-1.63523e+00	  0	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_1_F4_JPN_KK_SH(4)           
  -15.000	   15.000	 2.60705e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_2_F4_JPN_KK_SH(4)           
  -15.000	   15.000	-1.85740e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_3_F4_JPN_KK_SH(4)           
  -15.000	   15.000	 0.00000e+00	  4	50.000	0	 -4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_4_F4_JPN_KK_SH(4)           
    0.000	    5.000	 2.73030e-01	  4	50.000	0	  5	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_5_F4_JPN_KK_SH(4)           
   35.000	  250.000	 1.72299e+02	 50	 0.000	0	  2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F5_JPN_KK_DP(5)                 
  -15.000	   15.000	-6.00000e+00	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F5_JPN_KK_DP(5)                 
  -15.000	   15.000	 6.31504e+00	  0	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_3_F5_JPN_KK_DP(5)                 
  -15.000	   15.000	 6.18229e+00	  0	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_4_F5_JPN_KK_DP(5)                 
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_5_F5_JPN_KK_DP(5)                 
 -999.000	 -999.000	-9.99000e+02	  0	 5.000	0	 -2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_6_F5_JPN_KK_DP(5)                 
  -80.000	  200.000	-1.53097e+01	  0	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_1_F5_JPN_KK_DP(5)           
  -15.000	   15.000	-6.92640e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_2_F5_JPN_KK_DP(5)           
  -15.000	   15.000	 2.26300e-02	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_3_F5_JPN_KK_DP(5)           
  -15.000	   15.000	 0.00000e+00	  4	50.000	0	 -4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_4_F5_JPN_KK_DP(5)           
    0.000	    5.000	 1.47135e-01	  4	50.000	0	  5	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_5_F5_JPN_KK_DP(5)           
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F6_JPN_ENY_SHL(6)               
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F6_JPN_ENY_SHL(6)               
   35.000	  250.000	 1.73991e+02	 50	 0.000	0	  2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F7_JPN_ENY_DP(7)                
  -15.000	   15.000	-6.00000e+00	  0	 0.000	0	 -4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F7_JPN_ENY_DP(7)                
  -15.000	   15.000	 6.40144e+00	  0	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_3_F7_JPN_ENY_DP(7)                
  -15.000	   15.000	 6.73893e+00	  0	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_4_F7_JPN_ENY_DP(7)                
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_5_F7_JPN_ENY_DP(7)                
 -999.000	 -999.000	-9.99000e+02	  0	 5.000	0	 -2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_6_F7_JPN_ENY_DP(7)                
  -80.000	  200.000	-1.54564e+01	125	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_1_F7_JPN_ENY_DP(7)          
  -15.000	   15.000	-8.54843e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_2_F7_JPN_ENY_DP(7)          
  -15.000	   15.000	-5.83817e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_3_F7_JPN_ENY_DP(7)          
  -15.000	   15.000	 0.00000e+00	  4	50.000	0	 -4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_4_F7_JPN_ENY_DP(7)          
    0.000	    5.000	 5.03734e-01	  4	50.000	0	  5	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_5_F7_JPN_ENY_DP(7)          
   35.000	  250.000	 1.10416e+02	120	 0.000	0	 -2	0	0	   0	   0	0.5	2	2	#_SizeSel_P_1_F8_JPN_LG_MESH_EARLY(8)         
  -15.000	   15.000	-1.99692e+00	  0	 0.000	0	  3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F8_JPN_LG_MESH_EARLY(8)         
  -15.000	   15.000	 8.38994e+00	  5	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_3_F8_JPN_LG_MESH_EARLY(8)         
  -15.000	   15.000	 6.21242e+00	  5	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_4_F8_JPN_LG_MESH_EARLY(8)         
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_5_F8_JPN_LG_MESH_EARLY(8)         
 -999.000	 -999.000	-9.99000e+02	  0	 5.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_6_F8_JPN_LG_MESH_EARLY(8)         
  -80.000	  200.000	-4.83700e+00	125	50.000	0	 -4	0	0	   0	   0	0.0	2	2	#_SizeSel_PFemOff_1_F8_JPN_LG_MESH_EARLY(8)   
  -15.000	   15.000	-1.48476e+00	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_2_F8_JPN_LG_MESH_EARLY(8)   
  -15.000	   15.000	 6.34020e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_3_F8_JPN_LG_MESH_EARLY(8)   
  -15.000	   15.000	 0.00000e+00	  4	50.000	0	 -4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_4_F8_JPN_LG_MESH_EARLY(8)   
    0.000	    5.000	 5.45660e-01	  4	50.000	0	  5	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_5_F8_JPN_LG_MESH_EARLY(8)   
   35.000	  250.000	 1.32092e+02	120	 0.000	0	  2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F9_JPN_LG_MESH_LATE(9)          
  -15.000	   15.000	-6.00000e+00	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F9_JPN_LG_MESH_LATE(9)          
  -15.000	   15.000	 6.14067e+00	  5	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_3_F9_JPN_LG_MESH_LATE(9)          
  -15.000	   15.000	 6.85760e+00	  5	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_4_F9_JPN_LG_MESH_LATE(9)          
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_5_F9_JPN_LG_MESH_LATE(9)          
 -999.000	 -999.000	-9.99000e+02	  0	 5.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_6_F9_JPN_LG_MESH_LATE(9)          
  -80.000	  200.000	-4.83700e+00	125	50.000	0	 -4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_1_F9_JPN_LG_MESH_LATE(9)    
  -15.000	   15.000	-5.49376e-02	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_2_F9_JPN_LG_MESH_LATE(9)    
  -15.000	   15.000	-5.34531e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_3_F9_JPN_LG_MESH_LATE(9)    
  -15.000	   15.000	 0.00000e+00	  4	50.000	0	 -4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_4_F9_JPN_LG_MESH_LATE(9)    
    0.000	    5.000	 5.95938e-01	  4	50.000	0	  5	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_5_F9_JPN_LG_MESH_LATE(9)    
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F10_JPN_CST_Oth(10)             
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F10_JPN_CST_Oth(10)             
   25.000	  250.000	 2.78590e+01	120	 0.000	0	  2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F11_JPN_SM_MESH(11)             
  -15.000	   15.000	-6.00000e+00	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F11_JPN_SM_MESH(11)             
  -15.000	   15.000	-3.93768e+00	  5	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_3_F11_JPN_SM_MESH(11)             
  -15.000	   15.000	 7.77230e+00	  5	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_4_F11_JPN_SM_MESH(11)             
-1003.000	-1003.000	-1.00300e+03	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_5_F11_JPN_SM_MESH(11)             
 -999.000	 -999.000	-9.99000e+02	  0	 5.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_6_F11_JPN_SM_MESH(11)             
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F12_IATTC(12)                   
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F12_IATTC(12)                   
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F13_KOREA(13)                   
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F13_KOREA(13)                   
   35.000	  250.000	 1.67342e+02	 50	 0.000	0	  2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F14_NON_ISC(14)                 
  -15.000	   15.000	-6.00000e+00	  0	 0.000	0	 -4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F14_NON_ISC(14)                 
  -15.000	   15.000	 7.12531e+00	  0	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_3_F14_NON_ISC(14)                 
  -15.000	   15.000	 6.36943e+00	  0	 0.000	0	  4	0	0	   0	   0	0.5	0	0	#_SizeSel_P_4_F14_NON_ISC(14)                 
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_5_F14_NON_ISC(14)                 
 -999.000	 -999.000	-9.99000e+02	  0	 5.000	0	 -2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_6_F14_NON_ISC(14)                 
  -20.000	  200.000	-6.53552e+00	125	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_1_F14_NON_ISC(14)           
  -15.000	   15.000	-9.10712e-02	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_2_F14_NON_ISC(14)           
  -15.000	   15.000	-3.35981e-02	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_3_F14_NON_ISC(14)           
  -15.000	   15.000	 0.00000e+00	  4	50.000	0	 -4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_4_F14_NON_ISC(14)           
    0.000	    5.000	 7.99889e-01	  4	50.000	0	  5	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_5_F14_NON_ISC(14)           
   28.000	  250.000	 1.05621e+02	 50	 0.000	0	  3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F15_USA_GIILL(15)               
  -15.000	   15.000	-6.00000e+00	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F15_USA_GIILL(15)               
  -15.000	   15.000	 6.86387e+00	  0	 0.000	0	  3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_3_F15_USA_GIILL(15)               
  -15.000	   15.000	 8.24143e+00	  0	 0.000	0	  3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_4_F15_USA_GIILL(15)               
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_5_F15_USA_GIILL(15)               
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_6_F15_USA_GIILL(15)               
  -20.000	  200.000	 1.42848e+00	  0	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_1_F15_USA_GIILL(15)         
  -15.000	   15.000	 3.72869e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_2_F15_USA_GIILL(15)         
  -15.000	   15.000	-1.22617e+00	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_3_F15_USA_GIILL(15)         
  -15.000	   15.000	 0.00000e+00	  4	50.000	0	 -4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_4_F15_USA_GIILL(15)         
    0.000	    5.000	 6.56839e-01	  4	50.000	0	  5	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_5_F15_USA_GIILL(15)         
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F16_USA_SPORT(16)               
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F16_USA_SPORT(16)               
   35.000	  250.000	 1.73288e+02	 50	 0.000	0	  3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F17_USA_Lonline_DP(17)          
  -15.000	   15.000	-6.00000e+00	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F17_USA_Lonline_DP(17)          
  -15.000	   15.000	 6.51678e+00	  0	 0.000	0	  3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_3_F17_USA_Lonline_DP(17)          
  -15.000	   15.000	 8.92744e+00	  0	 0.000	0	  3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_4_F17_USA_Lonline_DP(17)          
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_5_F17_USA_Lonline_DP(17)          
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_6_F17_USA_Lonline_DP(17)          
  -80.000	  200.000	-1.86922e+01	  0	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_1_F17_USA_Lonline_DP(17)    
  -15.000	   15.000	-9.13869e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_2_F17_USA_Lonline_DP(17)    
  -15.000	   15.000	-2.03209e+00	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_3_F17_USA_Lonline_DP(17)    
  -15.000	   15.000	 0.00000e+00	  4	50.000	0	 -4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_4_F17_USA_Lonline_DP(17)    
    0.000	    5.000	 8.01517e-01	  4	50.000	0	  5	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_5_F17_USA_Lonline_DP(17)    
    0.000	    2.000	 2.00000e+00	  0	 1.000	0	-99	0	0	   0	   0	0.0	0	0	#_SizeSel_Spline_Code_F18_USA_Lonline_SH(18)  
   -0.001	    1.000	 2.35778e-01	  0	 0.001	1	  3	0	0	   0	   0	0.0	0	0	#_SizeSel_Spline_GradLo_F18_USA_Lonline_SH(18)
   -1.000	    0.001	-7.95877e-02	  0	 0.001	1	  3	0	0	   0	   0	0.0	0	0	#_SizeSel_Spline_GradHi_F18_USA_Lonline_SH(18)
    5.000	  300.000	 5.07096e+01	152	 1.000	0	-99	0	0	   0	   0	0.0	0	0	#_SizeSel_Spline_Knot_1_F18_USA_Lonline_SH(18)
    5.000	  300.000	 7.44584e+01	152	 1.000	0	-99	0	0	   0	   0	0.0	0	0	#_SizeSel_Spline_Knot_2_F18_USA_Lonline_SH(18)
    5.000	  300.000	 1.09805e+02	152	 1.000	0	-99	0	0	   0	   0	0.0	0	0	#_SizeSel_Spline_Knot_3_F18_USA_Lonline_SH(18)
    5.000	  300.000	 1.44765e+02	152	 1.000	0	-99	0	0	   0	   0	0.0	0	0	#_SizeSel_Spline_Knot_4_F18_USA_Lonline_SH(18)
    5.000	  300.000	 1.92065e+02	152	 1.000	0	-99	0	0	   0	   0	0.0	0	0	#_SizeSel_Spline_Knot_5_F18_USA_Lonline_SH(18)
   -9.000	    7.000	-1.82180e+00	  0	 0.001	1	  2	0	0	   0	   0	0.0	0	0	#_SizeSel_Spine_Val_1_F18_USA_Lonline_SH(18)  
   -9.000	    7.000	-1.21986e+00	  0	 0.001	1	  2	0	0	   0	   0	0.0	0	0	#_SizeSel_Spine_Val_2_F18_USA_Lonline_SH(18)  
   -9.000	    7.000	-1.00000e+00	  0	 1.000	0	-99	0	0	   0	   0	0.0	0	0	#_SizeSel_Spine_Val_3_F18_USA_Lonline_SH(18)  
   -9.000	    7.000	-6.03609e-01	  0	 0.001	1	  2	0	0	   0	   0	0.0	0	0	#_SizeSel_Spine_Val_4_F18_USA_Lonline_SH(18)  
   -9.000	    7.000	-1.46308e+00	  0	 0.001	1	  2	0	0	   0	   0	0.0	0	0	#_SizeSel_Spine_Val_5_F18_USA_Lonline_SH(18)  
  180.000	  190.000	 1.85000e+02	  0	 0.000	0	 -1	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_1_F18_USA_Lonline_SH(18)    
   -1.000	    1.000	 0.00000e+00	  0	 0.000	0	 -1	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_2_F18_USA_Lonline_SH(18)    
   -9.000	    9.000	 1.14940e-01	  0	 0.000	0	  3	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_3_F18_USA_Lonline_SH(18)    
   -9.000	    9.000	-5.43291e+00	  0	 0.000	0	  3	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_4_F18_USA_Lonline_SH(18)    
   35.000	  297.000	 2.13088e+02	 50	 0.000	0	  2	0	1	2004	2020	0.5	0	0	#_SizeSel_P_1_F19_TAIW_LG(19)                 
  -15.000	   15.000	-6.00000e+00	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F19_TAIW_LG(19)                 
  -15.000	   15.000	 7.55274e+00	  0	 0.000	0	  3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_3_F19_TAIW_LG(19)                 
  -15.000	   15.000	 6.83243e+00	  0	 0.000	0	  3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_4_F19_TAIW_LG(19)                 
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_5_F19_TAIW_LG(19)                 
 -999.000	 -999.000	-9.99000e+02	  0	 5.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_6_F19_TAIW_LG(19)                 
 -199.000	  200.000	-3.27503e+01	  9	50.000	0	  4	0	1	2004	2020	0.0	0	0	#_SizeSel_PFemOff_1_F19_TAIW_LG(19)           
  -15.000	   15.000	-7.13663e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_2_F19_TAIW_LG(19)           
  -15.000	   15.000	 5.44714e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_3_F19_TAIW_LG(19)           
  -15.000	   15.000	 0.00000e+00	  4	50.000	0	 -4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_4_F19_TAIW_LG(19)           
    0.000	    5.000	 3.21303e-01	  4	50.000	0	  5	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_5_F19_TAIW_LG(19)           
   35.000	  250.000	 1.75816e+02	 50	 0.000	0	  2	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_F20_TAIW_SM(20)                 
  -15.000	   15.000	-6.00000e+00	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_F20_TAIW_SM(20)                 
  -15.000	   15.000	 6.81310e+00	  0	 0.000	0	  3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_3_F20_TAIW_SM(20)                 
  -15.000	   15.000	 6.58646e+00	  0	 0.000	0	  3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_4_F20_TAIW_SM(20)                 
 -999.000	 -999.000	-9.99000e+02	  0	 0.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_5_F20_TAIW_SM(20)                 
 -999.000	 -999.000	-9.99000e+02	  0	 5.000	0	 -3	0	0	   0	   0	0.5	0	0	#_SizeSel_P_6_F20_TAIW_SM(20)                 
  -80.000	  200.000	-9.41270e+00	  9	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_1_F20_TAIW_SM(20)           
  -15.000	   15.000	-3.38038e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_2_F20_TAIW_SM(20)           
  -15.000	   15.000	 1.02736e-01	  4	50.000	0	  4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_3_F20_TAIW_SM(20)           
  -15.000	   15.000	 0.00000e+00	  4	50.000	0	 -4	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_4_F20_TAIW_SM(20)           
    0.000	    5.000	 3.36895e-01	  4	50.000	0	  5	0	0	   0	   0	0.0	0	0	#_SizeSel_PFemOff_5_F20_TAIW_SM(20)           
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_S1_HW_DP(21)                    
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_S1_HW_DP(21)                    
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_S2_HW_SH(22)                    
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_S2_HW_SH(22)                    
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_S3_TAIW_LG(23)                  
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_S3_TAIW_LG(23)                  
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_S4_TAIW_SM(24)                  
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_S4_TAIW_SM(24)                  
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_S5_JPN_EARLY(25)                
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_S5_JPN_EARLY(25)                
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_S6_JPN_LATE(26)                 
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_S6_JPN_LATE(26)                 
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_S7_JPN_RTV(27)                  
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_S7_JPN_RTV(27)                  
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_S8_SPC_OBS(28)                  
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_S8_SPC_OBS(28)                  
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_S9_SPC_OBS_TROPIC(29)           
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_S9_SPC_OBS_TROPIC(29)           
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_S10_MEX(30)                     
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_S10_MEX(30)                     
   -1.000	  200.000	 1.00000e+00	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_1_S11_DFA_LATE(31)                
   -1.000	  239.000	 6.00000e+01	 50	99.000	0	-99	0	0	   0	   0	0.5	0	0	#_SizeSel_P_2_S11_DFA_LATE(31)                
#_AgeSelex
#_No age_selex_parm
# timevary selex parameters 
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
 35.0000	250.00	161.0590	120.0	 0.0	0	 4	#_SizeSel_P_1_F8_JPN_LG_MESH_EARLY(8)_BLK2repl_1973      
 35.0000	250.00	185.1140	120.0	 0.0	0	 4	#_SizeSel_P_1_F8_JPN_LG_MESH_EARLY(8)_BLK2repl_1982      
-80.0000	200.00	-59.0849	125.0	50.0	0	 4	#_SizeSel_PFemOff_1_F8_JPN_LG_MESH_EARLY(8)_BLK2repl_1973
-80.0000	200.00	 11.6317	125.0	50.0	0	 4	#_SizeSel_PFemOff_1_F8_JPN_LG_MESH_EARLY(8)_BLK2repl_1982
  0.0001	  2.00	  0.6000	  0.6	 0.5	0	-5	#_SizeSel_P_1_F19_TAIW_LG(19)_dev_se                     
 -0.9900	  0.99	  0.0000	  0.0	 0.5	0	-6	#_SizeSel_P_1_F19_TAIW_LG(19)_dev_autocorr               
  0.0001	  2.00	  0.6000	  0.6	 0.5	0	-5	#_SizeSel_PFemOff_1_F19_TAIW_LG(19)_dev_se               
 -0.9900	  0.99	  0.0000	  0.0	 0.5	0	-6	#_SizeSel_PFemOff_1_F19_TAIW_LG(19)_dev_autocorr         
# info on dev vectors created for selex parms are reported with other devs after tag parameter section
#
0 #  use 2D_AR1 selectivity(0/1):  experimental feature
#_no 2D_AR1 selex offset used
# Tag loss and Tag reporting parameters go next
0 # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# Input variance adjustments factors: 
#_Factor	Fleet	Value
    1	21	0.000000	#_Variance_adjustment_list1 
    1	23	0.136000	#_Variance_adjustment_list2 
    1	25	0.188000	#_Variance_adjustment_list3 
    1	26	0.000000	#_Variance_adjustment_list4 
    1	27	0.000000	#_Variance_adjustment_list5 
    1	29	0.056000	#_Variance_adjustment_list6 
    1	30	0.134000	#_Variance_adjustment_list7 
    1	31	0.161100	#_Variance_adjustment_list8 
    4	 1	0.002161	#_Variance_adjustment_list9 
    4	 3	0.002161	#_Variance_adjustment_list10
    4	 4	0.000652	#_Variance_adjustment_list11
    4	 5	0.002161	#_Variance_adjustment_list12
    4	 7	0.002161	#_Variance_adjustment_list13
    4	 8	0.002161	#_Variance_adjustment_list14
    4	 9	0.002161	#_Variance_adjustment_list15
    4	14	0.002161	#_Variance_adjustment_list16
    4	15	0.002161	#_Variance_adjustment_list17
    4	17	0.002161	#_Variance_adjustment_list18
    4	18	0.002161	#_Variance_adjustment_list19
    4	19	0.002161	#_Variance_adjustment_list20
    4	20	0.002161	#_Variance_adjustment_list21
    7	11	0.002161	#_Variance_adjustment_list22
-9999	 0	0.000000	#_terminator                
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 67 changes to default Lambdas (default value is 1.0)
#_like_comp	fleet	phase	value	sizefreq_method
    1	 1	1	0	1	#_Surv_F1_MEX_Phz1                                  
    1	 2	1	0	1	#_Surv_F2_CAN_Phz1                                  
    1	 3	1	0	1	#_Surv_F3_CHINA_Phz1                                
    1	 4	1	0	1	#_Surv_F4_JPN_KK_SH_Phz1                            
    1	 5	1	0	1	#_Surv_F5_JPN_KK_DP_Phz1                            
    1	 6	1	0	1	#_Surv_F6_JPN_ENY_SHL_Phz1                          
    1	 7	1	0	1	#_Surv_F7_JPN_ENY_DP_Phz1                           
    1	 8	1	0	1	#_Surv_F8_JPN_LG_MESH_EARLY_Phz1                    
    1	 9	1	0	1	#_Surv_F9_JPN_LG_MESH_LATE_Phz1                     
    1	10	1	0	1	#_Surv_F10_JPN_CST_Oth_Phz1                         
    1	11	1	0	1	#_Surv_F11_JPN_SM_MESH_Phz1                         
    1	12	1	0	1	#_Surv_F12_IATTC_Phz1                               
    1	13	1	0	1	#_Surv_F13_KOREA_Phz1                               
    1	14	1	0	1	#_Surv_F14_NON_ISC_Phz1                             
    1	15	1	0	1	#_Surv_F15_USA_GIILL_Phz1                           
    1	16	1	0	1	#_Surv_F16_USA_SPORT_Phz1                           
    1	17	1	0	1	#_Surv_F17_USA_Lonline_DP_Phz1                      
    1	18	1	0	1	#_Surv_F18_USA_Lonline_SH_Phz1                      
    1	19	1	0	1	#_Surv_F19_TAIW_LG_Phz1                             
    1	20	1	0	1	#_Surv_F20_TAIW_SM_Phz1                             
    1	21	1	0	1	#_Surv_S1_HW_DP_Phz1                                
    1	22	1	0	1	#_Surv_S2_HW_SH_Phz1                                
    1	23	1	0	1	#_Surv_S3_TAIW_LG_Phz1                              
    1	24	1	0	1	#_Surv_S4_TAIW_SM_Phz1                              
    1	25	1	1	1	#_Surv_S5_JPN_EARLY_Phz1                            
    1	26	1	1	1	#_Surv_S6_JPN_LATE_Phz1                             
    1	27	1	0	1	#_Surv_S7_JPN_RTV_Phz1                              
    1	28	1	0	1	#_Surv_S8_SPC_OBS_Phz1                              
    1	29	1	0	1	#_Surv_S9_SPC_OBS_TROPIC_Phz1                       
    1	30	1	0	1	#_Surv_S10_MEX_Phz1                                 
    1	31	1	0	1	#_Surv_S11_DFA_LATE_Phz1                            
    4	 1	1	1	0	#_length_F1_MEX_sizefreq_method_0_Phz1              
    4	 2	1	0	0	#_length_F2_CAN_sizefreq_method_0_Phz1              
    4	 3	1	1	0	#_length_F3_CHINA_sizefreq_method_0_Phz1            
    4	 4	1	1	0	#_length_F4_JPN_KK_SH_sizefreq_method_0_Phz1        
    4	 5	1	1	0	#_length_F5_JPN_KK_DP_sizefreq_method_0_Phz1        
    4	 6	1	0	0	#_length_F6_JPN_ENY_SHL_sizefreq_method_0_Phz1      
    4	 7	1	1	0	#_length_F7_JPN_ENY_DP_sizefreq_method_0_Phz1       
    4	 8	1	1	0	#_length_F8_JPN_LG_MESH_EARLY_sizefreq_method_0_Phz1
    4	 9	1	1	0	#_length_F9_JPN_LG_MESH_LATE_sizefreq_method_0_Phz1 
    4	10	1	0	0	#_length_F10_JPN_CST_Oth_sizefreq_method_0_Phz1     
    4	11	1	0	0	#_length_F11_JPN_SM_MESH_sizefreq_method_0_Phz1     
    4	12	1	0	0	#_length_F12_IATTC_sizefreq_method_0_Phz1           
    4	13	1	0	0	#_length_F13_KOREA_sizefreq_method_0_Phz1           
    4	14	1	1	0	#_length_F14_NON_ISC_sizefreq_method_0_Phz1         
    4	15	1	1	0	#_length_F15_USA_GIILL_sizefreq_method_0_Phz1       
    4	16	1	0	0	#_length_F16_USA_SPORT_sizefreq_method_0_Phz1       
    4	17	1	1	0	#_length_F17_USA_Lonline_DP_sizefreq_method_0_Phz1  
    4	18	1	1	0	#_length_F18_USA_Lonline_SH_sizefreq_method_0_Phz1  
    4	19	1	1	0	#_length_F19_TAIW_LG_sizefreq_method_0_Phz1         
    4	20	1	1	0	#_length_F20_TAIW_SM_sizefreq_method_0_Phz1         
    4	21	1	0	0	#_length_S1_HW_DP_sizefreq_method_0_Phz1            
    4	22	1	0	0	#_length_S2_HW_SH_sizefreq_method_0_Phz1            
    4	23	1	0	0	#_length_S3_TAIW_LG_sizefreq_method_0_Phz1          
    4	24	1	0	0	#_length_S4_TAIW_SM_sizefreq_method_0_Phz1          
    4	25	1	0	0	#_length_S5_JPN_EARLY_sizefreq_method_0_Phz1        
    4	26	1	0	0	#_length_S6_JPN_LATE_sizefreq_method_0_Phz1         
    4	27	1	0	0	#_length_S7_JPN_RTV_sizefreq_method_0_Phz1          
    4	28	1	0	0	#_length_S8_SPC_OBS_sizefreq_method_0_Phz1          
    4	29	1	0	0	#_length_S9_SPC_OBS_TROPIC_sizefreq_method_0_Phz1   
    4	30	1	0	0	#_length_S10_MEX_sizefreq_method_0_Phz1             
    4	31	1	0	0	#_length_S11_DFA_LATE_sizefreq_method_0_Phz1        
    6	11	1	1	1	#_SizeFreq_F11_JPN_SM_MESH_sizefreq_method_1_Phz1   
    9	 1	1	1	0	#_init_equ_catch_F1_MEX_Phz1                        
   10	 1	1	1	0	#_recrdev_Phz1                                      
   11	 1	1	0	0	#_parm_prior_F1_MEX_Phz1                            
   12	 1	1	1	0	#_parm_dev_Phz1                                     
-9999	 0	0	0	0	#_terminator                                        
#
0 # 0/1 read specs for more stddev reporting
#
999
