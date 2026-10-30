#C file created using an r4ss function
#C file write time: 2026-10-30  02:32:36
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
3 #_Nblock_Patterns
1 2 2 #_blocks_per_pattern
#_begin and end years of blocks
16 16
17 76 77 144
17 100 101 144
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
3 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=Maunder_M;_6=Age-range_Lorenzen
#_ #_Age_natmort_by sex x growthpattern
#_Age_0	Age_1	Age_2	Age_3	Age_4	Age_5	Age_6	Age_7	Age_8	Age_9	Age_10	Age_11	Age_12	Age_13	Age_14	Age_15	Age_16	Age_17	Age_18	Age_19	Age_20	Age_21	Age_22	Age_23	Age_24	Age_25	Age_26	Age_27	Age_28	Age_29	Age_30	Age_31	Age_32	Age_33	Age_34	Age_35	Age_36	Age_37	Age_38	Age_39	Age_40
0.479	0.357	0.285847	0.238641	0.205259	0.180468	0.161387	0.146302	0.134126	0.124141	0.115851	0.108901	0.103034	0.1	0.10331	0.106621	0.109931	0.113242	0.116552	0.119863	0.123173	0.126484	0.129794	0.133104	0.136415	0.139725	0.143036	0.143036	0.143036	0.143036	0.143036	0.143036	0.143036	0.143036	0.143036	0.143036	0.143036	0.143036	0.143036	0.143036	0.143036	#_natM1
8 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr;5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0 #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
2 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
1 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
0 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
2 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn
    10	      45	    22.79	  26.7812	 10	0	 -3	0	0	0	0	0.5	0	0	#_L_at_Age0_Fem_GP_1 
    40	     250	  191.628	  192.446	 10	0	 -3	0	0	0	0	0.5	0	0	#_Linf_Fem_GP_1      
   -10	      50	  10.2127	     4.57	0.8	0	 -3	0	0	0	0	0.5	0	0	#_rmax_Fem_GP_1      
  0.05	    0.25	   0.1811	0.0676951	0.8	0	 -3	0	0	0	0	0.5	0	0	#_Cessation_Fem_GP_1 
   0.1	      10	   3.5308	0.0834877	0.8	0	  3	0	0	0	0	0.5	0	0	#_CV_young_Fem_GP_1  
   0.1	      50	  9.47985	  1.81575	0.8	0	 -3	0	0	0	0	0.5	0	0	#_CV_old_Fem_GP_1    
    -3	       3	3.661e-05	3.661e-05	0.8	0	 -3	0	0	0	0	0.5	0	0	#_Wtlen_1_Fem_GP_1   
    -3	       3	  2.90182	  2.90182	0.8	0	 -3	0	0	0	0	0.5	0	0	#_Wtlen_2_Fem_GP_1   
-10000	   10000	   133.76	       55	0.8	0	 -3	0	0	0	0	0.5	0	0	#_Mat50%_Fem_GP_1    
    -3	       3	-0.111862	    -0.25	0.8	0	 -3	0	0	0	0	0.5	0	0	#_Mat_slope_Fem_GP_1 
    -3	       3	        1	        1	0.8	6	 -3	0	0	0	0	0.5	0	0	#_Eggs_alpha_Fem_GP_1
    -3	       3	        0	        0	0.8	6	 -3	0	0	0	0	0.5	0	0	#_Eggs_beta_Fem_GP_1 
   0.1	      10	        1	        1	  1	0	 -1	0	0	0	0	  0	0	0	#_CohortGrowDev      
 1e-06	0.999999	      0.5	      0.5	0.5	0	-99	0	0	0	0	  0	0	0	#_FracFemale_GP_1    
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; 2=Ricker (2 parms); 3=std_B-H(2); 4=SCAA(2);5=Hockey(3); 6=B-H_flattop(2); 7=Survival(3);8=Shepard(3);9=Ricker_Power(3);10=B-H_a,b(4)
0 # 0/1 to use steepness in initial equ recruitment calculation
0 # future feature: 0/1 to make realized sigmaR a function of SR curvature
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn # parm_name
  3	20	9.44899	9.3	 10	0	  1	0	0	0	0	0	0	0	#_SR_LN(R0)  
0.2	 1	      1	  1	0.2	0	 -4	0	0	0	0	0	0	0	#_SR_BH_steep
  0	 2	    0.6	0.6	0.8	0	 -3	0	0	0	0	0	0	0	#_SR_sigmaR  
 -5	 5	      0	  0	  1	0	 -1	0	0	0	0	0	0	0	#_SR_regime  
  0	 0	      0	  0	 99	0	-99	0	0	0	0	0	0	0	#_SR_autocorr
#_no timevary SR parameters
2 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
17 # first year of main recr_devs; early devs can preceed this era
202 # last year of main recr_devs; forecast devs start in following year
2 #_recdev phase
1 # (0/1) to read 13 advanced options
-24 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
1 #_recdev_early_phase
0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
1 #_lambda for Fcast_recr_like occurring before endyr+1
6.8 #_last_yr_nobias_adj_in_MPD; begin of ramp
62.9 #_first_yr_fullbias_adj_in_MPD; begin of plateau
190 #_last_yr_fullbias_adj_in_MPD
206.5 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
0.5664 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
0 #_period of cycles in recruitment (N parms read below)
-15 #min rec_dev
15 #max rec_dev
0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
#Fishing Mortality info
0.1 # F ballpark
160 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
5 # max F or harvest rate, depends on F_Method
5 # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
0	1	 0.00327485	0	0	0	-1	#_InitF_seas_1_flt_2LL-n-A2
0	1	1.82458e-07	0	0	0	-1	#_InitF_seas_1_flt_19OBJ-A5
#
#_Q_setup for fleets with cpue or survey data
#_fleet	link	link_info	extra_se	biasadj	float  #  fleetname
   23	1	0	0	0	1	#_LL-Survey 
-9999	0	0	0	0	0	#_terminator
#_Q_parms(if_any);Qunits_are_ln(q)
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
-15	15	9.29847	0	-1	0	-1	0	0	0	0	0	0	0	#_LnQ_base_LL-Survey(23)
#_no timevary Q parameters
#
#_size_selex_patterns
#_Pattern	Discard	Male	Special
24	0	0	 0	#_1 LL-n-A1   
24	0	0	 0	#_2 LL-n-A2   
24	0	0	 0	#_3 LL-n-A3   
24	0	0	 0	#_4 LL-n-A4   
24	0	0	 0	#_5 LL-n-A5   
24	0	0	 0	#_6 LL-n-A6   
24	0	0	 0	#_7 LL-n-A7   
15	0	0	 1	#_8 LL-w-A1   
15	0	0	 2	#_9 LL-w-A2   
15	0	0	 3	#_10 LL-w-A3  
15	0	0	 4	#_11 LL-w-A4  
15	0	0	 5	#_12 LL-w-A5  
15	0	0	 6	#_13 LL-w-A6  
15	0	0	 7	#_14 LL-w-A7  
24	0	0	 0	#_15 OBJ-A1   
24	0	0	 0	#_16 OBJ-A2   
24	0	0	 0	#_17 OBJ-A3   
24	0	0	 0	#_18 OBJ-A4   
24	0	0	 0	#_19 OBJ-A5   
 0	0	0	 0	#_20 OBJ-Disc 
15	0	0	15	#_21 NOA-A1   
15	0	0	19	#_22 NOA-A2   
24	0	0	 0	#_23 LL-Survey
#
#_age_selex_patterns
#_Pattern	Discard	Male	Special
10	0	0	0	#_1 LL-n-A1   
10	0	0	0	#_2 LL-n-A2   
10	0	0	0	#_3 LL-n-A3   
10	0	0	0	#_4 LL-n-A4   
10	0	0	0	#_5 LL-n-A5   
10	0	0	0	#_6 LL-n-A6   
10	0	0	0	#_7 LL-n-A7   
10	0	0	0	#_8 LL-w-A1   
10	0	0	0	#_9 LL-w-A2   
10	0	0	0	#_10 LL-w-A3  
10	0	0	0	#_11 LL-w-A4  
10	0	0	0	#_12 LL-w-A5  
10	0	0	0	#_13 LL-w-A6  
10	0	0	0	#_14 LL-w-A7  
10	0	0	0	#_15 OBJ-A1   
10	0	0	0	#_16 OBJ-A2   
10	0	0	0	#_17 OBJ-A3   
10	0	0	0	#_18 OBJ-A4   
10	0	0	0	#_19 OBJ-A5   
11	0	0	0	#_20 OBJ-Disc 
10	0	0	0	#_21 NOA-A1   
10	0	0	0	#_22 NOA-A2   
10	0	0	0	#_23 LL-Survey
#
#_SizeSelex
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
 20	199	  132.914	50	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_LL-n-A1(1)   
-15	 15	  -11.624	 0	0	0	-3	0	0	0	0	0.5	0	0	#_SizeSel_P_2_LL-n-A1(1)   
-15	 15	  6.25251	 0	0	0	-3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_LL-n-A1(1)   
-15	 15	  5.87983	 0	0	0	-3	0	0	0	0	0.5	0	0	#_SizeSel_P_4_LL-n-A1(1)   
-15	 15	      -15	 0	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_5_LL-n-A1(1)   
-15	 15	     -999	 0	5	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_6_LL-n-A1(1)   
 20	199	  157.996	50	0	0	 2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_LL-n-A2(2)   
-15	 15	 -11.3249	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_2_LL-n-A2(2)   
-15	 15	   7.1142	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_LL-n-A2(2)   
-15	 15	  6.30691	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_4_LL-n-A2(2)   
-15	 15	      -15	 0	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_5_LL-n-A2(2)   
-15	 15	     -999	 0	5	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_6_LL-n-A2(2)   
 20	199	  119.395	50	0	0	 2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_LL-n-A3(3)   
-15	 15	 0.987523	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_2_LL-n-A3(3)   
-15	 15	  6.15076	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_LL-n-A3(3)   
-15	 15	-0.316846	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_4_LL-n-A3(3)   
-15	 15	      -15	 0	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_5_LL-n-A3(3)   
-15	 15	     -999	 0	5	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_6_LL-n-A3(3)   
 20	199	  155.921	50	0	0	 2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_LL-n-A4(4)   
-15	 15	 0.251241	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_2_LL-n-A4(4)   
-15	 15	  7.27387	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_LL-n-A4(4)   
-15	 15	-0.266777	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_4_LL-n-A4(4)   
-15	 15	      -15	 0	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_5_LL-n-A4(4)   
-15	 15	     -999	 0	5	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_6_LL-n-A4(4)   
 20	199	  146.047	50	0	0	 2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_LL-n-A5(5)   
-15	 15	 -9.77739	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_2_LL-n-A5(5)   
-15	 15	  6.84632	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_LL-n-A5(5)   
-15	 15	  6.53631	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_4_LL-n-A5(5)   
-15	 15	      -15	 0	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_5_LL-n-A5(5)   
-15	 15	     -999	 0	5	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_6_LL-n-A5(5)   
 20	199	  146.669	50	0	0	 2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_LL-n-A6(6)   
-15	 15	 -9.40124	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_2_LL-n-A6(6)   
-15	 15	  7.08678	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_LL-n-A6(6)   
-15	 15	  6.93032	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_4_LL-n-A6(6)   
-15	 15	      -15	 0	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_5_LL-n-A6(6)   
-15	 15	     -999	 0	5	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_6_LL-n-A6(6)   
 20	199	   126.16	50	0	0	 2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_LL-n-A7(7)   
-15	 15	 -10.7872	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_2_LL-n-A7(7)   
-15	 15	  6.31751	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_LL-n-A7(7)   
-15	 15	   6.8181	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_4_LL-n-A7(7)   
-15	 15	      -15	 0	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_5_LL-n-A7(7)   
-15	 15	     -999	 0	5	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_6_LL-n-A7(7)   
 20	199	  42.8728	50	0	0	 2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_OBJ-A1(15)   
-15	 15	  -11.192	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_2_OBJ-A1(15)   
-15	 15	  4.35639	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_OBJ-A1(15)   
-15	 15	  7.22551	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_4_OBJ-A1(15)   
-15	 15	      -15	 0	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_5_OBJ-A1(15)   
-15	 15	      -15	 0	5	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_6_OBJ-A1(15)   
 20	199	  45.6124	50	0	0	 2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_OBJ-A2(16)   
-15	 15	 -11.4558	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_2_OBJ-A2(16)   
-15	 15	  4.55505	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_OBJ-A2(16)   
-15	 15	  7.64285	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_4_OBJ-A2(16)   
-15	 15	      -15	 0	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_5_OBJ-A2(16)   
-15	 15	      -15	 0	5	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_6_OBJ-A2(16)   
 20	199	  48.5703	50	0	0	 2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_OBJ-A3(17)   
-15	 15	-0.430251	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_2_OBJ-A3(17)   
-15	 15	  5.10807	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_OBJ-A3(17)   
-15	 15	  5.96244	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_4_OBJ-A3(17)   
-15	 15	      -15	 0	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_5_OBJ-A3(17)   
-15	 15	      -15	 0	5	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_6_OBJ-A3(17)   
 20	199	  51.6743	50	0	0	 2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_OBJ-A4(18)   
-15	 15	 -1.40521	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_2_OBJ-A4(18)   
-15	 15	  4.83333	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_OBJ-A4(18)   
-15	 15	  7.33142	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_4_OBJ-A4(18)   
-15	 15	      -15	 0	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_5_OBJ-A4(18)   
-15	 15	      -15	 0	5	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_6_OBJ-A4(18)   
 20	199	  77.0964	50	0	0	 2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_OBJ-A5(19)   
-15	 15	 -1.23247	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_2_OBJ-A5(19)   
-15	 15	  6.60396	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_OBJ-A5(19)   
-15	 15	  6.84867	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_4_OBJ-A5(19)   
-15	 15	      -15	 0	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_5_OBJ-A5(19)   
-15	 15	     -999	 0	5	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_6_OBJ-A5(19)   
 20	199	  148.088	50	0	0	 2	0	0	0	0	0.5	0	0	#_SizeSel_P_1_LL-Survey(23)
-15	 15	  -7.8496	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_2_LL-Survey(23)
-15	 15	  6.96933	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_3_LL-Survey(23)
-15	 15	  6.84098	 0	0	0	 3	0	0	0	0	0.5	0	0	#_SizeSel_P_4_LL-Survey(23)
-15	 15	      -15	 0	0	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_5_LL-Survey(23)
-15	 15	     -999	 0	5	0	-2	0	0	0	0	0.5	0	0	#_SizeSel_P_6_LL-Survey(23)
#_AgeSelex
1	40	1	1	99	0	-99	0	0	0	0	0.5	0	0	#_AgeSel_P_1_OBJ-Disc(20)
1	40	3	3	99	0	-99	0	0	0	0	0.5	0	0	#_AgeSel_P_2_OBJ-Disc(20)
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
    4	15	  0.2424	#_Variance_adjustment_list1 
    4	16	    0.12	#_Variance_adjustment_list2 
    4	17	   0.028	#_Variance_adjustment_list3 
    4	18	    0.12	#_Variance_adjustment_list4 
    4	19	  0.0264	#_Variance_adjustment_list5 
    7	 2	     0.8	#_Variance_adjustment_list6 
    7	 3	 0.59316	#_Variance_adjustment_list7 
    7	 4	     0.4	#_Variance_adjustment_list8 
    7	 5	     0.4	#_Variance_adjustment_list9 
    7	 6	0.310067	#_Variance_adjustment_list10
    7	 7	0.417017	#_Variance_adjustment_list11
    7	23	     0.8	#_Variance_adjustment_list12
-9999	 0	       0	#_terminator                
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 8 changes to default Lambdas (default value is 1.0)
#_like_comp	fleet	phase	value	sizefreq_method
    5	16	1	0	1	#_age_OBJ-A2_Phz1                        
    6	 1	1	0	1	#_SizeFreq_LL-n-A1_sizefreq_method_1_Phz1
    4	21	1	0	1	#_length_NOA-A1_sizefreq_method_1_Phz1   
    4	22	1	0	1	#_length_NOA-A2_sizefreq_method_1_Phz1   
    9	 2	1	0	1	#_init_equ_catch_LL-n-A2_Phz1            
    9	19	1	0	1	#_init_equ_catch_OBJ-A5_Phz1             
   17	 1	1	0	1	#_F-ballpark_LL-n-A1_Phz1                
   18	 1	1	0	1	#_Regime-shift_LL-n-A1_Phz1              
-9999	 0	0	0	0	#_terminator                             
#
0 # 0/1 read specs for more stddev reporting
#
999
