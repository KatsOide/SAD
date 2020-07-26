      module wsbb
      use tfstk, only:sad_descriptor
      use macphys, only:elradi,elmass
      implicit none
      integer*4, parameter :: nslimax=500,nblist=1600
c      real*8, parameter:: re=2.81794092d-15, am_e=510999.06
      real*8, parameter:: re=elradi, am_e=elmass

      type (sad_descriptor) kvlum
      data kvlum%k /0/

      type sbeam 
      real*8 Benv(36),Benv5(25),v_cen(5),cod(6)
      real*8 xangle(8)  !(8)
      integer*8 iax
      real*8 zslice(nslimax*2)   !(nslimax*2+2)
!      real*8, pointer :: zslice(:)   !(nslimax*2+2)
      real*8 gamma,gambet,ce,Luminosity
      integer*4 nslice,bstrl
      end type sbeam
      end module wsbb

      module kyparam
      use wsbb, only:nblist
      integer*4, parameter ::
     $ky_L_DRFT=1,
     $ky_KIN_DRFT=2,
     $ky_COUPLE_DRFT=3,
     $ky_RADI_DRFT=4,
c
     $ky_AX_DRFT=5,
     $ky_AY_DRFT=6,
     $ky_RAD_DRFT=7,
c
     $ky_MAX_DRFT=8,
cc for STeering
c                         number of parameters
c  for bend
     $ky_L_BEND=1,
     $ky_ANGL_BEND=2,
     $ky_E1_BEND=3,
     $ky_E2_BEND=4,
     $ky_ROT_BEND=5,
c     ky_TILT_BEND=6,
     $ky_RAD_BEND=7,
     $ky_K1_BEND=8,
     $ky_DX_BEND=9,
     $ky_DY_BEND=10,
     $ky_K0_BEND=11,
     $ky_FRIN_BEND=12,
     $ky_EPS_BEND=13,
     $ky_RANK_BEND=14,
     $ky_F1_BEND=15,
     $ky_FRMD_BEND=16,
     $ky_COUPLE_BEND=17,
     $ky_DROT_BEND=18,
     $ky_AE1_BEND=19,
     $ky_AE2_BEND=20,
     $ky_FB1_BEND=21,
     $ky_FB2_BEND=22,
     $ky_INDX_BEND=23,
c
     $ky_AX_BEND=24,
     $ky_AY_BEND=25,
     $ky_LDEV_BEND=26,
     $ky_LRAD_BEND=27,
c
     $ky_MAX_BEND=28,
     $     p_L_BEND=ky_MAX_BEND+1,
     $     p_PSI1_BEND=p_L_BEND+1,
     $     p_PSI2_BEND=p_PSI1_BEND+1,
     $     p_COSPSI1_BEND=p_PSI2_BEND+1,
     $     p_SINPSI1_BEND=p_COSPSI1_BEND+1,
     $     p_COSPSI2_BEND=p_SINPSI1_BEND+1,
     $     p_SINPSI2_BEND=p_COSPSI2_BEND+1,
     $     p_COSTHETA_BEND=p_SINPSI2_BEND+1,
     $     p_SINTHETA_BEND=p_COSTHETA_BEND+1,
     $     p_COSW_BEND=p_SINTHETA_BEND+1,
     $     p_SINW_BEND=p_COSW_BEND+1,
     $     p_SQWH_BEND=p_SINW_BEND+1,
     $     p_SINWP1_BEND=p_SQWH_BEND+1,
     $     p_THETA_BEND=p_SINWP1_BEND+1,
     $     p_FB1_BEND=p_THETA_BEND+1,
     $     p_FB2_BEND=p_FB1_BEND+1,
     $     p_NPARAM_BEND=p_FB2_BEND-ky_MAX_BEND,
c  for quad
     $ky_L_QUAD=1,
     $ky_K1_QUAD=2,
c     ky_DK1_QUAD=3,
     $ky_ROT_QUAD=4,
     $ky_DX_QUAD=5,
     $ky_DY_QUAD=6,
     $ky_RAD_QUAD=7,
     $ky_CHRO_QUAD=8,
     $ky_FRIN_QUAD=9,
     $ky_F1_QUAD=10,
     $ky_F2_QUAD=11,
     $ky_FRMD_QUAD=12,
     $ky_EPS_QUAD=13,
     $ky_KIN_QUAD=14,
     $ky_COUPLE_QUAD=15,
     $ky_INDX_QUAD=16,
     $ky_AX_QUAD=17,
     $ky_AY_QUAD=18,
     $ky_LDEV_QUAD=19,
     $ky_LRAD_QUAD=20,
     $ky_F1K1F_QUAD=21,
     $ky_F2K1F_QUAD=22,
     $ky_F1K1B_QUAD=23,
     $ky_F2K1B_QUAD=24,
     $ky_MAX_QUAD=25,
c     $     p_SQRTK_QUAD=ky_MAX_QUAD+1,
c     $     p_COSTHETA_QUAD=p_SQRTK_QUAD+1,
c     $     p_SINTHETA_QUAD=p_COSTHETA_QUAD+1,
     $     p_THETA2_QUAD=ky_MAX_QUAD+1,
     $     p_AKF1F_QUAD=p_THETA2_QUAD+1,
     $     p_AKF2F_QUAD=p_AKF1F_QUAD+1,
     $     p_AKF1B_QUAD=p_AKF2F_QUAD+1,
     $     p_AKF2B_QUAD=p_AKF1B_QUAD+1,
     $     p_FRMD_QUAD=p_AKF2B_QUAD+1,
     $     p_NPARAM_QUAD=p_FRMD_QUAD-ky_MAX_QUAD,
c  for THIN
     $ky_L_THIN=1,
     $ky_K_THIN=2,
     $ky_ROT_THIN=4,
     $ky_DX_THIN=5,
     $ky_DY_THIN=6,
     $ky_RAD_THIN=7,
     $ky_FRIN_THIN=8,
     $ky_COUPLE_THIN=9,
     $ky_INDX_THIN=10,
c
     $ky_AX_THIN=11,
     $ky_AY_THIN=12,
     $ky_LDEV_THIN=13,
     $ky_LRAD_THIN=14,
c
     $ky_MAX_THIN=15,
     $     p_COSTHETA_THIN=ky_MAX_THIN+1,
     $     p_SINTHETA_THIN=p_COSTHETA_THIN+1,
     $     p_THETA_THIN=p_SINTHETA_THIN+1,
     $     p_NPARAM_THIN=p_THETA_THIN-ky_MAX_THIN,
c  for MULT
     $ky_L_MULT=1,
     $ky_DX_MULT=3,
     $ky_DY_MULT=4,
     $ky_DZ_MULT=5,
     $ky_CHI1_MULT=6,
     $ky_CHI2_MULT=7,
     $ky_CHI3_MULT=8,
     $ky_ROT_MULT=8,
     $ky_EPS_MULT=9,
     $ky_RAD_MULT=10,
     $ky_FRIN_MULT=11,
     $ky_F1_MULT=12,
     $ky_F2_MULT=13,
     $ky_FRMD_MULT=14,
     $ky_VOLT_MULT=15,
     $ky_HARM_MULT=16,
     $ky_PHI_MULT=17,
     $ky_FREQ_MULT=18,
     $ky_COUPLE_MULT=19,
     $ky_RADI_MULT=20,
     $ky_DPHI_MULT=21,
     $ky_W1_MULT=22,
c
     $ky_ANGL_MULT=23,
     $ky_E1_MULT=24,
     $ky_E2_MULT=25,
     $ky_DROT_MULT=26,
     $ky_K0FR_MULT=27,
     $ky_FB1_MULT=28,
     $ky_FB2_MULT=29,
c
     $ky_K0_MULT=30,
     $ky_SK0_MULT=31,
     $ky_K1_MULT=32,
     $ky_SK1_MULT=33,
     $ky_K2_MULT=34,
     $ky_SK2_MULT=35,
     $ky_K3_MULT=36,
     $ky_SK3_MULT=37,
     $ky_K4_MULT=38,
     $ky_SK4_MULT=39,
     $ky_K5_MULT=40,
     $ky_SK5_MULT=41,
     $ky_K6_MULT=42,
     $ky_SK6_MULT=43,
     $ky_K7_MULT=44,
     $ky_SK7_MULT=45,
     $ky_K8_MULT=46,
     $ky_SK8_MULT=47,
     $ky_K9_MULT=48,
     $ky_SK9_MULT=49,
     $ky_K10_MULT=50,
     $ky_SK10_MULT=51,
     $ky_K11_MULT=52,
     $ky_SK11_MULT=53,
     $ky_K12_MULT=54,
     $ky_SK12_MULT=55,
     $ky_K13_MULT=56,
     $ky_SK13_MULT=57,
     $ky_K14_MULT=58,
     $ky_SK14_MULT=59,
     $ky_K15_MULT=60,
     $ky_SK15_MULT=61,
     $ky_K16_MULT=62,
     $ky_SK16_MULT=63,
     $ky_K17_MULT=64,
     $ky_SK17_MULT=65,
     $ky_K18_MULT=66,
     $ky_SK18_MULT=67,
     $ky_K19_MULT=68,
     $ky_SK19_MULT=69,
     $ky_K20_MULT=70,
     $ky_SK20_MULT=71,
     $ky_K21_MULT=72,
     $ky_SK21_MULT=73,
     $ky_AE1_MULT=74,
     $ky_AE2_MULT=75,
c
     $ky_AX_MULT=76,
     $ky_AY_MULT=77,
     $ky_LDEV_MULT=78,
     $ky_LRAD_MULT=79,
     $ky_APHI_MULT=80,
     $ky_F1K1F_MULT=81,
     $ky_F2K1F_MULT=82,
     $ky_F1K1B_MULT=83,
     $ky_F2K1B_MULT=84,
     $ky_DVOLT_MULT=85,
     $     ky_PROF_MULT=86,
c
     $ky_MAX_MULT=ky_PROF_MULT+1,
     $     p_L_MULT=ky_MAX_MULT+1,
     $     p_ANGL_MULT=p_L_MULT+1,
     $     p_AKF1F_MULT=p_ANGL_MULT+1,
     $     p_AKF2F_MULT=p_AKF1F_MULT+1,
     $     p_AKF1B_MULT=p_AKF2F_MULT+1,
     $     p_AKF2B_MULT=p_AKF1B_MULT+1,
     $     p_PSI1_MULT=p_AKF2B_MULT+1,
     $     p_PSI2_MULT=p_PSI1_MULT+1,
     $     p_FB1_MULT=p_PSI2_MULT+1,
     $     p_FB2_MULT=p_FB1_MULT+1,
     $     p_CHI1_MULT=p_FB2_MULT+1,
     $     p_CHI2_MULT=p_CHI1_MULT+1,
     $     p_W_MULT=p_CHI2_MULT+1,
     $     p_VNOMINAL_MULT=p_W_MULT+1,
     $     p_FRMD_MULT=p_VNOMINAL_MULT+1,
     $     p_THETA2_MULT=p_FRMD_MULT+1,
     $     p_CR1_MULT=p_THETA2_MULT+1,
     $     p_CR1I_MULT=p_CR1_MULT+1,
     $     p_PROF_MULT=p_CR1I_MULT+1,
     $     p_NPARAM_MULT=p_PROF_MULT-ky_MAX_MULT,
c  for UNDULATOR
     $ky_L_UND=1,
     $ky_FBX_UND=2,
     $ky_FBY_UND=3,
     $ky_Kx_UND=4,
     $ky_Qy_UND=5,
     $ky_PHI_UND=6,
     $ky_SLI_UND=7,
     $ky_Pole_UND=8,
     $ky_COUPLE_UND=9,
     $ky_MAX_UND=10,
     $     p_PARAM_UND=ky_MAX_UND+1,
     $     p_NPARAM_UND=20,
c  for WIG
     $ky_L_WIG=1,
     $ky_BMAX_WIG=2,
     $ky_PRD_WIG=3,
     $ky_ROT_WIG=4,
     $ky_DX_WIG=5,
     $ky_DY_WIG=6,
     $ky_A3_WIG=7,
     $ky_A5_WIG=8,
     $ky_A7_WIG=9,
     $ky_A9_WIG=10,
     $ky_A11_WIG=11,
     $ky_A13_WIG=12,
     $ky_A15_WIG=13,
     $ky_A17_WIG=14,
     $ky_COUPLE_WIG=15,
     $ky_MAX_WIG=16,
     $     p_PARAM_WIG=ky_MAX_WIG+1,
c  for solenoid
     $ky_L_SOL=1,
     $ky_BZ_SOL=2,
     $ky_DX_SOL=3,
     $ky_DY_SOL=4,
     $ky_DZ_SOL=5,
     $ky_DPX_SOL=6,
     $ky_DPY_SOL=7,
     $ky_BND_SOL=8,
     $ky_CHI1_SOL=9,
     $ky_CHI2_SOL=10,
     $ky_CHI3_SOL=11,
     $ky_GEO_SOL=12,
     $ky_F1_SOL=13,
     $ky_DBZ_SOL=14,
     $ky_COUPLE_SOL=18,
     $ky_FRIN_SOL=19,
     $ky_RAD_SOL=20,
     $ky_FL_SOL=21,
     $ky_MAX_SOL=22,
     $     p_R11_SOL=ky_MAX_SOL+1,
     $     p_R12_SOL=p_R11_SOL+1,
     $     p_R13_SOL=p_R12_SOL+1,
     $     p_R21_SOL=p_R13_SOL+1,
     $     p_R22_SOL=p_R21_SOL+1,
     $     p_R23_SOL=p_R22_SOL+1,
     $     p_R31_SOL=p_R23_SOL+1,
     $     p_R32_SOL=p_R31_SOL+1,
     $     p_R33_SOL=p_R32_SOL+1,
     $     p_NPARAM_SOL=p_R33_SOL-ky_MAX_SOL,
cc for TEST
     $ky_L_TEST=1,
     $ky_ANGL_TEST=2,
     $ky_MAX_TEST=3,
cc for Cavity
     $ky_L_CAVI=1,
     $ky_VOLT_CAVI=2,
     $ky_HARM_CAVI=3,
     $ky_PHI_CAVI=4,
     $ky_FREQ_CAVI=5,
     $ky_RANV_CAVI=9,
     $ky_RANP_CAVI=10,
     $ky_LWAK_CAVI=11,
     $ky_TWAK_CAVI=12,
     $ky_DX_CAVI=13,
     $ky_DY_CAVI=14,
     $ky_ROT_CAVI=15,
     $ky_V1_CAVI=16,
     $ky_V20_CAVI=17,
     $ky_V11_CAVI=18,
     $ky_V02_CAVI=19,
     $ky_COUPLE_CAVI=20,
     $ky_DPHI_CAVI=21,
     $ky_FRIN_CAVI=22,
     $ky_FRMD_CAVI=23,
c
     $ky_AX_CAVI=24,
     $ky_AY_CAVI=25,
     $ky_LDEV_CAVI=26,
     $ky_APHI_CAVI=27,
     $ky_DVOLT_CAVI=28,
     $ky_MAX_CAVI=29,
     $     p_W_CAVI=ky_MAX_CAVI+1,
     $     p_VNOMINAL_CAVI=p_W_CAVI+1,
     $     p_FRMD_CAVI=p_VNOMINAL_CAVI+1,
     $     p_NPARAM_CAVI=p_FRMD_CAVI-ky_MAX_CAVI,
cc for t-Cavity
     $ky_L_TCAV=1,
     $ky_K0_TCAV=2,
     $ky_HARM_TCAV=3,
     $ky_PHI_TCAV=4,
     $ky_FREQ_TCAV=5,
     $ky_DX_TCAV=6,
     $ky_DY_TCAV=7,
     $ky_ROT_TCAV=8,
     $ky_RANK_TCAV=9,
     $ky_RANP_TCAV=10,
     $ky_LWAK_TCAV=11,
     $ky_TWAK_TCAV=12,
     $ky_COUPLE_TCAV=13,
c
     $ky_AX_TCAV=14,
     $ky_AY_TCAV=15,
     $ky_LDEV_TCAV=16,
     $ky_RAD_TCAV=17,
c
     $ky_MAX_TCAV=18,
cc for MAP
     $ky_L_MAP=1,
     $ky_MAX_MAP=2,
cc for INS
     $ky_AX_INS=1,
     $ky_BX_INS=2,
     $ky_PX_INS=3,
     $ky_AY_INS=4,
     $ky_BY_INS=5,
     $ky_PY_INS=6,
     $ky_EX_INS=7,
     $ky_EPX_INS=8,
     $ky_EY_INS=9,
     $ky_EPY_INS=10,
     $ky_R1_INS=11,
     $ky_R2_INS=12,
     $ky_R3_INS=13,
     $ky_R4_INS=14,
     $ky_DX_INS=15,
     $ky_DPX_INS=16,
     $ky_DY_INS=17,
     $ky_DPY_INS=18,
     $ky_DIR_INS=19,
     $ky_COUPLE_INS=62,
     $ky_MAX_INS=63,
cc for Coordinate transformation
     $ky_DX_COORD=1,
     $ky_DY_COORD=2,
     $ky_DZ_COORD=3,
     $ky_CHI1_COORD=4,
     $ky_CHI2_COORD=5,
     $ky_CHI3_COORD=6,
     $ky_DIR_COORD=7,
     $ky_COUPLE_COORD=8,
     $ky_MAX_COORD=9,
cc for BEAMBEAM
     $ky_AX_BEAM=1,
     $ky_BX_BEAM=2,
     $ky_AY_BEAM=3,
     $ky_BY_BEAM=4,
     $ky_R1_BEAM=5,
     $ky_R2_BEAM=6,
     $ky_R3_BEAM=7,
     $ky_R4_BEAM=8,
     $ky_EX_BEAM=9,
     $ky_EPX_BEAM=10,
     $ky_EY_BEAM=11,
     $ky_EPY_BEAM=12,
     $ky_ZX_BEAM=13,
     $ky_ZPX_BEAM=14,
     $ky_ZY_BEAM=15,
     $ky_ZPY_BEAM=16,
     $ky_DX_BEAM=17,
     $ky_DPX_BEAM=18,
     $ky_DY_BEAM=19,
     $ky_DPY_BEAM=20,
     $ky_XANGLE_BEAM=21,
     $ky_EMIX_BEAM=22,
     $ky_EMIY_BEAM=23,
     $ky_DP_BEAM=24,
     $ky_AZ_BEAM=25,
     $ky_DZ_BEAM=26,
     $ky_SIGZ_BEAM=27,
     $ky_SLI_BEAM=28,
     $ky_NP_BEAM=29,
     $ky_STURN_BEAM=30,
     $ky_R11_BEAM=31,
     $ky_R12_BEAM=32,
     $ky_R13_BEAM=33,
     $ky_R14_BEAM=34,
     $ky_R15_BEAM=35,
     $ky_R16_BEAM=36,
     $ky_R22_BEAM=37,
     $ky_R23_BEAM=38,
     $ky_R24_BEAM=39,
     $ky_R25_BEAM=40,
     $ky_R26_BEAM=41,
     $ky_R33_BEAM=42,
     $ky_R34_BEAM=43,
     $ky_R35_BEAM=44,
     $ky_R36_BEAM=45,
     $ky_R44_BEAM=46,
     $ky_R45_BEAM=47,
     $ky_R46_BEAM=48,
     $ky_R55_BEAM=49,
     $ky_R56_BEAM=50,
     $ky_R66_BEAM=51,
     $ky_COUPLE_BEAM=52,
     $ky_BSTRL_BEAM=53,
     $ky_MAX_BEAM=54,
     $     p_PARAM_BEAM=ky_MAX_BEAM+1,
     $     p_NPARAM_BEAM=nblist,
cc for PHSROT
     $ky_AX_Prot=1,
     $ky_BX_Prot=2,
     $ky_PX_Prot=3,
     $ky_AY_Prot=4,
     $ky_BY_Prot=5,
     $ky_PY_Prot=6,
     $ky_R1_Prot=7,
     $ky_R2_Prot=8,
     $ky_R3_Prot=9,
     $ky_R4_Prot=10,
     $ky_EX_Prot=11,
     $ky_EPX_Prot=12,
     $ky_EY_Prot=13,
     $ky_EPY_Prot=14,
     $ky_ZX_Prot=15,
     $ky_ZPX_Prot=16,
     $ky_ZY_Prot=17,
     $ky_ZPY_Prot=18,
     $ky_EMIX_Prot=19,
     $ky_EMIY_Prot=20,
     $ky_DP_Prot=21,
     $ky_AZ_Prot=22,
     $ky_DZ_Prot=23,
     $ky_SIGZ_Prot=24,
     $ky_PZ_Prot=25,
     $ky_EMIZ_Prot=26,
     $ky_BZ_Prot=27,
     $ky_JDY_Prot=28,
     $ky_JDPY_Prot=29,
     $ky_D11_Prot=30,
     $ky_D12_Prot=31,
     $ky_D13_Prot=32,
     $ky_D14_Prot=33,
     $ky_D15_Prot=34,
     $ky_D16_Prot=35,
     $ky_D21_Prot=36,
     $ky_D22_Prot=37,
     $ky_D23_Prot=38,
     $ky_D24_Prot=39,
     $ky_D25_Prot=40,
     $ky_D26_Prot=41,
     $ky_D31_Prot=42,
     $ky_D32_Prot=43,
     $ky_D33_Prot=44,
     $ky_D34_Prot=45,
     $ky_D35_Prot=46,
     $ky_D36_Prot=47,
     $ky_D41_Prot=48,
     $ky_D42_Prot=49,
     $ky_D43_Prot=50,
     $ky_D44_Prot=51,
     $ky_D45_Prot=52,
     $ky_D46_Prot=53,
     $ky_D51_Prot=54,
     $ky_D52_Prot=55,
     $ky_D53_Prot=56,
     $ky_D54_Prot=57,
     $ky_D55_Prot=58,
     $ky_D56_Prot=59,
     $ky_D61_Prot=60,
     $ky_D62_Prot=61,
     $ky_D63_Prot=62,
     $ky_D64_Prot=63,
     $ky_D65_Prot=64,
     $ky_D66_Prot=65,
     $ky_B11_Prot=66,
     $ky_B12_Prot=67,
     $ky_B13_Prot=68,
     $ky_B14_Prot=69,
     $ky_B15_Prot=70,
     $ky_B16_Prot=71,
     $ky_B22_Prot=72,
     $ky_B23_Prot=73,
     $ky_B24_Prot=74,
     $ky_B25_Prot=75,
     $ky_B26_Prot=76,
     $ky_B33_Prot=77,
     $ky_B34_Prot=78,
     $ky_B35_Prot=79,
     $ky_B36_Prot=80,
     $ky_B44_Prot=81,
     $ky_B45_Prot=82,
     $ky_B46_Prot=83,
     $ky_B55_Prot=84,
     $ky_B56_Prot=85,
     $ky_B66_Prot=86,
     $ky_COUPLE_Prot=87,
     $ky_MAX_Prot=88,
     $     p_PARAM_Prot=ky_MAX_Prot+1,
     $     p_NPARAM_Prot=240,
cc for MARK
     $ky_AX_Mark=1,
     $ky_BX_Mark=2,
     $ky_PX_Mark=3,
     $ky_AY_Mark=4,
     $ky_BY_Mark=5,
     $ky_PY_Mark=6,
     $ky_EX_Mark=7,
     $ky_EPX_Mark=8,
     $ky_EY_Mark=9,
     $ky_EPY_Mark=10,
     $ky_R1_Mark=11,
     $ky_R2_Mark=12,
     $ky_R3_Mark=13,
     $ky_R4_Mark=14,
     $ky_DETR_Mark=15,
     $ky_DX_Mark=16,
     $ky_DPX_Mark=17,
     $ky_DY_Mark=18,
     $ky_DPY_Mark=19,
     $ky_DZ_Mark=20,
     $ky_DDP_Mark=21,
     $ky_AZ_Mark=22,
     $ky_BZ_Mark=23,
     $ky_PZ_Mark=24,
     $ky_ZX_Mark=25,
     $ky_ZPX_Mark=26,
     $ky_ZY_Mark=27,
     $ky_ZPY_Mark=28,
     $ky_DP_Mark=29,
     $ky_OFFSET_Mark=30,
     $ky_SIGZ_Mark=31,
     $ky_SIGE_Mark=32,
     $ky_GEO_Mark=33,
     $ky_JDX_Mark=34,
     $ky_JDPX_Mark=35,
     $ky_JDY_Mark=36,
     $ky_JDPY_Mark=37,
     $ky_JDZ_Mark=38,
     $ky_JDPZ_Mark=39,
     $ky_EMIX_Mark=40,
     $ky_EMIY_Mark=41,
     $ky_EMIZ_Mark=42,
     $ky_COUPLE_Mark=43,
     $ky_MAX_Mark=44,
cc for apert
     $ky_DX1_Aprt=1,
     $ky_DX2_Aprt=2,
     $ky_DY1_Aprt=3,
     $ky_DY2_Aprt=4,
     $ky_JDPX_Aprt=5,
     $ky_JDPY_Aprt=6,
     $ky_DP_Aprt=7,
     $ky_COUPLE_Aprt=8,
     $ky_ROT_Aprt=9,
     $ky_AX_Aprt=10,
     $ky_AY_Aprt=11,
     $ky_DX_Aprt=12,
     $ky_DY_Aprt=13,      
     $ky_MAX_Aprt=14,
cc for mon
     $ky_DX_MONI=1,
     $ky_DY_MONI=2,
     $ky_OFFSET_MONI=3,
     $ky_COUPLE_MONI=4,
     $ky_ROT_MONI=5,
     $ky_MAX_MONI=6,
cc for spch
     $ky_AX_SPCH=1,
     $ky_BX_SPCH=2,
     $ky_PX_SPCH=3,
     $ky_AY_SPCH=4,
     $ky_BY_SPCH=5,
     $ky_PY_SPCH=6,
     $ky_R1_SPCH=7,
     $ky_R2_SPCH=8,
     $ky_R3_SPCH=9,
     $ky_R4_SPCH=10,
     $ky_EX_SPCH=11,
     $ky_EPX_SPCH=12,
     $ky_EY_SPCH=13,
     $ky_EPY_SPCH=14,
     $ky_ZX_SPCH=15,
     $ky_ZPX_SPCH=16,
     $ky_ZY_SPCH=17,
     $ky_ZPY_SPCH=18,
     $ky_MAX_SPCH=19

      contains
        logical*4 pure elemental function integv(k,ic)
        use maccode
        implicit none
        integer*4 ,intent(in):: k,ic
        select case (ic)
        case (icMULT)
          integv=k. eq. ky_L_MULT .or. k .eq. ky_ANGL_MULT .or.
     $         k .eq. ky_VOLT_MULT .or. k .eq. ky_DVOLT_MULT .or.
     $         k .ge. ky_K0_MULT .and. k .le. ky_SK21_MULT
        case default
          integv=.false.
        end select
        return
        end function
      end module

      subroutine initb1
      use kyparam
      use maccbk
      use mackw
      use macphys
      use macvar
      use macfile
      use tfmem, only:ktaloc
      use tfmem, only:maxstack
      use tfstk, only:dinfinity,dnotanumber
      implicit none
c
      integer*4 idummy,idummy1,hsrch,i
      integer*8 ktcaloc
c      character*132 stacksiz
c     external doline
      external doprin, doexpn, doread, dolist, docod, dostop, dotwis
      external dooffl, doonfl,dorvrs
      external ActLie,ActTra,ActPlt,ActGRA

      allocate(kytbl(0:kwMAX,0:icMXEL))
      kytbl=0

       call defglb('$PLOT$',icGLI,idummy)
       call IsetGL('$PLOT$',0,idummy)
       idummy=sethtb('PLOT   ',icACT,ktaloc(8))
       ilist(1,idval(idummy))=7
       call setfnp(klist(idval(idummy)+1),ActPlt)
       ilist(1,idval(idummy)+2)=hsrch('PTYPE')
       ilist(1,idval(idummy)+3)=hsrch('NPART')
       ilist(1,idval(idummy)+4)=hsrch('TURNS')
       ilist(1,idval(idummy)+5)=hsrch('SPAN    ')
       ilist(1,idval(idummy)+6)=hsrch('CENTER  ')
c
       idummy=sethtb('ID      ',icVAR,VarPt+VarLst)
       idummy1=sethtb('GTYPE   ',icVAR,VarStr)
c
c      idummy=sethtb('DRAW    ',icACT,mtaloc(8))
c      idummy=sethtb('graw    ',icACT,idval(idummy))
c      ilist(1,idval(idummy))=7
c      call setfnp(ilist(1,idval(idummy)+1),Act)
c
       idummy=sethtb('GRAPH   ',icACT,ktaloc(8))
       idummy1=sethtb('graph   ',icACT,idval(idummy))
       ilist(1,idval(idummy))=7
       call setfnp(klist(idval(idummy)+1),ActGra)
       ilist(1,idval(idummy)+2)=hsrch('ID')
       ilist(1,idval(idummy)+3)=hsrch('GTYPE')
       ilist(1,idval(idummy)+4)=hsrch('TURNS')
       ilist(1,idval(idummy)+5)=hsrch('SPAN    ')
       ilist(1,idval(idummy)+6)=hsrch('CENTER  ')
c
       kytbl(kwL,0)   =sethtb('L       ',icKWRD,kwL   )
       kytbl(kwKIN,0) =sethtb('DISKIN  ',icKWRD,kwKIN )
       kytbl(kwANGL,0)=sethtb('ANGLE   ',icKWRD,kwANGL)
       kytbl(kwROT ,0)=sethtb('ROT     ',icKWRD,kwROT)
       kytbl(kwROT ,0)=sethtb('TILT    ',icKWRD,kwROT)
       kytbl(kwROT ,0)=sethtb('ROTATE  ',icKWRD,kwROT)
c      kytbl(kwTILT,0)=sethtb('TILT    ',icKWRD,kwROT )
       kytbl(kwK0  ,0)=sethtb('K0      ',icKWRD,kwK0  )
       kytbl(kwK1  ,0)=sethtb('K1      ',icKWRD,kwK1  )
c      kytbl(kwDK1 ,0)=sethtb('DK1     ',icKWRD,kwDK1 )
       kytbl(kwK2  ,0)=sethtb('K2      ',icKWRD,kwK2  )
c      kytbl(kwDK2 ,0)=sethtb('DK2     ',icKWRD,kwDK2 )
       kytbl(kwK3  ,0)=sethtb('K3      ',icKWRD,kwK3  )
c      kytbl(kwDK3 ,0)=sethtb('DK3     ',icKWRD,kwDK3 )
       kytbl(kwK4  ,0)=sethtb('K4      ',icKWRD,kwK4  )
       kytbl(kwK5  ,0)=sethtb('K5      ',icKWRD,kwK5  )
       kytbl(kwK6  ,0)=sethtb('K6      ',icKWRD,kwK6  )
       kytbl(kwK7  ,0)=sethtb('K7      ',icKWRD,kwK7  )
       kytbl(kwK8  ,0)=sethtb('K8      ',icKWRD,kwK8  )
       kytbl(kwK9  ,0)=sethtb('K9      ',icKWRD,kwK9  )
       kytbl(kwK10 ,0)=sethtb('K10     ',icKWRD,kwK10  )
       kytbl(kwK10 ,0)=sethtb('K10    ',icKWRD,kwK10 )
       kytbl(kwK11 ,0)=sethtb('K11    ',icKWRD,kwK11 )
       kytbl(kwK12 ,0)=sethtb('K12    ',icKWRD,kwK12 )
       kytbl(kwK13 ,0)=sethtb('K13    ',icKWRD,kwK13 )
       kytbl(kwK14 ,0)=sethtb('K14    ',icKWRD,kwK14 )
       kytbl(kwK15 ,0)=sethtb('K15    ',icKWRD,kwK15 )
       kytbl(kwK16 ,0)=sethtb('K16    ',icKWRD,kwK16 )
       kytbl(kwK17 ,0)=sethtb('K17    ',icKWRD,kwK17 )
       kytbl(kwK18 ,0)=sethtb('K18    ',icKWRD,kwK18 )
       kytbl(kwK19 ,0)=sethtb('K19    ',icKWRD,kwK19 )
       kytbl(kwK20 ,0)=sethtb('K20    ',icKWRD,kwK20 )
       kytbl(kwK21 ,0)=sethtb('K21    ',icKWRD,kwK21 )
c
       kytbl(kwSK0 ,0)=sethtb('SK0     ',icKWRD,kwSK0 )
       kytbl(kwSK1 ,0)=sethtb('SK1     ',icKWRD,kwSK1 )
       kytbl(kwSK2 ,0)=sethtb('SK2     ',icKWRD,kwSK2 )
       kytbl(kwSK3 ,0)=sethtb('SK3     ',icKWRD,kwSK3 )
       kytbl(kwSK4 ,0)=sethtb('SK4     ',icKWRD,kwSK4 )
       kytbl(kwSK5 ,0)=sethtb('SK5     ',icKWRD,kwSK5 )
       kytbl(kwSK6 ,0)=sethtb('SK6     ',icKWRD,kwSK6 )
       kytbl(kwSK7 ,0)=sethtb('SK7     ',icKWRD,kwSK7  )
       kytbl(kwSK8 ,0)=sethtb('SK8     ',icKWRD,kwSK8  )
       kytbl(kwSK9 ,0)=sethtb('SK9     ',icKWRD,kwSK9  )
       kytbl(kwSK10,0)=sethtb('SK10    ',icKWRD,kwSK10 )
       kytbl(kwSK11,0)=sethtb('SK11    ',icKWRD,kwSK11 )
       kytbl(kwSK12,0)=sethtb('SK12    ',icKWRD,kwSK12 )
       kytbl(kwSK13,0)=sethtb('SK13    ',icKWRD,kwSK13 )
       kytbl(kwSK14,0)=sethtb('SK14    ',icKWRD,kwSK14 )
       kytbl(kwSK15 ,0)=sethtb('SK15    ',icKWRD,kwSK15 )
       kytbl(kwSK16 ,0)=sethtb('SK16    ',icKWRD,kwSK16 )
       kytbl(kwSK17 ,0)=sethtb('SK17    ',icKWRD,kwSK17 )
       kytbl(kwSK18 ,0)=sethtb('SK18    ',icKWRD,kwSK18 )
       kytbl(kwSK19 ,0)=sethtb('SK19    ',icKWRD,kwSK19 )
       kytbl(kwSK20 ,0)=sethtb('SK20    ',icKWRD,kwSK20 )
       kytbl(kwSK21 ,0)=sethtb('SK21    ',icKWRD,kwSK21 )
c
       kytbl(kwA3  ,0)=sethtb('A3      ',icKWRD,kwA3  )
       kytbl(kwA5  ,0)=sethtb('A5      ',icKWRD,kwA5  )
       kytbl(kwA7  ,0)=sethtb('A7      ',icKWRD,kwA7  )
       kytbl(kwA9  ,0)=sethtb('A9      ',icKWRD,kwA9  )
       kytbl(kwA11 ,0)=sethtb('A11     ',icKWRD,kwA11 )
       kytbl(kwA13 ,0)=sethtb('A13     ',icKWRD,kwA13 )
       kytbl(kwA15 ,0)=sethtb('A15     ',icKWRD,kwA15 )
       kytbl(kwA17 ,0)=sethtb('A17     ',icKWRD,kwA17 )
       kytbl(kwE1  ,0)=sethtb('E1      ',icKWRD,kwE1  )
       kytbl(kwE2  ,0)=sethtb('E2      ',icKWRD,kwE2  )
       kytbl(kwDX  ,0)=sethtb('DX      ',icKWRD,kwDX  )
       kytbl(kwDY  ,0)=sethtb('DY      ',icKWRD,kwDY  )
       kytbl(kwKICK,0)=sethtb('KICK    ',icKWRD,kwKICK)
       kytbl(kwVOLT,0)=sethtb('VOLT    ',icKWRD,kwVOLT)
       kytbl(kwDVOLT,0)=sethtb('DVOLT   ',icKWRD,kwDVOLT)
       kytbl(kwPHI ,0)=sethtb('PHI     ',icKWRD,kwPHI )
       kytbl(kwDPHI,0)=sethtb('DPHI    ',icKWRD,kwDPHI)
       kytbl(kwFREQ,0)=sethtb('FREQ    ',icKWRD,kwFREQ )
       kytbl(kwHARM,0)=sethtb('HARM    ',icKWRD,kwHARM)
       kytbl(kwLWAK,0)=sethtb('LWAKE   ',icKWRD,kwLWAK)
       kytbl(kwTWAK,0)=sethtb('TWAKE   ',icKWRD,kwTWAK)
       kytbl(kwAX  ,0)=sethtb('ALPHAX  ',icKWRD,kwAX  )
       kytbl(kwAX  ,0)=sethtb('AX      ',icKWRD,kwAX  )
       kytbl(kwAY  ,0)=sethtb('ALPHAY  ',icKWRD,kwAY  )
       kytbl(kwAY  ,0)=sethtb('AY      ',icKWRD,kwAY  )
       kytbl(kwBX  ,0)=sethtb('BETAX   ',icKWRD,kwBX  )
       kytbl(kwBX  ,0)=sethtb('BX      ',icKWRD,kwBX  )
       kytbl(kwBY  ,0)=sethtb('BETAY   ',icKWRD,kwBY  )
       kytbl(kwBY  ,0)=sethtb('BY      ',icKWRD,kwBY  )
       kytbl(kwEMIX,0)=sethtb('EMIX    ',icKWRD,kwEMIX)
       kytbl(kwEMIX,0)=sethtb('EMITX   ',icKWRD,kwEMIX)
       kytbl(kwEMIY,0)=sethtb('EMIY    ',icKWRD,kwEMIY)
       kytbl(kwEMIY,0)=sethtb('EMITY   ',icKWRD,kwEMIY)
       kytbl(kwEMIZ,0)=sethtb('EMIZ    ',icKWRD,kwEMIZ)
       kytbl(kwEMIZ,0)=sethtb('EMITZ   ',icKWRD,kwEMIZ)
       kytbl(kwPX  ,0)=sethtb('PSIX    ',icKWRD,kwPX  )
       kytbl(kwPY  ,0)=sethtb('PSIY    ',icKWRD,kwPY  )
       kytbl(kwPZ  ,0)=sethtb('PSIZ    ',icKWRD,kwPZ  )
       kytbl(kwDP  ,0)=sethtb('DP      ',icKWRD,kwDP  )
       kytbl(kwSIGZ,0)=sethtb('SIGZ    ',icKWRD,kwSIGZ)
       kytbl(kwSIGZ,0)=sethtb('SIGMAZ  ',icKWRD,kwSIGZ)
       kytbl(kwSIGE,0)=sethtb('SIGE    ',icKWRD,kwSIGZ)
       kytbl(kwGEO ,0)=sethtb('GEO     ',icKWRD,kwGEO )
       kytbl(kwR1  ,0)=sethtb('R1      ',icKWRD,kwR1  )
       kytbl(kwR2  ,0)=sethtb('R2      ',icKWRD,kwR2  )
       kytbl(kwR3  ,0)=sethtb('R3      ',icKWRD,kwR3  )
       kytbl(kwR4  ,0)=sethtb('R4      ',icKWRD,kwR4  )
       kytbl(kwDETR,0)=sethtb('DETR    ',icKWRD,kwR4  )
       kytbl(kwEX  ,0)=sethtb('EX      ',icKWRD,kwEX  )
       kytbl(kwEPX ,0)=sethtb('EPX     ',icKWRD,kwEPX )
       kytbl(kwEY  ,0)=sethtb('EY      ',icKWRD,kwEY  )
       kytbl(kwEPY ,0)=sethtb('EPY     ',icKWRD,kwEPY )
       kytbl(kwZX  ,0)=sethtb('ZX      ',icKWRD,kwZX  )
       kytbl(kwZPX ,0)=sethtb('ZPX     ',icKWRD,kwZPX )
       kytbl(kwZY  ,0)=sethtb('ZY      ',icKWRD,kwZY  )
       kytbl(kwZPY ,0)=sethtb('ZPY     ',icKWRD,kwZPY )
       kytbl(kwDPX ,0)=sethtb('DPX     ',icKWRD,kwDPX )
       kytbl(kwDPY ,0)=sethtb('DPY     ',icKWRD,kwDPY )
       kytbl(kwAZ  ,0)=sethtb('AZ      ',icKWRD,kwAZ  )
       kytbl(kwD11  ,0)=sethtb('D11      ',icKWRD,kwD11 )
       kytbl(kwD12  ,0)=sethtb('D12      ',icKWRD,kwD12 )
       kytbl(kwD13  ,0)=sethtb('D13      ',icKWRD,kwD13 )
       kytbl(kwD14  ,0)=sethtb('D14      ',icKWRD,kwD14 )
       kytbl(kwD15  ,0)=sethtb('D15      ',icKWRD,kwD15 )
       kytbl(kwD16  ,0)=sethtb('D16      ',icKWRD,kwD16 )
       kytbl(kwD21  ,0)=sethtb('D21      ',icKWRD,kwD21 )
       kytbl(kwD22  ,0)=sethtb('D22      ',icKWRD,kwD22 )
       kytbl(kwD23  ,0)=sethtb('D23      ',icKWRD,kwD23 )
       kytbl(kwD24  ,0)=sethtb('D24      ',icKWRD,kwD24 )
       kytbl(kwD25  ,0)=sethtb('D25      ',icKWRD,kwD25 )
       kytbl(kwD26  ,0)=sethtb('D26      ',icKWRD,kwD26 )
       kytbl(kwD31  ,0)=sethtb('D31      ',icKWRD,kwD31 )
       kytbl(kwD32  ,0)=sethtb('D32      ',icKWRD,kwD32 )
       kytbl(kwD33  ,0)=sethtb('D33      ',icKWRD,kwD33 )
       kytbl(kwD34  ,0)=sethtb('D34      ',icKWRD,kwD34 )
       kytbl(kwD35  ,0)=sethtb('D35      ',icKWRD,kwD35 )
       kytbl(kwD36  ,0)=sethtb('D36      ',icKWRD,kwD36 )
       kytbl(kwD41  ,0)=sethtb('D41      ',icKWRD,kwD41 )
       kytbl(kwD42  ,0)=sethtb('D42      ',icKWRD,kwD42 )
       kytbl(kwD43  ,0)=sethtb('D43      ',icKWRD,kwD43 )
       kytbl(kwD44  ,0)=sethtb('D44      ',icKWRD,kwD44 )
       kytbl(kwD45  ,0)=sethtb('D45      ',icKWRD,kwD45 )
       kytbl(kwD46  ,0)=sethtb('D46      ',icKWRD,kwD46 )
       kytbl(kwD51  ,0)=sethtb('D51      ',icKWRD,kwD51 )
       kytbl(kwD52  ,0)=sethtb('D52      ',icKWRD,kwD52 )
       kytbl(kwD53  ,0)=sethtb('D53      ',icKWRD,kwD53 )
       kytbl(kwD54  ,0)=sethtb('D54      ',icKWRD,kwD54 )
       kytbl(kwD55  ,0)=sethtb('D55      ',icKWRD,kwD55 )
       kytbl(kwD56  ,0)=sethtb('D56      ',icKWRD,kwD56 )
       kytbl(kwD61  ,0)=sethtb('D61      ',icKWRD,kwD61 )
       kytbl(kwD62  ,0)=sethtb('D62      ',icKWRD,kwD62 )
       kytbl(kwD63  ,0)=sethtb('D63      ',icKWRD,kwD63 )
       kytbl(kwD64  ,0)=sethtb('D64      ',icKWRD,kwD64 )
       kytbl(kwD65  ,0)=sethtb('D65      ',icKWRD,kwD65 )
       kytbl(kwD66  ,0)=sethtb('D66      ',icKWRD,kwD66 )
       kytbl(kwB11  ,0)=sethtb('B11      ',icKWRD,kwB11 )
       kytbl(kwB12  ,0)=sethtb('B12      ',icKWRD,kwB12 )
       kytbl(kwB13  ,0)=sethtb('B13      ',icKWRD,kwB13 )
       kytbl(kwB14  ,0)=sethtb('B14      ',icKWRD,kwB14 )
       kytbl(kwB15  ,0)=sethtb('B15      ',icKWRD,kwB15 )
       kytbl(kwB16  ,0)=sethtb('B16      ',icKWRD,kwB16 )
       kytbl(kwB22  ,0)=sethtb('B22      ',icKWRD,kwB22 )
       kytbl(kwB23  ,0)=sethtb('B23      ',icKWRD,kwB23 )
       kytbl(kwB24  ,0)=sethtb('B24      ',icKWRD,kwB24 )
       kytbl(kwB25  ,0)=sethtb('B25      ',icKWRD,kwB25 )
       kytbl(kwB26  ,0)=sethtb('B26      ',icKWRD,kwB26 )
       kytbl(kwB33  ,0)=sethtb('B33      ',icKWRD,kwB33 )
       kytbl(kwB34  ,0)=sethtb('B34      ',icKWRD,kwB34 )
       kytbl(kwB35  ,0)=sethtb('B35      ',icKWRD,kwB35 )
       kytbl(kwB36  ,0)=sethtb('B36      ',icKWRD,kwB36 )
       kytbl(kwB44  ,0)=sethtb('B44      ',icKWRD,kwB44 )
       kytbl(kwB45  ,0)=sethtb('B45      ',icKWRD,kwB45 )
       kytbl(kwB46  ,0)=sethtb('B46      ',icKWRD,kwB46 )
       kytbl(kwB55  ,0)=sethtb('B55      ',icKWRD,kwB55 )
       kytbl(kwB56  ,0)=sethtb('B56      ',icKWRD,kwB56 )
       kytbl(kwB66  ,0)=sethtb('B66      ',icKWRD,kwB66 )
       kytbl(kwR11  ,0)=sethtb('R11      ',icKWRD,kwR11 )
       kytbl(kwR12  ,0)=sethtb('R12      ',icKWRD,kwR12 )
       kytbl(kwR13  ,0)=sethtb('R13      ',icKWRD,kwR13 )
       kytbl(kwR14  ,0)=sethtb('R14      ',icKWRD,kwR14 )
       kytbl(kwR15  ,0)=sethtb('R15      ',icKWRD,kwR15 )
       kytbl(kwR16  ,0)=sethtb('R16      ',icKWRD,kwR16 )
       kytbl(kwR22  ,0)=sethtb('R22      ',icKWRD,kwR22 )
       kytbl(kwR23  ,0)=sethtb('R23      ',icKWRD,kwR23 )
       kytbl(kwR24  ,0)=sethtb('R24      ',icKWRD,kwR24 )
       kytbl(kwR25  ,0)=sethtb('R25      ',icKWRD,kwR25 )
       kytbl(kwR26  ,0)=sethtb('R26      ',icKWRD,kwR26 )
       kytbl(kwR33  ,0)=sethtb('R33      ',icKWRD,kwR33 )
       kytbl(kwR34  ,0)=sethtb('R34      ',icKWRD,kwR34 )
       kytbl(kwR35  ,0)=sethtb('R35      ',icKWRD,kwR35 )
       kytbl(kwR36  ,0)=sethtb('R36      ',icKWRD,kwR36 )
       kytbl(kwR44  ,0)=sethtb('R44      ',icKWRD,kwR44 )
       kytbl(kwR45  ,0)=sethtb('R45      ',icKWRD,kwR45 )
       kytbl(kwR46  ,0)=sethtb('R46      ',icKWRD,kwR46 )
       kytbl(kwR55  ,0)=sethtb('R55      ',icKWRD,kwR55 )
       kytbl(kwR56  ,0)=sethtb('R56      ',icKWRD,kwR56 )
       kytbl(kwR66  ,0)=sethtb('R66      ',icKWRD,kwR66 )
       kytbl(kwRAD ,0)=sethtb('DISRAD  ',icKWRD,kwRAD )
       kytbl(kwCHRO,0)=sethtb('ACHROMA ',icKWRD,kwCHRO)
       kytbl(kwFRIN,0)=sethtb('DISFRIN ',icKWRD,kwFRIN)
       kytbl(kwF1  ,0)=sethtb('F1      ',icKWRD,kwF1  )
       kytbl(kwF2  ,0)=sethtb('F2      ',icKWRD,kwF2  )
       kytbl(kwF1K1F,0)=sethtb('F1K1F  ',icKWRD,kwF1K1F)
       kytbl(kwF2K1F,0)=sethtb('F2K1F  ',icKWRD,kwF2K1F)
       kytbl(kwF1K1B,0)=sethtb('F1K1B  ',icKWRD,kwF1K1B)
       kytbl(kwF2K1B,0)=sethtb('F2K1B  ',icKWRD,kwF2K1B)
       kytbl(kwF2  ,0)=sethtb('F2      ',icKWRD,kwF2  )
       kytbl(kwFRMD,0)=sethtb('FRINGE  ',icKWRD,kwFRMD)
       kytbl(kwK0FR,0)=sethtb('DISK0FR ',icKWRD,kwK0FR)
       kytbl(kwEPS ,0)=sethtb('EPS     ',icKWRD,kwEPS )
       kytbl(kwRANK,0)=sethtb('RANKICK ',icKWRD,kwRANK)
       kytbl(kwRANV,0)=sethtb('RANVOLT ',icKWRD,kwRANV)
       kytbl(kwRANP,0)=sethtb('RANPHASE',icKWRD,kwRANP)
       kytbl(kwDZ  ,0)=sethtb('DZ      ',icKWRD,kwDZ  )
       kytbl(kwDDP ,0)=sethtb('DDP     ',icKWRD,kwDZ  )
       kytbl(kwINDX,0)=sethtb('INDEX   ',icKWRD,kwINDX)
       kytbl(kwBMAX,0)=sethtb('BMAX    ',icKWRD,kwBMAX)
       kytbl(kwBND ,0)=sethtb('BOUND   ',icKWRD,kwBND )
       kytbl(kwPRD ,0)=sethtb('PERIOD  ',icKWRD,kwPRD )
       kytbl(kwBZ  ,0)=sethtb('BZ      ',icKWRD,kwBZ  )
       kytbl(kwDBZ ,0)=sethtb('DBZ     ',icKWRD,kwDBZ )
       kytbl(kwCHI1,0)=sethtb('CHI1    ',icKWRD,kwCHI1)
       kytbl(kwCHI2,0)=sethtb('CHI2    ',icKWRD,kwCHI2)
       kytbl(kwCHI3,0)=sethtb('CHI3    ',icKWRD,kwCHI3)
       kytbl(kwDIR ,0)=sethtb('DIR     ',icKWRD,kwDIR )
       kytbl(kwSLI ,0)=sethtb('SLICE   ',icKWRD,kwSLI )
       kytbl(kwSTURN ,0)=sethtb('STURN   ',icKWRD,kwSTURN )
       kytbl(kwBSTRL ,0)=sethtb('BSTRL   ',icKWRD,kwBSTRL )
       kytbl(kwXANGLE ,0)=sethtb('XANGLE  ',icKWRD,kwXANGLE )
       kytbl(kwNP  ,0)=sethtb('NP      ',icKWRD,kwNP  )
       kytbl(kwKx  ,0)=sethtb('KX      ',icKWRD,kwKx  )
       kytbl(kwQy  ,0)=sethtb('QY      ',icKWRD,kwQy  )
       kytbl(kwFBx ,0)=sethtb('B0X     ',icKWRD,kwFBx  )
       kytbl(kwFBy ,0)=sethtb('B0Y     ',icKWRD,kwFBy  )
       kytbl(kwPole,0)=sethtb('POLE    ',icKWRD,kwPole  )
       kytbl(kwJDX ,0)=sethtb('JDX     ',icKWRD,kwJDX  )
       kytbl(kwJDY  ,0)=sethtb('JDY     ',icKWRD,kwJDY  )
       kytbl(kwJDZ  ,0)=sethtb('JDZ     ',icKWRD,kwJDZ  )
       kytbl(kwJDPX ,0)=sethtb('JDPX    ',icKWRD,kwJDPX  )
       kytbl(kwJDPY ,0)=sethtb('JDPY    ',icKWRD,kwJDPY  )
       kytbl(kwJDPZ ,0)=sethtb('JDPZ    ',icKWRD,kwJDPZ  )
       kytbl(kwOFFSET ,0)=sethtb('OFFSET    ',icKWRD,kwOFFSET  )
       kytbl(kwCOUPLE ,0)=sethtb('COUPLE    ',icKWRD,kwCOUPLE  )
       kytbl(kwV1   ,0)=sethtb('V1      ',icKWRD,kwV1  )
       kytbl(kwV20  ,0)=sethtb('V20     ',icKWRD,kwV20  )
       kytbl(kwV11  ,0)=sethtb('V11     ',icKWRD,kwV11  )
       kytbl(kwV02  ,0)=sethtb('V02     ',icKWRD,kwV02  )
       kytbl(kwDX1  ,0)=sethtb('DX1     ',icKWRD,kwDX1  )
       kytbl(kwDX2  ,0)=sethtb('DX2     ',icKWRD,kwDX2  )
       kytbl(kwDY1  ,0)=sethtb('DY1     ',icKWRD,kwDY1  )
       kytbl(kwDY2  ,0)=sethtb('DY2     ',icKWRD,kwDY2  )
       kytbl(kwRADI ,0)=sethtb('RADIUS  ',icKWRD,kwRADI )
       kytbl(kwW1   ,0)=sethtb('W1      ',icKWRD,kwW1 )
       kytbl(kwDROT ,0)=sethtb('DROT    ',icKWRD,kwDROT)
       kytbl(kwDROT ,0)=sethtb('DROTATE ',icKWRD,kwDROT)
       kytbl(kwAE1 ,0) =sethtb('AE1     ',icKWRD,kwAE1)
       kytbl(kwAE2 ,0) =sethtb('AE2     ',icKWRD,kwAE2)
       kytbl(kwFB1 ,0) =sethtb('FB1     ',icKWRD,kwFB1)
       kytbl(kwFB2 ,0) =sethtb('FB2     ',icKWRD,kwFB2)
       kytbl(kwLDEV,0) =sethtb('LDEV    ',icKWRD,kwLDEV)
       kytbl(kwLRAD,0) =sethtb('LRAD    ',icKWRD,kwLRAD)
       kytbl(kwFL,0)   =sethtb('FLAT    ',icKWRD,kwLRAD)
       kytbl(kwAPHI,0) =sethtb('AUTOPHI ',icKWRD,kwAPHI)
       kytbl(kwPROF,0) =sethtb('PROFILE ',icKWRD,kwPROF)
c  for drift 
      idummy=sethtb('drift   ',icDEF,icDRFT)
      kytbl(0,icDRFT)=sethtb('DRIFT   ',icDEF,icDRFT)
      kytbl(kwL,icDRFT)=ky_L_DRFT
      kytbl(kwKIN,icDRFT)=ky_KIN_DRFT
      kytbl(kwCOUPLE,icDRFT)=ky_COUPLE_DRFT
      kytbl(kwRADI,icDRFT)=ky_RADI_DRFT
c
      kytbl(kwAX,icDRFT)=ky_AX_DRFT
      kytbl(kwAY,icDRFT)=ky_AY_DRFT
      kytbl(kwRAD,icDRFT)=ky_RAD_DRFT
c
      kytbl(kwMAX,icDRFT)=ky_MAX_DRFT
cc for STeering
      idummy=sethtb('st      ',icDEF,icBEND)
      idummy=sethtb('ST      ',icDEF,icBEND)
c                         number of parameters
c  for bend
      kytbl(0,icBEND)=sethtb('bend    ',icDEF,icBEND)
      kytbl(0,icBEND)=sethtb('BEND    ',icDEF,icBEND)
      kytbl(kwL,icBEND)=ky_L_BEND
      kytbl(kwANGL,icBEND)=ky_ANGL_BEND
      kytbl(kwE1,icBEND)=ky_E1_BEND
      kytbl(kwE2,icBEND)=ky_E2_BEND
      kytbl(kwROT,icBEND)=ky_ROT_BEND
c     kytbl(kwTILT,icBEND)=ky_TILT_BEND
      kytbl(kwRAD,icBEND)=ky_RAD_BEND
      kytbl(kwK1,icBEND)=ky_K1_BEND
      kytbl(kwDX,icBEND)=ky_DX_BEND
      kytbl(kwDY,icBEND)=ky_DY_BEND
      kytbl(kwK0,icBEND)=ky_K0_BEND
      kytbl(kwFRIN,icBEND)=ky_FRIN_BEND
      kytbl(kwEPS,icBEND)=ky_EPS_BEND
      kytbl(kwRANK,icBEND)=ky_RANK_BEND
      kytbl(kwF1,icBEND)=ky_F1_BEND
      kytbl(kwFRMD,icBEND)=ky_FRMD_BEND
      kytbl(kwCOUPLE,icBEND)=ky_COUPLE_BEND
      kytbl(kwDROT,icBEND)=ky_DROT_BEND
      kytbl(kwAE1,icBEND)=ky_AE1_BEND
      kytbl(kwAE2,icBEND)=ky_AE2_BEND
      kytbl(kwFB1,icBEND)=ky_FB1_BEND
      kytbl(kwFB2,icBEND)=ky_FB2_BEND
      kytbl(kwINDX,icBEND)=ky_INDX_BEND
c
      kytbl(kwAX,icBEND)=ky_AX_BEND
      kytbl(kwAY,icBEND)=ky_AY_BEND
      kytbl(kwLDEV,icBEND)=ky_LDEV_BEND
      kytbl(kwLRAD,icBEND)=ky_LRAD_BEND
c
      kytbl(kwMAX,icBEND)=ky_MAX_BEND
      kytbl(kwNPARAM,icBEND)=p_NPARAM_BEND
c  for quad
      kytbl(0,icQUAD)=sethtb('quad    ',icDEF,icQUAD)
      kytbl(0,icQUAD)=sethtb('QUAD    ',icDEF,icQUAD)
      kytbl(kwL,icQUAD)=ky_L_QUAD
      kytbl(kwK1,icQUAD)=ky_K1_QUAD
c     kytbl(kwDK1,icQUAD)=ky_DK1_QUAD
      kytbl(kwROT,icQUAD)=ky_ROT_QUAD
      kytbl(kwDX,icQUAD)=ky_DX_QUAD
      kytbl(kwDY,icQUAD)=ky_DY_QUAD
      kytbl(kwRAD,icQUAD)=ky_RAD_QUAD
      kytbl(kwCHRO,icQUAD)=ky_CHRO_QUAD
      kytbl(kwFRIN,icQUAD)=ky_FRIN_QUAD
      kytbl(kwF1,icQUAD)=ky_F1_QUAD
      kytbl(kwF2,icQUAD)=ky_F2_QUAD
      kytbl(kwFRMD,icQUAD)=ky_FRMD_QUAD
      kytbl(kwEPS,icQUAD)=ky_EPS_QUAD
      kytbl(kwKIN,icQUAD)=ky_KIN_QUAD
      kytbl(kwCOUPLE,icQUAD)=ky_COUPLE_QUAD
      kytbl(kwINDX,icQUAD)=ky_INDX_QUAD
c
      kytbl(kwAX,icQUAD)=ky_AX_QUAD
      kytbl(kwAY,icQUAD)=ky_AY_QUAD
      kytbl(kwLDEV,icQUAD)=ky_LDEV_QUAD
      kytbl(kwLRAD,icQUAD)=ky_LRAD_QUAD
      kytbl(kwF1K1F,icQUAD)=ky_F1K1F_QUAD
      kytbl(kwF2K1F,icQUAD)=ky_F2K1F_QUAD
      kytbl(kwF1K1B,icQUAD)=ky_F1K1B_QUAD
      kytbl(kwF2K1B,icQUAD)=ky_F2K1B_QUAD
c
      kytbl(kwMAX,icQUAD)=ky_MAX_QUAD
      kytbl(kwNPARAM,icQUAD)=p_NPARAM_QUAD
c  for THIN
      kytbl(0,icSEXT)=sethtb('sext    ',icDEF,icsext)
      kytbl(0,icSEXT)=sethtb('SEXT    ',icDEF,icsext)
      kytbl(0,icOCTU)=sethtb('oct    ',icDEF,icOCTU)
      kytbl(0,icOCTU)=sethtb('OCT    ',icDEF,icOCTU)
      kytbl(0,icDECA)=sethtb('deca    ',icDEF,icDECA)
      kytbl(0,icDECA)=sethtb('DECA    ',icDEF,icDECA)
      kytbl(0,icDODECA)=sethtb('dodeca    ',icDEF,icdodeca)
      kytbl(0,icDODECA)=sethtb('DODECA    ',icDEF,icdodeca)

      do i=icSEXT,icDODECA,2
        kytbl(kwL,i)=ky_L_THIN
c     kytbl(kwDK2,i)=ky_DK2_THIN
        kytbl(kwROT,i)=ky_ROT_THIN
        kytbl(kwDX,i)=ky_DX_THIN
        kytbl(kwDY,i)=ky_DY_THIN
        kytbl(kwRAD,i)=ky_RAD_THIN
        kytbl(kwFRIN,i)=ky_FRIN_THIN
        kytbl(kwCOUPLE,i)=ky_COUPLE_THIN
        kytbl(kwINDX,i)=ky_INDX_THIN
c     
        kytbl(kwAX,i)=ky_AX_THIN
        kytbl(kwAY,i)=ky_AY_THIN
        kytbl(kwLDEV,i)=ky_LDEV_THIN
        kytbl(kwLRAD,i)=ky_LRAD_THIN
c     
        kytbl(kwMAX,i)=ky_MAX_THIN
        kytbl(kwNPARAM,i)=p_NPARAM_THIN
      enddo
      kytbl(kwK2,icSEXT)=ky_K_THIN
      kytbl(kwK3,icOCTU)=ky_K_THIN
      kytbl(kwK4,icDECA)=ky_K_THIN
      kytbl(kwK5,icDODECA)=ky_K_THIN
c  for octu
c  for MULT
      kytbl(0,icMULT)=sethtb('mult    ',icDEF,icMULT)
      kytbl(0,icMULT)=sethtb('MULT    ',icDEF,icMULT)
      kytbl(kwL,icMULT)=ky_L_MULT
      kytbl(kwDX,icMULT)=ky_DX_MULT
      kytbl(kwDY,icMULT)=ky_DY_MULT
      kytbl(kwDZ,icMULT)=ky_DZ_MULT
      kytbl(kwCHI1,icMULT)=ky_CHI1_MULT
      kytbl(kwCHI2,icMULT)=ky_CHI2_MULT
      kytbl(kwCHI3,icMULT)=ky_CHI3_MULT
      kytbl(kwROT,icMULT)=ky_ROT_MULT
      kytbl(kwEPS,icMULT)=ky_EPS_MULT
      kytbl(kwRAD,icMULT)=ky_RAD_MULT
      kytbl(kwFRIN,icMULT)=ky_FRIN_MULT
      kytbl(kwF1,icMULT)=ky_F1_MULT
      kytbl(kwF2,icMULT)=ky_F2_MULT
      kytbl(kwFRMD,icMULT)=ky_FRMD_MULT
      kytbl(kwVOLT,icMULT)=ky_VOLT_MULT
      kytbl(kwDVOLT,icMULT)=ky_DVOLT_MULT
      kytbl(kwHARM,icMULT)=ky_HARM_MULT
      kytbl(kwPHI,icMULT)=ky_PHI_MULT
      kytbl(kwFREQ,icMULT)=ky_FREQ_MULT
      kytbl(kwCOUPLE,icMULT)=ky_COUPLE_MULT
      kytbl(kwRADI,icMULT)=ky_RADI_MULT
      kytbl(kwDPHI,icMULT)=ky_DPHI_MULT
      kytbl(kwW1,icMULT)=ky_W1_MULT
c
      kytbl(kwANGL,icMULT)=ky_ANGL_MULT
      kytbl(kwE1,icMULT)=ky_E1_MULT
      kytbl(kwE2,icMULT)=ky_E2_MULT
      kytbl(kwDROT,icMULT)=ky_DROT_MULT
      kytbl(kwK0FR,icMULT)=ky_K0FR_MULT
      kytbl(kwFB1,icMULT)=ky_FB1_MULT
      kytbl(kwFB2,icMULT)=ky_FB2_MULT
c
      kytbl(kwK0,icMULT)=ky_K0_MULT
      kytbl(kwSK0,icMULT)=ky_SK0_MULT
      kytbl(kwK1,icMULT)=ky_K1_MULT
      kytbl(kwSK1,icMULT)=ky_SK1_MULT
      kytbl(kwK2,icMULT)=ky_K2_MULT
      kytbl(kwSK2,icMULT)=ky_SK2_MULT
      kytbl(kwK3,icMULT)=ky_K3_MULT
      kytbl(kwSK3,icMULT)=ky_SK3_MULT
      kytbl(kwK4,icMULT)=ky_K4_MULT
      kytbl(kwSK4,icMULT)=ky_SK4_MULT
      kytbl(kwK5,icMULT)=ky_K5_MULT
      kytbl(kwSK5,icMULT)=ky_SK5_MULT
      kytbl(kwK6,icMULT)=ky_K6_MULT
      kytbl(kwSK6,icMULT)=ky_SK6_MULT
      kytbl(kwK7,icMULT)=ky_K7_MULT
      kytbl(kwSK7,icMULT)=ky_SK7_MULT
      kytbl(kwK8,icMULT)=ky_K8_MULT
      kytbl(kwSK8,icMULT)=ky_SK8_MULT
      kytbl(kwK9,icMULT)=ky_K9_MULT
      kytbl(kwSK9,icMULT)=ky_SK9_MULT
      kytbl(kwK10,icMULT)=ky_K10_MULT
      kytbl(kwSK10,icMULT)=ky_SK10_MULT
      kytbl(kwK11,icMULT)=ky_K11_MULT
      kytbl(kwSK11,icMULT)=ky_SK11_MULT
      kytbl(kwK12,icMULT)=ky_K12_MULT
      kytbl(kwSK12,icMULT)=ky_SK12_MULT
      kytbl(kwK13,icMULT)=ky_K13_MULT
      kytbl(kwSK13,icMULT)=ky_SK13_MULT
      kytbl(kwK14,icMULT)=ky_K14_MULT
      kytbl(kwSK14,icMULT)=ky_SK14_MULT
      kytbl(kwK15,icMULT)=ky_K15_MULT
      kytbl(kwSK15,icMULT)=ky_SK15_MULT
      kytbl(kwK16,icMULT)=ky_K16_MULT
      kytbl(kwSK16,icMULT)=ky_SK16_MULT
      kytbl(kwK17,icMULT)=ky_K17_MULT
      kytbl(kwSK17,icMULT)=ky_SK17_MULT
      kytbl(kwK18,icMULT)=ky_K18_MULT
      kytbl(kwSK18,icMULT)=ky_SK18_MULT
      kytbl(kwK19,icMULT)=ky_K19_MULT
      kytbl(kwSK19,icMULT)=ky_SK19_MULT
      kytbl(kwK20,icMULT)=ky_K20_MULT
      kytbl(kwSK20,icMULT)=ky_SK20_MULT
      kytbl(kwK21,icMULT)=ky_K21_MULT
      kytbl(kwSK21,icMULT)=ky_SK21_MULT
      kytbl(kwAE1,icMULT)=ky_AE1_MULT
      kytbl(kwAE2,icMULT)=ky_AE2_MULT
c
      kytbl(kwAX,icMULT)=ky_AX_MULT
      kytbl(kwAY,icMULT)=ky_AY_MULT
      kytbl(kwLDEV,icMULT)=ky_LDEV_MULT
      kytbl(kwLRAD,icMULT)=ky_LRAD_MULT
      kytbl(kwAPHI,icMULT)=ky_APHI_MULT
      kytbl(kwF1K1F,icMULT)=ky_F1K1F_MULT
      kytbl(kwF2K1F,icMULT)=ky_F2K1F_MULT
      kytbl(kwF1K1B,icMULT)=ky_F1K1B_MULT
      kytbl(kwF2K1B,icMULT)=ky_F2K1B_MULT
      kytbl(kwPROF,icMULT)=ky_PROF_MULT
c
      kytbl(kwMAX,icMULT)=ky_MAX_MULT
      kytbl(kwNPARAM,icMULT)=p_NPARAM_MULT
c
c  for UNDULATOR
      kytbl(0,icUND)=sethtb('und    ',icDEF,icUND)
      kytbl(0,icUND)=sethtb('UND    ',icDEF,icUND)
      kytbl(kwL,icUND)=ky_L_UND
      kytbl(kwFBX,icUND)=ky_FBX_UND
      kytbl(kwFBY,icUND)=ky_FBY_UND
      kytbl(kwKx,icUND)=ky_Kx_UND
      kytbl(kwQy,icUND)=ky_Qy_UND
      kytbl(kwPHI,icUND)=ky_PHI_UND
      kytbl(kwSLI,icUND)=ky_SLI_UND
      kytbl(kwPole,icUND)=ky_Pole_UND
      kytbl(kwCOUPLE,icUND)=ky_COUPLE_UND
      kytbl(kwMAX,icUND)=ky_MAX_UND
      kytbl(kwNPARAM,icUND)=p_NPARAM_UND
c  for WIG
      kytbl(0,icWIG)=sethtb('wig    ',icDEF,icWIG)
      kytbl(0,icWIG)=sethtb('WIG    ',icDEF,icWIG)
      kytbl(kwL,icWIG)=ky_L_WIG
      kytbl(kwBMAX,icWIG)=ky_BMAX_WIG
      kytbl(kwPRD,icWIG)=ky_PRD_WIG
      kytbl(kwROT,icWIG)=ky_ROT_WIG
      kytbl(kwDX,icWIG)=ky_DX_WIG
      kytbl(kwDY,icWIG)=ky_DY_WIG
      kytbl(kwA3,icWIG)=ky_A3_WIG
      kytbl(kwA5,icWIG)=ky_A5_WIG
      kytbl(kwA7,icWIG)=ky_A7_WIG
      kytbl(kwA9,icWIG)=ky_A9_WIG
      kytbl(kwA11,icWIG)=ky_A11_WIG
      kytbl(kwA13,icWIG)=ky_A13_WIG
      kytbl(kwA15,icWIG)=ky_A15_WIG
      kytbl(kwA17,icWIG)=ky_A17_WIG
      kytbl(kwCOUPLE,icWIG)=ky_COUPLE_WIG
      kytbl(kwMAX,icWIG)=ky_MAX_WIG
c  for solenoid
      kytbl(0,icSOL)=sethtb('sol     ',icDEF,icSOL)
      kytbl(0,icSOL)=sethtb('SOL     ',icDEF,icSOL)
      kytbl(kwL,icSOL)=ky_L_SOL
      kytbl(kwBZ,icSOL)=ky_BZ_SOL
      kytbl(kwDX,icSOL)=ky_DX_SOL
      kytbl(kwDY,icSOL)=ky_DY_SOL
      kytbl(kwDZ,icSOL)=ky_DZ_SOL
      kytbl(kwDPX,icSOL)=ky_DPX_SOL
      kytbl(kwDPY,icSOL)=ky_DPY_SOL
      kytbl(kwBND,icSOL)=ky_BND_SOL
      kytbl(kwCHI1,icSOL)=ky_CHI1_SOL
      kytbl(kwCHI2,icSOL)=ky_CHI2_SOL
      kytbl(kwCHI3,icSOL)=ky_CHI3_SOL
      kytbl(kwGEO,icSOL)=ky_GEO_SOL
      kytbl(kwF1,icSOL)=ky_F1_SOL
      kytbl(kwDBZ,icSOL)=ky_DBZ_SOL
      kytbl(kwCOUPLE,icSOL)=ky_COUPLE_SOL
      kytbl(kwFRIN,icSOL)=ky_FRIN_SOL
      kytbl(kwRAD,icSOL)=ky_RAD_SOL
      kytbl(kwFL,icSOL)=ky_FL_SOL
      kytbl(kwMAX,icSOL)=ky_MAX_SOL
      kytbl(kwNPARAM,icSOL)=p_NPARAM_SOL
cc for TEST
      idummy=sethtb('test    ',icDEF,icTEST)
      kytbl(0,icTEST)=sethtb('TEST    ',icDEF,icTEST)
      kytbl(kwL,icTEST)=ky_L_TEST
      kytbl(kwANGL,icTEST)=ky_ANGL_TEST
      kytbl(kwMAX,icTEST)=ky_MAX_TEST
cc for Cavity
      kytbl(0,icCAVI)=sethtb('cavi    ',icDEF,icCavi)
      kytbl(0,icCAVI)=sethtb('CAVI    ',icDEF,icCavi)
      kytbl(kwL,icCAVI)=ky_L_CAVI
      kytbl(kwVOLT,icCAVI)=ky_VOLT_CAVI
      kytbl(kwDVOLT,icCAVI)=ky_DVOLT_CAVI
      kytbl(kwHARM,icCAVI)=ky_HARM_CAVI
      kytbl(kwPHI,icCAVI)=ky_PHI_CAVI
      kytbl(kwFREQ,icCAVI)=ky_FREQ_CAVI
      kytbl(kwRANV,icCAVI)=ky_RANV_CAVI
      kytbl(kwRANP,icCAVI)=ky_RANP_CAVI
      kytbl(kwLWAK,icCAVI)=ky_LWAK_CAVI
      kytbl(kwTWAK,icCAVI)=ky_TWAK_CAVI
      kytbl(kwDX,icCAVI)=ky_DX_CAVI
      kytbl(kwDY,icCAVI)=ky_DY_CAVI
      kytbl(kwROT,icCAVI)=ky_ROT_CAVI
      kytbl(kwV1,icCAVI)=ky_V1_CAVI
      kytbl(kwV20,icCAVI)=ky_V20_CAVI
      kytbl(kwV11,icCAVI)=ky_V11_CAVI
      kytbl(kwV02,icCAVI)=ky_V02_CAVI
      kytbl(kwCOUPLE,icCAVI)=ky_COUPLE_CAVI
      kytbl(kwDPHI,icCAVI)=ky_DPHI_CAVI
      kytbl(kwFRIN,icCAVI)=ky_FRIN_CAVI
      kytbl(kwFRMD,icCAVI)=ky_FRMD_CAVI
c
      kytbl(kwAX,icCAVI)=ky_AX_CAVI
      kytbl(kwAY,icCAVI)=ky_AY_CAVI
      kytbl(kwLDEV,icCAVI)=ky_LDEV_CAVI
      kytbl(kwAPHI,icCAVI)=ky_APHI_CAVI
      kytbl(kwMAX,icCAVI)=ky_MAX_CAVI
      kytbl(kwNPARAM,icCAVI)=p_NPARAM_CAVI
cc for t-Cavity
      kytbl(0,icTCAV)=sethtb('tcavi   ',icDEF,icTCAV)
      kytbl(0,icTCAV)=sethtb('TCAVI   ',icDEF,icTCAV)
      kytbl(kwL,icTCAV)=ky_L_TCAV
      kytbl(kwK0,icTCAV)=ky_K0_TCAV
      kytbl(kwHARM,icTCAV)=ky_HARM_TCAV
      kytbl(kwPHI,icTCAV)=ky_PHI_TCAV
      kytbl(kwFREQ,icTCAV)=ky_FREQ_TCAV
      kytbl(kwDX,icTCAV)=ky_DX_TCAV
      kytbl(kwDY,icTCAV)=ky_DY_TCAV
      kytbl(kwROT,icTCAV)=ky_ROT_TCAV
      kytbl(kwRANK,icTCAV)=ky_RANK_TCAV
      kytbl(kwRANP,icTCAV)=ky_RANP_TCAV
      kytbl(kwLWAK,icTCAV)=ky_LWAK_TCAV
      kytbl(kwTWAK,icTCAV)=ky_TWAK_TCAV
      kytbl(kwCOUPLE,icTCAV)=ky_COUPLE_TCAV
c
      kytbl(kwAX,icTCAV)=ky_AX_TCAV
      kytbl(kwAY,icTCAV)=ky_AY_TCAV
      kytbl(kwLDEV,icTCAV)=ky_LDEV_TCAV
c     kytbl(kwLRAD,icTCAV)=ky_LRAD_TCAV
c
      kytbl(kwMAX,icTCAV)=ky_MAX_TCAV
cc for MAP
      kytbl(0,icMAP)=sethtb('map    ',icDEF,icMAP)
      kytbl(0,icMAP)=sethtb('MAP    ',icDEF,icMAP)
      kytbl(kwL,icMAP)=ky_L_MAP
      kytbl(kwMAX,icMAP)=ky_MAX_MAP
cc for INS
      kytbl(0,icINS)=sethtb('ins    ',icDEF,icINS)
      kytbl(0,icINS)=sethtb('INS    ',icDEF,icINS)
      kytbl(kwAX,icINS)=ky_AX_INS
      kytbl(kwBX,icINS)=ky_BX_INS
      kytbl(kwPX,icINS)=ky_PX_INS
      kytbl(kwAY,icINS)=ky_AY_INS
      kytbl(kwBY,icINS)=ky_BY_INS
      kytbl(kwPY,icINS)=ky_PY_INS
      kytbl(kwEX,icINS)=ky_EX_INS
      kytbl(kwEPX,icINS)=ky_EPX_INS
      kytbl(kwEY,icINS)=ky_EY_INS
      kytbl(kwEPY,icINS)=ky_EPY_INS
      kytbl(kwR1,icINS)=ky_R1_INS
      kytbl(kwR2,icINS)=ky_R2_INS
      kytbl(kwR3,icINS)=ky_R3_INS
      kytbl(kwR4,icINS)=ky_R4_INS
      kytbl(kwDX,icINS)=ky_DX_INS
      kytbl(kwDPX,icINS)=ky_DPX_INS
      kytbl(kwDY,icINS)=ky_DY_INS
      kytbl(kwDPY,icINS)=ky_DPY_INS
      kytbl(kwDIR,icINS)=ky_DIR_INS
      kytbl(kwCOUPLE,icINS)=ky_COUPLE_INS
      kytbl(kwMAX,icINS)=ky_MAX_INS
cc for Coordinate transformation
      kytbl(0,     icCOORD  )=sethtb('COORD   ',icDEF,icCOORD)
      kytbl(kwDX,icCOORD)=ky_DX_COORD
      kytbl(kwDY,icCOORD)=ky_DY_COORD
      kytbl(kwDZ,icCOORD)=ky_DZ_COORD
      kytbl(kwCHI1,icCOORD)=ky_CHI1_COORD
      kytbl(kwCHI2,icCOORD)=ky_CHI2_COORD
      kytbl(kwCHI3,icCOORD)=ky_CHI3_COORD
      kytbl(kwDIR,icCOORD)=ky_DIR_COORD
      kytbl(kwCOUPLE,icCOORD)=ky_COUPLE_COORD
      kytbl(kwMAX,icCOORD)=ky_MAX_COORD
cc for BEAMBEAM
      kytbl(0,icBEAM)=sethtb('beambeam',icDEF,icBEAM)
      kytbl(0,icBEAM)=sethtb('BEAMBEAM',icDEF,icBEAM)
      kytbl(kwAX,icBEAM)=ky_AX_BEAM
      kytbl(kwBX,icBEAM)=ky_BX_BEAM
      kytbl(kwAY,icBEAM)=ky_AY_BEAM
      kytbl(kwBY,icBEAM)=ky_BY_BEAM
      kytbl(kwR1,icBEAM)=ky_R1_BEAM
      kytbl(kwR2,icBEAM)=ky_R2_BEAM
      kytbl(kwR3,icBEAM)=ky_R3_BEAM
      kytbl(kwR4,icBEAM)=ky_R4_BEAM
      kytbl(kwEX,icBEAM)=ky_EX_BEAM
      kytbl(kwEPX,icBEAM)=ky_EPX_BEAM
      kytbl(kwEY,icBEAM)=ky_EY_BEAM
      kytbl(kwEPY,icBEAM)=ky_EPY_BEAM
      kytbl(kwZX,icBEAM)=ky_ZX_BEAM
      kytbl(kwZPX,icBEAM)=ky_ZPX_BEAM
      kytbl(kwZY,icBEAM)=ky_ZY_BEAM
      kytbl(kwZPY,icBEAM)=ky_ZPY_BEAM
      kytbl(kwDX,icBEAM)=ky_DX_BEAM
      kytbl(kwDPX,icBEAM)=ky_DPX_BEAM
      kytbl(kwDY,icBEAM)=ky_DY_BEAM
      kytbl(kwDPY,icBEAM)=ky_DPY_BEAM
      kytbl(kwXANGLE,icBEAM)=ky_XANGLE_BEAM
      kytbl(kwEMIX,icBEAM)=ky_EMIX_BEAM
      kytbl(kwEMIY,icBEAM)=ky_EMIY_BEAM
      kytbl(kwDP,icBEAM)=ky_DP_BEAM
      kytbl(kwAZ,icBEAM)=ky_AZ_BEAM
      kytbl(kwDZ,icBEAM)=ky_DZ_BEAM
      kytbl(kwSIGZ,icBEAM)=ky_SIGZ_BEAM
      kytbl(kwSLI,icBEAM)=ky_SLI_BEAM
      kytbl(kwNP,icBEAM)=ky_NP_BEAM
      kytbl(kwSTURN,icBEAM)=ky_STURN_BEAM
      kytbl(kwR11,icBEAM)=ky_R11_BEAM
      kytbl(kwR12,icBEAM)=ky_R12_BEAM
      kytbl(kwR13,icBEAM)=ky_R13_BEAM
      kytbl(kwR14,icBEAM)=ky_R14_BEAM
      kytbl(kwR15,icBEAM)=ky_R15_BEAM
      kytbl(kwR16,icBEAM)=ky_R16_BEAM
      kytbl(kwR22,icBEAM)=ky_R22_BEAM
      kytbl(kwR23,icBEAM)=ky_R23_BEAM
      kytbl(kwR24,icBEAM)=ky_R24_BEAM
      kytbl(kwR25,icBEAM)=ky_R25_BEAM
      kytbl(kwR26,icBEAM)=ky_R26_BEAM
      kytbl(kwR33,icBEAM)=ky_R33_BEAM
      kytbl(kwR34,icBEAM)=ky_R34_BEAM
      kytbl(kwR35,icBEAM)=ky_R35_BEAM
      kytbl(kwR36,icBEAM)=ky_R36_BEAM
      kytbl(kwR44,icBEAM)=ky_R44_BEAM
      kytbl(kwR45,icBEAM)=ky_R45_BEAM
      kytbl(kwR46,icBEAM)=ky_R46_BEAM
      kytbl(kwR55,icBEAM)=ky_R55_BEAM
      kytbl(kwR56,icBEAM)=ky_R56_BEAM
      kytbl(kwR66,icBEAM)=ky_R66_BEAM
      kytbl(kwCOUPLE,icBEAM)=ky_COUPLE_BEAM
      kytbl(kwBSTRL,icBEAM)=ky_BSTRL_BEAM
      kytbl(kwMAX,icBEAM)=ky_MAX_BEAM
      kytbl(kwNPARAM,icBEAM)=p_NPARAM_BEAM
cc for PHSROT
      idummy=sethtb('phsrot  ',icDEF,icProt)
      kytbl(0,icProt)=sethtb('PHSROT  ',icDEF,icProt)
      kytbl(kwAX,icProt)=ky_AX_Prot
      kytbl(kwBX,icProt)=ky_BX_Prot
      kytbl(kwPX,icProt)=ky_PX_Prot
      kytbl(kwAY,icProt)=ky_AY_Prot
      kytbl(kwBY,icProt)=ky_BY_Prot
      kytbl(kwPY,icProt)=ky_PY_Prot
      kytbl(kwR1,icProt)=ky_R1_Prot
      kytbl(kwR2,icProt)=ky_R2_Prot
      kytbl(kwR3,icProt)=ky_R3_Prot
      kytbl(kwR4,icProt)=ky_R4_Prot
      kytbl(kwEX,icProt)=ky_EX_Prot
      kytbl(kwEPX,icProt)=ky_EPX_Prot
      kytbl(kwEY,icProt)=ky_EY_Prot
      kytbl(kwEPY,icProt)=ky_EPY_Prot
      kytbl(kwZX,icProt)=ky_ZX_Prot
      kytbl(kwZPX,icProt)=ky_ZPX_Prot
      kytbl(kwZY,icProt)=ky_ZY_Prot
      kytbl(kwZPY,icProt)=ky_ZPY_Prot
      kytbl(kwEMIX,icProt)=ky_EMIX_Prot
      kytbl(kwEMIY,icProt)=ky_EMIY_Prot
      kytbl(kwDP,icProt)=ky_DP_Prot
      kytbl(kwAZ,icProt)=ky_AZ_Prot
      kytbl(kwDZ,icProt)=ky_DZ_Prot
      kytbl(kwSIGZ,icProt)=ky_SIGZ_Prot
      kytbl(kwPZ,icProt)=ky_PZ_Prot
      kytbl(kwEMIZ,icProt)=ky_EMIZ_Prot
      kytbl(kwBZ,icProt)=ky_BZ_Prot
      kytbl(kwJDY,icProt)=ky_JDY_Prot
      kytbl(kwJDPY,icProt)=ky_JDPY_Prot
      kytbl(kwD11,icProt)=ky_D11_Prot
      kytbl(kwD12,icProt)=ky_D12_Prot
      kytbl(kwD13,icProt)=ky_D13_Prot
      kytbl(kwD14,icProt)=ky_D14_Prot
      kytbl(kwD15,icProt)=ky_D15_Prot
      kytbl(kwD16,icProt)=ky_D16_Prot
      kytbl(kwD21,icProt)=ky_D21_Prot
      kytbl(kwD22,icProt)=ky_D22_Prot
      kytbl(kwD23,icProt)=ky_D23_Prot
      kytbl(kwD24,icProt)=ky_D24_Prot
      kytbl(kwD25,icProt)=ky_D25_Prot
      kytbl(kwD26,icProt)=ky_D26_Prot
      kytbl(kwD31,icProt)=ky_D31_Prot
      kytbl(kwD32,icProt)=ky_D32_Prot
      kytbl(kwD33,icProt)=ky_D33_Prot
      kytbl(kwD34,icProt)=ky_D34_Prot
      kytbl(kwD35,icProt)=ky_D35_Prot
      kytbl(kwD36,icProt)=ky_D36_Prot
      kytbl(kwD41,icProt)=ky_D41_Prot
      kytbl(kwD42,icProt)=ky_D42_Prot
      kytbl(kwD43,icProt)=ky_D43_Prot
      kytbl(kwD44,icProt)=ky_D44_Prot
      kytbl(kwD45,icProt)=ky_D45_Prot
      kytbl(kwD46,icProt)=ky_D46_Prot
      kytbl(kwD51,icProt)=ky_D51_Prot
      kytbl(kwD52,icProt)=ky_D52_Prot
      kytbl(kwD53,icProt)=ky_D53_Prot
      kytbl(kwD54,icProt)=ky_D54_Prot
      kytbl(kwD55,icProt)=ky_D55_Prot
      kytbl(kwD56,icProt)=ky_D56_Prot
      kytbl(kwD61,icProt)=ky_D61_Prot
      kytbl(kwD62,icProt)=ky_D62_Prot
      kytbl(kwD63,icProt)=ky_D63_Prot
      kytbl(kwD64,icProt)=ky_D64_Prot
      kytbl(kwD65,icProt)=ky_D65_Prot
      kytbl(kwD66,icProt)=ky_D66_Prot
      kytbl(kwB11,icProt)=ky_B11_Prot
      kytbl(kwB12,icProt)=ky_B12_Prot
      kytbl(kwB13,icProt)=ky_B13_Prot
      kytbl(kwB14,icProt)=ky_B14_Prot
      kytbl(kwB15,icProt)=ky_B15_Prot
      kytbl(kwB16,icProt)=ky_B16_Prot
      kytbl(kwB22,icProt)=ky_B22_Prot
      kytbl(kwB23,icProt)=ky_B23_Prot
      kytbl(kwB24,icProt)=ky_B24_Prot
      kytbl(kwB25,icProt)=ky_B25_Prot
      kytbl(kwB26,icProt)=ky_B26_Prot
      kytbl(kwB33,icProt)=ky_B33_Prot
      kytbl(kwB34,icProt)=ky_B34_Prot
      kytbl(kwB35,icProt)=ky_B35_Prot
      kytbl(kwB36,icProt)=ky_B36_Prot
      kytbl(kwB44,icProt)=ky_B44_Prot
      kytbl(kwB45,icProt)=ky_B45_Prot
      kytbl(kwB46,icProt)=ky_B46_Prot
      kytbl(kwB55,icProt)=ky_B55_Prot
      kytbl(kwB56,icProt)=ky_B56_Prot
      kytbl(kwB66,icProt)=ky_B66_Prot
      kytbl(kwCOUPLE,icProt)=ky_COUPLE_Prot
      kytbl(kwMAX,icProt)=ky_MAX_Prot
      kytbl(kwNPARAM,icProt)=p_NPARAM_Prot
cc for MARK
      kytbl(0,icMARK)=sethtb('mark    ',icDEF,icMark)
      kytbl(0,icMARK)=sethtb('MARK    ',icDEF,icMark)
      kytbl(kwAX,icMark)=ky_AX_Mark
      kytbl(kwBX,icMark)=ky_BX_Mark
      kytbl(kwPX,icMark)=ky_PX_Mark
      kytbl(kwAY,icMark)=ky_AY_Mark
      kytbl(kwBY,icMark)=ky_BY_Mark
      kytbl(kwPY,icMark)=ky_PY_Mark
      kytbl(kwEX,icMark)=ky_EX_Mark
      kytbl(kwEPX,icMark)=ky_EPX_Mark
      kytbl(kwEY,icMark)=ky_EY_Mark
      kytbl(kwEPY,icMark)=ky_EPY_Mark
      kytbl(kwR1,icMark)=ky_R1_Mark
      kytbl(kwR2,icMark)=ky_R2_Mark
      kytbl(kwR3,icMark)=ky_R3_Mark
      kytbl(kwR4,icMark)=ky_R4_Mark
      kytbl(kwDETR,icMark)=ky_DETR_Mark
      kytbl(kwDX,icMark)=ky_DX_Mark
      kytbl(kwDPX,icMark)=ky_DPX_Mark
      kytbl(kwDY,icMark)=ky_DY_Mark
      kytbl(kwDPY,icMark)=ky_DPY_Mark
      kytbl(kwDZ,icMark)=ky_DZ_Mark
      kytbl(kwDDP,icMark)=ky_DDP_Mark
      kytbl(kwAZ,icMark)=ky_AZ_Mark
      kytbl(kwBZ,icMark)=ky_BZ_Mark
      kytbl(kwPZ,icMark)=ky_PZ_Mark
      kytbl(kwZX,icMark)=ky_ZX_Mark
      kytbl(kwZPX,icMark)=ky_ZPX_Mark
      kytbl(kwZY,icMark)=ky_ZY_Mark
      kytbl(kwZPY,icMark)=ky_ZPY_Mark
      kytbl(kwDP,icMark)=ky_DP_Mark
      kytbl(kwOFFSET,icMark)=ky_OFFSET_Mark
      kytbl(kwSIGZ,icMark)=ky_SIGZ_Mark
      kytbl(kwSIGE,icMark)=ky_SIGE_Mark
      kytbl(kwGEO,icMark)=ky_GEO_Mark
      kytbl(kwJDX,icMark)=ky_JDX_Mark
      kytbl(kwJDPX,icMark)=ky_JDPX_Mark
      kytbl(kwJDY,icMark)=ky_JDY_Mark
      kytbl(kwJDPY,icMark)=ky_JDPY_Mark
      kytbl(kwJDZ,icMark)=ky_JDZ_Mark
      kytbl(kwJDPZ,icMark)=ky_JDPZ_Mark
      kytbl(kwEMIX,icMark)=ky_EMIX_Mark
      kytbl(kwEMIY,icMark)=ky_EMIY_Mark
      kytbl(kwEMIZ,icMark)=ky_EMIZ_Mark
      kytbl(kwCOUPLE,icMark)=ky_COUPLE_Mark
      kytbl(kwMAX,icMark)=ky_MAX_Mark
cc for apert
      kytbl(0,icAprt)=sethtb('apert   ',icDEF,icAprt)
      kytbl(0,icAprt)=sethtb('APERT   ',icDEF,icAprt)
      kytbl(kwDX1,icAprt)=ky_DX1_Aprt
      kytbl(kwDX2,icAprt)=ky_DX2_Aprt
      kytbl(kwDY1,icAprt)=ky_DY1_Aprt
      kytbl(kwDY2,icAprt)=ky_DY2_Aprt
      kytbl(kwJDPX,icAprt)=ky_JDPX_Aprt
      kytbl(kwJDPY,icAprt)=ky_JDPY_Aprt
      kytbl(kwDP,icAprt)=ky_DP_Aprt
      kytbl(kwCOUPLE,icAprt)=ky_COUPLE_Aprt
      kytbl(kwROT,icAprt)=ky_ROT_Aprt
      kytbl(kwAX,icAprt)=ky_AX_Aprt
      kytbl(kwAY,icAprt)=ky_AY_Aprt
      kytbl(kwDX,icAprt)=ky_DX_Aprt
      kytbl(kwDY,icAprt)=ky_DY_Aprt      
      kytbl(kwMAX,icAprt)=ky_MAX_Aprt
cc for mon
      kytbl(0,icMONI)=sethtb('moni    ',icDEF,icMONI)
      kytbl(0,icMONI)=sethtb('MONI    ',icDEF,icMONI)
      kytbl(kwDX,icMONI)=ky_DX_MONI
      kytbl(kwDY,icMONI)=ky_DY_MONI
      kytbl(kwOFFSET,icMONI)=ky_OFFSET_MONI
      kytbl(kwCOUPLE,icMONI)=ky_COUPLE_MONI
      kytbl(kwROT,icMONI)=ky_ROT_MONI
      kytbl(kwMAX,icMONI)=ky_MAX_MONI
cc for spch
      idummy=sethtb('spch    ',icDEF,icSPCH)
      kytbl(0,icSPCH)=sethtb('SPCH    ',icDEF,icSPCH)
      kytbl(kwAX,icSPCH)=ky_AX_SPCH
      kytbl(kwBX,icSPCH)=ky_BX_SPCH
      kytbl(kwPX,icSPCH)=ky_PX_SPCH
      kytbl(kwAY,icSPCH)=ky_AY_SPCH
      kytbl(kwBY,icSPCH)=ky_BY_SPCH
      kytbl(kwPY,icSPCH)=ky_PY_SPCH
      kytbl(kwR1,icSPCH)=ky_R1_SPCH
      kytbl(kwR2,icSPCH)=ky_R2_SPCH
      kytbl(kwR3,icSPCH)=ky_R3_SPCH
      kytbl(kwR4,icSPCH)=ky_R4_SPCH
      kytbl(kwEX,icSPCH)=ky_EX_SPCH
      kytbl(kwEPX,icSPCH)=ky_EPX_SPCH
      kytbl(kwEY,icSPCH)=ky_EY_SPCH
      kytbl(kwEPY,icSPCH)=ky_EPY_SPCH
      kytbl(kwZX,icSPCH)=ky_ZX_SPCH
      kytbl(kwZPX,icSPCH)=ky_ZPX_SPCH
      kytbl(kwZY,icSPCH)=ky_ZY_SPCH
      kytbl(kwZPY,icSPCH)=ky_ZPY_SPCH
      kytbl(kwMAX,icSPCH)=ky_MAX_SPCH

      call initkyindex

c.....for debug
c     do 9999 i=1,kwMAX
c       print *,i,kytbl(i,0),pname(kytbl(i,0))
c9999 continue
c     do 9998 i=1,icQUAD
c       do 9998 j=0,kwMAX
c       print *,j,i,kytbl(j,i),pname(kytbl(0,i)),pname(kytbl(j,0))
c9998 continue
c.....end debug
c
       rlist(idval(sethtb('%       ',icUNIT,ktcaloc(3))))=0.01d0
       rlist(idval(sethtb('rad     ',icUNIT,ktcaloc(3))))=1.00d0
       rlist(idval(sethtb('RAD     ',icUNIT,ktcaloc(3))))=1.00d0
       rlist(idval(sethtb('mrad    ',icUNIT,ktcaloc(3))))=1.00d-3
       rlist(idval(sethtb('MRAD    ',icUNIT,ktcaloc(3))))=1.00d-3
       rlist(idval(sethtb('DEG     ',icUNIT,ktcaloc(3))))=pi/180.d0
       rlist(idval(sethtb('deg     ',icUNIT,ktcaloc(3))))=pi/180.d0
       rlist(idval(sethtb('M       ',icUNIT,ktcaloc(3))))=1.00d0
       rlist(idval(sethtb('m       ',icUNIT,ktcaloc(3))))=1.00d0
       rlist(idval(sethtb('cm      ',icUNIT,ktcaloc(3))))=1.00d-2
       rlist(idval(sethtb('CM      ',icUNIT,ktcaloc(3))))=1.0d-2
       rlist(idval(sethtb('mm      ',icUNIT,ktcaloc(3))))=1.0d-3
       rlist(idval(sethtb('MM      ',icUNIT,ktcaloc(3))))=1.0d-3
       rlist(idval(sethtb('T       ',icUNIT,ktcaloc(3))))=1.0d0
       rlist(idval(sethtb('GAUSS   ',icUNIT,ktcaloc(3))))=1.0d-4
       rlist(idval(sethtb('gauss   ',icUNIT,ktcaloc(3))))=1.0d-4
c      Standard energy unit wasJoule.(MKSA unit system)
c      Standard energy unit is Now eV.
       rlist(idval(sethtb('JOULE   ',icUNIT,ktcaloc(3))))
     &                 = 1.0d0/elemch
       rlist(idval(sethtb('eV      ',icUNIT,ktcaloc(3))))
     &                 = 1.0d0
       rlist(idval(sethtb('KeV     ',icUNIT,ktcaloc(3))))
     &                 = 1.d3
       rlist(idval(sethtb('MeV     ',icUNIT,ktcaloc(3))))
     &                 = 1.d6
       rlist(idval(sethtb('GeV     ',icUNIT,ktcaloc(3))))
     &                 = 1.d9
       rlist(idval(sethtb('EV      ',icUNIT,ktcaloc(3))))
     &                 = 1.0d0
       rlist(idval(sethtb('KEV     ',icUNIT,ktcaloc(3))))
     &                 = 1.d3
       rlist(idval(sethtb('MEV     ',icUNIT,ktcaloc(3))))
     &                 = 1.d6
       rlist(idval(sethtb('GEV     ',icUNIT,ktcaloc(3))))
     &                 = 1.d9
c
       rlist(idval(sethtb('V       ',icUNIT,ktcaloc(3))))
     &                 = 1.d0
       rlist(idval(sethtb('KV      ',icUNIT,ktcaloc(3))))
     &                 = 1.d3
       rlist(idval(sethtb('MV      ',icUNIT,ktcaloc(3))))
     &                 = 1.d6
       rlist(idval(sethtb('GV      ',icUNIT,ktcaloc(3))))
     &                 = 1.d9
c
       rlist(idval(sethtb('HZ      ',icUNIT,ktcaloc(3))))
     &                 = 1.d0
       rlist(idval(sethtb('KHZ     ',icUNIT,ktcaloc(3))))
     &                 = 1.d3
       rlist(idval(sethtb('MHZ     ',icUNIT,ktcaloc(3))))
     &                 = 1.d6
       rlist(idval(sethtb('GHZ     ',icUNIT,ktcaloc(3))))
     &                 = 1.d9
c
       idummy=sethtb('normal  ',icRAND,1 )
       idummy=sethtb('NORMAL  ',icRAND,1 )
       idummy=sethtb('uniform ',icRAND,2 )
       idummy=sethtb('UNIFORM ',icRAND,2 )
c
       call defglb('STACKSIZ',icGLR,idummy)
c       call get_environment_variable('SAD_STACKSIZ',stacksiz)
c       read(stacksiz,'(i20)')iss
c       call RsetGL('STACKSIZ',max(2d0**18,dble(iss)),idummy)
       call rsetGL('STACKSIZ',dble(maxstack),idummy)
       call defglb('$SEED',icGLR,idummy)
       call RsetGL('$SEED',17.d0,idummy)
       call defglb('SEED',icGLR,idummy)
       call RsetGL('SEED',17.d0,idummy)
       call defglb('PI',icGLR,idummy)
       call RsetGL('PI',asin(1.d0)*2.d0,idummy)
       call defglb('INF',icGLR,idummy)
       call RsetGL('INF',dinfinity,idummy)
       call defglb('INF.',icGLR,idummy)
       call RsetGL('INF.',dinfinity,idummy)
       call defglb('NaN',icGLR,idummy)
       call RsetGL('NaN',dnotanumber,idummy)
       call defglb('NAN',icGLR,idummy)
       call RsetGL('NAN',dnotanumber,idummy)
       call defglb('MOMENTUM',icGLR,idummy)
       call RsetGL('MOMENTUM',30.d9,idummy)
       call defglb('PBUNCH',icGLR,idummy)
       call RsetGL('PBUNCH',1.d10,idummy)
       call defglb('NBUNCH',icGLR,idummy)
       call RsetGL('NBUNCH',1.d0,idummy)
       call defglb('EMITDIV',icGLR,idummy)
       call RsetGL('EMITDIV',1.d0,idummy)
       call defglb('EMITDIVB',icGLR,idummy)
       call RsetGL('EMITDIVB',1.d0,idummy)
       call defglb('EMITDIVQ',icGLR,idummy)
       call RsetGL('EMITDIVQ',1.d0,idummy)
       call defglb('EMITDIVS',icGLR,idummy)
       call RsetGL('EMITDIVS',1.d0,idummy)
       call defglb('FRINGDIV',icGLR,idummy)
       call RsetGL('FRINGDIV',1.d0,idummy)
       call defglb('MINCOUP',icGLR,idummy)
       call RsetGL('MINCOUP',1.d-2,idummy)
       call defglb('LOSSAMPL',icGLR,idummy)
       call RsetGL('LOSSAMPL',1.d0,idummy)
       call defglb('LOSSDZ',icGLR,idummy)
       call RsetGL('LOSSDZ',100.d0,idummy)
       call defglb('BBCUT',icGLR,idummy)
       call RsetGL('BBCUT',0.d0,idummy)
       call defglb('EMITX',icGLR,idummy)
       call RsetGL('EMITX',0.0d0,idummy)
       call defglb('EMITY',icGLR,idummy)
       call RsetGL('EMITY',0.0d0,idummy)
       call defglb('EMITZ',icGLR,idummy)
       call RsetGL('EMITZ',0.0d0,idummy)
       call defglb('EMITXE',icGLR,idummy)
       call RsetGL('EMITXE',0.0d0,idummy)
       call defglb('EMITYE',icGLR,idummy)
       call RsetGL('EMITYE',0.0d0,idummy)
       call defglb('EMITZE',icGLR,idummy)
       call RsetGL('EMITZE',0.0d0,idummy)
       call defglb('EMITXC',icGLR,idummy)
       call RsetGL('EMITXC',0.0d0,idummy)
       call defglb('EMITYC',icGLR,idummy)
       call RsetGL('EMITYC',0.0d0,idummy)
       call defglb('EMITZC',icGLR,idummy)
       call RsetGL('EMITZC',0.0d0,idummy)
       call defglb('SIGE',icGLR,idummy)
       call RsetGL('SIGE',0.0d0,idummy)
       call defglb('SIGZ',icGLR,idummy)
       call RsetGL('SIGZ',0.0d0,idummy)
       call defglb('GCUT',icGLR,idummy)
       call RsetGL('GCUT',1.d300,idummy)
       call defglb('TDXI',icGLR,idummy)
       call RsetGL('TDXI',0.0d0,idummy)
       call defglb('TDYI',icGLR,idummy)
       call RsetGL('TDYI',0.0d0,idummy)
       call defglb('TDZI',icGLR,idummy)
       call RsetGL('TDZI',0.0d0,idummy)
       call defglb('$MASS$',icGLR,idummy)
       call RsetGL('$MASS$',elmass,idummy)
       call defglb('MASS',icGLR,idummy)
       call RsetGL('MASS',elmass,idummy)
       call defglb('NP',icGLR,idummy)
       call RsetGL('NP',0.d0,idummy)
       call defglb('TURNS',icGLR,idummy)
       call RsetGL('TURNS',0.d0,idummy)
       call defglb('CHARGE',icGLR,idummy)
       call RsetGL('CHARGE',0.d0,idummy)
c
       call defglb('OMEGA0',icGLR,idummy)
       call RsetGL('OMEGA0',0.0d0,idummy)
       call defglb('DTSYNCH',icGLR,idummy)
       call RsetGL('DTSYNCH',0.0d0,idummy)
       call defglb('PHICAV',icGLR,idummy)
       call RsetGL('PHICAV',0.0d0,idummy)
       call defglb('EFFVCRATIO',icGLR,idummy)
       call RsetGL('EFFVCRATIO',1.0d0,idummy)
       call defglb('EFFVC',icGLR,idummy)
       call RsetGL('EFFVC',1.0d0,idummy)
       call defglb('EFFRFFREQ',icGLR,idummy)
       call RsetGL('EFFRFFREQ',1.0d0,idummy)
       call defglb('FSHIFT',icGLR,idummy)
       call RsetGL('FSHIFT',0.0d0,idummy)
       call defglb('PSPAN',icGLR,idummy)
       call RsetGL('PSPAN',0.0d0,idummy)
       call defglb('ESHIFT',icGLR,idummy)
       call RsetGL('ESHIFT',0.0d0,idummy)
       call defglb('PHIS',icGLR,idummy)
       call RsetGL('PHIS',0.D00,idummy)
       call defglb('MAXORDER',icGLR,idummy)
       call RsetGL('MAXORDER',6.D00,idummy)
       call defglb('ADDDENSE',icGLR,idummy)
       call RsetGL('ADDDENSE',6.D00,idummy)
       call defglb('ADDTERMS',icGLR,idummy)
       call RsetGL('ADDTERMS',0.D00,idummy)
       call defglb('DJPLOT',icGLR,idummy)
       call RsetGL('DJPLOT',0.D00,idummy)
       call defglb('PHSPLOTS',icGLR,idummy)
       call RsetGL('PHSPLOTS',0.D00,idummy)
       call defglb('DAPWIDTH',icGLR,idummy)
       call RsetGL('DAPWIDTH',7.D00,idummy)
c       call defglb('RADIAT',icGLI,idummy)
c       call IsetGL('RADIAT',-1,idummy)
       call defglb('PSPACNX',icGLR,idummy)
       call RsetGL('PSPACNX',128.0D0,idummy)
       call defglb('PSPACNY',icGLR,idummy)
       call RsetGL('PSPACNY',128.0D0,idummy)
       call defglb('PSPACNZ',icGLR,idummy)
       call RsetGL('PSPACNZ',1.0D0,idummy)
       call defglb('PSPACDX',icGLR,idummy)
       call RsetGL('PSPACDX',1.0d-3,idummy)
       call defglb('PSPACDY',icGLR,idummy)
       call RsetGL('PSPACDY',1.0d-3,idummy)
       call defglb('PSPACDZ',icGLR,idummy)
       call RsetGL('PSPACDZ',1.0d-3,idummy)
       call defglb('PSPACNTURN',icGLR,idummy)
       call RsetGL('PSPACNTURN',1.0D0,idummy)
       call defglb('PSPACNTURNCALC',icGLR,idummy)
       call RsetGL('PSPACNTURNCALC',0.0D0,idummy)
       call defglb('BSTRL',icGLR,idummy)
       call RsetGL('BSTRL',0.d0,idummy)
       call defglb('NPARA',icGLR,idummy)
       call RsetGL('NPARA',1.d0,idummy)
c
c      call defglb('$RADIZZZ',icGLI,idummy)
c      call IsetGL('$RADIZZZ',-1,idummy)
c
       call defflg('RAD',FLAGOF)
       call defflg('RFSW',FLAGOF)
       call defflg('RADCOD',FLAGOF)
       call defflg('COD',FLAGON)

       call defflg('INTRA',FLAGOF)
       call defflg('TRPT',FLAGOF)
       call defflg('EMIOUT',FLAGOF)
       call defflg('GAUSS',FLAGOF)

       call defflg('PHOTONS',FLAGOF)
       call defflg('EMIT',FLAGON)
       call defflg('POL'  ,FLAGOF)
       call defflg('RADPOL',FLAGON)
       call defflg('DAPERT',FLAGOF)
       call defflg('FLUC',FLAGON)
       call defflg('K64' ,FLAGON)
       call defflg('FOURIE',FLAGOF)
       call defflg('SMEAR',FLAGON)
       call defflg('GEOCAL',FLAGON)
       call defflg('CALC6D',FLAGOF)
       call defflg('DEBUG',FLAGOF)
       call defflg('ECHO',FLAGOF)
       call defflg('LOG',FLAGOF)
       call defflg('CTIME',FLAGOF)
       call defflg('INTRES',FLAGON)
       call defflg('HALFRES',FLAGON)
       call defflg('SUMRES',FLAGON)
       call defflg('DIFFRES',FLAGOF)
       msglvl=0
       return
       end
      
