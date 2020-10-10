# column spec:
# parameter number (index)  Phase 1                 Phase 2                 Phase 1                 Phase 2                 P1 fixed/   P2 fixed/   enabled/
#   parameter name          initial value           initial value           constraint value        constraint value        notfixed    notfixed    disabled
#                                       initial error       initial error               constraint error        contraint error
#
#   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38
#   1                       2           3           4           5           6           7           8           9           10          11
#   name                    Phase 1                 Phase 2                 Phase 1                 Phase2                  P1 fixed    P2 fixed    enabled
#                           initial     initial     initial     initial     constraint  constraint  constraint  constraint  spec        spec        spec
#                           value       error       value       error       value       error       value       error
BREAK
# TODO: a 3rd value is required, fixed, free and float
# or                             hard, soft, free
# fixed is hard fixed value, cannot change
# float is soft fixed (constrained value)
# free is completely free, no fixed value or penalty term applied
# none is a placeholder for 0.0
# doesn't do anything just signifies that this value is N/A or ignored
BREAK
# first col is the parameter index/number
# second col is a list of mc files which are combined in the same fit
BREAK
# notes regarding where values are used:
# code makes use of initial and constrained values in 2 places, these are
# 1) when defining the parameters
# 2) when computing the loglikelihood value
#
# in the case of 1):
# parameter is initialized using initial value & error, if constrain mode
# is either free or soft
# if constrain mode is hard, then parameter is intilized using constrain
# value & error
#
# in the case of 2):
# if parameter is hard constrained, then Minuit does not vary that parameter
# in the fit, therefore it is essentially ignored
# if parameter is free, Minuit varies the paramter internally, no penalty term
# is applied
# if parameter is soft constrained, then penalty terms are applied in the
# logliklihood function, using the values and errors from the constrained
# inputs
# while at the same time, Minuit varies the parameters that were initialized
# using the initial parameter value & error
#
# simplified instructions:
#
# 1: If parameter is "free", set phase 1/2 initial value & error, set phase 1/2
#    constraint value & error to "none"
# 2: If parameter is "soft" constrained, set phase 1/2 initial value & error,
#    set phase 1/2 constraint value & error
# 3: If parameter is "hard" constrained, set phase 1/2 initial value & error to
#    the values to be used for the constraining value & error, then set the
#    constraint value & error to "useinit"
# this is counter intuitive, but is just how the code works at the moment, it
# requires an update
# TODO
BREAK
#   enabled/    phase 1                 phase 2                 phase 1                             phase 2                             list
#   disabled    initial                 initial                 constraint                          constraint                          of
#               value       error       value       error       mode        value       error       mode        value       error       MC files
#   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39  40
#example:
#0	enabled     3.5350e-04	5.9777e-06  same        same        free        none        none        free        none        none        nd150_rot_2n2b_m4       END
BREAK
#   list options:
#   "enabled"   [value]     [value]     [value]     [value]     [value]     [value]     [value]     [value]     "free"  "free"  [MC datafile name]              "END"
#   "disabled"                          "same"      "same"      "none"      "none"      "none"      "none"      "soft"  "soft"  list separated by space/tab
#                                                               "useinit"   "useinit"   "same"      "same"      "hard"  "hard"
#
#   ---------------------------------------------------------------------------------------------------------------------------------------------------------------
#   enabled/    phase 1                 phase 2                 phase 1                 phase 2                 P1      P2      list
#   disabled    initial                 initial                 constraint              constraint              cnstrt  cnstrt  of
#               value       error       value       error       value       error       value       error       mode    mode    MC files
#   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31  32  33  34  35  36  37  38  39
#
# TODO: check these
#0	enabled     3.65e-04   	6.0e-06     same        same        none        none        none        none        free    free    nd150_rot_2n2b_m4               END
#0	enabled     4.10996e-04 6.19840e-06 same        same        none        none        none        none        free    free    nd150_rot_2n2b_m4               END
0	enabled     3.45e-04    3.45e-05    same        same        none        none        none        none        free    free    nd150_rot_2n2b_m4               END
#0	enabled     1.459       0.1         same        same        none        none        none        none        hard    hard    nd150_rot_2n2b_m4               END
BREAK
# TODO: bug: line with whitespace crashes
BREAK
1	enabled     1.70e-3   	0.02e-3     same        same        1.70e-3   	0.02e-3     same        same        soft    soft    ac228_int_rot   bi212_int_rot   END
2   enabled     0.62e-3     0.01e-3     same        same        0.62e-3     0.01e-3     same        same        soft    soft    tl208_int_rot                   END 
3   enabled     5.66e-3     0.15e-3     same        same        5.66e-3     0.15e-3     same        same        soft    soft    bi207_int_rot  	                END
4   enabled     0.29e-3     0.01e-3     same        same        0.29e-3     0.01e-3     same        same        soft    soft    bi214_int_rot   pb214_int_rot   END
5   enabled  	2.98e-3     0.21e-3     same        same        2.98e-3     0.21e-3     same        same        soft    soft    eu152_int_rot                   END
6   enabled     1.08e-3     0.25e-3     same        same        1.08e-3     0.25e-3     same        same        soft    soft    eu154_int_rot  	                END
7   enabled     10.14e-3    0.08e-3     same        same        10.14e-3    0.08e-3     same        same        soft    soft    k40_int_rot	                    END
8   enabled     1.54e-3     0.04e-3     same        same        1.54e-3     0.04e-3     same        same        soft    soft    pa234m_int_rot                  END
BREAK
9   enabled     1.6e-3      0.7e-3      same        same        1.6e-3      0.7e-3      same        same        soft    soft    bi214_mylar     pb214_mylar     END
BREAK
10  enabled     4.86e-03    0.04e-05    same        same        4.86e-03    0.04e-05    same        same        soft    soft    mo100_99_rot_2n2b_m14           END
11  enabled     0.134e-03   0.002e-03   same        same        0.134e-03   0.002e-03   same        same        soft    soft    mo100_99_rot_bi214	            END
12  enabled     0.73e-03    0.03e-03    same        same        0.73e-03    0.03e-03    same        same        soft    soft    mo100_99_rot_pa234m	            END
13  enabled     4.27e-03    0.07e-03    same        same        4.27e-03    0.07e-03    same        same        soft    soft    mo100_99_rot_k40	            END
BREAK
14  enabled     0.041e-03	0.003e-03   same        same        0.041e-03	0.003e-03   same        same        soft    soft    zr96_rot_2n2b_m4	            END
15  enabled     0.19e-03    0.02e-03    same        same        0.19e-03    0.02e-03    same        same        soft    soft    zr96_rot_bi214	                END
16  enabled     0.49e-03    0.01e-03    same        same        0.49e-03    0.01e-03    same        same        soft    soft    zr96_rot_pa234m	                END
17  enabled     19.7e-03	0.1e-03     same        same        19.7e-03	0.1e-03     same        same        soft    soft    zr96_rot_k40	                END 
BREAK
18  enabled     0.031e-03	0.003e-03   same        same        0.031e-03	0.003e-03   same        same        soft    soft    ca48_63_rot_2n2b_m4	            END
19  enabled     0.08e-03    0.01e-03    same        same        0.08e-03    0.01e-03    same        same        soft    soft    ca48_63_rot_bi214	            END
20  enabled     0.3e-03     0.1e-03     same        same        0.3e-03     0.1e-03     same        same        soft    soft    ca48_63_rot_pa234m	            END
21  enabled     29.6e-03    0.1e-03     same        same        29.6e-03    0.1e-03     same        same        soft    soft    ca48_63_rot_y90	                END
BREAK
#10  enabled     4.86e-03    0.04e-03    same        same        4.86e-03    0.04e-03    same        same        hard    hard    mo100_99_rot_2n2b_m14           END
#11  disabled    0.134e-03   0.002e-03   same        same        0.134e-03   0.002e-03   same        same        hard    hard    mo100_99_rot_bi214	            END
#12  disabled    0.73e-03    0.03e-03    same        same        0.73e-03    0.03e-03    same        same        hard    hard    mo100_99_rot_pa234m	            END
#13  disabled    4.27e-03    0.07e-03    same        same        4.27e-03    0.07e-03    same        same        hard    hard    mo100_99_rot_k40	            END
BREAK
#14  disabled    0.041e-03	0.003e-03   same        same        0.041e-03	0.003e-03   same        same        hard    hard    zr96_rot_2n2b_m4	            END
#15  disabled    0.19e-03    0.02e-03    same        same        0.19e-03    0.02e-03    same        same        hard    hard    zr96_rot_bi214	                END
#16  disabled    0.49e-03    0.01e-03    same        same        0.49e-03    0.01e-03    same        same        hard    hard    zr96_rot_pa234m	                END
#17  disabled    19.7e-03	0.1e-03     same        same        19.7e-03	0.1e-03     same        same        hard    hard    zr96_rot_k40	                END 
BREAK
#18  disabled    0.031e-03	0.003e-03   same        same        0.031e-03	0.003e-03   same        same        hard    hard    ca48_63_rot_2n2b_m4	            END
#19  disabled    0.08e-03    0.01e-03    same        same        0.08e-03    0.01e-03    same        same        hard    hard    ca48_63_rot_bi214	            END
#20  disabled    0.3e-03     0.1e-03     same        same        0.3e-03     0.1e-03     same        same        hard    hard    ca48_63_rot_pa234m	            END
#21  disabled    29.6e-03    0.1e-03     same        same        29.6e-03    0.1e-03     same        same        hard    hard    ca48_63_rot_y90	                END
BREAK
# MC is not split corrently TODO
# TODO: make soft after initial debug run
# TODO what numbers here - disabled?
22  disabled    1.0e-3      0.6e-3      same        same        useinit     useinit     same        same        hard    hard    bi214_sfoil_rot pb214_sfoil	    END
#11  disabled    0.0091	    0.0048      0.00341     0.00096     useinit     useinit     same        same        soft    soft    bi214_sfoil_rot                 END
#11  disabled    0.0091	    0.0048      0.00341	    0.00096     useinit     useinit     same        same        soft    soft    pb214_sfoil	                    END
BREAK
# unsure about these numbers, and MC is not split correctly TODO
# TODO what numbers here - disabled?
# TODO: changed error to smaller value, MINUIT complaining about larger error
#23  enabled     1.3e-3      1.6e-3      same        same        useinit     useinit     same        same        soft    soft    bi214_swire     pb214_swire 	END
23  enabled     1.3e-3      1.0e-3      same        same        useinit     useinit     same        same        soft    soft    bi214_swire     pb214_swire 	END
#12  disabled    0.556	    0.058       0.085	    0.012       useinit     useinit     same        same        soft    soft    bi214_swire 	                END
#12  disabled    0.556	    0.058       0.085       0.012       useinit     useinit     same        same        soft    soft    pb214_swire 	                END
BREAK
# TODO think this is neglegiable?
#24  enabled     0.0035	    0.0004      0.00290     0.00040     useinit     useinit     same        same        hard    hard    tl208_swire 	                END
# TODO what numbers here - disabled?
24  disabled    0.0         0.0         same        same        useinit     useinit     same        same        hard    hard    tl208_swire 	                END
BREAK
# TODO what numbers here
# not included in background model?
25  disabled     0.380	    0.038       0.380       0.0380      useinit     useinit     same        same        hard    hard    bi214_sscin	    pb214_sscin	    END
#14  disabled    0.380	    0.038       0.380       0.0380      useinit     useinit     same        same        hard    hard    bi214_sscin	                    END
#14  disabled    0.380	    0.380       0.380       0.0380      useinit     useinit     same        same        hard    hard    pb214_sscin	                    END
BREAK
# TODO what numbers here - disabled?
26  enabled     1.3e-3      0.2e-3      same        same        useinit     useinit     same        same        soft    soft    bi210_sfoil 	                END
# TODO duplicate measurement?
# TODO what numbers here - disabled?
27  enabled     11.7e-3     1.2e-3      same        same        useinit     useinit     same        same        soft    soft    bi210_swire 	                END
#27  enabled     9.0         0.0         same        same        useinit     useinit     same        same        hard    hard    bi210_swire 	                END
# TODO what numbers here - disabled?
28  enabled     34.4        4.0         same        same        useinit     useinit     same        same        soft    soft    bi210_sscin 	                END
BREAK
# these are disabled due to low numbers of events / zero events
# if backgrounds are not listed in the header files, it does not matter if
# they are in this parameters file as they will be ignored
29  disabled    557.0	    38.0        same        same        557.0	    38.0        0.0         0.0         hard    hard    bi214_air	                    END
30  disabled    566.0	    56.6        same        same        566.0	    56.6        0.0         0.0         hard    hard    pb214_air	                    END
31  disabled    11.5  	    1.150       same        same        11.5  	    1.150       0.0         0.0         hard    hard    tl208_air   	                END
BREAK
# TODO: re-enable
32  enabled     186.0	    62.0        same        same        186.0	    62.0        same        same        soft    soft    bi214_pmt	                    END
33  enabled     36.0	    5.0         same        same        36.0	    5.0         same        same        soft    soft    ac228_pmt	                    END
34  enabled     36.0	    5.0         same        same        36.0	    5.0         same        same        soft    soft    tl208_pmt	                    END
35  enabled     9295.0	    2691        same        same        9295.0	    2691        same        same        soft    soft    bi214_feShield	                END
36  enabled     1090.0      529.0       same        same        1090.0      529.0       same        same        soft    soft    ac228_feShield	                END
37  enabled     1090.0      529.0       same        same        1090.0      529.0       same        same        soft    soft    tl208_feShield	                END
BREAK
# TODO what numbers here - disabled?
38  enabled     21.5        0.9         same        same        useinit     useinit     same        same        soft    soft    k40_sscin	                    END
BREAK
# TODO what numbers here - disabled?
39  enabled     3.1         0.5         same        same        3.1         0.5         same        same        soft    soft    pa234m_sscin                    END
40  enabled     60.0    	22.0        same        same        60.0    	22.0        same        same        soft    soft    co60_cuTower	                END
41  enabled     1124.0	    93.0        same        same        1124.0	    93.0        same        same        soft    soft    k40_pmt	                        END
# TODO: should be scint not scintIn ?
# TODO what numbers here - disabled?
42  enabled     16.0        4.0         same        same        16.0        4.0         same        same        soft    soft    k40_scintIN	                    END
BREAK
BREAK
BREAK
#1	ac228_int_rot  	        1.363   	0.056       same        same        useinit     useinit     same        same        soft    soft    enabled
#30  bi212_int_rot  	    1.363	    0.056       same        same        useinit     useinit     same        same        soft    soft    enabled
#31  tl208_int_rot 	        0.490	    0.020       same        same        useinit     useinit     same        same        soft    soft    enabled
#2   bi207_int_rot  	    6.666	    0.606       same        same        useinit     useinit     same        same        soft    soft    enabled
#3   bi214_int_rot  	    0.172	    0.020       same        same        useinit     useinit     same        same        soft    soft    enabled
#32  pb214_int_rot 	        0.172	    0.020       same        same        useinit     useinit     same        same        soft    soft    enabled
#4   eu152_int_rot  	    4.545	    0.303       same        same        useinit     useinit     same        same        soft    soft    enabled
#5   eu154_int_rot  	    0.758	    0.101       same        same        useinit     useinit     same        same        soft    soft    enabled
#6   k40_int_rot	        9.696	    0.101       same        same        useinit     useinit     same        same        soft    soft    enabled
#7   pa234m_int_rot 	    2.778	    0.035       same        same        useinit     useinit     same        same        soft    soft    enabled
BREAK
#18  enabled     0.18e-03    0.02e-03    same        same        soft        useinit     useinit     soft        same        same        bi214_mylar pb214_mylar END
#18  pb214_mylar	            3.24e-03    0.74e-03    same        same        useinit     useinit     same        same        soft    soft    enabled
BREAK
#18  bi214_mylar	            0.18e-03    0.02e-03    same        same        useinit     useinit     same        same        soft    soft    enabled
#18  pb214_mylar	            3.24e-03    0.74e-03    same        same        useinit     useinit     same        same        soft    soft    enabled
BREAK
#9   mo100_99_rot_2n2b_m14	4.86e-03    0.04e-03    same        same        useinit     useinit     same        same        soft    soft    enabled
#36  mo100_99_rot_bi214	    0.135e-03   0.002e-03   same        same        useinit     useinit     same        same        hard    hard    disabled
#37  mo100_99_rot_pa234m	    0.73e-03    0.03e-03    same        same        useinit     useinit     same        same        hard    hard    disabled
#38  mo100_99_rot_k40	    5.27e-03    0.07e-03    same        same        useinit     useinit     same        same        hard    hard    disabled
BREAK
#10  zr96_rot_2n2b_m4	    0.051e-03	0.003e-03   same        same        useinit     useinit     same        same        soft    soft    enabled
#39  zr96_rot_bi214	        0.19e-03    0.02e-03    same        same        useinit     useinit     same        same        hard    hard    disabled
#40  zr96_rot_pa234m	        0.49e-03    0.01e-03    same        same        useinit     useinit     same        same        hard    hard    disabled
#41  zr96_rot_k40	        19.7e-03	0.1e-03     same        same        useinit     useinit     same        same        hard    hard    disabled
BREAK
#8   ca48_63_rot_2n2b_m4	    0.031e-03	0.003e-03   same        same        useinit     useinit     same        same        soft    soft    enabled
#33  ca48_63_rot_bi214	    0.08e-03    0.01e-03    same        same        useinit     useinit     same        same        hard    hard    disabled
#34  ca48_63_rot_pa234m	    0.3e-03     0.1e-03     same        same        useinit     useinit     same        same        hard    hard    disabled
#35  ca48_63_rot_y90	        29.6e-03    0.1e-03     same        same        useinit     useinit     same        same        hard    hard    disabled
BREAK
#11  bi214_sfoil_rot         0.0091	    0.0048      0.00341     0.00096     useinit     useinit     same        same        soft    soft    enabled
#11  pb214_sfoil	            0.0091	    0.0048      0.00341	    0.00096     useinit     useinit     same        same        soft    soft    enabled
BREAK
#12  bi214_swire 	        0.556	    0.058       0.085	    0.012       useinit     useinit     same        same        soft    soft    enabled
#12  pb214_swire 	        0.556	    0.058       0.085       0.012       useinit     useinit     same        same        soft    soft    enabled
BREAK
#13  tl208_swire 	        0.0035	    0.0004      0.00290     0.00040     useinit     useinit     same        same        soft    soft    enabled
BREAK
#14  bi214_sscin	            0.380	    0.038       0.380       0.0380      useinit     useinit     same        same        hard    hard    disabled
#14  pb214_sscin	            0.380	    0.380       0.380       0.0380      useinit     useinit     same        same        hard    hard    disabled
BREAK
#15  bi210_sfoil 	        4.244	    0.023       same        same        useinit     useinit     same        same        hard    hard    disabled
#16  bi210_swire 	        9.0         0.0         same        same        useinit     useinit     same        same        hard    hard    disabled
#17  bi210_sscin 	        30.400	    3.040       same        same        useinit     useinit     same        same        hard    hard    disabled
BREAK
#19  bi214_air	            557.0	    38. 0       0.0         0.0         useinit     useinit     same        same        free    hard    disabled
#19  pb214_air	            566.000	    56.600      0.0         0.0         useinit     useinit     same        same        free    hard    disabled
#20  tl208_air   	        11.500  	1.150       0.0         0.0         useinit     useinit     same        same        free    hard    disabled
BREAK
#21  bi214_pmt	            324.0	    1.0         same        same        useinit     useinit     same        same        soft    soft    disabled
#22  ac228_pmt	            27.0	    0.6         same        same        useinit     useinit     same        same        soft    soft    disabled
#22  tl208_pmt	            27.0	    0.6         same        same        useinit     useinit     same        same        soft    soft    disabled
#23  bi214_feShield	        7360.0	    200.0       same        same        useinit     useinit     same        same        soft    soft    disabled
#24  ac228_feShield	        484.0	    24.0        same        same        useinit     useinit     same        same        soft    soft    disabled
#24  tl208_feShield	        484.0	    24.0        same        same        useinit     useinit     same        same        soft    soft    disabled
#25  k40_sscin	            0.0	        1.0         same        same        useinit     useinit     same        same        free    free    disabled
#26  pa234m_sscin            0.0         1.0         same        same        useinit     useinit     same        same        free    free    disabled
#27  co60_cuTower	        18.4	    0.8         same        same        useinit     useinit     same        same        soft    soft    disabled
#28  k40_pmt	                1078.0	    32.0        same        same        useinit     useinit     same        same        soft    soft    disabled
#29  k40_scintIN	            21.5        0.9         same        same        useinit     useinit     same        same        soft    soft    disabled
BREAK
# TODO: there are some missing externals, see Page40 of position paper, Table 4.4
# NOTE: these are reported from NEMO3 background paper, so not necessarily totally relevant here
