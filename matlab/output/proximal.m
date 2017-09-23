
% ************************************************************************************

%   OUTPUT FROM MODELS WITH Y VARIABLES BY SPECIES, PROXIMAL SEGMENT ONLY

% ************************************************************************************




% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%       log_intestinal_diameter ~ taxon + log_mass

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% bGLS_approx_conf_int_for_b =
%
%    -1.5155   -2.0354   -1.4689   -0.9023
%     0.0667   -0.1421    0.0658    0.2737
%     0.3035    0.1500    0.2926    0.4352
%
% Parameters with (lb,mean,ub) from simulation
%
% b0 (intercept) = -1.5155 (-2.7234, -1.4357, -0.26421)
%             b1 = 0.066684 (-1.038, 0.08375, 1.1851)
%             b2 = 0.30348 (0.075105, 0.30722, 0.5402)
%
%         sigma2 = 0.025478 (0.024458, 0.06611, 0.13027)
%
%              d = 6.7485e-67 (2.1238e-73, 0.045139, 0.41371)
%
% Independent variable means, variances, and correlations
% aX1 = 0.56286
% aX2 = 3.5064
% s2X1 = -0.11655
% s2X2 = -0.022316
% r(1,2) = -0.15702
%
% LnLikelihood = -0.66777
% -2LL = 1.3355
% AIC(par=6) = 13.3355
%
%
% ans =
%
%    -0.6678






% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%       villus_height ~ taxon + log_mass

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



% bGLS_approx_conf_int_for_b =
%
%    -0.0741   -0.3856   -0.0386    0.3084
%     0.0822   -0.0402    0.0872    0.2145
%     0.1663    0.0718    0.1596    0.2475
%
% Parameters with (lb,mean,ub) from simulation
%
% b0 (intercept) = -0.074079 (-0.76569, -0.0050897, 0.52831)
%             b1 = 0.082169 (-0.4835, 0.078652, 0.64603)
%             b2 = 0.16634 (0.056159, 0.16856, 0.29521)
%
%         sigma2 = 0.0087501 (0.0050542, 0.01635, 0.032039)
%
%              d = 4.8708e-57 (5.7673e-73, 0.057463, 0.47884)
%
% Independent variable means, variances, and correlations
% aX1 = 0.56233
% aX2 = 3.5032
% s2X1 = 0.11645
% s2X2 = 0.022879
% r(1,2) = -0.16112
%
% LnLikelihood = 8.1156
% -2LL = -16.2312
% AIC(par=6) = -4.2312
%
%
% ans =
%
%     8.1156





% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%       villus_width ~ taxon

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


% bGLS_approx_conf_int_for_b =
%
%     0.0995    0.0815    0.0990    0.1165
%    -0.0304   -0.0536   -0.0302   -0.0068
%
% Parameters with (lb,mean,ub) from simulation
%
% b0 (intercept) = 0.099467 (0.047321, 0.098759, 0.14576)
%             b1 = -0.030353 (-0.10964, -0.028849, 0.057096)
%
%         sigma2 = 0.00028166 (0.00014608, 0.00038704, 0.0007367)
%
%              d = 0.053721 (7.2563e-70, 0.087587, 0.51546)
%
% Independent variable means, variances, and correlations
% aX1 = 0.56235
% s2X1 = 0.11866
%
% LnLikelihood = 60.6175
% -2LL = -121.235
% AIC(par=3) = -115.235
%
%
% ans =
%
%    60.6175








% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%       crypt_width ~ taxon

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



% bGLS_approx_conf_int_for_b =
%
%     0.0404    0.0309    0.0409    0.0509
%    -0.0108   -0.0238   -0.0108    0.0021
%
% Parameters with (lb,mean,ub) from simulation
%
% b0 (intercept) = 0.040445 (0.023254, 0.04091, 0.058467)
%             b1 = -0.0108 (-0.040055, -0.011532, 0.017313)
%
%         sigma2 = 2.797e-05 (7.8912e-06, 2.7944e-05, 5.5915e-05)
%
%              d = 0.61621 (1.5718e-66, 0.34795, 1.1601)
%
% Independent variable means, variances, and correlations
% aX1 = 0.56235
% s2X1 = 0.11866
%
% LnLikelihood = 80.6525
% -2LL = -161.305
% AIC(par=3) = -155.305
%
%
% ans =
%
%    80.6525






% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%       sef ~ taxon

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


% bGLS_approx_conf_int_for_b =
%
%    12.2213    8.9993   12.2031   15.4070
%     4.5537    0.2676    4.5290    8.7905
%
% Parameters with (lb,mean,ub) from simulation
%
% b0 (intercept) = 12.2213 (3.4387, 12.1924, 21.003)
%             b1 = 4.5537 (-10.0795, 4.5775, 19.6817)
%
%         sigma2 = 9.6227 (4.4399, 11.8119, 22.1731)
%
%              d = 0.09318 (2.8847e-69, 0.097494, 0.51169)
%
% Independent variable means, variances, and correlations
% aX1 = 0.56235
% s2X1 = 0.11866
%
% LnLikelihood = -33.7159
% -2LL = 67.4319
% AIC(par=3) = 73.4319
%
%
% ans =
%
%   -33.7159






% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%       enterocyte_diameter ~ taxon

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


% bGLS_approx_conf_int_for_b =
%
%     0.0071    0.0060    0.0071    0.0083
%    -0.0007   -0.0022   -0.0007    0.0009
%
% Parameters with (lb,mean,ub) from simulation
%
% b0 (intercept) = 0.007143 (0.0035258, 0.0071221, 0.010468)
%             b1 = -0.00066108 (-0.0065495, -0.00061737, 0.0054594)
%
%         sigma2 = 1.5917e-06 (7.3699e-07, 1.9193e-06, 3.6611e-06)
%
%              d = 0.11354 (7.1852e-69, 0.10522, 0.54802)
%
% Independent variable means, variances, and correlations
% aX1 = 0.56235
% s2X1 = 0.11866
%
% LnLikelihood = 107.5112
% -2LL = -215.0223
% AIC(par=3) = -209.0223
%
%
% ans =
%
%   107.5112





% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% log_enterocyte_density ~ taxon

% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


% bGLS_approx_conf_int_for_b =
%
%    17.0791   16.7000   17.1185   17.5371
%     0.5845    0.0051    0.5925    1.1799
%
% Parameters with (lb,mean,ub) from simulation
%
% b0 (intercept) = 17.0791 (15.6435, 17.091, 18.6166)
%             b1 = 0.58453 (-1.9234, 0.57243, 3.0527)
%
%         sigma2 = 0.24882 (0.14798, 0.37223, 0.68682)
%
%              d = 0.00097176 (1.4481e-71, 0.054235, 0.398)
%
% Independent variable means, variances, and correlations
% aX1 = 0.56235
% s2X1 = 0.11866
%
% LnLikelihood = -0.44839
% -2LL = 0.89677
% AIC(par=3) = 6.8968
%
%
% ans =
%
%    -0.4484