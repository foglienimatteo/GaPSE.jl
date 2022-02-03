# -*- encoding: utf-8 -*-
#
# This file is part of GaPSE
# Copyright (C) 2022 Matteo Foglieni
#
# GaPSE is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GaPSE is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GaPSE. If not, see <http://www.gnu.org/licenses/>.
#


const ZS = [0.04247970795382, 0.04981568656961, 0.05720310657579, 0.06464232613797,
     0.07213370591695, 0.07967760908732, 0.08727440135637, 0.09492445098303,
     0.102628128797, 0.1103858082182, 0.1181978652761, 0.1260646786292,
     0.1339866295853, 0.1419641021209, 0.1499974829017, 0.1580871613026,
     0.1662335294283, 0.1744369821336, 0.1826979170446, 0.1910167345789,
     0.1993938379675, 0.2078296332749]

const CONF_TIME = [9600.77243363, 9579.242853425, 9557.637777233, 9535.958012965999,
     9514.204385993999, 9492.377738984, 9470.478931732, 9448.508840974,
     9426.468360169, 9404.358399338, 9382.17988477, 9359.933758846999,
     9337.62097975, 9315.242521228, 9292.799372296999, 9270.292536953,
     9247.723033885, 9225.091896146, 9202.400170845, 9179.648918789999,
     9156.839214152, 9133.972144094]

const COM_DIST = [126.12245883194998, 147.65203903268, 169.2571152265,
     190.93687949574, 212.69050646823, 234.51715347585,
     256.41596072631, 278.38605148773996, 300.42653228769996,
     322.53649312415996, 344.71500769082996, 366.96113361466996,
     389.27391270725997, 411.65237122806997, 434.09552016158995,
     456.60235550608, 479.17185857534, 501.80299631216997,
     524.49472161387, 547.2459736699899, 570.05567831037,
     592.9227483653599]

const ANG_DIST = [120.98313076947998, 140.64567801913998, 160.09895749812998,
     179.34368642694, 198.38058004743996, 217.21035196252, 235.83371447572998,
     254.25137893105997, 272.46405605076, 290.4724562733, 308.27729008916,
     325.87926837507996, 343.27910272593, 360.47750578455, 377.47519156844004,
     394.2728757933, 410.87127619306995, 427.27111283612, 443.47310843717,
     459.47798866445, 475.28648244217, 490.89932224774]

const LUM_DIST = [131.48010404956, 155.00742673052, 178.93914802756998, 203.27948353187998,
     228.03266091311997, 253.20291955482998, 278.79451019686, 304.81169458655995,
     331.25874513732, 358.13994459753997, 385.45958572849, 413.22197099323,
     441.43141225638, 470.09223049536996, 499.20875552470994, 528.78532573216,
     558.82628782904, 589.33599661447, 620.31881475361, 651.77911257193, 683.72126786383,
     716.1496657181999]

const GROWTH_FACTOR_D = [0.978380534605, 0.9746812066927, 0.9709667334511, 0.9672373848039,
     0.9634934355829, 0.9597351654362, 0.9559628587335, 0.9521768044659,
     0.9483772961421, 0.94456463168, 0.9407391132946, 0.9369010473811,
     0.9330507443948, 0.929188518727, 0.9253146885768, 0.9214295758196,
     0.917533505872, 0.9136268075535, 0.9097098129447, 0.9057828572426,
     0.9018462786131, 0.8979004180403]

const GROWTH_FACTOR_F = [0.5380746269755, 0.5423655052041, 0.5466585304504, 0.5509530414522,
     0.5552483721526, 0.5595438521047, 0.563838806887, 0.568132558527,
     0.5724244259342, 0.5767137253412, 0.5809997707523, 0.5852818743989,
     0.5895593472022, 0.5938314992407, 0.5980976402232, 0.6023570799664,
     0.6066091288758, 0.61085309843, 0.6150883016669, 0.6193140536721,
     0.6235296720668, 0.6277344774972]

const COM_H = [0.0003262902207941897, 0.00032513406832466497, 0.00032399793165289537,
     0.00032288173953739976, 0.0003217854203739545, 0.0003207089021754946,
     0.0003196521125537181, 0.0003186149786999385, 0.00031759742736752936,
     0.0003165993848549706, 0.0003156207769894872, 0.00031466152911081526,
     0.000313721566057328, 0.0003128008121510082, 0.00031189919118453015,
     0.00031101662640907766, 0.00031015304052184887, 0.00030930835565644126,
     0.0003084824933720576, 0.00030767537464485695, 0.0003068869198595986,
     0.0003061170488029462]


const COM_H_P = [5.4262324798057274e-8, 5.3141774159118e-8, 5.203416611112423e-8,
     5.093945428912454e-8, 4.985739959051383e-8, 4.8787817979581974e-8,
     4.773051394735852e-8, 4.668529828121139e-8, 4.5651983333984464e-8,
     4.4630384164018666e-8, 4.362031821820588e-8, 4.2621605403525496e-8,
     4.1634067956486114e-8, 4.0657530528081654e-8, 3.96918199795169e-8,
     3.8736765570664605e-8, 3.7792198595880137e-8, 3.685795324651405e-8,
     3.593386334008363e-8, 3.501977445649015e-8, 3.411549830517371e-8,
     3.322098316855463e-8]