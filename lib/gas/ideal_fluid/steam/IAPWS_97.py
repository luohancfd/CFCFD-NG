# -*- coding: utf-8 -*-
"""
IAPWS-97.py
===========
Description
-----------
Created on Fri Oct 11 11:39:10 2013
@author: uqbwheat

"""
# Specific gas constant for ordinary water
R = 461.526 # J/kg/K    Wagner (2000) equation 1
Tc = 647.096 # K Critical Temperature Wagner (2000) equation 4
Pc = 22.064 # MPa Critical pressure. Wagner (2000) equation 5
rho_c = 322 # kg/m3  Critical density. Wagner (2000) equation 6
Tt = 273.16 #K Triple point. Wagner (2000) equation 7
Pt = 611.657 # Pa Triple point pressure. Wagner (2000) equation 8

properties = ["Specific Volume", "Specific Enthalpy", "Specific Internal Energy", "Specific Entropy", "Specific Isobaric Heat Capacity", "Speed of Sound"]

units = ["m3/kg", "kJ/kg", "kJ/kg", "kJ/kg/K", "kJ/kg/K", "m/s"]

import scipy as sp
from copy import *

class IAPWS:
    """
    A class to implement the International Association for the Properties of Water and Steam (IAPWS) Industrial Formulation 1997 for the Thermodynamic Properties of water and Steam.

    Reference:
        Wagner et al., 2000, The IAPWS Industrial Formulation 1997 for the Thermodynamic Properties of Water and Steam, ASME J. Eng. Gas Turbines and Power, Vol 122, 150-182.

    Description
    -----------
    Wagner (2000) lists the numerical information for the equations of IAPWS-IF97. This class implements those equations.

    Credit
    -----
    This code is built on Matlab code developed by David Buttsworth.

    Range of Validity
    -----------------

    This class is valid for:
        Temperatures from 0 C to 800 C
        Pressures less than 100 MPa.

    For high temperatures, the class is valid for:
        Temperatures from 800 C to 2000 C
        Pressures less than 10 MPa

    Accuracy
    --------
    Refer to Wagner (2000) for accuracy details.
    """

    def __init__(self):


        self.tableA1 = [0.34805185628969e3, -0.11671859879975e1, 0.10192970039326e-2, 0.57254459862746e3, 0.13918839778870e2]

        self.tableA2 = []
        I = [0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,4,4,5,8,8,21,23,29,30,31,32]
        J = [-2,-1,0,1,2,3,4,5,-9,-7,-1,0,1,3,-3,0,1,6,17,-4,0,6,-5,-2,10,-8,-11,-6,-29,-31,-38,-39,-40,-41]
        n = [0.14632971213167e+00, -0.84548187169114e+00, -0.37563603672040e+01, 0.33855169168385e+01, -0.95791963387872e+00, 0.15772038513228e+00, -0.16616417199501e-01, 0.81214629983568e-03, 0.28319080123804e-03, -0.60706301565874e-03, -0.18990068218419e-01, -0.32529748770505e-01, -0.21841717175414e-01, -0.52838357969930e-04, -0.47184321073267e-03, -0.30001780793026e-03, 0.47661393906987e-04, -0.44141845330846e-05, -0.72694996297594e-15, -0.31679644845054e-04, -0.28270797985312e-05, -0.85205128120103e-09, -0.22425281908000e-05, -0.65171222895601e-06, -0.14341729937924e-12, -0.40516996860117e-06, -0.12734301741641e-08, -0.17424871230634e-09, -0.68762131295531e-18, 0.14478307828521e-19, 0.26335781662795e-22, -0.11947622640071e-22, 0.18228094581404e-23, -0.93537087292458e-25]
        for i in range(len(I)):
            self.tableA2.append([I[i], J[i], n[i]])

        Condition_1 = [0.100215168e-2, 0.115331273e3, 0.112324818e3, 0.392294792, 0.417301218e1, 0.150773921e4]
        Condition_2 = [0.971180894e-3, 0.184142828e3, 0.106448356e3, 0.368563852, 0.401008987e1, 0.1563469054e4]
        Condition_3 = [0.120241800e-2, 0.975542239e3, 0.971934985e3, 0.258041912e1, 0.465580682e1, 0.124071337e4]
        self.tableA3 = [Condition_1, Condition_2, Condition_3]

        self.tableA4 = []
        J = [0, 1, -5, -4, -3, -2, -1, 2, 3]
        n = [-0.96927686500217e1, 0.10086655968018e2, -0.56087911283020e-2, 0.71452738081455e-1, -0.40710498223928, 0.14240819171444e1, -0.43839511319450e1, -0.28408632460772, 0.21268463753307e-1]
        for i in range(len(J)):
            self.tableA4.append([J[i], n[i]])

        self.tableA4a = copy(self.tableA4)
        self.tableA4a[0][1] = -0.96937268393049e1
        self.tableA4a[1][1] = 0.10087275970006e2

        self.tableA5 = []
        I = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24]
        J = [0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58]
        n = [-0.0017731742473213, -0.017834862292358, -0.045996013696365, -0.057581259083432, -0.05032527872793, -3.3032641670203e-05, -0.00018948987516315, -0.0039392777243355, -0.043797295650573, -2.6674547914087e-05, 2.0481737692309e-08, 4.3870667284435e-07, -3.227767723857e-05, -0.0015033924542148, -0.040668253562649, -7.8847309559367e-10, 1.2790717852285e-08, 4.8225372718507e-07, 2.2922076337661e-06, -1.6714766451061e-11, -0.0021171472321355, -23.895741934104, -5.905956432427e-18, -1.2621808899101e-06, -0.038946842435739, 1.1256211360459e-11, -8.2311340897998, 1.9809712802088e-08, 1.0406965210174e-19, -1.0234747095929e-13, -1.0018179379511e-09, -8.0882908646985e-11, 0.10693031879409, -0.33662250574171, 8.9185845355421e-25, 3.0629316876232e-13, -4.2002467698208e-06, -5.9056029685639e-26, 3.7826947613457e-06, -1.2768608934681e-15, 7.3087610595061e-29, 5.5414715350778e-17, -9.436970724121e-07]
        for i in range(len(I)):
            self.tableA5.append([I[i], J[i], n[i]])



        Condition_1 = [0.394913866e2, 0.254991145e4, 0.241169160e1, 0.852238967e1, 0.191300162e1, 0.427924172e3]
        Condition_2 = [0.923015898e2, 0.333568375e4, 0.301262819e4, 0.101749996e2, 0.208141274e1, 0.644289068e3]
        Condition_3 = [0.542946619e-2, 0.263149474e4, 0.246861076e4, 0.517540298e1, 0.103505092e2, 0.480386523e3]
        self.tableA6 = [Condition_1, Condition_2, Condition_3]

        self.tableA7=[]
        I = [1,1,1,1,2,2,2,3,3,4,4,5,5]
        J = [0,2,5,11,1,7,16,4,16,7,10,9,10]
        n = [-0.73362260186506e-02, -0.88223831943146e-01, -0.72334555213245e-01, -0.40813178534455e-02, 0.20097803380207e-02, -0.53045921898642e-01, -0.76190409086970e-02, -0.63498037657313e-02, -0.86043093028588e-01,  0.75321581522770e-02, -0.79238375446139e-02, -0.22888160778447e-03, -0.26456501482810e-02]
        for i in range(len(I)):
            self.tableA7.append([I[i], J[i], n[i]])

        Condition_1 = [0.192516540, 0.276881115e4, 0.257629461e4, 0.656660377e1, 0.276349265e1, 0.498408101e3]
        Condition_2 = [0.186212297, 0.274015123e4, 0.255393894e4, 0.650218759e1, 0.298166443e1, 0.489363295e3]
        Condition_3 = [0.121685206, 0.272134539e4, 0.253881758e4, 0.629170440e1, 0.362795578e1, 0.481941819e3]
        self.tableA8 = [Condition_1, Condition_2, Condition_3]

        self.tableA9 = []
        I = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11]
        J = [0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26]
        n = [1.0658070028513, -15.732845290239, 20.944396974307, -7.6867707878716, 2.6185947787954, -2.808078114862, 1.2053369696517, -0.0084566812812502, -1.2654315477714, -1.1524407806681, 0.88521043984318, -0.64207765181607, 0.38493460186671, -0.85214708824206, 4.8972281541877, -3.0502617256965, 0.039420536879154, 0.12558408424308, -0.2799932969871, 1.389979956946, -2.018991502357, -0.0082147637173963, -0.47596035734923, 0.0439840744735, -0.44476435428739, 0.90572070719733, 0.70522450087967, 0.10770512626332, -0.32913623258954, -0.50871062041158, -0.022175400873096, 0.094260751665092, 0.16436278447961, -0.013503372241348, -0.014834345352472, 0.00057922953628084, 0.0032308904703711, 8.0964802996215e-05, -0.00016557679795037, -4.4923899061815e-05]
        for i in range(len(I)):
            self.tableA9.append([I[i], J[i], n[i]])

        Condition_1 = [0.255837018e2, 0.186343019e4, 0.181226279e4, 0.405427273e1, 0.138935717e2, 0.502005554e3]
        Condition_2 = [0.222930643e2, 0.237512401e4, 0.226365868e4, 0.485438792e1, 0.446579342e2, 0.383444594e3]
        Condition_3 = [0.783095639e2, 0.225868845e4, 0.210206932e4, 0.446971906e1, 0.634165359e1, 0.760696041e3]
        self.tableA10 = [Condition_1, Condition_2, Condition_3]


        self.tableA11=[0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849, 0.65017534844798e3]

        self.tableA12 = {300: 0.353658941e-2, 500: 0.263889776e1, 600: 0.123443146e2}

        self.tableA13=[]
        J = [0, 1, -3, -2, -1, 2]
        n = [-0.13179983674201e2,  0.68540841634434e1, -0.24805148933466e-1, 0.36901534980333, -0.31161318213925e1, -0.32961626538917]
        for i in range(len(J)):
            self.tableA13.append([J[i], n[i]])

        self.tableA14=[]
        I = [1,1,1,2,3]
        J = [0,1,3,9,3]
        n = [-0.12563183589592e-3, 0.21774678714571e-2, -0.459428208999e-2, -0.39724828359569e-5, 0.12919228289784e-6]
        for i in range(len(I)):
            self.tableA14.append([I[i], J[i], n[i]])

        Condition_1 = [0.138455354e1, 0.521976332e4, 0.452748654e4, 0.965408431e1, 0.261610228e1, 0.917071933e3]
        Condition_2 = [0.865156616e-1, 0.520609634e4, 0.451397105e4, 0.836546724e1, 0.264453866e1, 0.919708859e3]
        Condition_3 = [0.115743146, 0.658380291e4, 0.565785774e4, 0.915671044e1, 0.285306750e1, 0.105435806e4]
        self.tableA15 = [Condition_1, Condition_2, Condition_3]

        self.tableA16 = []
        I = [0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,3,3,4,5,6]
        J = [0,1,2,6,22,32,0,1,2,3,4,10,32,10,32,10,32,32,32,32]
        n = [-0.23872488924521e3, 0.40421188637945e3, 0.11349746881718e3, -0.58457616048039e1, -0.15285482413140e-3, -0.10866707695377e-5, -0.13391744872602e2, 0.43211039183559e2, -0.54010067170506e2, 0.30535892203916e2, -0.65964749423638e1, 0.93965400878363e-2, 0.11573647505340e-6, -0.25858641282073e-4, -0.40644363084799e-8, 0.66456186191635e-7, 0.80670734103027e-10, -0.93477771213947e-12, 0.58265442020601e-14, -0.15020185953503e-16]
        for i in range(len(I)):
            self.tableA16.append([I[i], J[i], n[i]])

        self.tableA17 = []
        Pressure = [3.0, 80.0 ,80.0]
        Enthalpy = [500.0, 500.0, 1500.0]
        Temperature = [0.391798509e3, 0.37810862e3, 0.611041229e3]
        for i in range(len(Pressure)):
            self.tableA17.append([Pressure[i], Enthalpy[i], Temperature[i]])

        self.tableA18 = []
        I = [0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,4]
        J = [0,1,2,3,11,31,0,1,2,3,12,31,0,1,2,9,31,10,32,32]
        n = [0.17478268058307e3, 0.34806930892873e2, 0.65292584978455e1, 0.33039981775489, -0.19281382923196e-6, -0.24909197244573e-22, -0.26107636489332, 0.22592965981586, -0.64256463395226e-1, 0.78876289270526e-2, 0.35672110607366e-9, 0.17332496994895e-23, 0.56608900654837e-3, -0.32635483139717e-3, 0.44778286690632e-4, -0.51322156908507e-9, -0.42522657042207e-25, 0.26400441360689e-12, 0.78124600459723e-28, -0.30732199903668e-30]
        for i in range(len(I)):
            self.tableA18.append([I[i], J[i], n[i]])

        self.tableA19 = []
        Pressure = [3.0, 80.0 ,80.0]
        Entropy = [0.5, 0.5, 3.0]
        Temperature = [0.307842258e3, 0.309979785e3, 0.565899909e3]
        for i in range(len(Pressure)):
            self.tableA19.append([Pressure[i], Entropy[i], Temperature[i]])

        self.tableA20 = [0.9058478514723e3, -0.67955786399241, 0.12809002730136e-36, 0.26526571908428e4, 0.45257578905948e1]

        self.tableA21 = []
        I = [0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,4,4,4,5,5,5,6,6,7]
        J = [0,1,2,3,7,20,0,1,2,3,7,9,11,18,44,0,2,7,36,38,40,42,44,24,44,12,32,44,32,36,42,34,44,28]
        n = [0.10898952318288e4, 0.84951654495535e3, -0.10781748091826e3, 0.33153654801263e2, -0.74232016790248e1, 0.11765048724356e2, 0.18445749355790e1, -0.41792700549624e1, 0.62478196935812e1, -0.17344563108114e2, -0.20058176862096e3, 0.27196065473796e3, -0.45511318285818e3, 0.30919688604755e4, 0.25226640357872e6, -0.61707422868339e-2, -0.31078046629583, 0.11670873077107e2, 0.12812798404046e9, -0.98554909623276e9, 0.28224546973002e10, -0.35948971410703e10, 0.17227349913198e10, -0.1355133440775e5, 0.12848734664650e8, 0.13865724283226e1, 0.23598832556514e6, -0.13105236545054e8, 0.73999835474766e4, -0.55196697030060e6, 0.37154085996233e7, 0.19127729239660e5, -0.41535164835634e6, -0.62459855192507e2]
        for i in range(len(I)):
            self.tableA21.append([I[i], J[i],n[i]])

        self.tableA22 = []
        I = [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,4,4,5,5,5,6,7,7,9,9]
        J = [0,1,2,12,18,24,28,40,0,2,6,12,18,24,28,40,2,8,18,40,1,2,12,24,2,12,18,24, 28,40,18,24,40,28,2,28,1,40]
        n = [0.14895041079516e4, 0.74307798314034e3, -0.97708318797837e2, 0.24742464705674e1, -0.63281320016026, 0.11385952129658e1, -0.4781186348625, 0.85208123431544e-2, 0.93747147377932, 0.33593118604916e1, 0.33809355601454e1, 0.16844539671904, 0.73875745236695, -0.47128737436186, 0.15020273139707, -0.21764114219750e-2, -0.21810755324761e-1, -0.10829784477, -0.46333324635812e-1, 0.71280351951e-4, 0.11032831789999e-3, 0.18955248387902e-3, 0.30891541160537e-2, 0.13555504554949e-2, 0.28640237477456e-6, -0.10779857357512e-4, -0.76462712454814e-4, 0.14052392818316e-4, -0.31083814331434e-4, -0.10302738212103e-5, 0.28217281635040e-6, 0.12704902271945e-5, 0.73803353468292e-7, -0.11030139238909e-7, -0.81456365207833e-13, -0.25180545682962e-10, -0.17565233969407e-17, 0.86934156344163e-14]
        for i in range(len(I)):
            self.tableA22.append([I[i], J[i],n[i]])

        self.tableA23 = []
        I = [-7,-7,-6,-6,-5,-5,-2,-2,-1,-1,0,0,1,1,2,6,6,6,6,6,6,6,6]
        J = [0,4,0,2,0,2,0,1,0,2,0,1,4,8,4,0,1,4,10,12,16,20,22]
        n = [-0.3268398555242e13, 0.73263350902181e13, 0.35825089945447e12, -0.58340131851590e12, -0.10783068217470e11, 0.20825544563171e11, 0.61074783564516e6, 0.85977722535580e6, -0.25745723604170e5, 0.31081088422714e5, 0.12082315865936e4, 0.48219755109255e3, 0.37966001272486e1, -0.10842984880077e2, -0.45364172676660e-1, 0.14559115658698e-12, 0.11261597407230e-11, -0.17804982240686e-10, 0.12324579690832e-6, -0.11606921130984e-5, 0.27846367088554e-4, -0.59270038474176e-3, 0.12918582991878e-2]
        for i in range(len(I)):
            self.tableA23.append([I[i], J[i],n[i]])

        self.tableA24 = {"43":[[0.001,3000.0,0.534433241e3],[3.0, 3000.0,0.575373370e3],[3.0, 4000.0,0.101077577e4]],\
        "44":[[5.0,3500.0,0.801299102e3],[5.0,4000.0, 0.101531583e4],[25.0,3500.0, 0.875279054e3]],\
        "45":[[40.0,2700.0,0.743056411e3],[60.0,2700.0,0.791137067e3],[60.0,3200.0, 0.882756860e3]]}

        self.tableA28 = {"49":[[0.1,7.5,0.399517097e3],[0.1, 8.0,0.514127081e3],[2.5,8.0,0.103984917e4]],\
        "50":[[8.0,6.0,0.600484040e3],[8.0,7.5,0.106495556e4],[90.0, 6.0, 0.103801126e4]],\
        "52":[[20.0,5.75,0.697992849e3],[80.0,5.25,0.854011484e3],[80.0,5.75, 0.949017998e3]]}


        self.tableA29 = {0.1:0.372755919e3, 1.0:0.453035632e3, 10:0.584149488e3}

        self.region = "1"

    def findRegion(self, P=3, T = 300, printRegion = False):
        """

        """
        if P<=0 or P>100:
            print "\nPressure out of bounds. Pressure must be greater than zero and less than 100 MPa."

        if T<273.15 or T > 2273.15:
            print "\nTemperature out of bounds. Temperature must be greater than 273.15 K and less than 2273.15 K."
        if T>=1073.15 and T<=2273.15:
            if P>10:
                print "\nPressure out of bounds. When temperature is above 1073.15 K, pressure must be less than 10 MPa."
            else:
                self.region = "5"

        """Calculate Saturation temperature and pressure to determine if the point is on the saturation line (Region 4)"""
        if T>=273.15 and T<=623.15:                              # Critical Point
            self.region = "1"                                    # Default is region 1
            if P>611.213e-6:                                # Possible Saturation condition
                if P<22.064:
                    Ps, Ts = self.calcSatLine(P,T)   # Saturation Pressure in MPa
                    if T < Ts:
                        if P < Ps:
                            self.region = "2.1"
                        if P > Ps and P < 10:
                            self.region = "1"
                    else:
                        if P < Ps:
                            self.region = "2"
                        if P > Ps and P < 10:
                            self.region = "2.1"

        if T >= 623.15 and T <=1073.15:
            self.region = "2"
            if T<= 863.15:
                if P >= 16.5292 and P<= 100:
                    Pi, Theta = self.Eqn_B23(P,T)
                    if P >= Pi and T <= Theta:
                        self.region = "3"

        if printRegion:
            print "Temperature = ", T, "K"
            print "Pressure = ", P, "MPa"
            print "Region = ", self.region
            if self.region == "2.1":
                print "Saturation Pressure = ", Ps, "MPa"
                print "Saturation Temperature = ", Ts, "K"
            if self.region == "2" or self.region == "3":
                if T>= 623.15 and T <= 863.15 and P>=16.5292 and P<100:
                    print "Pressure at boundary between 2 and 3 is ", Pi, "MPa"


    def calcSatLine(self, P = 611.213e-6, T = 273.15):
        """
        Saturation pressure for a given saturation temperature. Region 4

        Region of validity
        ------------------
        Region 4:
            273.15 K < T < 647.096 K

        Input
        ------
        T:
            Temperature in K. Default is 273.15 K.

        Returns
        -------
        beta:
            Saturation pressure ratio defined as Ps/P* where P* = 1 MPa and Ps is the saturation pressure.

        Parameters
        ----------
        n[i]:
            Wagner (2000) Table A11 coefficients

        Method
        ------
        Refer to Wagner (2000) equations 28 and 28a.

        See Also
        --------
        calcSatTemp
        """

        n = self.tableA11
        theta = T/1.0 + n[8]/((T/1.0) - n[9])
        A = n[1] + theta * (n[0] + theta)
        B = n[4] + theta * (n[3] + n[2] * theta)
        C = n[7] + theta * (n[6] + n[5] * theta)
        beta = 2.0*C/(sp.sqrt(B*B-4.0*A*C)-B)

        Ps = sp.power(beta, 4)

        beta = sp.power((P/1.0), 0.25)

        """Equation 55a"""
        E = n[5] + beta * (n[2] + beta)
        F = n[6] + beta * (n[3] + n[0] * beta)
        G = n[7] + beta * (n[4] + n[1] *beta)
        D = -2.0*G/(F+sp.sqrt(F*F-4.0*E*G))

        """Equation 55"""
        theta = n[9] + D
        theta -= sp.sqrt(sp.power((n[9]+D),2) - 4.0*(n[8] + n[9] * D))
        theta /= 2.0

        Ts = theta

        return Ps, Ts



    def Eqn_B23(self, P=0.165291643e2, T = 0.62315e3):
        """
        Auxiliary equation for the boundary between regions 2 and 3
        ------------------------------------------------------------

        This function defines the boundary between regions 2 and 3. For a given temperature ratio, this function returns the pressure ratio

        Input
        -----
        theta - Temperature ratio. Default value is 623.15 K.

        Returns
        -------
        Pi - Pressure ratio

        Parameters
        ----------

        Pi:
            p/p* - Pressure ratio
        Theta:
            T/T* - Temperature ratio
        p* = 1 MPa
        T* = 1K

        Notes
        -----
        Describes roughly an isentropic line. The entropy values along this boundary are between s = 5.047 kJ/kg/K and s = 5.261 kJ/kg/K

        Equation 10 of Wagner (2000).

        For verification, Equations 10 and 11 must meet the following T-p point.
            T= 623.150000K
            p = 16.5291643 MPa

        Example
        -------
        >>> from IAPWS_97 import *
        >>> i=IAPWS()
        >>> pi, theta = i.Eqn_B23()
        >>> pi
        16.529164252621626
        >>> theta
        >>> 623.15
        >>>

        """
        n = self.tableA1
        if T>= 623.15 and T<=863.15:
            if P >= 16.5292 and P <= 100:
                Pi = n[0] + T * (n[1] + n[2]*T)
                Theta = n[3] + sp.power((P - n[4])/n[2], 0.5)
            else:
                print "The boundary between regions 2 and 3 can only be calculated when the\
                pressure is between 16.5292 MPa and 100 MPa."
        else:
            print "The boundary between regions 2 and 3 can only be calculated when the\
            temperature is between 623.15 K and 863.15 K."

        return Pi, Theta

    def calcProperties(self, P = 16.53, T=1386):
        """
        This function is used by the other functions to set up the inputs for gamma and its derivatives.

        Inputs
        ------
        P:
            Pressure in MPa. Default is 16.53 MPa
        T:
            Temperature in K. Default is 1386 K

        Returns
        -------
        table:
            Link to the appropriate table
        Pi:
            Non-dimensional pressure ratio. Pi = P/16.53
        Tau:
            Non-dimensional temperature ratio. Tau = 1386/T

        """
#        self.findRegion(P,T, printRegion = False)
        gamma = 0.0
        gamma0 = 0.0
        gamma1 = 0.0
        gamma_pi = 0.0
        gamma_pi_0 = 0.0
        gamma_pi_1 = 0.0
        gamma_pi2 = 0.0
        gamma_pi2_0 = 0.0
        gamma_pi2_1 = 0.0
        gamma_tau = 0.0
        gamma_tau_0 = 0.0
        gamma_tau_1 = 0.0
        gamma_tau2 = 0.0
        gamma_tau2_0 = 0.0
        gamma_tau2_1 = 0.0
        gamma_pi_tau = 0.0
        gamma_pi_tau_0 = 0.0
        gamma_pi_tau_1 = 0.0

        if self.region == "1":
            Pi = P/16.53        # Non-dimensional pressure ratio
            Tau = 1386.0/T    # Non-dimensional temperature ratio
            n1 = self.tableA2

            for i in range(len(n1)):
                gamma += n1[i][2] *sp.power((7.1 - Pi), n1[i][0]) * sp.power((Tau - 1.222), n1[i][1])
                gamma_pi += -1*n1[i][2]*n1[i][0] *sp.power((7.1 - Pi), (n1[i][0]-1)) * sp.power((Tau - 1.222), n1[i][1])
                gamma_pi2 += -1*n1[i][2]*n1[i][0]*(n1[i][0]-1)*sp.power((7.1-Pi),(n1[i][0]-2))*sp.power((Tau - 1.222), n1[i][1])
                gamma_tau += n1[i][2]*n1[i][1]*sp.power((7.1-Pi),(n1[i][0]))*sp.power((Tau-1.222), n1[i][1]-1)
                gamma_tau2 += n1[i][2]*n1[i][1]*(n1[i][1]-1)*sp.power((7.1-Pi),(n1[i][0]))*sp.power((Tau-1.222), n1[i][1]-2)
                gamma_pi_tau += n1[i][2]*n1[i][0]*n1[i][1]*sp.power((7.1-Pi),(n1[i][0]-1))*sp.power((Tau-1.222), n1[i][1]-1)

            nu = Pi*gamma_pi * R * T/P/1e6
            u = (Tau*gamma_tau-Pi*gamma_pi) * R * T
            s = (Tau*gamma_tau-gamma)*R
            h = (Tau*gamma_tau)*R*T
            Cp = (-1*Tau*Tau*gamma_tau2)*R
            Cv = -1*Tau*Tau*gamma_tau2
            Cv += sp.power((gamma_pi-Tau*gamma_pi_tau),2)/gamma_pi2
            a2 = gamma_pi*gamma_pi
            c1 = sp.power((gamma_pi-Tau*gamma_pi_tau),2)
            c2 = Tau*Tau*gamma_tau2
            c3 = c1/c2-gamma_pi2
            a2 /= c3
            a2 *= R * T
            a = sp.sqrt(abs(a2))
            rho = 1/nu

        if self.region == "2" or self.region == "2.1":
            Pi =P/1.0
            Tau = 540.0/T
            if self.region == "2":
                n1 = self.tableA4
                n2 = self.tableA5
            else:
                n1 = self.tableA4a
                n2 = self.tableA7


            gamma0 = sp.log(Pi)
            gamma_pi_0 = 1.0/Pi
            gamma_pi2_0 = -1.0/Pi/Pi

            for i in range(len(n1)):
                gamma0 += n1[i][1] * sp.power(Tau, n1[i][0])
                gamma_tau_0 += n1[i][1]*n1[i][0]*sp.power(Tau,n1[i][0]-1)
                gamma_tau2_0 += n1[i][1]*n1[i][0]*(n1[i][0]-1)*sp.power(Tau,n1[i][0]-2)
            for i in range(len(n2)):
                gamma1 += n2[i][2] * sp.power(Pi, n2[i][0]) * sp.power((Tau - 0.5),n2[i][1])
                gamma_pi_1 += n2[i][2] * n2[i][0]* sp.power(Pi, n2[i][0]-1) * sp.power((Tau - 0.5),n2[i][1])
                gamma_pi2_1 += n2[i][2]*n2[i][0]*(n2[i][0]-1)*sp.power(Pi, n2[i][0]-2)*sp.power((Tau-0.5),n2[i][1])
                gamma_tau_1 += n2[i][2] * n2[i][1]* sp.power(Pi, n2[i][0]) * sp.power((Tau - 0.5),n2[i][1]-1)
                gamma_tau2_1 += n2[i][2]*n2[i][1]*(n2[i][1]-1)*sp.power(Pi, n2[i][0]) * sp.power((Tau - 0.5),n2[i][1]-2)
                gamma_pi_tau_1 += n2[i][2]*n2[i][0]*n2[i][1]*sp.power(Pi, n2[i][0]-1)*sp.power((Tau-0.5),n2[i][1]-1)
            gamma = gamma0 + gamma1
            gamma_pi = gamma_pi_0 + gamma_pi_1
            gamma_pi2 = gamma_pi2_0 + gamma_pi2_1
            gamma_tau = gamma_tau_0 + gamma_tau_1
            gamma_tau2 = gamma_tau2_0 + gamma_tau2_1
            gamma_pi_tau = gamma_pi_tau_0 + gamma_pi_tau_1

            nu = Pi*gamma_pi * R * T/P/1e6
            u = (Tau*gamma_tau-Pi*gamma_pi) * R * T
            s = (Tau*gamma_tau-gamma)*R
            h = (Tau*gamma_tau)*R*T
            Cp = (-1*Tau*Tau*gamma_tau2)*R
            Cv = -1*Tau*Tau*gamma_tau2
            Cv -= sp.power((1+Pi*gamma_pi_1-Tau*gamma_pi_tau_1),2)/(1-Pi*Pi*gamma_pi2_1)
            Cv *= R
            c1 = 1.0+2.0*Pi*gamma_pi_1+Pi*Pi*gamma_pi_1*gamma_pi_1
            c2 = 1.0-Pi*Pi*gamma_pi2_1
            c3 = 1+Pi*gamma_pi_1-Tau*Pi*gamma_pi_tau_1
            c3 = sp.power(c3,2.0)
            c4 = Tau*Tau*gamma_tau2
            a2 = c1/(c2+c3/c4)
            a2 *= R * T
            a = sp.sqrt(abs(a2))
            rho = 1/nu

        if self.region == "3":
            rho = P             # Density and Temperature are used for region 3
            delta = rho/rho_c
            Tau = Tc/T
            n1 = self.tableA9

            phi = n1[0][2] * sp.log(Pi)
            phi_delta = n1[0][2] / delta
            phi_delta2 = -n1[0][2]/delta/delta
            phi_tau = 0.0
            phi_tau2 = 0.0
            phi_delta_tau = 0.0


            for i in range(1, len(n1)):
                phi += n1[i][2] * sp.power(delta, n1[i][0]) * sp.power(Tau, n1[i][1])
                phi_delta += n1[i][2] * n1[i][0] * sp.power(delta, n1[i][0]-1) * sp.power(Tau, n1[i][1])
                phi_delta2 += n1[i][2]*n1[i][0]*(n1[i][0]-1)*sp.power(delta, n1[i][0]-2)*sp.power(Tau, n1[i][0])
                phi_tau += n1[i][2]*n1[i][1]*sp.power(delta,n1[i][0])*sp.power(Tau,n1[i][1]-1)
                phi_tau2 += n1[i][2]*n1[i][1]*(n1[i][1]-1)*sp.power(delta,n1[i][0])*sp.power(Tau,n1[i][1]-2)
                phi_delta_tau += n1[i][2]*n1[i][0]*n1[i][1]*sp.power(delta,n1[i][0]-1)*sp.power(Tau,n1[i][1]-1)

            f = phi * R * T
            Pressure = delta * phi_delta * rho * R * T
            u = Tau * phi_tau * R * T
            s = (Tau * phi_tau - phi) * R
            h = (Tau * phi_tau + delta * phi_delta) * R * T
            Cv = -1 * Tau * Tau * phi_tau2 * R

            c1 = sp.power(delta * phi_delta - delta*tau*phi(delta_tau), 2)
            c2 = delta * (2 * phi_delta +delta * delta * phi_delta2)
            Cp = Cv + c1/c2*R

            a2 = c2 - c1 / Tau/Tau/phi_tau2
            a2 *= R * T
            a = sp.sqrt(abs(a2))


        if self.region == "4":
            n1 = self.tableA11
            self.calcSatLine(P,T)


        if self.region == "5":
            Pi = P/1.0
            Tau = 1000.0/T
            n1 = self.tableA13
            n2 = self.tableA14
            gamma0 = sp.log(Pi)
            gamma_pi_0 = 1.0/Pi
            gamma_pi2_0 = -1.0/Pi/Pi
            gamma_tau_0 = 1.0/Pi

            for i in range(len(n1)):
                gamma0 += n1[i][1] * sp.power(Tau, n1[i][0])
                gamma_tau_0 += n1[i][1]*n1[i][0]*sp.power(Tau, n1[i][0]-1)
                gamma_tau2_0 += n1[i][2]*n1[i][1]*(n1[i][1]-1)*sp.power(Tau, n1[i][1]-2)

            for i in range(len(n2)):
                gamma1 += n2[i][2] * sp.power(Pi, n2[i][0]) * sp.power(Tau, n2[i][1])
                gamma_pi_1 += n2[i][2]*n2[i][0]*sp.power(Pi, n2[i][0]-1) * sp.power(Tau, n2[i][1])
                gamma_pi2_1 += n2[i][2]*n2[i][0]*(n2[i][0]-1)* sp.power(Pi, n2[i][0]-2)*sp.power(Tau, n2[i][1])
                gamma_tau_1 += n2[i][2]*n2[i][1]*sp.power(Pi, n2[i][0]) * sp.power(Tau, n2[i][1]-1)
                gamma_tau2_1 += n2[i][2]*n2[i][1]*(n2[i][1]-1)*sp.power(Pi, n2[i][0])*sp.power(Tau, n2[i][1]-2)
                gamma_pi_tau_1 += n2[i][2]*n2[i][0]*n2[i][1]*sp.power(Pi, n2[i][0]-1) * sp.power(Tau, n2[i][1]-1)

            gamma = gamma0 +gamma1
            gamma_pi = gamma_pi_0 +gamma_pi_1
            gamma_pi2 = gamma_pi2_0 +gamma_pi2_1
            gamma_tau = gamma_tau_0 +gamma_tau_1
            gamma_tau2 = gamma_tau2_0 +gamma_tau2_1
            gamma_pi_tau = gamma_pi_tau_0 +gamma_pi_tau_1

            nu = Pi*gamma_pi * R * T/P/1e6
            u = (Tau*gamma_tau-Pi*gamma_pi) * R * T
            s = (Tau*gamma_tau-gamma)*R
            h = (Tau*gamma_tau)*R*T
            Cp = (-1*Tau*Tau*gamma_tau2)*R
            Cv = -1*Tau*Tau*gamma_tau2
            Cv -= sp.power((1+Pi*gamma_pi_1-Tau*gamma_pi_tau_1),2)/(1-Pi*Pi*gamma_pi2_1)
            Cv *= R
            c1 = 1.0+2.0*Pi*gamma_pi_1+Pi*Pi*gamma_pi_1*gamma_pi_1
            c2 = 1.0-Pi*Pi*gamma_pi2_1
            c3 = 1+Pi*gamma_pi_1-Tau*Pi*gamma_pi_tau_1
            c3 = sp.power(c3,2.0)
            c4 = Tau*Tau*gamma_tau2

            a2 = c1/(c2+c3/c4)
            a2 *= R * T
            a = sp.sqrt(abs(a2))
            rho = 1/nu

        if self.region == "3":
            Properties = [Pressure, u, s, h, Cv, Cp, a, P, T, rho]
        Properties = [nu, u, s, h, Cp, Cv, a, P, T, rho]

        return Properties


    def verifyAll(self):
        for region in ["1", "2", "2.1", "3", "5"]:
            print "\n\nRegion ", Region
            Conditions, table, table_no = self.selectConditions()
            for (P,T) in Conditions:
                Properties = self.calcProperties(P,T)
                errors = self.calcErrors(Properties, table)

            self.printResults(Conditions, Properties, errors, table_no)

    def selectConditions(self):
        """A function to set the conditions and applicable table to be used when verifying results.
        Input
        -----
        self.region:
            The region in which you want to verify the results. Refer to figure 2 of Wagner (2000) for an explanation of the regions. This variable is set using the findRegion function.
        Outputs
        -------
        Conditions:
            List-like. Contains the verification conditions in the order [Pressure (MPa), Temperature (K)]. Each line represents a different verification condition.
        table:
            index. Links to a table of verification data listed in Wagner (2000).
        Notes
        -----
        The following conditions and tables are used.
        Region 1:
            3 Mpa, 300 K; 80 MPa, 300 K; 3 MPa, 500 K
            Table A3
        Region 2:
            0.0035 MPa, 300 K; 0.0035 MPa, 700 K; 30 MPa, 700 K
            Table A6
        Region 2 Metastable:
            1 MPa, 450 K; 1 MPa, 440 K; 1.5 MPa 450 K
            Table A8
        Region 3:


        """

        if self.region == "1":
            Conditions=[[3.0, 300.0], [80.0, 300.0], [3.0, 500.0]]
            table = self.tableA3
            table_no = "A3"
        if self.region == "2":
            Conditions = [[0.0035, 300.0], [0.0035, 700.0], [30.0, 700.0]]
            table = self.tableA6
            table_no = "A6"
        if self.region == "2.1":
            Conditions = [[1.0, 450.0], [1.0, 440.0], [1.5, 450.0]]
            table = self.tableA8
            table_no = "A8"
        if self.region == "3":
            Conditions = [[500.0, 650.0], [200.0, 650.0], [500.0, 750.0]]
            table = self.tableA10
            table_no = "A10"
        if self.region == "5":
            Conditions = [[0.5, 1500.0], [8.0, 1500.0], [8.0, 2000.0]]
            table = self.tableA15
            table_no = "A15"

        return Conditions, table, table_no

    def calcErrors(self, values, table):
        errors = []
        for i in range(len(values)):
            for j in range(len(values[i])):
                errors[i][j] = (values[i][j] - table[i][j])/table[i][j]*100
        return errors

    def printResults(self, Conditions, values, errors, table_no):
        s = "Table "+ table_no+ " Conditions and Properties"
        print s
        for i in range(len(Conditions)):
            print "\nCondition ", i+1
            print "Pressure = ", Conditions[i][0], " MPa"
            print "Temperature = ", Conditions[i][1], " K"
            for j in range(len(properties)):
                print properties[j], ": ", values[i][j], units[j]

        print "\nRelative Errors in calculated values versus Tabulated Values"
        for i in range(len(Conditions)):
            print "\nCondition ", i+1
            for j in range(len(properties)):
                print properties[j], ": ", errors[i][j], "%"

