#! /usr/bin/env python
# spectral_functions.py
#
# A collection of functions required for spectral processing in the X-labs
#
# Fabian Zander
# Mon 6th February 2012
# Dungeon, St Lucia
#
# The file contains a collection of different functions for processing spectra in regards 
# to grating efficiency's, pmt quantum efficiencies etc. etc.

import matplotlib.pyplot as mplt

"""
All the grating functions are taken from the data sheets that come with the spectrometer.
BE CAREFUL as the data is only supplied as graphs and these numbers are FZ's interpretation
of these plots - if you need high accuracy you need to do more calibration
"""

def ir_150(wavelength, intensity):
    """
    IR grating 150g/mm BLZ 500nm
    This function converts a 'read' intensity into an actual intensity going into the
    spectrometer based on the gratings efficiency.
    Linear interpolation used between points.
                             ir_150(wavelength, intensity) = 
    This gives pre-grating efficiency ie. ir_150(500, 100) = 136.99
    """
    efficiency = [[300,	0.18], [305, 0.202], [310, 0.224], [315, 0.246], [320, 0.268], [325, 0.29], [330, 0.312], [335,	0.334], [340, 0.356], [345, 0.378], [350, 0.4],
                  [355,	0.42], [360, 0.44], [365, 0.46], [370, 0.48], [375, 0.5], [380,	0.52], [385, 0.54], [390, 0.56], [395, 0.58], [400, 0.6], [405,	0.61],
                  [410,	0.62], [415, 0.63], [420, 0.64], [425, 0.65], [430, 0.66], [435, 0.67], [440, 0.68], [445, 0.69], [450,	0.7], [455, 0.703], [460, 0.706],
                  [465,	0.709], [470, 0.712], [475, 0.715], [480, 0.718], [485,	0.721], [490, 0.724], [495, 0.727], [500, 0.73], [505, 0.727], [510, 0.724], [515, 0.721],
                  [520,	0.718], [525, 0.715], [530, 0.712], [535, 0.709], [540,	0.706], [545, 0.703], [550, 0.7], [555,	0.697], [560, 0.694], [565, 0.691], [570, 0.688],
                  [575,	0.685], [580, 0.682], [585, 0.679], [590, 0.676], [595,	0.673], [600, 0.67], [605, 0.663], [610, 0.656], [615, 0.649], [620, 0.642], [625, 0.635],
                  [630,	0.628], [635, 0.621], [640, 0.614], [645, 0.607], [650,	0.6], [655, 0.594], [660, 0.588], [665,	0.582], [670, 0.576], [675, 0.57], [680, 0.564],
                  [685,	0.558], [690, 0.552], [695, 0.546], [700, 0.54], [705,	0.534], [710, 0.528], [715, 0.522], [720, 0.516], [725,	0.51], [730, 0.504], [735, 0.498],
                  [740,	0.492], [745, 0.486], [750, 0.48], [755, 0.475], [760, 0.47], [765, 0.465], [770, 0.46], [775, 0.455], [780, 0.45], [785, 0.445], [790,	0.44],
                  [795,	0.435], [800, 0.43], [850, 0.4], [900, 0.37], [950, 0.34], [1000, 0.31]]
    signal_perc = 0
    for i in range(len(efficiency)-1):
        if (efficiency[i][0] == wavelength):
            signal_perc = efficiency[i][1]
        elif (efficiency[i][0] < wavelength and efficiency[i+1][0] > wavelength):
            signal_perc = ( (efficiency[i][1] - efficiency[i+1][1]) / (efficiency[i][0] - efficiency[i+1][0]) ) * ( wavelength - efficiency[i][0] ) + efficiency[i][1] 
            #             { (      y1         -       y2          ) / (        x1       -          x2       ) } * {     x      -       x1         } +        y1
    adj_int = intensity / signal_perc
    return adj_int

def ir_150_inv(wavelength, intensity):
    """
    IR grating 150g/mm BLZ 500nm
    This function calculates the intensity coming out of the spectrometer for a given input.
    Linear interpolation used between points.
                             ir_150(wavelength, intensity) = 
    This gives pre-grating efficiency ie. ir_150(500, 100) = 73
    """
    efficiency = [[300,	0.18], [305, 0.202], [310, 0.224], [315, 0.246], [320, 0.268], [325, 0.29], [330, 0.312], [335,	0.334], [340, 0.356], [345, 0.378], [350, 0.4],
                  [355,	0.42], [360, 0.44], [365, 0.46], [370, 0.48], [375, 0.5], [380,	0.52], [385, 0.54], [390, 0.56], [395, 0.58], [400, 0.6], [405,	0.61],
                  [410,	0.62], [415, 0.63], [420, 0.64], [425, 0.65], [430, 0.66], [435, 0.67], [440, 0.68], [445, 0.69], [450,	0.7], [455, 0.703], [460, 0.706],
                  [465,	0.709], [470, 0.712], [475, 0.715], [480, 0.718], [485,	0.721], [490, 0.724], [495, 0.727], [500, 0.73], [505, 0.727], [510, 0.724], [515, 0.721],
                  [520,	0.718], [525, 0.715], [530, 0.712], [535, 0.709], [540,	0.706], [545, 0.703], [550, 0.7], [555,	0.697], [560, 0.694], [565, 0.691], [570, 0.688],
                  [575,	0.685], [580, 0.682], [585, 0.679], [590, 0.676], [595,	0.673], [600, 0.67], [605, 0.663], [610, 0.656], [615, 0.649], [620, 0.642], [625, 0.635],
                  [630,	0.628], [635, 0.621], [640, 0.614], [645, 0.607], [650,	0.6], [655, 0.594], [660, 0.588], [665,	0.582], [670, 0.576], [675, 0.57], [680, 0.564],
                  [685,	0.558], [690, 0.552], [695, 0.546], [700, 0.54], [705,	0.534], [710, 0.528], [715, 0.522], [720, 0.516], [725,	0.51], [730, 0.504], [735, 0.498],
                  [740,	0.492], [745, 0.486], [750, 0.48], [755, 0.475], [760, 0.47], [765, 0.465], [770, 0.46], [775, 0.455], [780, 0.45], [785, 0.445], [790,	0.44],
                  [795,	0.435], [800, 0.43], [850, 0.4], [900, 0.37], [950, 0.34], [1000, 0.31]]
    signal_perc = 0
    for i in range(len(efficiency)-1):
        if (efficiency[i][0] == wavelength):
            signal_perc = efficiency[i][1]
        elif (efficiency[i][0] < wavelength and efficiency[i+1][0] > wavelength):
            signal_perc = ( (efficiency[i][1] - efficiency[i+1][1]) / (efficiency[i][0] - efficiency[i+1][0]) ) * ( wavelength - efficiency[i][0] ) + efficiency[i][1] 
            #             { (      y1         -       y2          ) / (        x1       -          x2       ) } * {     x      -       x1         } +        y1
    adj_int = intensity * signal_perc
    return adj_int

def uv_150(wavelength, intensity):
    """
    UV grating 150g/mm BLZ 300nm
    This function converts a 'read' intensity into an actual intensity going into the
    spectrometer based on the gratings efficiency.
    Linear interpolation used between points. 
    NOTE: The [0,0] point at the start has been added to avoid errors but it should 
    but there is no calibration under 200nm.
                             uv_150(wavelength, intensity) = 
    This gives pre-grating efficiency ie. uv_150(500, 100) = 200
    """
    efficiency = [[0.0, 0.0], [200, 0.28], [250, 0.45], [300, 0.58], [350, 0.72], [400, 0.72], [450, 0.57], [500, 0.50], [550, 0.43], [600, 0.38], [650, 0.33], [700, 0.29], 
                  [750, 0.25], [800, 0.22], [850, 0.2], [900, 0.17], [950, 0.15], [1000, 0.13]]
    signal_perc = 0
    for i in range(len(efficiency)-1):
        if (efficiency[i][0] == wavelength):
            signal_perc = efficiency[i][1]
        elif (efficiency[i][0] < wavelength and efficiency[i+1][0] > wavelength):
            signal_perc = ( (efficiency[i][1] - efficiency[i+1][1]) / (efficiency[i][0] - efficiency[i+1][0]) ) * ( wavelength - efficiency[i][0] ) + efficiency[i][1] 
            #             { (      y1         -       y2          ) / (        x1       -          x2       ) } * {     x      -       x1         } +        y1
    adj_int = intensity / signal_perc
    return adj_int

def uv_600(wavelength, intensity):
    """
    UV Grating 600g/mm BLZ 300nm
    This function account for the efficiency of the UV spectrometer
    grating. A given input intensity will be returned as a scaled
    value based on the Grating Efficiency Curve
    www.princetoninstruments.com/Uploads/Princeton/Images/Grating%20curves/1-060-300-1336-Sept2010.gif
    Linear interpolation is used between data points.
                             uv_600(wavelength, intensity) = 
    This gives pre-grating efficiency ie. uv_600(500, 100) = 200
    """
    # Using the S-Plane data at this stage - I have no idea what S & P plane are.....I do now know but the data will stay the same for the moment
    efficiency = [[200, 0.28], [220, 0.38], [240, 0.46], [260, 0.54], [280, 0.60], [300, 0.63], [320, 0.65], [340, 0.67], [360, 0.66], [380, 0.64], [400, 0.62],
                  [420, 0.59], [440, 0.56], [460, 0.52], [480, 0.50], [500, 0.50], [520, 0.46], [540, 0.42], [560, 0.38], [580, 0.36], [600, 0.34], [620, 0.32],
                  [640, 0.30], [660, 0.29], [680, 0.28], [700, 0.31], [720, 0.27], [740, 0.24], [760, 0.22], [780, 0.20], [800, 0.19], [820, 0.17], [840, 0.17],
                  [860, 0.16], [880, 0.15], [900, 0.15], [901, 0.15]]
    signal_perc = 0
    for i in range(len(efficiency)-1):
        if (efficiency[i][0] == wavelength):
            signal_perc = efficiency[i][1]
        elif (efficiency[i][0] < wavelength and efficiency[i+1][0] > wavelength):
            signal_perc = ( (efficiency[i][1] - efficiency[i+1][1]) / (efficiency[i][0] - efficiency[i+1][0]) ) * ( wavelength - efficiency[i][0] ) + efficiency[i][1] 
            #             { (      y1         -       y2          ) / (        x1       -          x2       ) } * {     x      -       x1         } +        y1
    adj_int = intensity / signal_perc
    return adj_int

def tungsten_lamp_output():
    """
    This function simply returns a list giving the tungsten calibration lamp output.
    Linear interpolation has been used between some data points. This hasn't been done
    for the longer wavelengths as it doesn't interest FZ at the moment.
    """
    lamp_output = [[300, 2.08E-007], [305, 2.5176265E-007], [310, 2.95E-007], [315, 3.4854415E-007], [320, 4.02E-007], [325, 4.699143E-007], [330, 5.38E-007], [335, 6.205823E-007],
                   [340, 7.03E-007], [345, 8.029409E-007], [350, 9.03E-007], [355, 1.01833785E-006], [360, 1.13E-006], [365, 1.2676995E-006], [370, 1.40E-006], [375, 1.5538285E-006],
                   [380, 1.71E-006], [385, 1.882886E-006], [390, 2.06E-006], [395, 2.250851E-006], [400, 2.44E-006], [405, 2.6902056E-006], [410, 2.9381022E-006], [415, 0.000003186],
                   [420, 3.4338954E-006], [425,	3.681792E-006], [430, 3.9296886E-006], [435, 4.1775852E-006], [440, 4.4254818E-006], [445, 4.6733784E-006], [450, 4.92E-006],
                   [455, 5.2446827E-006], [460,	5.5680904E-006], [465, 5.8914981E-006], [470, 6.2149058E-006], [475, 6.5383135E-006], [480, 6.8617212E-006], [485, 7.1851289E-006],
                   [490, 7.5085366E-006], [495,	7.8319443E-006], [500, 8.16E-006], [505, 8.51454636363637E-006], [510, 8.87374072727273E-006], [515, 9.23293509090909E-006], 
                   [520, 9.59212945454545E-006], [525, 9.95132381818182E-006], [530, 1.03105181818182E-005], [535, 1.06697125454545E-005], [540, 1.10289069090909E-005], 
                   [545, 1.13881012727273E-005], [550, 1.17472956363636E-005], [555, 1.21E-005], [560, 1.24558611111111E-005], [565, 1.28052322222222E-005], [570, 1.31546033333333E-005],
                   [575, 0.000013504], [580, 1.38533455555556E-005], [585, 1.42027166666667E-005], [590, 1.45520877777778E-005], [595, 1.49014588888889E-005], [600, 1.53E-005],
                   [605, 1.55684254545455E-005], [610, 0.000015886], [615, 1.62036163636364E-005], [620, 1.65212118181818E-005], [625, 1.68388072727273E-005], [630, 1.71564027272727E-005],
                   [635, 0.000017474], [640, 1.77915936363636E-005], [645, 1.81091890909091E-005], [650, 1.84267845454545E-005], [655, 1.87E-005], [660, 1.90111611111111E-005],
                   [665, 1.92779422222222E-005], [670, 1.95447233333333E-005], [675, 1.98115044444444E-005], [680, 2.00782855555556E-005], [685, 2.03450666666667E-005],
                   [690, 2.06118477777778E-005], [695, 2.08786288888889E-005], [700, 2.11E-005], [705, 0.000021315], [710, 2.1484506E-005], [715, 2.1654054E-005], [720, 2.1823602E-005],
                   [725, 2.199315E-005], [730, 2.2162698E-005], [735, 2.2332246E-005], [740, 2.2501794E-005], [745, 2.2671342E-005], [750, 2.284089E-005], [755, 2.3010438E-005],
                   [760, 0.00002318], [765, 2.3349534E-005], [770, 2.3519082E-005], [775, 2.368863E-005], [780,	2.3858178E-005], [785, 2.4027726E-005], [790, 2.4197274E-005],
                   [795, 2.4366822E-005], [800,	2.45E-005], [900, 2.563871E-05], [1050, 2.433556E-05], [1100, 2.344202E-05]]
    return lamp_output

def uv_iccd_qe(wavelength, intensity):
    """
    This function applies the inverse quantum efficiency to a signal and return the input intensity.
    This data is the quantum efficiency for a PI-MAX ICCD 'SB Slow Gate' as quoted by TE and US.
    NOTE: the (0,0) point has been added to avoid errors but there is no calibration under 185nm.
    """
    qe_iccd = [[0, 0], [205, 0.2757], [241, 0.2757], [264, 0.2854], [271.5, 0.2854], [279.2, 0.28125], [296.5, 0.2653], [300, 0.263], [311, 0.263], [315.3, 0.2646], [331, 0.2653],
              [344.4, 0.2465], [352, 0.244], [354, 0.24375], [356, 0.2452], [357, 0.249], [359, 0.254], [361, 0.254], [379.8, 0.243], [400, 0.2284], [450, 0.1854],
              [481.25, 0.15], [500, 0.134], [550, 0.093], [600, 0.0576], [650, 0.0305], [700, 0.0097], [750, 0.0]]

    signal_perc = 0
    for i in range(len(qe_iccd)-1):
        if (qe_iccd[i][0] == wavelength):
            signal_perc = qe_iccd[i][1]
        elif (qe_iccd[i][0] < wavelength and qe_iccd[i+1][0] > wavelength):
            signal_perc = ( (qe_iccd[i][1] - qe_iccd[i+1][1]) / (qe_iccd[i][0] - qe_iccd[i+1][0]) ) * ( wavelength - qe_iccd[i][0] ) + qe_iccd[i][1] 
            #             { (    y1        -       y2       ) / (     x1       -       x2       ) } * {     x      -        x1     } +     y1
    return intensity / signal_perc

def ir_iccd_qe(wavelength, intensity):
    """
    This function applies the inverse quantum efficiency to a signal and return the input intensity.
    This data is the quantum efficiency for a PI-MAX ICCD 'Unigen' as supplied by Tim - including
    graph data sheet
    NOTE: these values have been read from a graph by FZ
    NOTE: the (0,0) point has been added to avoid errors but there is no calibration under 155nm.
    """
    qe_iccd = [[0, 0], [150, 0.19], [200, 0.20], [250, 0.18], [300, 0.14], [350, 0.10], [400, 0.17], [450, 0.35], [470, 0.38], [500, 0.39], [550, 0.40],
              [600, 0.38], [700, 0.36], [750, 0.32], [800, 0.35], [825, 0.325], [850, 0.255], [900, 0.0]]

    signal_perc = 0
    for i in range(len(qe_iccd)-1):
        if (qe_iccd[i][0] == wavelength):
            signal_perc = qe_iccd[i][1]
        elif (qe_iccd[i][0] < wavelength and qe_iccd[i+1][0] > wavelength):
            signal_perc = ( (qe_iccd[i][1] - qe_iccd[i+1][1]) / (qe_iccd[i][0] - qe_iccd[i+1][0]) ) * ( wavelength - qe_iccd[i][0] ) + qe_iccd[i][1] 
            #             { (     y1       -       y2       ) / (     x1       -       x2       ) } * {     x      -        x1     } +     y1
    return intensity / signal_perc

def pmt_qe(wavelength, intensity):
    """
    This function applies the inverse quantum efficiency to a signal and return the input intensity.
    This data is the quantum efficiency for a Hamamatsu R928 Photomultiplier tube.
    NOTE: the (0,0) point has been added to avoid errors but there is no calibration under 185nm.
    """
    qe_pmt = [[0, 0], [185, 0.0804], [190, 0.0979], [200, 0.158], [210, 0.218], [220, 0.242], [230, 0.253], [240, 0.253], [250, 0.254], [260, 0.254], [270, 0.248], [280, 0.246], 
              [290, 0.244], [300, 0.238], [310, 0.235], [320, 0.233], [330, 0.233], [340, 0.233], [350, 0.236], [360, 0.237], [370, 0.237], [380, 0.238], [390, 0.234], 
              [400, 0.229], [410, 0.223], [420, 0.216], [430, 0.208], [440, 0.200], [450, 0.193], [460, 0.185], [470, 0.178], [480, 0.171], [490, 0.163], [500, 0.154], 
              [510, 0.147], [520, 0.140], [530, 0.133], [540, 0.126], [550, 0.120], [560, 0.113], [570, 0.107], [580, 0.102], [590, 0.0975], [600, 0.0930], [610, 0.0888], 
              [620, 0.0848], [630, 0.0817], [640, 0.0771], [650, 0.0734], [660, 0.0699], [670, 0.0664], [680, 0.0631], [690, 0.0598], [700, 0.0567], [710, 0.0536],
              [720, 0.0506], [730, 0.0477], [740, 0.0449], [750, 0.0422], [760, 0.0388], [770, 0.0358], [780, 0.0326], [790, 0.0297], [800, 0.0267], [810, 0.0219], 
              [820, 0.0175], [830, 0.0134], [840, 0.0094], [850, 0.0058], [860, 0.0026], [870, 0.0011], [880, 0.0005], [890, 0.0002], [900, 0.0001], [901, 0.0001]]
    signal_perc = 0
    for i in range(len(qe_pmt)-1):
        if (qe_pmt[i][0] == wavelength):
            signal_perc = qe_pmt[i][1]
        elif (qe_pmt[i][0] < wavelength and qe_pmt[i+1][0] > wavelength):
            signal_perc = ( (qe_pmt[i][1] - qe_pmt[i+1][1]) / (qe_pmt[i][0] - qe_pmt[i+1][0]) ) * ( wavelength - qe_pmt[i][0] ) + qe_pmt[i][1] 
            #             { (    y1       -      y2       ) / (    x1       -      x2       ) } * {     x      -       x1     } +     y1
    return intensity / signal_perc

def plot_2d_spectra(wave, spec_image_data, x_size=14, y_size=9, im_aspect=1.0, max_int=50000):
    """
    This function creates a 2d spectral plot with wavelength on the x axis and pixel position on the y axis.
    NOTE: that as of writing this plots the spectra 'upside down' compared to winspec.
    Inputs: wave - the wavelengths, as taken from spe file using pyspec
            spec_image_data - the spectral data as taken from spe file using pyspec
            x_size - figure size in the x direction
            y_size - figure size in the y direction
            im_aspect - is the aspect for the plot part of the image, changes depending on grating
                                   recommend 1.0 for 150g/mm grating
                                             0.2 for 600g/mm grating
            max_int - the maximum value for the contour colour bar of the plot
    """
    mplt.figure(figsize=(x_size,y_size))
    im = mplt.imshow(spec_image_data[0], interpolation='bilinear', origin='lower', cmap=mplt.cm.hot, extent=(wave[0],wave[-1],0,256), aspect=im_aspect, vmax=max_int)
    mplt.title('Spectral plot', fontsize=28)
    mplt.xlabel('Wavelength (nm)', fontsize=16)
    mplt.ylabel('Position (pixel)', fontsize=16)
    CBI = mplt.colorbar(im, orientation='horizontal', shrink=1.0, pad=0.08)
    CBI.set_label('Intensity', fontsize=16)
    mplt.show()

def wave_index(desired_wave, wavelengths):
    """
    Find the index in wavelengths which best represents the desired wavelength.
    This is used in other functions within this file but can be used elsewhere as well.
    """
    wavelength_index = 0
    for i in range(len(wavelengths)):
        if (wavelengths[i] == desired_wave): wavelength_index = i
        elif (wavelengths[i] < desired_wave and wavelengths[i+1] > desired_wave):
            wavelength_index = i
            break
    return wavelength_index

def plot_wave(plot_wave, wave_data, spec_image_data):
    """
    Plot a single wavelength
    """
    index = wave_index(plot_wave, wave_data)
    plot_data = []
    for item in spec_image_data[0]: plot_data.append(item[index])
    mplt.figure()
    mplt.plot(plot_data)
    mplt.title('%.2fnm plot v. Pixel position' % plot_wave, fontsize=28)
    mplt.xlabel('Pixel position', fontsize=16)
    mplt.ylabel('Intensity', fontsize=16)
    mplt.xlim(0,256)
    mplt.show()

def plot_slit_pos(pixel_index, wave_data, spec_image_data):
    """
    Plot a single slit position v. wavelength
    """
    mplt.figure()
    mplt.plot(wave_data, spec_image_data[0][pixel_index])
    mplt.title('Pixel %i  v. Wavelength position' % pixel_index, fontsize=28)
    mplt.xlabel('Wavelength', fontsize=16)
    mplt.ylabel('Intensity', fontsize=16)
    mplt.show()
