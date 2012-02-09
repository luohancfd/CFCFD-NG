;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Marc Haeming inserted x1und_plot_edges to show where the several
;; edges are.  This should ease to find appropreate gap sizes

PRO x1und_plot_edges

OPlot, [285,285,300,300], [8.*(10.)^15,0,0,8.*(10.)^15], Color=150
OPlot, [540,540,560,560], [8.*(10.)^15,0,0,8.*(10.)^15], Color=150
OPlot, [400,400,420,420], [8.*(10.)^15,0,0,8.*(10.)^15], Color=150
OPlot, [700,700,730,730], [8.*(10.)^15,0,0,8.*(10.)^15], Color=150
OPlot, [840,840,860,860], [8.*(10.)^15,0,0,8.*(10.)^15], Color=150
OPlot, [920,920,960,960], [8.*(10.)^15,0,0,8.*(10.)^15], Color=150

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO X1und, k, gap_mm, lambda_nm, inten, $
           setup_onaxis=setup_onaxis, build_xdr=build_xdr, $
           psfile=psfile, epsfile=epsfile, $
           nodisplay=nodisplay, nocontour=nocontour, $
           plotgap=plotgap, ev=ev, nm=nm, xrange=xrange, $
           high_energy=high_energy, low_energy=low_energy, $
           help=help
  
  IF keyword_set(help) THEN BEGIN
     print, 'x1und,k,gap,lambda,inten'
     print, '  The values for k, gap (in mm), lambda (in nm), and'
     print, '    inten (photons/sec/250 mA/0.1%BW/mrad^2)'
     print, '    are returned to you, and are plotted as well.'
     print, '  Options:'
     print, '    /LOW_ENERGY   For 2.584 GeV (otherwise 2.800 GeV)'
     print, '    /NODISPLAY     Just returns numbers without displaying plot'
     print, '    /NOCONTOUR     No contour plot'
     print, '    /EPSFILE       '+$
            'Write plot to an EPS file (default is "x1und.eps")'
     print, '    /PSFILE       '+$
            'Write plot to a PostScript file (default is "x1und.ps")'
     print, '    /SETUP_ONAXIS  Make the file "x1_many_onaxis.com" for SMUT'
     print, '    /BUILD_XDR     Read SMUT MAPPER files and build ' + $
            '"x1und.xdr" or "x1und_28gev.xdr"'
     print, '    PLOTGAP=mm     Plots the spectrum for a gap '+$
            'in "mm" millimeters.'
     print, '      With PLOTGAP, use /EV or /NM to pick the X axis,'
     print, '      and XRANGE=[min,max] if you wish.'
     return
  ENDIF
  
  IF keyword_set(low_energy) THEN BEGIN
     xdr_filename = 'x1und.xdr'
     energy_text = '2.584'
  ENDIF ELSE BEGIN
     xdr_filename = 'x1und_28.xdr'
     energy_text = '2.80'
  ENDELSE
  
  EMASS = 511.003414e3          ; electron rest energy
  C = 2.99792458e8              ; Velocity of light in meters/sec
  B_REMNANT = 0.534             ; Remnant field of X1 permanent
                                ; magnets (Tesla).
                                ; From p. 198 of X1 Log IV.
  LAMBDA0 = 0.08                ; X1 undulator magnet period length in
                                ; meters.
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  IF keyword_set(setup_onaxis) THEN BEGIN
     filename = 'x1_many_onaxis.com'
     print, 'Creating file "'+filename+'"'
     
     k_step = 0.025
     k_min = 0.75
     k_max = 2.325
     n_k = fix(1.1+(k_max-k_min)/k_step)
     k = k_min+k_step*findgen(n_k)
     
     get_lun, lun
     openw, lun, filename
     printf, lun, '$ LOAD_UND'
     printf, lun, '$ SET DEFAULT UND_SPECTRUM_DIR:'
     printf, lun, '$ SET VERIFY'
     
     FOR i = 0, (n_k-1) DO BEGIN
        k_name = strtrim(string(1000.*k(i), format = '(i4)'), 2)
        IF (strlen(k_name) EQ 2) THEN BEGIN
           k_name = '00'+k_name
        ENDIF ELSE IF (strlen(k_name) EQ 3) THEN BEGIN
           k_name = '0'+k_name
        ENDIF
        
        printf, lun, '$ RUN ONAXIS_SPECTRUM'
        printf, lun, 'x1.par'
        printf, lun, '13.0'
        printf, lun, strtrim(string(k(i), format = '(f10.3)'), 2)
        printf, lun, '10.'
        printf, lun, '50.'
        printf, lun, '0.004'
        printf, lun, 'x1_k'+k_name+'_onaxis'
        printf, lun, '$ RUN ONAXIS_SPECTRUM'
        printf, lun, 'x1.par'
        printf, lun, '13.0'
        printf, lun, strtrim(string(k(i), format = '(f10.3)'), 2)
        printf, lun, '6.2'
        printf, lun, '10.'
        printf, lun, '0.002'
        printf, lun, 'x1_k'+k_name+'_highe_onaxis'
     ENDFOR
     
     printf, lun, '$ set noverify'
     
     close, lun
     free_lun, lun
     
     return
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ENDIF ELSE IF keyword_set(build_xdr) THEN BEGIN
     IF keyword_set(high_energy) THEN BEGIN
        print, 'Building file "'+xdr_filename+'"'
        
        k_step = 0.025
        k_min = 0.80
        k_max = 2.35
        n_k = fix(1.1+(k_max-k_min)/k_step)
        k = k_min+k_step*findgen(n_k)
        print, 'n_k = ', n_k, ', k_max = ', max(k)
        gap_mm = fltarr(n_k)

        test_gap = 30.+(80.-30.)*findgen(2001)/float(2000)

        test_gap_meters = 1.e-3*test_gap
        test_b0 = B_REMNANT/sinh(!Pi*test_gap_meters/LAMBDA0)
        test_k = C*LAMBDA0*test_b0/(EMASS*2.0*!Pi)

        FOR i = 0, (n_k-1) DO BEGIN
           dummy = min(abs(k(i)-test_k), test_index)
           gap_mm(i) = test_gap(test_index)
        ENDFOR
        
        k_name = strtrim(string(1000.*k(0), format = '(i4)'), 2)
        IF (strlen(k_name) EQ 2) THEN BEGIN
           k_name = '00'+k_name
        ENDIF ELSE IF (strlen(k_name) EQ 3) THEN BEGIN
           k_name = '0'+k_name
        ENDIF
        
        lowe_mapfile = 'x1_28gev_k'+k_name+'_l100_l500_onaxis.map'
        print, 'Reading "'+lowe_mapfile+'"'
        read_mapper, lowe_mapfile, lowe_lambda, lowe_inten, /quiet
        n_lambda = n_elements(lowe_lambda)
        lambda_nm = 0.1*lowe_lambda
        
        inten = fltarr(n_lambda, n_k)
        
        FOR i = 0, (n_k-1) DO BEGIN
           k_name = strtrim(string(1000.*k(i), format = '(i4)'), 2)
           IF (strlen(k_name) EQ 2) THEN BEGIN
              k_name = '00'+k_name
           ENDIF ELSE IF (strlen(k_name) EQ 3) THEN BEGIN
              k_name = '0'+k_name
           ENDIF
           
           lowe_mapfile = 'x1_28gev_k'+k_name+'_l100_l500_onaxis.map'
           print, 'Reading "'+lowe_mapfile+'"'
           read_mapper, lowe_mapfile, lowe_lambda, lowe_inten, /quiet
           inten(0:(n_lambda-1), i) = lowe_inten
        ENDFOR
        
        get_lun, lun
        openw, lun, xdr_filename, /xdr
        writeu, lun, fix(n_k), fix(n_lambda)
        writeu, lun, float(k), float(gap_mm), float(lambda_nm), float(inten)
        close, lun
        free_lun, lun
        
        print, 'Built file "'+xdr_filename+'" with '+$
               strtrim(string(n_k, format = '(i5)'), 2)+' K values.'
        print,'Move the file to an appropriate directory'
        return
     ENDIF ELSE BEGIN
        print, 'Building file "'+xdr_filename+'"'
        
        k_step = 0.025
        k_min = 0.75
        k_max = 2.325
        n_k = fix(1.1+(k_max-k_min)/k_step)
        k = k_min+k_step*findgen(n_k)
        print, 'n_k = ', n_k, ', k_max = ', max(k)
        gap_mm = fltarr(n_k)
        
        test_gap = 30.+(80.-30.)*findgen(2001)/float(2000)
        
        test_gap_meters = 1.e-3*test_gap
        test_b0 = B_REMNANT/sinh(!Pi*test_gap_meters/LAMBDA0)
        test_k = C*LAMBDA0*test_b0/(EMASS*2.0*!Pi)
        
        FOR i = 0, (n_k-1) DO BEGIN
           dummy = min(abs(k(i)-test_k), test_index)
           gap_mm(i) = test_gap(test_index)
        ENDFOR
        
        k_name = strtrim(string(1000.*k(0), format = '(i4)'), 2)
        IF (strlen(k_name) EQ 2) THEN BEGIN
           k_name = '00'+k_name
        ENDIF ELSE IF (strlen(k_name) EQ 3) THEN BEGIN
           k_name = '0'+k_name
        ENDIF
        
        lowe_mapfile = 'x1_k'+k_name+'_onaxis.map'
        print, 'Reading "'+lowe_mapfile+'"'
        read_mapper, lowe_mapfile, lowe_lambda, lowe_inten, /quiet
        highe_mapfile = 'x1_k'+k_name+'_highe_onaxis.map'
        print, 'Reading "'+highe_mapfile+'"'
        read_mapper, highe_mapfile, highe_lambda, highe_inten, /quiet
        n_lambda = n_elements(highe_lambda)+n_elements(lowe_lambda)
        n_highe_lambda = n_elements(highe_lambda)
        lambda_nm = fltarr(n_lambda)
        lambda_nm(0:(n_highe_lambda-1)) = 0.1*highe_lambda
        lambda_nm(n_highe_lambda:(n_lambda-1)) = 0.1*lowe_lambda
        
        inten = fltarr(n_lambda, n_k)
        
        FOR i = 0, (n_k-1) DO BEGIN
           k_name = strtrim(string(1000.*k(i), format = '(i4)'), 2)
           IF (strlen(k_name) EQ 2) THEN BEGIN
              k_name = '00'+k_name
           ENDIF ELSE IF (strlen(k_name) EQ 3) THEN BEGIN
              k_name = '0'+k_name
           ENDIF
           
           lowe_mapfile = 'x1_k'+k_name+'_onaxis.map'
           print, 'Reading "'+lowe_mapfile+'"'
           read_mapper, lowe_mapfile, lowe_lambda, lowe_inten, /quiet
           highe_mapfile = 'x1_k'+k_name+'_highe_onaxis.map'
           print, 'Reading "'+highe_mapfile+'"'
           read_mapper, highe_mapfile, highe_lambda, highe_inten, /quiet
           inten(0:(n_highe_lambda-1), i) = highe_inten
           inten(n_highe_lambda:(n_lambda-1), i) = lowe_inten
        ENDFOR
        
        get_lun, lun
        openw, lun, xdr_filename, /xdr
        writeu, lun, fix(n_k), fix(n_lambda)
        writeu, lun, float(k), float(gap_mm), float(lambda_nm), float(inten)
        close, lun
        free_lun, lun
        
        print,'Built file "'+xdr_filename+'" with '+$
              strtrim(string(n_k, format = '(i5)'), 2)+' K values.'
        print,'Move the file to an appropriate directory'
        return
     ENDELSE
     
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ENDIF ELSE BEGIN
     ;; Read the data from file. 
     this_xdr_filename = file_which(xdr_filename)
     get_lun,lun
     openr, lun, this_xdr_filename, /xdr
     n_k = fix(0)
     n_lambda = fix(0)
     readu, lun, n_k, n_lambda
     k = fltarr(n_k)
     gap_mm = fltarr(n_k)
     lambda_nm = fltarr(n_lambda)
     inten = fltarr(n_lambda, n_k)
     readu, lun, k, gap_mm, lambda_nm, inten
     close, lun
     free_lun, lun
     
     IF keyword_set(nodisplay) THEN BEGIN
        ;; If no plot was asked for, we're done!
        return
     ENDIF ELSE IF keyword_set(plotgap) THEN BEGIN
        ;; We were asked for a simple plot
        IF keyword_set(epsfile) OR keyword_set(psfile) THEN BEGIN
           lambda_nm_title = '!Ml!X (nm)'
           old_plot = !D.name
           old_font = !P.font
           set_plot, 'ps'
           !P.font = 0
           xinches = 6.
           yinches = 4.
           
           ;; Did the user specify a file name?
           IF keyword_set(epsfile) THEN BEGIN
              svec = size(epsfile)
              IF (svec(1) EQ 7) THEN BEGIN
                 filename = epsfile
              ENDIF ELSE BEGIN
                 filename = 'x1undgap.eps'
              ENDELSE
              device, file = filename, /encap, /inch, $
                      xsize = xinches, ysize = yinches
           ENDIF ELSE BEGIN
              svec = size(psfile)
              IF (svec(1) EQ 7) THEN BEGIN
                 filename = psfile
              ENDIF ELSE BEGIN
                 filename = 'x1undgap.ps'
              ENDELSE
              device, file = filename, /inch, $
                      xsize = xinches, ysize = yinches, xoffset = 1., $
                      yoffset = 1.
           ENDELSE
           
           print, 'Writing file "'+filename+'" of size '+$
                  strtrim(string(xinches, format = '(f10.2)'), 2)+'x'+$
                  strtrim(string(yinches, format = '(f10.2)'), 2)+' inches.'
        ENDIF ELSE BEGIN
           lambda_nm_title = 'Wavelength (nm)'
        ENDELSE
        
        dummy = min(abs(gap_mm-plotgap), gap_index)
        title = 'For gap='+$
                strtrim(string(gap_mm(gap_index), format = '(f10.2)'), 2)+$
                ', k='+$
                strtrim(string(k(gap_index), format = '(f10.2)'), 2)
        IF keyword_set(ev) THEN BEGIN
           IF keyword_set(xrange) THEN BEGIN
              plot, 1239.85/lambda_nm, inten(*, gap_index), $
                    xtitle = 'Energy (eV)', $
                    ytitle = 'I (ph/mrad!E2!N/250 mA)@'+energy_text+' GeV', $
                    title = title, $
                    xrange = xrange, xstyle = 1
              
              ;; plot_edges
           ENDIF ELSE BEGIN
              print,'plotting'
              plot, 1239.85/lambda_nm, inten(*, gap_index), $
                    xtitle = 'Energy (eV)', $
                    ytitle = 'I (ph/mrad!E2!N/250 mA)@'+energy_text+' GeV', $
                    title = title
           ENDELSE
        ENDIF ELSE BEGIN
           IF keyword_set(xrange) THEN BEGIN
              plot, lambda_nm, inten(*, gap_index), $
                    xtitle = lambda_nm_title, $
                    ytitle = 'I (ph/mrad!E2!N/250 mA)@'+energy_text+' GeV', $
                    title = title, $
                    xrange = xrange, xstyle = 1
           ENDIF ELSE BEGIN
              plot, lambda_nm, inten(*, gap_index), $
                    xtitle = lambda_nm_title, $
                    ytitle = 'I (ph/mrad!E2!N/250 mA)@'+energy_text+' GeV', $
                    title = title
           ENDELSE
        ENDELSE
        
        IF keyword_set(epsfile) OR keyword_set(psfile) THEN BEGIN
           device, /close
           set_plot, old_plot
           !P.font = old_font
        ENDIF
        return
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ENDIF ELSE BEGIN
        ;; In this case we want the contour plot
        IF keyword_set(epsfile) OR keyword_set(psfile) THEN BEGIN
           ;; we could open an EPS file temporarily to get these
           ;; numbers if we were really picky...
           xcharinches = 0.0874016
           ycharinches = 0.138583
           thick = 1.5
           linecolor = 0
           lambda_nm_title = '!Ml!X (nm)'
        ENDIF ELSE BEGIN
           xcharinches = !D.x_ch_size/(2.54*!D.x_px_cm)
           ycharinches = !D.y_ch_size/(2.54*!D.y_px_cm)
           thick = 1
           linecolor = 255
           lambda_nm_title = 'lambda (nm)'
        ENDELSE
        energy_ev_title = 'Energy (eV)'
        k_title = 'K'
        gap_mm_title = 'Gap (mm)'
        
        bigcharsize = 1.3
        tickinches = 0.08
        ximinch = 4.
        yiminch = 4.
        xlborder = 2.*tickinches+5.*xcharinches+$
                   1.5*bigcharsize*ycharinches
        xrborder = xlborder
        ybborder = 2.*tickinches+ycharinches+$
                   1.5*bigcharsize*ycharinches
        ytborder = 2.*tickinches+ycharinches+$
                   2.*bigcharsize*ycharinches+$
                   1.5*bigcharsize*ycharinches
        xinches = ximinch+xlborder+xrborder
        yinches = yiminch+ybborder+ytborder
        
        IF keyword_set(epsfile) OR keyword_set(psfile) THEN BEGIN
           xpixels = fix(xinches*75.)
           ypixels = fix(yinches*75.)
           bits = 8
           old_plot = !D.name
           old_font = !P.font
           set_plot, 'ps'
           !P.font = 0
           
           ;; Did the user specify a file name?
           IF keyword_set(epsfile) THEN BEGIN
              svec = size(epsfile)
              IF (svec(1) EQ 7) THEN BEGIN
                 filename = epsfile
              ENDIF ELSE BEGIN
                 filename = 'x1und.eps'
              ENDELSE
              device, file = filename, /encap, bits = bits, /inch, $
                      xsize = xinches, ysize = yinches
           ENDIF ELSE BEGIN
              svec = size(psfile)
              IF (svec(1) EQ 7) THEN BEGIN
                 filename = psfile
              ENDIF ELSE BEGIN
                 filename = 'x1und.ps'
              ENDELSE
              device, file = filename, bits = bits, /inch, $
                      xsize = xinches, ysize = yinches, xoffset = 1., $
                      yoffset = 1.
           ENDELSE
           
           print, 'Writing file "'+filename+'" of size '+$
                  strtrim(string(xinches, format = '(f10.2)'), 2)+'x'+$
                  strtrim(string(yinches, format = '(f10.2)'), 2)+' inches.'
        ENDIF ELSE BEGIN
           xpixels = fix(xinches*2.54*!D.x_px_cm)
           ypixels = fix(yinches*2.54*!D.y_px_cm)
           print, 'Opening a '+$
                  strtrim(string(xpixels), 2)+'x'+$
                  strtrim(string(ypixels), 2)+' pixel screen as window 1.'
           window, 1, xsize = xpixels, ysize = ypixels
           wset, 1
        ENDELSE
        
        IF keyword_set(nocontour) THEN BEGIN
           image_maxinten = 1.
        ENDIF ELSE BEGIN
           image_maxinten = 1.5
        ENDELSE
        image_gamma = 0.3
        
        inten = inten
        energy_ev = 1239.85/lambda_nm
        
        ;; Because lambda is not on a straight linear scale,
        ;; we must remap it.
        lambda_nm_min = min(lambda_nm, max = lambda_nm_max)
        n_linear = xpixels
        lambda_linear = lambda_nm_min+(lambda_nm_max-lambda_nm_min)*$
                        findgen(n_linear)/float(n_linear-1)
        inten_linear = fltarr(n_linear, n_k)
        FOR i_k = 0, n_elements(k)-1 DO BEGIN
           inten_linear(0:(n_linear-1), i_k) = $
              interpol(inten(0:(n_lambda-1), i_k), lambda_nm, lambda_linear)
        ENDFOR
        
        byteimage = bytscl( (inten_linear/max(inten_linear))^image_gamma, $
                            min = 0., max = image_maxinten, $
                            top = !D.table_size-1 )
        IF keyword_set(epsfile) OR keyword_set(psfile) THEN BEGIN
           tv, rebin(byteimage, n_linear, 4*n_k), $
               xlborder, ybborder, $
               xsize = ximinch, ysize = yiminch, /inch
        ENDIF ELSE BEGIN
           xpix = fix(xpixels*ximinch/xinches)
           ypix = fix(ypixels*yiminch/yinches)
           IF (strmid(!Version.release, 0, 1) EQ '3') THEN BEGIN
              tv, congrid(byteimage, xpix, ypix, /interp), $
                  fix(xpixels*xlborder/xinches), $
                  fix(ypixels*ybborder/yinches), /device
           ENDIF ELSE BEGIN
              tv, congrid(byteimage, xpix, ypix, /cubic), $
                  fix(xpixels*xlborder/xinches), $
                  fix(ypixels*ybborder/yinches), /device
           ENDELSE
        ENDELSE
        
        IF (NOT keyword_set(nocontour)) THEN BEGIN
           c_levels = [0.5, 1., 2., 4., 6., 8., 10.]
           c_annotation = ['0.5', '1', '2', '4', '6', '8', '10']
           c_color = 255+bytarr(n_elements(c_levels))
           c_charsize = 1.
           contour, 1.e-15*rebin(inten, n_lambda, 4*n_k), $
                    lambda_nm, rebin(k, 4*n_k), $
                    xrange = [min(lambda_nm), max(lambda_nm)], xstyle = 5, $
                    yrange = [min(k), max(k)], ystyle = 5, $
                    levels = c_levels, c_annot = c_annotation, $
                    c_color = c_color, $
                    c_charsize = c_charsize, /noerase, color = c_color(0), $
                    position = [(xlborder/xinches), $
                                (ybborder/yinches), $
                                (xlborder+ximinch)/xinches, $
                                (ybborder+yiminch)/yinches], /norm
        ENDIF
        
        ;; Borders
        plots, (xlborder+[0., 1.]*ximinch)/xinches, $
               (ybborder+[0., 0.])/yinches, $
               /norm, color = linecolor, thick = thick
        plots, (xlborder+[0., 1.]*ximinch)/xinches, $
               (ybborder+yiminch+[0., 0.])/yinches, $
               /norm, color = linecolor, thick = thick
        plots, (xlborder+[0., 0.])/xinches, $
               (ybborder+[0., 1.]*yiminch)/yinches, $
               /norm, color = linecolor, thick = thick
        plots, (xlborder+ximinch+[0., 0.])/xinches, $
               (ybborder+[0., 1.]*yiminch)/yinches, $
               /norm, color = linecolor, thick = thick
        
        ;; Lambda axis
        x = lambda_nm
        xstep = 0.5
        xtickformat = '(f10.1)'
        xtitle = lambda_nm_title
        xrange = [x(0), x(n_elements(x)-1)]
        xtick_minindex = fix(0.02+min(xrange)/xstep)
        xtick_maxindex = fix(0.02+max(xrange)/xstep)
        xticks = (xtick_maxindex-xtick_minindex)+1
        xtickpos = fltarr(xticks)
        FOR i_xtick = xtick_minindex, xtick_maxindex DO BEGIN
           xtickpos(i_xtick-xtick_minindex) = float(i_xtick)*xstep
        ENDFOR
        xtickstr = strarr(xticks)
        tickdir = -1
        IF (tickdir EQ 1) THEN BEGIN
           ;; At top of image
           yoffsetinches = ybborder+yiminch
           ytextposinches = yoffsetinches+1.5*tickinches
           ytitleposinches = ytextposinches+1.5*ycharinches+$
                             0.5*bigcharsize*ycharinches
        ENDIF ELSE BEGIN
           ;; At bottom of image
           yoffsetinches = ybborder
           ytextposinches = yoffsetinches-1.5*tickinches-ycharinches
           ytitleposinches = ytextposinches-$
                             1.5*bigcharsize*ycharinches
        ENDELSE
        FOR i = 0, (xticks-1) DO BEGIN
           xdatanormpos = (xtickpos(i)-xrange(0))/(xrange(1)-xrange(0))
           IF (xdatanormpos GE -0.02) AND (xdatanormpos LE 1.02) THEN BEGIN
              xtickstr(i) = strtrim(string(xtickpos(i), $
                                           format = xtickformat), 2)
              plots, ([0, 0]+xlborder+ximinch*xdatanormpos)/xinches, $
                     (yoffsetinches+tickinches*tickdir*[0., 1.])/yinches, $
                     /norm, color = linecolor, thick = thick
              xyouts, (xlborder+ximinch*xdatanormpos)/xinches, $
                      (ytextposinches)/yinches, $
                      /norm, xtickstr(i), color = linecolor, align = 0.5
           ENDIF
        ENDFOR
        xyouts, (xlborder+0.5*ximinch)/xinches, $
                (ytitleposinches)/yinches, /norm, $
                xtitle, charsize = bigcharsize, color = linecolor, align = 0.5
        
        ;; energy axis
        x = energy_ev
        xtickformat = '(i4)'
        xtitle = energy_ev_title
        xtickpos = [250., 300., 350., 400., 500., 600., 750., 1000., 1500.]
        xticks = n_elements(xtickpos)
        xtickstr = strarr(xticks)
        tickdir = 1
        IF (tickdir EQ 1) THEN BEGIN
           ;; At top of image
           yoffsetinches = ybborder+yiminch
           ytextposinches = yoffsetinches+1.5*tickinches
           ytitleposinches = ytextposinches+ycharinches+$
                             0.5*bigcharsize*ycharinches
        ENDIF ELSE BEGIN
           ;; At bottom of image
           yoffsetinches = ybborder
           ytextposinches = yoffsetinches-1.5*tickinches-ycharinches
           ytitleposinches = ytextposinches-$
                             1.5*bigcharsize*ycharinches
        ENDELSE
        FOR i = 0, (xticks-1) DO BEGIN
           this_lambda_nm = 1239.85/xtickpos(i)
           xdatanormpos = (this_lambda_nm-lambda_nm_min)/$
                          (lambda_nm_max-lambda_nm_min)
           IF (xdatanormpos GE -0.02) AND (xdatanormpos LE 1.02) THEN BEGIN
              xtickstr(i) = strtrim(string(xtickpos(i), $
                                           format = xtickformat), 2)
              plots, ([0, 0]+xlborder+ximinch*xdatanormpos)/xinches, $
                     (yoffsetinches+tickinches*tickdir*[0., 1.])/yinches, $
                     /norm, color = linecolor, thick = thick
              xyouts, (xlborder+ximinch*xdatanormpos)/xinches, $
                      (ytextposinches)/yinches, $
                      /norm, xtickstr(i), color = linecolor, align = 0.5
           ENDIF
        ENDFOR
        xyouts, (xlborder+0.5*ximinch)/xinches, $
                (ytitleposinches)/yinches, /norm, $
                xtitle, charsize = bigcharsize, color = linecolor, align = 0.5
        
        ;; K axis
        y = k
        ystep = 0.2
        ytickformat = '(f10.1)'
        ytitle = k_title
        yrange = [y(0), y(n_elements(y)-1)]
        ytick_minindex = fix(0.02+min(yrange)/ystep)
        ytick_maxindex = fix(0.02+max(yrange)/ystep)
        yticks = (ytick_maxindex-ytick_minindex)+1
        ytickpos = fltarr(yticks)
        FOR i_ytick = ytick_minindex, ytick_maxindex DO BEGIN
           ytickpos(i_ytick-ytick_minindex) = float(i_ytick)*ystep
        ENDFOR
        ytickstr = strarr(yticks)
        tickdir = -1
        IF (tickdir EQ 1) THEN BEGIN
           ;; At right of image
           xoffsetinches = xlborder+ximinch
           xtextposinches = xoffsetinches+1.5*tickinches
           xalign = 0.
           xorient = 0.
           xtitleposinches = xtextposinches+5.*xcharinches+$
                             0.5*bigcharsize*ycharinches
        ENDIF ELSE BEGIN
           ;; At left of image
           xoffsetinches = xlborder
           xtextposinches = xoffsetinches-1.5*tickinches
           xalign = 1.
           xorient = 0.
           xtitleposinches = xtextposinches-5.*xcharinches-$
                             0.5*bigcharsize*ycharinches
        ENDELSE
        FOR i = 0, (yticks-1) DO BEGIN
           ydatanormpos = (ytickpos(i)-yrange(0))/(yrange(1)-yrange(0))
           IF (ydatanormpos GE -0.02) AND (ydatanormpos LE 1.02) THEN BEGIN
              ytickstr(i) = strtrim(string(ytickpos(i), $
                                           format = ytickformat), 2)
              plots, (xoffsetinches+tickinches*tickdir*[0., 1.])/xinches, $
                     ([0., 0.]+ybborder+yiminch*ydatanormpos)/yinches, $
                     /norm, color = linecolor, thick = thick
              xyouts, (xtextposinches)/xinches, $
                      (ybborder+yiminch*ydatanormpos-0.5*ycharinches)/yinches, $
                      /norm, ytickstr(i), color = linecolor, align = xalign, $
                      orientation = xorient
           ENDIF
        ENDFOR
        xyouts, (xtitleposinches)/xinches, $
                (ybborder+0.5*yiminch)/yinches, /norm, $
                ytitle, charsize = bigcharsize, color = linecolor, $
                align = 0.5, orientation = 90
        
        ;; gap axis
        y = gap_mm
        ystep = 2.
        ytickformat = '(f10.1)'
        ytitle = gap_mm_title
        yrange = [y(0), y(n_elements(y)-1)]
        ytick_minindex = fix(0.02+min(yrange)/ystep)
        ytick_maxindex = fix(0.02+max(yrange)/ystep)
        yticks = (ytick_maxindex-ytick_minindex)+1
        ytickpos = fltarr(yticks)
        FOR i_ytick = ytick_minindex, ytick_maxindex DO BEGIN
           ytickpos(i_ytick-ytick_minindex) = float(i_ytick)*ystep
        ENDFOR
        ytickstr = strarr(yticks)
        tickdir = 1
        IF (tickdir EQ 1) THEN BEGIN
           ;; At right of image
           xoffsetinches = xlborder+ximinch
           xtextposinches = xoffsetinches+1.5*tickinches
           xalign = 0.
           xorient = 0.
           xtitleposinches = xtextposinches+5.*xcharinches+$
                             0.5*bigcharsize*ycharinches
        ENDIF ELSE BEGIN
           ;; At left of image
           xoffsetinches = xlborder
           xtextposinches = xoffsetinches-1.5*tickinches
           xalign = 1.
           xorient = 0.
           xtitleposinches = xtextposinches-5.*xcharinches-$
                             0.5*bigcharsize*ycharinches
        ENDELSE
        min_k = min(k, max = max_k)
        FOR i = 0, (yticks-1) DO BEGIN
           this_b0 = B_REMNANT/sinh(!Pi*1.e-3*ytickpos(i)/LAMBDA0)
           this_k = C*LAMBDA0*this_b0/(EMASS*2.0*!Pi)
           ydatanormpos = (this_k-min_k)/(max_k-min_k)
           IF (ydatanormpos GE -0.02) AND (ydatanormpos LE 1.02) THEN BEGIN
              ytickstr(i) = strtrim(string(ytickpos(i), $
                                           format = ytickformat), 2)
              plots, (xoffsetinches+tickinches*tickdir*[0., 1.])/xinches, $
                     ([0., 0.]+ybborder+yiminch*ydatanormpos)/yinches, $
                     /norm, color = linecolor, thick = thick
              xyouts, (xtextposinches)/xinches, $
                      (ybborder+yiminch*ydatanormpos-0.5*ycharinches)/yinches, $
                      /norm, ytickstr(i), color = linecolor, align = xalign, $
                      orientation = xorient
           ENDIF
        ENDFOR
        xyouts, (xtitleposinches)/xinches, $
                (ybborder+0.5*yiminch)/yinches, /norm, $
                ytitle, charsize = bigcharsize, color = linecolor, $
                align = 0.5, orientation = -90.
        
        title = 'X1 undulator intensity (10!E15!N/250 mA)@'+energy_text+' GeV'
        xyouts, $
           (xlborder+0.5*ximinch)/xinches, $
           (yinches-0.5*tickinches-1.5*bigcharsize*ycharinches)/yinches, $
           /norm, title, color = linecolor, align = 0.5, $
           charsize = bigcharsize
        
        IF keyword_set(epsfile) OR keyword_set(psfile) THEN BEGIN
           device, /close
           set_plot, old_plot
           !P.font = old_font
        ENDIF
     ENDELSE
     
  ENDELSE

END
