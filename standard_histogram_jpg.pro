pro standard_histogram_jpg_placeholder,indir=indir,outdir=outdir,region=region,x_x=x_x,y_y=y_y,peak_location=peak_location,inputVals=inputVals,inputHist=inputHist,prob_exp=prob_exp,yfit=yfit,param=param
	time_of_start=systime()
	; Set directories from what is created in reduction code
	if not (keyword_set(indir)) then begin
		indir = 'H:\QIS1-1\conversiongain_sp\conversiongain_sp.13Jan21_4\row_43_col_19\'
	endif
	if not (keyword_set(clusindir)) then begin
		clusindir =  indir       ;"H:\QIS1-1\conversiongain_sp\set1\"
	endif
	if not (keyword_set(outdir)) then begin
		outdir = indir + 'output\' ;test the new code on Oct 31
	endif
	; Sets the files to be used by the program for analysis
	if not (keyword_set(list)) then begin
		list = 'conversiongain_sp.lst\'
	endif
	if not (keyword_set(cluslist)) then begin
		cluslist = 'conversiongain_sp.lst\'
	endif
	if not (keyword_set(region)) then begin
		region = [4,1019,4,1019]
	endif

	histdir=outdir+'histograms\'
	mapdir=outdir+'maps\'

	conversiongain = param[0]
	RN = param[1]
	H = param[2]
	offset = param[3]
	chisq = param[4]
	est_rn = param[5]


;	################### Plotting #################### Plotting #################### PLotting #################;


				; Ploting the average readnoise vs integration time
    set_plot, 'z'
    ; set resolution and font
    device, set_resolution=[8000,6000]
    device, set_font='Courier'
    ; Setting Colors
    white=255
    black=0
    !P.MULTI = 0    ; Set up for one plot per page
    jpgfile = outdir + 'plots\'
    fsub = jpgfile  ; set text to print at bottom of plot
    ; set plot margin & tick marks
    !Y.MINOR=0
    !y.margin=[6,4]
    !x.margin=[16,4]
    ; set some size values
    asize = 8.0     ; charsize of annotations
    athick = 15.0   ; charthick of annotations
    lthick = 20.0   ; line thickness
    ; Define a special plotting symbol to plot big, fat points
    _radius = asize  ; Offset radius is about size of a character
    _X = fltarr(10)  ; Place for x verticies, offset from point
    _Y = fltarr(10)  ; Place for y verticies, offset from point
    _theta = 0
    _delta_theta = 2 * !PI / 10
    for i=0, 9 do begin
        _X[i] = _radius * cos(_theta)
        _Y[i] = _radius * sin(_theta)
        _theta = _theta + _delta_theta
    endfor
    ; set the symbol created above as IDL symbol #8
    usersym, _X, _Y, /FILL


    xtitle='ADU'
    ytitle='Probability'
    ; Plot the result

	if (keyword_set(prob_exp)) then begin
		plot, inputVals, prob_exp, TITLE='Single Pixel Histogram', XTITLE=xtitle, $
     		YTITLE=ytitle, thick=15, XRANGE=[MIN(inputVals)*0.975,MAX(inputVals)*1.1], XSTYLE=1, YRANGE=[MIN(prob_exp)*0.975,MAX(prob_exp)*1.025], YSTYLE=1, BACKGROUND=white, $
     		charsize=asize*1.25, charthick=athick, xthick=lthick, ythick=lthick, COLOR=black
	endif


	if not (keyword_set(prob_exp)) then begin
		plot, inputVals, inputHist, TITLE='Single Pixel Histogram', XTITLE=xtitle, $
     		YTITLE=ytitle, thick=15, XRANGE=[MIN(inputVals)*0.975,MAX(inputVals)*1.1], XSTYLE=1, YRANGE=[MIN(inputHist)*0.975,MAX(inputHist)*1.025], YSTYLE=1, BACKGROUND=white, $
     		charsize=asize*1.25, charthick=athick, xthick=lthick, ythick=lthick, COLOR=black
	endif

	if (keyword_set(peak_location)) then begin
		oplot,inputVals,yfit,linestyle=2,thick=10,COLOR=black
		for j =0, n_elements(peak_location)-1 do begin
			oplot, fltarr(n_elements(inputHist))+peak_location[j], inputHist/max(inputHist), thick=5, color=black
			;wait, 0.25
		endfor
	endif


    timeStamp= systime()
    legenditems=[$
                'Detector: '+'QIS1-1',$
                'Analysis Region: '+'Single Pixel',$
                'Position: '+ strtrim(string(x_x),2) + '_' + strtrim(string(y_y),2), $
                'Date: '+ strtrim(string(timeStamp),2)+' UT',$
                'Detector Temperature: '+strtrim(string(297,format='(G10.4)'),2)+' K',$
                'Conversion Gain: ' +strtrim(string(conversiongain, format='(G10.4)'),2)+' ADU/e!E-',$
                'Read Noise: ' + strtrim(string(RN,format='(G10.4)'),2)+ ' e!E-!N rms', $
                'Est Read Noise: ' + strtrim(string(est_rn,format='(G10.4)'),2)+ ' e!E-!N rms', $
                'Quanta Exposure: ' + strtrim(string(H,format='(G10.4)'),2)+ ' e!E-', $
                'Chi-Squared: ' + strtrim(string(chisq,format='(G10.4)'),2),$
                'Offset: ' +strtrim(string(offset,format='(G10.4)'),2)+ 'e!E-']

	;stop
    ; make legend
    legend,legenditems,/top,/left,box=0,textcolors=replicate(0,n_elements(legenditems)),charsize=8.,charthick=12.
    ; write the file name outside the margin of the plot
    xyouts, 1.0, 40, fsub, charsize = asize, charthick = athick, $
    color = black, /device
    ; For the JPEG plot, the plot image must be read in from memory and
    ; written out to a file
    jpgfile = outdir + 'plots\'
    time_of_end=systime()
    x_frame = x_x + region[0]
    y_frame = y_y + region[2]
    jpg_name = 'Pixel_'+ strtrim(string(x_frame),2)+'_'+ strtrim(string(y_frame),2)+'.jpg'
    jpgfile = jpgfile + jpg_name
    jpgimg = tvrd()
    write_jpeg, jpgfile, congrid(jpgimg, 1600, 1200, /center, /interp), quality=100
end
