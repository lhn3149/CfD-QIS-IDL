; NAME:
;      makehistogram.pro (v1.3)
;
; LOCATION:
;      \\wolf\RIDL\internal\software\IDL\makehistogram.pro
;
; PURPOSE:
;      This procedure reads from a list of fits images (containing a single 2D plane) and generates a histogram
;      for its values.
;
; CALLING SEQUENCE:
;      makehistogram, testdir, infile, intitle, inxtitle [, region [, outdir]]
;
; INPUTS:
;      testdir  - the full path name of the directory containing an input data file
;      infile   - name of the file containing a list of input images
;      intitle  - title of the plot--the word "Histogram" is appended to this string
;      inxtitle	- x-axis title for the plot
;
; KEYWORD PARAMETERS:
;      outdir - the full path name of a directory in which the output plots will be written.  If unspecified,
;               plots will be placed in the directory indicated by the testdir input variable
;      region - 4-vector specifying the outer boundaries of the analysis region to be considered when making
;               the histogram
;		ADU   - plot in ADUs
; OUTPUT:
;
;
; EXAMPLE
;  IDL> makehistogram, '.', 'slopefiles.txt', 'Dark Current', 'Dark Current (ADU/sec)', region=[4,511,4,2043]
;
;
; REFERENCE:
;
;
; MODIFICATION HISTORY:
;      Created by:  Ernie Morse, IDTL, April 29, 2003 (v1.1)
;
;      Modified by: Ernie Morse, IDTL, May 14, 2003 (v1.2)
;                    Changed # of significant digits on plot output notation
;                    Corrected creation of jpgfile and psfile so that the base output filenames are the input file names up until the
;                    last occurrence of '.'
;
;                   Ernie Morse, IDTL, May 22, 2003 (v1.3)
;                    Increased width of initial histogram to median +/- 1.5 SD instead of median +/- 0.5 SD.  The old limits don't work to well
;                    if the distribution is not very Gaussian.
;
;                   Ernie Morse, IDTL, August 12, 2004 (v1.4)
;                    Added colinc, rowinc keywords
;
;                   Nick Cox, RIDL, January 5, 2009
;			         Added header information
;
;					Joong Lee, RIDL, August 2, 2017
;					 changed code to plot in electrons unless spcified to plot in ADUs explicitly
;
;					Valerie Fleischauer, June 3, 2019
;					 Added option for QIS-type detector conversion gain data reduction with keyword /QIScg
;
;					Don Figer, June 6, 2019
;					 Removed one instance of two in which the conversion gain was applied when e was not set and /ADU was not set
;					 This change assumes that the units of the input data for CDS are electrons when those keywords are not set
;----------------------------------------------------------------------------------------------------------------------------------------------
Pro makehistogram_conversiongain, testdir, infile, intitle, inxtitle, region=region, rowinc=rowinc, colinc=colinc, outdir=outdir, e=e, ADU=ADU, QIScg=QIScg

	; set directory according to whether using a PC or workstation
	pathdelim=path_sep()

	; begin by setting defaults for keyword values not provided by the calling context
	if (not keyword_set(outdir)) then outdir=testdir

	; check to see if input and outdirectory names end in path delimiter--if not, add it to the end
	if (strmid(testdir, strlen(testdir)-1, 1) ne pathdelim) then testdir=testdir+pathdelim
	if (strmid(outdir, strlen(outdir)-1, 1) ne pathdelim) then outdir=outdir+pathdelim

    if (not keyword_set(rowinc)) then begin
        rowinc=1
    endif
    if (not keyword_set(colinc)) then begin
        colinc=1
    endif
	xlength = (region[1]-region[0]+1)/10
	ylength = (region[3]-region[2]+1)/10
	total_pixel = xlength * ylength ; per 100

	; read all of the input image names into an array
	if file_test(testdir+infile) eq 0 then return
	readcol, testdir+infile, inputimgs, format='(A)'

	; loop through all of the images
	for imgcount=0, n_elements(inputimgs)-1 do begin

		; read the input image and its header
		fits_read, testdir+inputimgs[imgcount], inimage, imghead, /noscale

		; if the keyword "e" is set, that implies that the input data are in ADU and need to be converted into electrons
		if keyword_set(e) then inimage=inimage*e else begin

			; the following code will ensure that the output is in electrons as a default unless /ADU is specified
			if not keyword_set(ADU) then begin
				; electronics name
				electronics = strcompress(sxpar(imghead, 'DAQ_ELEC',count=count_elec),/remove_all)
				if strpos(electronics,'LEACH') ne -1 then electronics='LEACH'
				; detector analysis profile
				detname = sxpar(imghead, 'detname')
				analysisprof='QIS1-1-QISPF'+'.cfg'

				searchresults=file_search('C:\RIDL\software\DetectorAnalysisProfile\'+analysisprof,count=profcount)
				; if analysis profile exists

				if profcount gt 0 then begin
					egainexp = sxpar(imghead, 'E_GAIN',count=count_egain)
					if count_egain gt 0 then begin
						CDSpos=strpos(inxtitle,'CDS Noise')
						Darkpos=strpos(inxtitle,'Dark Current')

						; to get to this point in the code, neither e nor /ADU has been set
						; the code assumes that in this case that if the data are CDS noise data, they are alredy in electrons
						; this will be a false assumption in the case that the input data are in ADU but the user did not specify /ADU and did not set e,
						; as can be the case when calling this routine from within batch_darkcurrent in an autoreduce script
						; in this case, the x-axis will say electrons but the units will really be ADU
						if CDSpos eq 0 then inxtitle='CDS Noise (e!E-!N)'

						; the following applies the conversion gain in the case that e was not specified and ADU output was not chosen
						; this only applies to dark current data because the code assumes that CDS noise data are already in electrons
						if Darkpos eq 0 then begin
							inxtitle='Dark Current (e!E-!N/s)'
							determineanalysisparameters, '\\hawk\C\RIDL\internal\software\DetectorAnalysisProfile\', analysisprof, ConvGainHigh,$
		 						ConvGainLow, nonlinearitycorrfile, eGainHigh, eGainLow,pixelpitch, lowgain=lowgain
							inimage=inimage*ConvGainHigh*egainexp/eGainHigh
						endif
					endif
				endif
			endif
		endelse

		; IDTL FITS headers are non-standard, so the following function calls are necessary to extract
		; time and date info from the header date and time keywords
;		timeStamp='' 		; An empty string to hold the result
;		dateStamp='' 		; Likewise
;		gettime, imghead, timeStamp
;		getdate, imghead, dateStamp
		timestamp=sxpar(imghead,'TIME')
		timestamp=' '+timestamp
		MONTHS=['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG', 'SEP','OCT','NOV','DEC']
		MONTH=long(where(strupcase(strmid(timestamp,5,3)) eq MONTHS))+1
		day=strmid(timestamp,9,2)
		year=strmid(timestamp,21,4)
		plotdate=string(month)+'/'+string(day)+'/'+string(year)
		timestamp=sxpar(imghead,'TIME')
		dettemp=sxpar(imghead, 'DETTEMP')

		; make the standard labels from the header
        getdetname,imghead,scaname
	    ; extract SCA type from SCA name if there is a dash in the name
	    scaName = "QIS-1'
	    dummy=strpos(scaName, '-')
	    if dummy ne -1 then scatype=strmid(scaName, 0, strpos(scaName, '-'))
        ; set sca_type to the MOSISN designation in the case that the detname is for a MOSIS device
        if strmid(scaname,0,3) eq 'MOS' then scatype=strmid(scaname,0,6)
		; MOSIS3 headers use DETTEMP instead of MOLYTABL for detector temperature
		if dettemp eq 0 and scatype eq 'MOSIS3' then dettemp=float(sxpar(imghead, 'dettemp'))
		if dettemp eq 0 and scatype eq 'MOSIS4' then dettemp=float(sxpar(imghead, 'dettemp'))
		if dettemp eq 0 and scatype eq 'MOSIS5' then dettemp=float(sxpar(imghead, 'dettemp'))

		; define default regions based on detector type
		scaType=(strsplit(scaName, '-', /extract))[0]

		if (not keyword_set(region)) then begin
			case scaType OF
	            'HRG' : region=[132, 4091, 4, 4091]
	            'H4RG' : region=[4, 4091, 4, 4091]
				'H1RG' : region=[4, 1019, 4, 1019]
				'H2RG' : region=[4, 2043, 4, 2043]
				'SB304': region=[4, 2051, 0, 2047]
                'MOSIS3': region=[1,126,1,126]
                'VIRGO' : region=[0,2047,1,2048]    ;jylpci 9/13/2012
                'QIS' : region=[4, 204, 4, 204]	;vefcfd 6/3/2019
			ELSE   : begin
				print, 'Detector type '+scaType+' not supported.  Exiting procedure.'
				return
						 end
			endcase
		endif

		; define the output plot names
		jpgfile=outdir+strmid(inputimgs[imgcount],0, strpos(inputimgs[imgcount], '.', /reverse_search))+'_hist.jpg'
		psfile=outdir+strmid(inputimgs[imgcount],0, strpos(inputimgs[imgcount], '.', /reverse_search))+'_hist.ps'

		; define forground and background colors
		white=255
		black=0

		ptype=['ps','jpg']
		fg=0
		;imgregion=inimage[region[0]:region[1]:colinc, region[2]:region[3]:rowinc]
		imgregion=inimage[0:region[1]-region[0]-1:colinc, 0:region[3]-region[2]-1:rowinc]
		theMedian=median(imgregion)
		valid_pixels=where(finite(imgregion),nvalid_pixels)

		sdev=stddev(imgregion, /nan)

		; Extend the max/min values for QIS detector type conversion gain analysis - vef 6/3/2019
		if keyword_set(QIScg) then begin
			_min=(theMedian - (sdev*3))
			_max=(theMedian + (sdev*3))
		endif else begin
			_min=(theMedian - (sdev*3))
			_max=(theMedian+(sdev*3))
		endelse

		_binsize=(_max - _min)/128.0

		if _binsize le 0 then _binsize=1.
		if _max lt _min then _max=_min+1

		theHist=histogram(imgregion, binsize=_binsize, min=_min, max=_max, /nan)
		theVals=findgen(n_elements(theHist))* _binsize+_min+_binsize/2.0

		; For QIS type detector conversion gain experiment, smooth the histogrammed data
		; and produce an output text file for analysis instead of a plot. - vef 6/3/2019
		if keyword_set(QIScg) then begin
			smHist = DOUBLE(SMOOTH(theHist, 4))
			vals = transpose([[theHist], [theVals], [smHist]])
			CD, outdir
			valsfile = 'cghistvals.txt'
			openw, 1, valsfile
			printf, 1, 'theHist			theVals			smHist'
			printf, 1, vals
			close, 1
			GOTO, endQIScg

		endif

		; Fit a Gaussian to the binned data
		if (n_elements(thehist) gt 1) then begin
			gaussparms=gaussfit(theVals, theHist, A,nterms=3)
			;sdev=A[2]
			;_min=(A[1] - (sdev*10))
			;_max=(A[1] + (sdev*10))
			_binsize=(_max - _min)/100.0  ; 256
			if _binsize le 0 then _binsize=1.
			if _max lt _min then _max=_min+1
			theHist=histogram(imgregion, binsize=_binsize , min=_min, max=_max, /nan)
			theVals=findgen(n_elements(theHist))* _binsize+_min+_binsize/2.0

			; calculate cumulative counts
			cumulative=fltarr(n_elements(theHist))
			for j=0,n_elements(theHist)-1 do begin
				cumulative(j)=float(total(theHist(0:j)))/nvalid_pixels
			endfor
			percentile99=where(cumulative ge .99)
			if percentile99[0] ne -1 then begin

				theHist=histogram(imgregion, binsize=_binsize , min=_min, max=thevals[percentile99[0]]+_binsize/2.0, /nan)
				theVals=findgen(n_elements(theHist))* _binsize+_min+_binsize/2.0
			endif

			; loop through plot types
			; this has been changed to select only JPG. change back to i=0,1 to get PS and JPG files
			for i=1, 1 do Begin
				plottype=ptype[i]
				!x.minor=10
				!y.minor=10
				if (plottype eq 'ps') then begin
					;;set up for landscape Postscript plotting
					set_plot, 'ps'
					!p.font=-1
					device, filename= psfile, /landscape
					fsub=psfile           ;filename to be written at bottom of plot

					;; set line and text characteristics for publication quality
					!P.CHARSIZE=1.2
					!P.CHARTHICK=4.
					!X.THICK=8.
					!Y.THICK=8.
					!P.THICK=8.
					!P.MULTI=0
					!P.TICKLEN=0.03
					symsize=1.5
					filenamecharsize=1.1

				endif

				if (plottype eq 'jpg') then begin
					;set up to plot to memory buffer.  This plot will later be written from memory to a JPEG file.
					set_plot, 'z'
					device, set_resolution=[8000,6000]
					device, set_font='Courier'
					fsub=jpgfile          ;filename to be written at bottom of plot

					; set line and text characteristics for publication quality
					!P.CHARSIZE=10
					!P.CHARTHICK=15.
					!X.THICK=20.
					!Y.THICK=20.
					!P.THICK=20.
					symsize=10.
					filenamecharsize=8.

				endif

				; draw plot
				erase
				plot,theVals,theHist,psym=10,title=intitle,xtitle=inxtitle,charsize=1.25*!P.CHARSIZE,ytitle='N', $
					color=black,background=white,position=[0.15,0.15,0.8,0.9],ystyle=8,$
					;xrange=[(A[1] - (sdev*7)),(A[1]+(sdev*7))],xstyle=1, /ylog, yrange=[1,max(theHist)]
					xrange=[_min,_max],xstyle=1
;					xrange=[1e-3,_max],xstyle=1;,/xlog
				; add vertical line through zero
				;oplot,[0,0],[!Y.CRANGE(0),!Y.CRANGE(1)],color=black

				; add file name to bottom of plot
				xyouts, 1.0, 40.0, fsub, color=0, charsize=filenamecharsize, /device

				; overplot labels
				; Add SCA, date & time stamps, read noise and gain to plot


				timeStamp= systime()


				legenditems=[$
					'Detector: '+scaName,$
					'Analysis Region: '+'['+strtrim(string(region[0]),2)+':'+strtrim(string(region[1]),2)+', '+ $
						strtrim(string(region[2]),2)+':'+strtrim(string(region[3]),2)+']',$
						'Date: '+ strtrim(string(timeStamp),2)+' UT',$
					;'Date: '+timeStamp+' UT',$
					;'Detector Temperature: '+strtrim(string(dettemp,format='(G10.4)'),2)+' K',$
					'X!ICenter!N='+strtrim(string(A[1],format='(G10.4)'),2),$
					'sigma='+strtrim(string(A[2],format='(G10.4)'),2), $
					'Operability: ' + strtrim(string(double(nvalid_pixels)/double(total_pixel),format='(G10.4)'),2)+'%']

					;'Refsub Version: '+strtrim(string(sxpar(imghead, 'refver'), format='(F4.1)'),2),$
					;'Refsub Type: '+strtrim(sxpar(imghead, 'reftype'),2),$
				; make legend
				legend,legenditems,/top,/left,box=0,textcolors=[0,0,0,0,0,0,0,0],charsize=8.,charthick=12.


				; plot vertical line through median value
				xpeak=[A[1],A[1]]
				ypeak=[0,!Y.CRANGE(1)]
				oplot,xpeak,ypeak,color=black,thick=!P.THICK*0.7

				; draw curve for cumulative plot and second y-axis
				axis,yaxis=1,yrange=[0,100],/save,ytitle='Cumulative %',charsize=1.25*!P.CHARSIZE,ytickv=indgen(11)*10.,yticks=10,yminor=5,color=black

				; calculate cumulative counts
				cumulative=fltarr(n_elements(theHist))
				for j=0,n_elements(theHist)-1 do cumulative(j)=float(total(theHist(0:j)))/nvalid_pixels

				; overplot cumulative counts
				oplot,theVals,cumulative*100.,linestyle=0,color=black,thick=!P.THICK*0.8

				;; write JPG image file
				if (plottype eq 'jpg') then begin
					jpgimg=tvrd()
					write_jpeg, jpgfile, congrid(jpgimg, 1280, 1024, /interp, /center), quality=100
				endif

				device, /close
			endfor
		endif
		endQIScg:

	endfor	;end imgcount loop

End


