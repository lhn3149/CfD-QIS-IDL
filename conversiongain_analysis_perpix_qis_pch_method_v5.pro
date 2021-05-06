;======================================================================================
	;================This is a prototype of the conversiongain analysis procedure================================
	;======================================================================================


pro conversiongain_analysis_perpix_qis_pch_method_v5,indir=indir,outdir=outdir,list=list,region=region,jpgfile=jpgfile,histdir=histdir,mapdir=mapdir,test_pix_flag=test_pix_flag

	time_of_start=systime()

	if not (keyword_set(indir)) then begin
		indir = "H:\QIS1-1\warm1\QISPF\conversiongain\conversiongain.03Mar21_6\"
	endif


	if not (keyword_set(outdir)) then begin
		outdir = indir + 'results\'
	endif

	; Sets the files to be used by the program for analysis
	if not (keyword_set(list)) then begin
		list = 'conversiongain.lst'
	endif

	list_rn = 'readnoise.lst'
	list_h = 'quantaexposure.lst'
	list_os = 'offset.lst'
	list_est_rn = 'est_rn.lst'

	if not (keyword_set(list)) then begin
		framelist = list
	endif

	if not (keyword_set(region)) then begin
		;region = [100,125,100,125]
		region = [350,450,350,450]
	endif

	if not (keyword_set(jpgfile)) then begin
		jpgfile =  outdir + 'plots\'
	endif

	if not (keyword_set(histdir)) then begin
		histdir=outdir+'histograms\'
	endif

	if not (keyword_set(mapdir)) then begin
		mapdir=outdir+'maps\'
	endif

	;TESTING: set flag to 1 to select a specific pixel [x,y]
	; [x,y] are coordinates referenced to the region (0,0), not the orginal pixel location
	if keyword_set(test_pix_flag) then begin
		test_pix_flag = 1
		;folder to write all plots of selected pixel to
		test_outdir=indir+'\results\plots\testing\'

		FILE_MKDIR,test_outdir

		; These are default set x and y locations of interesting pixel which can loop over if you comment out the stop found below. Feel free to change this to your hearts content.
		test_pix_all_x=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,464,464,464,464,464,464,464,464,464]
		test_pix_all_y=[8,38,214,14,28,40,59,131,140,143,146,302,305,306,375,126,162,158,165,192,191,198,197,293,297,310,330,336,338,340,346,348,351,352,357]

	endif else begin
		test_pix_flag = 0
	endelse

	; read input frames file names
	readcol,indir+list,inputframe,format='(A)'

	; create final data arrays
	x_length = region[1]-region[0]
	y_length = region[3]-region[2]
	finalarray_cg_f = make_array(x_length+1,y_length+1)
	finalarray_rn_f = make_array(x_length+1,y_length+1)
	finalarray_H_f = make_array(x_length+1,y_length+1)
	finalarray_os_f = make_array(x_length+1,y_length+1)
	estimate_rn_f = make_array(x_length+1,y_length+1)

	; reading header to determine size and information from raw datacubes
	time_at_start_of_read=systime()
	FITS_READ,indir+inputframe[0],framedata,header,/noscale,/header_only
	NAXIS1=FIX(STRCOMPRESS(sxpar(header,'NAXIS1'),/REMOVE_ALL))
	NAXIS2=FIX(STRCOMPRESS(sxpar(header,'NAXIS2'),/REMOVE_ALL))
	NAXIS3=FIX(STRCOMPRESS(sxpar(header,'NAXIS3'),/REMOVE_ALL))
	xmax=NAXIS1-1
	ymax=NAXIS2-1
	zmax=NAXIS3-1
	number_of_datacubes=long(n_elements(inputframe))
	frames_per_cube=long(NAXIS3)
	pix_data=make_array(x_length+1,y_length+1,number_of_datacubes*frames_per_cube)

	; read in all data as extract data within region
	time_at_start_of_read=systime()
	for framecount=0,n_elements(inputframe)-1 do begin
		FITS_READ,indir+inputframe[framecount],framedata,header,/noscale
		start=long(framecount)*long(frames_per_cube)
		stop=(long(framecount+1)*long(frames_per_cube))-1L
		pix_data[*,*,start:stop]=float(framedata[region[0]:region[1],region[2]:region[3],*])
	endfor

	time_of_concatenation=systime()
	count_nan=0
	count_high_CG=0
	count_std_CG=0
	count_low_CG=0
	count_cg_nan=0
	count_nan_worst=0
	count_H=0
	count_RN=0
	count_fit=0
	x_range=x_length+1
	y_range=y_length+1

	if test_pix_flag eq 1 then begin
		; Loop location, loop counter, and extracting specific pixel to test (test_pix[x,y])
		test_loop_num=-1
		test_loop:
		test_loop_num=test_loop_num+1
		test_pix=[test_pix_all_x[test_loop_num],test_pix_all_y[test_loop_num]]

		stop ;at this stop, you may now select another pixel: enter "test_pix= [x,y]" into the IDL command line and then continue (.c)
	endif




	; loop over each pixel in the region datacube
	for x_x= 0, x_range-1  do begin ;ymaxclus
		for y_y= 0, y_range-1 do begin	;xmaxclus

			; forcing pixel to be a specific pixel for testing purposes
			if test_pix_flag eq 1 then begin
				x_x=test_pix[0]
				y_y=test_pix[1]
			endif

			; Extract pixel data and make float for Histogram to use the input binsize as a float
			outputcg = 0.0
			H = 0.0
			RN = 0.0
			offset = 0.0
			chisq=  0.0
			un_est = 0.0

			;take all values of a pixel at x,y, array 1x1x(all frames)
			pix_data_single= float(pix_data[x_x,y_y,*])

			; computing paramters for histogram max, min, and binsize
			theMedian=median(pix_data_single)
			sdev=stddev(pix_data_single, /nan)
			_min=floor((theMedian - (sdev*4)))
			_max=ceil((theMedian + (sdev*4)))
			_binsize = 0.01 ;divide points on the graph into small bins
			;_binsize=1 ;This Binsize considers the QIS Pathfinder #4 which has an ADC that outputs in integer values
			if _binsize le 0 then _binsize=1.;less than or equal to
			if _max lt _min then _max=_min+1 ;less than

			;create histogram arrays that represent output measurements (ADU) of a single pixel
			;theHist (y-axis) includes the number of pixel measurements in each ADU bin. infinite and /NaN data is treated as missing data
			;theVals (x-axis) represents the bin location (the location of the element in the middle of the bin)
			theHist=histogram(pix_data_single, binsize=_binsize, min=_min, max=_max, /nan)
			theVals=findgen(n_elements(theHist))* _binsize+_min+_binsize/2.0 ;bin location (the location of the element in the middle of the bin)

			;remove any bins that do not contain data
			tempidx=where(theHist gt 0,count) ;count: number of element greater than 0,
			temphist= theHist[tempidx] ; tempidx: subscripts of elements that is greater than 0
			tempvals= round(theVals[tempidx]) ;rounded values of the bin's location (that is greater than 0)
			inputVals = double(tempvals) ;tempval: bin's location
			inputHist = double(temphist) ;temphist value of element that is greater than 0

	;----------------------------------------------------------------------------------------------------------------------------------------------


			; FIND MINIMA!!!!
			; Smooth out inputhist to migitate effect of noise, work very well;
			; However,a smoothed histogram's minima will be slightly shifted from original histogram
			; However, minima from smoothed histogram prove to be more reliable and accurate, so as the distance between minima;

			; This condition removes pixels that did not contain enough measurements to have more than 12 unempty bins
			; because we slide first 4 and last 5 of 70% data, so we need at least 12.8 pixels
			if n_elements(inputHist) lt 12 then begin
				outputcg = -1
				goto, NOPEAKS
			endif

			; Hardcoded smooth parameters
			sm_param = 5
			sm_param_normal = 3

			; Smooth the histogram
			input_hist_smooth = smooth(inputHist,sm_param, /EDGE_TRUNCATE)
			input_hist_smooth_normal = smooth(inputHist,sm_param_normal, /EDGE_TRUNCATE)

			; This is used to ignore the last 10% of the histogram where bin shot noise dominates peaks in histogram
			start_array = ceil(1*n_elements(input_hist_smooth)/20)
			end_array = ceil(16*n_elements(input_hist_smooth)/20)
			inputHist1 = input_hist_smooth[start_array:end_array] ; ignoring the last 10% data

			; Find all peak locations in the histogram that are higher then neighboring 4 points and exist above meanlimit * mean(inputHist)
			; This is to prevent small histogram shotnoise peaks from being recognized.
			meanlimit= 0.8
			peak_loc   = WHERE (inputHist1(2:n_elements(inputHist1)-4) ge inputHist1(3:n_elements(inputHist1)-3) $ // +1
							and inputHist1(2:n_elements(inputHist1)-4) ge inputHist1(1:n_elements(inputHist1)-5) $ // -1
							and inputHist1(2:n_elements(inputHist1)-4) ge inputHist1(4:n_elements(inputHist1)-2) $ // +2
							and inputHist1(2:n_elements(inputHist1)-4) ge inputHist1(0:n_elements(inputHist1)-6) $ // -2
							and inputHist1(2:n_elements(inputHist1)-4) gt meanlimit*mean(inputHist1) )
			peak_loc =  peak_loc +2+ start_array ;TRICKY, WHERE function will give your subscript of element in the input(2:n_elements(inputHist)-5)not inputHist

			; Pixel is assigned a NAN conversion gain if not enough peaks are found
			if n_elements(peak_loc) lt 3 then begin
				outputcg = -2
				goto, NOPEAKS
			endif

	;------------------------------------------------------------------------------------------------------------------------------------------------

			; This section removes
			; This part to choose 1 peaks out of multiple close minima/peaks (in case of all equal points)
			;-> 4 is not arbitrary, is the region we compare (3 to the left, 3 to the right) plus 1
			min_distance = make_array(n_elements(peak_loc)-1)
			fix_peak_loc = make_array(n_elements(peak_loc)) ;fix_peak_loc: index of peak that needs to be remove

			; This sets the minimum separation distance for peaks before they are considered too close to actually be neighboring photon number peaks
			; This value is in ADU and is based upon the CG for pixels being 14-16 ADU/e-.
			peak_minimum_separation_limit = 6

			; This loop determines when peak locations j and j+1 are too close
			; If peaks j and j+1 are too close, peak j is shifted to be the center point between j and j+1
			for j = 0, n_elements(peak_loc)-2 do begin
				if peak_loc[j+1] - peak_loc [j] lt peak_minimum_separation_limit then begin  ;
					peak_loc[j] = round((peak_loc[j+1]+peak_loc[j])/2)
					fix_peak_loc [j] = j+1
				endif else begin
					peak_loc[j] = peak_loc[j]
				endelse
			endfor

			; Determine how many peaks where shifted
			number_fix_minima=where(fix_peak_loc gt 0 )

			; If peak j and j+1 are too close then j gets shifted.
			; As a result, j+1 is then removed as a peak location
			if number_fix_minima[0] ge 0 then begin
				for j=0,n_elements(fix_peak_loc)-2 do begin
					if fix_peak_loc[j] gt 0 then peak_loc[fix_peak_loc[j]]=0	;remove the j+1 peak
				endfor
			endif

			; Remove peaks determined not to correspond to photon number peaks
			peak_loc = peak_loc (where(peak_loc gt 0)) ;filter out the latter of the "bad" minima pair
			peak_loc = peak_loc (where(peak_loc gt 0.01*n_elements(inputVals)))

			; If there is only one peak, it is not possible to measure conversion gain and CG for pixel is assigned as NAN
			if n_elements(peak_loc) lt 3 then begin
					outputcg = -3
					goto, NOPEAKS
			endif


	;-----------------------------------------------------------------------------------------------------------------------------------;



	;___________________________________________________________________________________________________________________________________________;

			; Linear fit
			; divide histogram (smooth) by average distance and round it (it will be electron number), if 2 repeated number -> nopeaks or try to save that pixel

			;stop

			peak_location = inputVals (peak_loc)
			peak_val = input_hist_smooth(peak_loc)

			; Check to ensure at least one peak location is greater than 0 ADU
			check_all_zeros_array = 0
			for j =0, n_elements(peak_location)-1 do begin
				if peak_location [j] gt 0 then check_all_zeros_array += 1
			endfor

			; If there is no non-zero ADU peak location, this pixel either did not receive light or only outputs 0 photons
			if check_all_zeros_array eq 0 then begin
				outputcg = -7
				goto, NOPEAKS
			endif

			; compute the distance (delta) between each peak
			deltas = double(fltarr(n_elements(peak_location)-1))
			for m = 0, n_elements(peak_location)-2 do begin
				delta = peak_location[m+1] - peak_location[m]
				deltas[m] = delta
			endfor
			peak_loc_est = peak_location

			; determine if distances (delta) ended up being less than zero
			; to deal with peak_loc all 0 -> delta all 0
			deltas_zeros = where (deltas le 0, count)
			if count eq n_elements(deltas)  then begin
				outputcg = -9
				GOTO, NOPEAKS
			endif

			; remove any deltas equal to or less than zero
			deltas = deltas(where (deltas gt 0))

			; compute pixel's estimate conversion gain based upon average distance (delta) between peaks
			est_CG = double(mean(deltas))

			; use pixel's estimate conversiongain to convert location from ADU to units of electron
			; Important: Since estimate CG will not always be the exact CG, during conversion
			; from ADU to electron, peak locations may not correspond to intereger numbers electrons.
			; Therefore rounding is used to force each peak to the closes intereger electron number and then is linar fit.
			peak_location_e = floor(double(peak_location/est_CG))

			if test_pix_flag eq 1 then begin
				print,'x',x_x
				print,'y',y_y
				print,'adu_location',peak_location
				print,'estimate_CG',est_CG
				print,'e_location',peak_location/est_CG
				print,'e_location_round',peak_location_e
			endif

			; If there is only one peak, it is not possible to measure conversion gain and CG for pixel is assigned as NAN
			num_peak = n_elements(peak_location)
			if num_peak lt 3 then begin
				outputcg = -5
				goto,NOPEAKS
			endif

			; Rounding may cause two peaks to correspond to a single electron number.
			; This code here determines if too many peaks are alighted.
			; In the future, this will instead shift the two peaks correctly to improve linear fitting
			num_error = 0
			for j =0, n_elements(peak_location_e)-2 do begin
				if peak_location_e [j] eq peak_location_e[j+1] then begin
					num_error = num_error + 1
					if num_error eq (num_peak - 2) then begin
						outputcg = -4
						goto,NOPEAKS
					endif
				endif
			endfor

			; Compute pixel's CG from the linear fit of peak location in ADU as a function of peak electron number
			CG_extract = linfit(peak_location_e, peak_location, yfit = ylinfit)
			outputcg = CG_extract[1]

			; Create a plot of the linear fit for a subset of all pixels in region
			peak_location_fit=peak_location
			if (count_fit lt 20) then begin

				outdir_lin = outdir+ "plots\linfit\"

				;
				if test_pix_flag eq 1 then begin
					outdir_lin=test_outdir
				endif


				standard_linfit_hist_jpg_placeholder,outdir=outdir_lin,region=region,x_x=x_x,y_y=y_y,$
						CG_extract=CG_extract,peak_location_fit=peak_location_fit,peak_location_e=peak_location_e,ylinfit=ylinfit
				count_fit=count_fit+1
			endif

	;--------------------------------------PREPARE INTIAL GUESS------------------------------------------------------------------;
	;				From Dartmouth's paper H and u_n is estimated:
	;				H = (kp+1) P[kp+1]/P[kp] where 	P [kp] is probability of highest peak
	;				un = sqrt (1/(8*ln(2/(1-VPM)))  where VPM = 1 -Pv/Pp  where Pp=1/2 (P[kp]+P[kp+1])where	Pv: prob of valley between peak kp, kp+1
	;
	;			##### we know how to cal prob of specific bin, based on y_val
	;			##### question is how many bins that are considered as valley/ peak
	; 			----> solution can have a parameter, default will be 1 (easiest because just divide y_val) if other than 1, will have to add up y_val of neighboring peaks
	;


			peak_val_sort = peak_val[REVERSE(sort(peak_val))]  ; sort peak descending order
			peak_loc_sort = peak_location[reverse(sort(peak_val))]
			peak_index_sort = peak_loc[reverse(sort(peak_val))]


			; In 2 highest peaks, if highest peak comes first, range is from index of highest to 2nd highest
			; if highest peak comes after, range is from index of 2nd highest to highest
			if peak_index_sort[0] gt peak_index_sort[1] then begin
				m = peak_index_sort[1]
				n = peak_index_sort[0]
			endif

			if peak_index_sort[0] lt peak_index_sort[1] then begin
				m = peak_index_sort[0]
				n = peak_index_sort[1]
			endif

			; determine locations of local minimum or valleys
			valley   = WHERE (input_hist_smooth(m:n) lt input_hist_smooth(m+1:n+1) $ // +1
								and input_hist_smooth(m:n) lt input_hist_smooth(m-1:n-1)$ // -1
								and input_hist_smooth(m:n) lt input_hist_smooth(m+2:n+2)$ // +2
								and input_hist_smooth(m:n) lt input_hist_smooth(m-2:n-2)$// -2
								and inputHist1(*) gt 0.03*mean(inputHist) )
			;valley1 = min(input_hist_smooth(m:n), Min_Subscript)
			;valley = Min_Subscript
			valley = valley + m
			min_val = input_hist_smooth[valley]
			min_val_sort = min_val [REVERSE(sort(peak_val))]

			; Compute value for average quanta exposure and readnoise
			H_est = (peak_loc_sort[0]+1)*peak_val_sort[0]/peak_val_sort[1]
			Pp = (peak_val_sort[0]+peak_val_sort[1])/2
			VPM = 1 - min_val_sort[0]/Pp
			H_est = H_est/outputcg ; in e-
			un_est = sqrt(1/(8*alog(2/(1-VPM))))

			; append estimate readnoise (un_est) and average photons (H_est) to array for use as input estimates to curvefit
			A = double([un_est, H_est,0.0])


	;----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------;
	;PCH METHOD TO CALCULATE CG;

			; use CURVEFIT ~> can fit a set of data to a user-defined function (without unlimited parameters) and output fit parameters ----> exactly what we need
			; first create a procedure (pro) to define the function with partial derivatives
			; provide X, Y, weights (?!), and intial guess!!!

			; convert bin locations in units of ADU to electroncs for curvefit
			X = double(inputVals/outputcg)

			; convert bin values into a probability by normalizing bin values to the total number of elements
			prob_exp = outputcg * input_hist_smooth_normal/total(input_hist_smooth_normal)


			; Perform curvefit with user generated function 'gfunct' and hardcoded weight and maximum interations
			weights = 1.0
			itmax=450
			function_name='PCH_funct' ;"The user must define a function which returns the model value."- mpcurvefit documentation
			yfit = mpcurvefit(X,prob_exp,weights,A,function_name=function_name,chisq=chisq,itmax=itmax,quiet=1)


			; determine peak locations not in the raw data but now from the poisson-gaussian output from curvefit
			peak_fit   = WHERE (yfit(2:n_elements(yfit)-4) ge yfit(3:n_elements(yfit)-3) $ // +1
								and yfit(2:n_elements(yfit)-4) ge yfit(1:n_elements(yfit)-5) $ // -1
								and yfit(2:n_elements(yfit)-4) ge yfit(4:n_elements(yfit)-2) $ // +2
								and yfit(2:n_elements(yfit)-4) ge yfit(0:n_elements(yfit)-6) $ // -2
								and yfit(*) gt meanlimit*mean(yfit) )

			; set final computed CG,RN,and H to write to final fits file and maps
			conversiongain = outputcg ;This CG is not from CURVEFIT
			RN = A[0]
			H = A[1]
			offset = A[2]

			; outputs specific to MPCURVEFIT
			if test_pix_flag eq 1 then begin
				testing_dir='D:\QIS1-1\warm1\QISPF\conversiongain\conversiongain.12Mar21\results\plots\testing\'

				fitname='inputVals.fits'
				writefits,testing_dir+fitname,inputVals,header

				fitname='inputVals_e.fits'
				writefits,testing_dir+fitname,X,header

				fitname='input_hist_smooth_normal.fits'
				writefits,testing_dir+fitname,input_hist_smooth_normal,header

				fitname='prob_exp.fits'
				writefits,testing_dir+fitname,prob_exp,header

				fitname='yfit.fits'
				writefits,testing_dir+fitname,yfit,header

				openw,unit,testing_dir+'parameters.txt',/get_lun
				printf,unit,'RN estimate',un_est
				printf,unit,'H estimate',H_est
				printf,unit,'RN fit',RN
				printf,unit,'H fit',H
				printf,unit,'offset fit',offset
				printf,unit,'chisq',chisq
				printf,unit,'peak_fit',peak_fit
				free_lun,unit

				stop
			endif






			; Write CG of pixel (x_x,y_y) in region to final CG array
			finalarray_cg_f[x_x,y_y]= conversiongain

			; If Curvefit output a poor fit, then the number of peaks in peak_fit will be significantly different from the number in peak_lcoation
			if(n_elements(peak_fit) lt n_elements(peak_location)-1)  then begin
				if (n_elements(peak_location) lt 4) then begin
					outputcg = -11
					goto, NOPEAKS
				endif
				offset = !VALUES.F_NaN
				RN = !VALUES.F_NaN
				H = !VALUES.F_NaN
			endif

			; Write parameters to final arrays
			finalarray_rn_f[x_x,y_y]= RN		;Computed RN from curvefit
			finalarray_H_f[x_x,y_y]= H			;Computed H from curvefit
			finalarray_os_f[x_x,y_y] = offset	;Computed horizontal shift/offset of of all peaks from curvefit
			estimate_rn_f [x_x,y_y] = un_est	;Computed RN from dartmouth valley/peak approximation equation (before curvefit)

			; Combine all calculated parameters of a pixel into an array
			param = [conversiongain,RN,H,offset,chisq,un_est]


	;________________________________################ Plotting cases ##############________________________________;

			if test_pix_flag eq 1 then begin
				standard_histogram_jpg_placeholder,outdir=test_outdir,region=region,x_x=x_x,y_y=y_y,peak_location=peak_location,$
										prob_exp=prob_exp,yfit=yfit,param=param,inputvals=inputvals,inputHist=input_hist_smooth_normal
				goto,test_loop
			endif


			if conversiongain eq 12 then begin
				if count_low_CG le 20 then begin
					outdir_low = outdir + "plots\low_CG\"

					count_low_CG = count_low_CG + 1
					standard_histogram_jpg_placeholder,outdir=outdir_low,region=region,x_x=x_x,y_y=y_y,peak_location=peak_location,$
											prob_exp=prob_exp,yfit=yfit,param=param,inputvals=inputvals,inputHist=input_hist_smooth_normal
				endif
			endif


			if ((conversiongain ge 13.5) and conversiongain le 14.5 ) then begin
				if count_std_CG le 20 then begin
					outdir_std = outdir+ "plots\std_CG\"
					count_std_CG = count_std_CG + 1
					standard_histogram_jpg_placeholder,outdir=outdir_std,region=region,x_x=x_x,y_y=y_y,peak_location=peak_location,$
											prob_exp=prob_exp,yfit=yfit,param=param,inputvals=inputvals,inputHist=input_hist_smooth_normal

				endif
			endif

			if conversiongain eq 17 then begin
				if count_high_CG le 20 then begin
					outdir_high = outdir+ "plots\high_CG\"
					count_high_CG = count_high_CG +1
					standard_histogram_jpg_placeholder,outdir=outdir_high,region=region,x=x_x,y_y=y_y,peak_location=peak_location,$
											prob_exp=prob_exp,yfit=yfit,param=param,inputvals=inputvals, inputHist=input_hist_smooth_normal
				endif
			endif

			if H gt 10 then begin
				if count_H le 20 then begin
					outdir_high = outdir+ "plots\high_H\"
					count_H = count_H +1
					standard_histogram_jpg_placeholder,outdir=outdir_high,region=region,x=x_x,y_y=y_y,peak_location=peak_location,$
											prob_exp=prob_exp,yfit=yfit,param=param,inputvals=inputvals, inputHist=input_hist_smooth_normal
				endif
			endif


			if RN gt 0.5 then begin
				if count_RN le 20 then begin
					outdir_high = outdir+ "plots\high_RN\"
					count_RN = count_RN +1
					standard_histogram_jpg_placeholder,outdir=outdir_high,region=region,x=x_x,y_y=y_y,peak_location=peak_location,$
											prob_exp=prob_exp,yfit=yfit,param=param,inputvals=inputvals, inputHist=input_hist_smooth_normal
				endif
			endif

				NOPEAKS:

					if test_pix_flag eq 1 then begin
						param = [outputcg,RN,H,offset,chisq,un_est]
						standard_histogram_jpg_placeholder,outdir=test_outdir,region=region,x=x_x,y_y=y_y,peak_location=peak_location,$
											prob_exp=prob_exp,yfit=yfit,param=param,inputvals=inputvals, inputHist=input_hist_smooth_normal
						goto,test_loop
					endif


					if ((outputcg le 0) and (RN ne 0) and (count_cg_nan le 50) and (outputcg ne -1)) then begin
						param = [outputcg,RN,H,offset,chisq,un_est]
						count_cg_nan = count_cg_nan + 1
						outdir_cg_nan = outdir + "plots\nan\plots\CG_nan_only\"
						standard_histogram_jpg_placeholder,outdir=outdir_cg_nan,region=region,x=x_x,y_y=y_y,peak_location=peak_location,$
											prob_exp=prob_exp,yfit=yfit,param=param,inputvals=inputvals, inputHist=input_hist_smooth_normal
					endif

					if outputcg lt 0 then begin
						param = [outputcg,RN,H,offset,chisq,un_est]
						;prob_exp = outputcg * input_hist_smooth/total(input_hist_smooth)

						if ((count_nan le 50) and (RN eq 0) and (outputcg ne -1)) then begin
							outdir_nan = outdir + "plots\nan\plots\raw\"
							count_nan = count_nan + 1
							standard_histogram_jpg_placeholder,	outdir=outdir_nan,region=region,x=x_x,y_y=y_y,$
												param=param,inputvals=inputvals, inputHist=inputHist
						endif

						if ((count_nan le 50) and (RN eq 0)and(outputcg ne -1)) then begin
							outdir_nan_sm = outdir + "plots\nan\plots\smooth\"
							count_nan = count_nan + 1
							standard_histogram_jpg_placeholder,	outdir=outdir_nan_sm,region=region,x=x_x,y_y=y_y,$
												param=param,inputvals=inputvals, inputHist=input_hist_smooth
						endif

						if ((outputcg eq -1) and (count_nan_worst le 50)) then begin
							param = [0, 0, 0, 0, 0, 0]
							outdir_nan_worst = outdir + "plots\nan\plots\worst\"
							count_nan_worst = count_nan_worst + 1
							standard_histogram_jpg_placeholder,	outdir=outdir_nan_worst,region=region,x=x_x,y_y=y_y,$
												param=param,inputvals=inputvals, inputHist=inputHist
						endif
						if (keyword_set(test_flag)and(outputcg ne -1)) then begin
							out_test= outdir+ "plots\test\"
							standard_histogram_jpg_placeholder,outdir=outdir_high,region=region,x=x_x,y_y=y_y,peak_location=peak_location,$
											prob_exp=prob_exp,yfit=yfit,param=param,inputvals=inputvals, inputHist=input_hist_smooth_normal
						endif
						conversiongain = !VALUES.F_NaN
						finalarray_cg_f[x_x,y_y]= conversiongain
						RN =!VALUES.F_NaN
						finalarray_rn_f[x_x,y_y]= RN
						estimate_rn_f [x_x,y_y] = RN
						H =!VALUES.F_NaN
						finalarray_H_f[x_x,y_y]= H
						offset =!VALUES.F_NaN
						finalarray_os_f[x_x,y_y] = offset
						chisq=  !VALUES.F_NaN
						un_est = !VALUES.F_NaN

					endif

		endfor
	endfor

	;Saving the final 2D CG array to a fits file and a jpg image
	writefits,outdir+'conversiongain.fits',finalarray_cg_f,header
	writefits,outdir+'readnoise.fits',finalarray_rn_f,header
	writefits,outdir+'quantaexposure.fits',finalarray_H_f,header
	writefits,outdir+'offset.fits',finalarray_os_f,header
	writefits,outdir+'est_rn.fits',estimate_rn_f,header

	;Save the final fits file name to a text file
	openw,unit,outdir+list,/get_lun
	printf,unit,'conversiongain.fits'
	free_lun,unit

	; Plot the image map
	plot_image_jpg,outdir,finalarray_cg_f,outfilename='conversiongain_map.jpg',outdir=mapdir

	openw,unit,outdir+list_rn,/get_lun
	printf,unit,'readnoise.fits'
	free_lun,unit

	; Plot the image mapo
	plot_image_jpg,outdir,finalarray_rn_f,outfilename='readnoise_map.jpg',outdir=mapdir

	openw,unit,outdir+list_h,/get_lun
	printf,unit,'quantaexposure.fits'
	free_lun,unit

	; Plot the image map
	plot_image_jpg,outdir,finalarray_H_f,outfilename='quantaexposure_map.jpg',outdir=mapdir


	openw,unit,outdir+list_os,/get_lun
	printf,unit,'offset.fits'
	free_lun,unit

	; Plot the image map
	;plot_image_jpg,outdir,finalarray_os_f,outfilename='offset_map.jpg',outdir=mapdir


	openw,unit,outdir+list_est_rn,/get_lun
	printf,unit,'est_rn.fits'
	free_lun,unit

	; Plot the image map
	plot_image_jpg,outdir,estimate_rn_f,outfilename='est_rn_map.jpg',outdir=mapdir


	; Make histogram of the output fits file
	makehistogram_conversiongain, outdir, list,'Conversion Gain','Conversion Gain (ADU/e!E-!N)',region=region, outdir=histdir

	; Make histogram of the output fits file
	makehistogram_conversiongain, outdir, list_rn,'Read Noise','Read Noise (e!E-!N rms)',region=region, outdir=histdir

	; Make histogram of the output fits file
	makehistogram_conversiongain, outdir, list_h,'Quanta Exposure','H (e!E-!N)',region=region, outdir=histdir

	; Make histogram of the output fits file
	makehistogram_conversiongain, outdir, list_os,'Offset','Offset (e!E-!N)',region=region, outdir=histdir

	; Make histogram of the output fits file
	makehistogram_conversiongain, outdir, list_est_rn,'Estimated Read Noise','Estimated Read Noise (e!E-!N rms)',region=region, outdir=histdir



	; Print the final time of running conversion gain analysis
	time_of_end=systime()
	print,'Time of Start:',time_of_start
	print, 'Time of concatenation end', time_of_concatenation
	print,'Time of end:',time_of_end
	print, 'Conversion Gain Analysis Complete'
	stop


	;goto,superloop




END