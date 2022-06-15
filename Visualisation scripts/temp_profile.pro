pro temp_profile, ds

;; DESCRIPTION

; IDL script used to produce the right panel of Figure S3
; in the supplementary material for the paper,
; Torsional oscillations within a magnetic pore in the solar
; photosphere, M.Stangalini, R.Erdelyi, C.Boocock et al.

; A coloured contour of the temperature profile is plotted across the 
; vertical midplane of the domain for the loaded simulation output.

; Figures were produced from the outputs of the original simulation.

;; USAGE

; To use this script data must first be loaded from the chosen 
; output using the getdata function provided for Lare3d in IDL.
; The loaded output must then be specified when calling this script.

;; SCRIPT

; Load colour table
ct = COLORTABLE(13)

; Set up window
w = WINDOW(DIMENSIONS=[800,800])

; Plot temperature contour
c = CONTOUR(reform(ds.temperature[*,100,*])*11556,$
            ds.x*0.2,ds.z*0.2,/fill, $
            POSITION=[0.1,0.22,0.95,0.9],RGB_TABLE=ct , $
            TITLE='Temperature profile',N_LEVELS=100,$
            XTITLE='X / Mm', YTITLE='Height / Mm',FONT_SIZE=15,$
            XSTYLE=1,YSTYLE=1,/CURRENT)

; Add a colorbar.
cb = COLORBAR(target=c, TICKFORMAT='(I4)',$
 POSITION=[0.15,0.1,0.95,0.125], $
 TITLE='Temperature / K',FONT_SIZE=15)

; Save as either a png or eps file
c.Save, 'temp_prof.png'
;c.Save, 'temp_prof.eps', BITMAP=1, width=10 

end
