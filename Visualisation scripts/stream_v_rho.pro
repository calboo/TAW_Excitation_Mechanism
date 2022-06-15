pro stream_v_rho, ds, z

;; DESCRIPTION

; IDL script used to produce Figure S4 in the 
; supplementary material for the paper,
; Torsional oscillations within a magnetic pore in the solar
; photosphere, M.Stangalini, R.Erdelyi, C.Boocock et al.

; Velocity streamlines from a Lare3d simulation output, ds,
; are plotted across a plane at a chosen height, z.
; The streamlines are plotted over a coloured contour of the
; density across the plane at height z.

; Figures were produced from the outputs of the original simulation.

;; USAGE

; To use this script data must first be loaded from the chosen 
; output using the getdata function provided for Lare3d in IDL.

; The loaded output must then be specified when calling this script,
; additionally the height at which to plot the contour must be given.
; ds is the chosen output and z is the height in simulation cells.

;; SCRIPT

; Load positions and velocities, positions are scaled into Mm,
; velocities are normalised into kms-1 and then scaled for aesthetics.
x = ds.x*0.2
y = ds.y*0.2
xdot = ds.vx[1:*,1:*,z]*8.9206/250
ydot = ds.vy[1:*,1:*,z]*8.9206/250

; Set contour levels
clevels =  findgen(100)*0.4/100

; Set up window
w = WINDOW(DIMENSIONS=[800,800])

; Plot density contour
loadct, 1
cn = CONTOUR(ds.rho[*,*,z]*0.1,x,y,/fill, $
     POSITION=[0.1,0.22,0.95,0.9], $ 
     XTITLE='X / Mm', YTITLE='Y / Mm',FONT_SIZE=15, $         
     N_LEVELS=100, ASPECT_RATIO=1,c_value=clevels,/current)

; Plot partial streamlines
loadct, 3
s = STREAMLINE(xdot, ydot, x, y, $
     STREAMLINE_STEPSIZE=10.0, $
     POSITION=[0.1,0.22,0.95,0.9], $
     X_STREAMPARTICLES=20, Y_STREAMPARTICLES=20, $
     XTITLE='X / Mm', YTITLE='Y / Mm', $
     ;TITLE='Title of Plot',$
     THICK=3,AUTO_COLOR = 1,AUTO_RANGE=[0.0,0.02],$
     RGB_TABLE=3,/overplot)

; Add a colorbar.
c = COLORBAR(TARGET=s, MINOR=0, $
             POSITION=[0.15,0.1,0.95,0.125],$
             FONT_SIZE=15, $
             TITLE='Velocity / kms!U-1',$
             TICKNAME=['0','0.1','0.2','0.3','0.4','0.5'])

; Save as either a png or eps file
c.save, 'stream_v_rho.png'
;c.Save, 'stream_v_rho.eps';, BITMAP=1, width=10

end

