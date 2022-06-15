pro stream_diff, ds, z

;; DESCRIPTION

; IDL script used to produce Figure 4 in the main body 
; and the first two panels in Figure S6
; in the supplementary material for the paper
; Torsional oscillations within a magnetic pore in the solar
; photosphere, M.Stangalini, R.Erdelyi, C.Boocock et al.

; Velocity streamlines from a Lare3d simulation output, ds,
; are plotted across a plane at a chosen height, z.
; The streamlines are plotted over a coloured contour of the
; difference between the density of the chosen output and the  
; initial density across the plane at height z.

; Figures were produced from the outputs of the original simulation.

;; USAGE

; To use this script data must first be loaded from the chosen 
; output using the getdata function provided for Lare3d in IDL.

; The loaded output must then be specified when calling this script,
; additionally the height at which to plot the contour must be given.
; ds is the chosen output and z is the height in simulation cells.

;; SCRIPT

; Load initial density
ds0 = getdata(0, /rho)

; Load positions and velocities, positions are scaled into Mm,
; velocities are normalised into kms-1 and then scaled for aesthetics.
x = ds.x*0.2
y = ds.y*0.2
xdot = ds.vx[1:*,1:*,z]*8.9206/250
ydot = ds.vy[1:*,1:*,z]*8.9206/250

; Set levels for density contour
clevels =  (findgen(100)-50.0)/250.0

; Plot density contour
loadct, 1
cn = CONTOUR(((ds.rho[*,*,z]-ds0.rho[*,*,z])/ds0.rho[*,*,z]),x,y,/fill, $
     POSITION=[0.1,0.22,0.95,0.9], $ 
     XTITLE='X / Mm', YTITLE='Y / Mm', $         
     ASPECT_RATIO=1,c_value=clevels,RGB_TABLE=5)

; Plot partial streamlines
loadct, 3
s = STREAMLINE(xdot, ydot, x, y, $
     STREAMLINE_STEPSIZE=10.0, $
     POSITION=[0.1,0.22,0.95,0.9], $
     X_STREAMPARTICLES=20, Y_STREAMPARTICLES=20, $
     XTITLE='X / Mm', YTITLE='Y / Mm', $
     ;TITLE='Title of Plot',$
     THICK=2,/overplot)

; Add a colorbar.
c = COLORBAR(TARGET=cn, $
             POSITION=[0.25,0.1,0.8,0.125], $
             TITLE='Relative Difference in Density')

; Save as either a png or eps file
c.save, 'stream_diff.png'
;c.Save, 'stream_diff.eps';, BITMAP=1, width=10

end

