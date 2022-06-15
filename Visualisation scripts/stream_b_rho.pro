pro stream_b_rho, ds, z

;; DESCRIPTION

; IDL script used to produce Figure S4 in the 
; supplementary material for the paper,
; Torsional oscillations within a magnetic pore in the solar
; photosphere, M.Stangalini, R.Erdelyi, C.Boocock et al.

; Partial magnetic fieldlines from a Lare3d simulation output, ds,
; are plotted across a plane at a chosen height, z.
; The field lines are plotted over a coloured contour of the
; density across the plane at height z.

; Figures were produced from the outputs of the original simulation.

;; USAGE

; To use this script data must first be loaded from the chosen 
; output using the getdata function provided for Lare3d in IDL.

; The loaded output must then be specified when calling this script,
; additionally the height at which to plot the contour must be given.
; ds is the chosen output and z is the height in simulation cells.

;; SCRIPT

; Load initial magnetic field data
ds0 = getdata(0, /bx,/by,/grid)

; Select perpendicular magnetic fields
; in the horizontal plane at height z
bx0 = ds0.bx[1:*,*,z]
by0 = ds0.by[*,1:*,z]

; Load positions and scale into Mm
x = ds.x*0.2
y = ds.y*0.2

; Deduct the background magnetic field from the 
; perpendicular magnetic field data and scale
xdot = (ds.bx[0:199,0:199,z]-bx0)/10.0/250
ydot = (ds.by[0:199,0:199,z]-by0)/10.0/250

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

; Plot partial magnetic field lines
loadct, 3
s = STREAMLINE(xdot, ydot, x, y, $
     STREAMLINE_STEPSIZE=500.0, $
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
             TITLE='Magnetic Field Strength / G',$
             TICKNAME=['0','20','40','60','80','100'])

; Save as either a png or eps file
c.save, 'stream_b_rho.png'
;c.Save, 'stream_b_rho.eps';, BITMAP=1, width=10

end

