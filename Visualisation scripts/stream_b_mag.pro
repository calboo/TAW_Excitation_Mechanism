pro stream_b_mag, ds, z

;; DESCRIPTION

; IDL script used to produce Figure 1C in the main body 
; Figure S5 in the supplementary material for the paper,
; Torsional oscillations within a magnetic pore in the solar
; photosphere, M.Stangalini, R.Erdelyi, C.Boocock et al.

; Partial magnetic fieldlines from a Lare3d simulation output, ds,
; are plotted across a plane at a chosen height, z.
; The field lines are plotted over a coloured contour of the
; total magnetic field strength across the plane at height z.

; Figures were produced from the outputs of the wider simulation variant.

;; USAGE

; To use this script data must first be loaded from the chosen 
; output using the getdata function provided for Lare3d in IDL.

; The loaded output must then be specified when calling this script,
; additionally the height at which to plot the contour must be given.
; ds is the chosen output and z is the height in simulation cells.

;; SCRIPT

; Set up arrays dimensions
nx=200
ny=200

; Set up arrays for calculating field strengths
arrx      =dblarr(nx,ny)
arry      =dblarr(nx,ny)
arrz      =dblarr(nx,ny)
magfield  =dblarr(nx,ny)

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

; Calculate the magnetic fields in each direction across the plane
; Averages must be taken to acount for magnetic grid staggering
for i= 0,nx-1 do begin
   for j = 0,ny-1 do begin
      arrx[i,j]=0.5*ds.bx[i,j,z]+ds.bx[i+1,j,z]
      arry[i,j]=0.5*ds.by[i,j,z]+ds.by[i,j+1,z]
      arrz[i,j]=0.5*ds.bz[i,j,z]+ds.bz[i,j,z+1]
   endfor
endfor

; Calculate the magnetic field strength across the plane
magfield  = sqrt(arrx^2+arry^2+arrz^2)

; Define colour table
ct = colortable([[0,0,0],[0,0,10],[0,0,20],[0,0,40],[0,0,50],[0,0,60],[0,0,100],[0,100,255],[255,255,255]], /TRANSPOSE)

; Set contour levels
clevels =  ((findgen(100)*0.30)/100)+1.10

; Set up window
w = WINDOW(DIMENSIONS=[800,800])

; Plot magnetic field strength contour
cn = CONTOUR(magfield,x,y,/fill, $
             POSITION=[0.1,0.22,0.95,0.9], $ 
             XTITLE='X / km', YTITLE='Y / km', $         
             N_LEVELS=100,c_value=clevels,$;/buffer,$
             RGB_TABLE=ct,  ASPECT_RATIO=1, FONT_SIZE=15, $
             XRANGE = [-8,8], YRANGE = [-8,8], /current)

; Plot partial magnetic field lines
s = STREAMLINE(xdot, ydot, x, y, $
     STREAMLINE_STEPSIZE=500.0, STREAMLINE_NSTEPS=200, $
     POSITION=[0.1,0.22,0.95,0.9], $
     X_STREAMPARTICLES=20, Y_STREAMPARTICLES=20, $
     XTITLE='X / Mm', YTITLE='Y / Mm', $
     ;TITLE='Title of Plot',$
     THICK=3,AUTO_COLOR = 1,AUTO_RANGE=[0.0,0.02],$
     RGB_TABLE=3,/overplot))

; Add a colorbar.
c = COLORBAR(TARGET=s, MINOR=0, $
             POSITION=[0.15,0.1,0.95,0.125],$
             FONT_SIZE=15, $
             TITLE='Magnetic Field Strength / G',$
             TICKNAME=['0','20','40','60','80','100'])

; Save as either a png or eps file
c.save, 'stream_b_mag.png'
;c.Save, 'stream_b_mag.eps';, BITMAP=1, width=10

end

