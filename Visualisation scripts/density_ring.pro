pro density_ring

;; DESCRIPTION

; IDL script used to produce third panel in Figure S6
; in the supplementary material for the paper
; Torsional oscillations within a magnetic pore in the solar
; photosphere, M.Stangalini, R.Erdelyi, C.Boocock et al.

; The density around a ring of radius 2 Mm, at a chosen height z
; is read at each from each Lare3d output and the difference
; in density around  this ring, relative to the initial density
; around the ring, is plotted as a colored contour. 

;; USAGE

; This script must be run from the directory containing the outputs
; from the Lare3d simulation. There must be 1501 output files.
; This script has only been used on the original simulation
; not the wider variant.

;; SCRIPT

; Set height in simulation cells
z = 25

; Create arrays  
r = dblarr(200)
yinds = dblarr(200)
xinds= dblarr(200)
rho_ring = dblarr(200)
angle = dblarr(200)
density_ring = dblarr(1501,200)
difference = dblarr(1501,200)

; Load grid and density data at first timestep
ds0 = getdata(0,/grid,/rho)

; Calculate the x & y indices of points around the ring r = 2 Mm
for i = 50,149 do begin
   for j = 0,199 do begin
      ; Calculate r^2 for a particular y value
      r[j] = ds0.y[i]^2+ds0.x[j]^2 
   endfor
   ; Locate where the radius is 2 Mm either side of x = 0
   loc1 = VALUE_LOCATE(r[0:100], 100)
   loc2 = VALUE_LOCATE(r[100:*], 100)
   ; Store indices in periodic arrays xinds and yinds
   yinds[i-50] = i
   xinds[i-50] = loc1+1 
   yinds[-i+49] = i
   xinds[-i+49] = 100+loc2
endfor

; Calculate density and angle of each point around the ring
; NB the index of rho_ring and angle has been shifted 
; to coincide with the x axis
for i = 0, 199 do begin
   rho_ring[i] =  ds0.rho[xinds[i-50],yinds[i-50],z]
   angle[i] = atan(ds0.y[yinds[i-50]],-ds0.x[xinds[i-50]])*(180/!pi)+180
endfor

; Plot the ring at r = 2 Mm

; Quick plot
window, 1
plot, ds0.x[xinds]*0.2, ds0.y[yinds]*0.2,psym=4, $
xrange=[min(ds0.x)*0.2,max(ds0.x)*0.2], $
yrange=[min(ds0.y)*0.2,max(ds0.y)*0.2], /iso

; Fancy rainbow ring plot
; Takes a while to load, comment out if unecessary
p = PLOT(ds0.x[xinds]*0.2,ds0.y[yinds]*0.2, $
               xrange=[min(ds0.x)*0.2,max(ds0.x)*0.2], $
               yrange=[min(ds0.y)*0.2,max(ds0.y)*0.2], $
               DIM=[500,500],LINESTYLE=6,$
               XTITLE='X / Mm', YTITLE='Y / Mm')
for c = 0,199 do begin
   s = SYMBOL(ds0.x[xinds[c-50]]*0.2, ds0.y[yinds[c-50]]*0.2, $
              'D',SYM_COLOR=[50+c,250-2.0*sqrt((100-c)^2),200-c],/SYM_FILLED,/DATA)
endfor
; Save rainbow ring plot as a png 
p.Save, 'ring.png'

; Load the density at each timestep and store the densities
; around the ring in the array density_ring.
;
; This process takes a long time, 
; it is recommended to save density_ring in .sav file upon first use
; then reload from .sav file instead of running this very long loop each time.

for t = 0,1500 do begin
   ds = getdata(t,/rho)
   for i = 0, 199 do begin
      density_ring[t,i] =  ds.rho[xinds[i-50],yinds[i-50],z]
   endfor
   print, t
endfor
save, density_ring, FILENAME = 'density_ring_z25.sav'

restore, 'density_ring_z25.sav'

; Calculates the difference between the density around the ring at a
; given time in the simulation and the initial density around the ring

for t = 0,1500 do begin
   for i = 0, 199 do begin
      difference[t,i] = (density_ring[t,i]-rho_ring[i])/rho_ring[i]
   endfor
endfor

; Plot a coloured contour of the difference in density around the ring
; as it varies with the position around the ring measured in degrees
; and the simulation time given in units of the characteristic Alfven time.

; Create simple time index for x-axis
time = indgen(1501)

; Simple contour of density difference around ring

window, 2
contour, difference, time ,angle, nlevels=100, /fill

; Fancy contour of density difference around ring
; Used in Figure S6 in supplementary material for the paper
; Torsional oscillations within a magnetic pore in the solar
; photosphere, M.Stangalini, R.Erdelyi, C.Boocock et al.

clevels =  (findgen(101)-50.0)/100
cn = CONTOUR(difference, time, angle, /fill, $
            POSITION=[0.1,0.22,0.95,0.9], $ 
            XTITLE='Time / $\tau_A$', YTITLE='Angle / degrees', $
            RGB_TABLE=5,c_value=clevels, YRANGE=[0,360])

c = COLORBAR(TARGET=cn, POSITION=[0.2,0.1,0.8,0.125], $
            FONT_SIZE=10, $
            TITLE='Relative Difference in Density')
; Save Fancy contour as a png 
cn.Save, 'density_ring.png'

end
