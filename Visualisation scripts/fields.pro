pro fields, ds

;; DESCRIPTION

; IDL script used to produce Figure S2 in the 
; supplementary material for the paper,
; Torsional oscillations within a magnetic pore in the solar
; photosphere, M.Stangalini, R.Erdelyi, C.Boocock et al.

; Magnetic field lines are plotted over a coloured contour of
; the magnetic potential across the vertical midplane of the domain.
; This is done both using the initial simulation output and using the
; analytical formulae for the magnetic field coordinates,
; in order to demonstrate their equivalence.

; Figures were produced from the outputs of the original simulation.

;; USAGE

; To use this script data must first be loaded from the chosen 
; output using the getdata function provided for Lare3d in IDL.
; The loaded output must then be specified when calling this script.

;; SCRIPT

; Set up arrays for x and z
; Scale to Mm using length normalisation
x1 = ds.x*0.2
x2 = [-4.0,x1+0.02]
z1 = ds.z*0.2
z2 = [0.0,z1+0.01]

; Enter parameter for the magnetic scale height in Mm
; This should be the same as in the simulation
; but note the length normalisation that is applied

H = 10.0 ; This is equivalent to 50 in simulation units

; Read the magnetic fields in x and z directions 
; across the vertical midplane of the domain
bx = reform(ds.bx[*,100,*])
bz = reform(ds.bz[*,100,*])
sx = size(bx)
sz = size(bz)

; Set up the colouring for plots

TVLCT, 255, 255, 255, 254 ; White color
TVLCT, 0, 0, 0, 253       ; Black color
!P.Color = 253
!P.Background = 254


; Numerical calculation from B-fields to demonstrate it is the same as
; the analytic expression below

; Calculate phi as the integral of Bz dz

phi = dblarr(sz(1),sz(2))
psi = dblarr(sz(1),sz(2))

for i = 0, sz(1)-1 do begin
      phi(i,0) = -H*beselj(abs(x2[i])/H,0)
endfor
for i = 0, sz(1)-1 do begin
   for j = 1, sz(2)-1 do begin
      phi(i,j) = -H*beselj(abs(x2[i])/H,0) + int_tabulated(z2[0:j],bz[i,0:j])
   endfor
endfor

;Calculate psi as the integral of r*Bz dr / H

for i = ((sz(1)-1)/2)+1, sz(1)-1 do begin
   for j = 0, sz(2)-1 do begin
      psi(i,j) = int_tabulated(x2[(sz(1)-1)/2:i],abs(x2((sz(1)-1)/2:i))*bz[(sz(1)-1)/2:i,j])
   endfor
endfor
for i = 0, (sz(1)-1)/2 do begin
   for j = 0, sz(2)-1 do begin
      psi(i,j) = psi(sz(1)-1-i,j)
   endfor
endfor
psi = psi/H

; For plotting the central contour
psi(99,*) = 0.0
psi(100,*) = 0.0

; Shaded surface of phi calculation
window,1
shade_surf, phi
; Shaded surface of psi calculation
window,2
shade_surf, psi
; Coloured contour of phi with contours of psi overplotted
; This effectively shows the field lines of the setup
window,3
contour, phi,x1,z2, XTITLE= 'x (Mm)', YTITLE= 'Height (Mm)', CHARSIZE=2, CHARTHICK=2,/fill, nlevels=1000
contour, psi,x1,z2,/overplot, THICK=2, levels = (((findgen(10))*0.1)^2.0)+2.0e-7


; Analytic expression of B-fields to demonstrate it is the same as
; the numerical calculation above 

phi = dblarr(sx(1),sz(2))
psi = dblarr(sx(1),sz(2))

for i = 0,sx(1)-1 do begin
   for j = 0, sz(2)-1 do begin
      phi(i,j) = -H*exp(-z2[j]/H)*beselj(abs(x2[i])/H,0)
      psi(i,j) = abs(x2[i])*exp(-z2[j]/H)*beselj(abs(x2[i])/H,1)
   endfor
endfor

; Shaded surface of phi calculation
window,11
shade_surf, phi
; Shaded surface of psi calculation
window,12
shade_surf, psi
; Coloured contour of phi with contours of psi overplotted
; This effectively shows the field lines of the setup
window,13
contour, phi,x2,z2, XTITLE= 'x (Mm)', YTITLE= 'Height (Mm)', CHARSIZE=2, CHARTHICK=2,/fill, nlevels=1000
contour, psi,x2,z2,/overplot, THICK=2, levels = (((findgen(10))*0.1)^2.0)+2.0e-7

; Write out the final plot as a png file

WRITE_PNG, 'fields2d.png', TVRD(/TRUE)

end

        
