pro video_b_rho

;; DESCRIPTION

; IDL script used to produce supplemenary videos 4,5 and 6
; given with the paper,
; Torsional oscillations within a magnetic pore in the solar
; photosphere, M.Stangalini, R.Erdelyi, C.Boocock et al.

; This script produces a video that shows partial magnetic fieldlines
; from a Lare3d simulation plotted across a plane at a chosen height,
; z, at different times throughout the simulation.
; The field lines are plotted over a coloured contour of the
; density across the plane at height z.

; Figures were produced from the outputs of the original simulation.

;; USAGE

; This script must be run from the directory containing the outputs
; from the Lare3d simulation. 

;; SCRIPT

; Set the height of the horizontal plane in units of simulation cells
  z = 25

; Set up mp4 file
  compile_opt idl2
  video_file = 'B_rho_500km.mp4'
  video = idlffvideowrite(video_file)
  framerate = 50
  framedims = [640,512]
  stream = video.addvideostream(framedims[0], framedims[1], framerate)
  set_plot, 'z', /copy
  device, set_resolution=framedims, Z_Buffer=0, set_pixel_depth=24, decomposed=0

; Load initial magnetic field data
  ds0 = getdata(0, /bx,/by,/grid)

; Select perpendicular magnetic fields
; in the horizontal plane at height z
  bx0 = ds0.bx[1:*,*,z]
  by0 = ds0.by[*,1:*,z]

; Select number of frames to be included in the video
  nframes = 1250

; Set contour levels
  clevels =  findgen(100)*0.3/100

; Begin loop for adding frames to video
for t=250,250+nframes do begin

   ; Load output at time t
   ds = getdata(t, /bx, /by, /rho, /grid)

   ; Load positions and scale into Mm
   x = ds.x*0.2
   y = ds.y*0.2

   ; Deduct the background magnetic field from the 
   ; perpendicular magnetic field data and scale
   xdot = (ds.bx[1:*,*,z]-bx0)/10/250
   ydot = (ds.by[*,1:*,z]-by0)/10/250
   
   ; Plot density contour
   loadct, 1
   cn = CONTOUR(ds.rho[*,*,z]*0.1,x,y,/fill, $
                POSITION=[0.1,0.22,0.95,0.9], $ 
                XTITLE='X / km', YTITLE='Y / km', $         
                N_LEVELS=100,c_value=clevels,/buffer)

   ; Plot partial magnetic field lines
   loadct, 3
   s = STREAMLINE(xdot, ydot, x, y, $
                  STREAMLINE_STEPSIZE=500.0, $
                  POSITION=[0.1,0.22,0.95,0.9], $
                  X_STREAMPARTICLES=20, Y_STREAMPARTICLES=20, $
                  XTITLE='X / Mm', YTITLE='Y / Mm', $
                  TITLE='Title of Plot',$
                  THICK=2,AUTO_COLOR = 1,AUTO_RANGE=[0.0,0.02],$
                  RGB_TABLE=3,/overplot,/buffer)

   ; Add a colorbar.
   c = COLORBAR(TARGET=s, MINOR=0, $
                POSITION=[0.35,0.1,0.7,0.125], $
                TITLE='Magnetic Field Strength / G',$
                TICKNAME=['0','20','40','60','80','100'])

   ; Add timer  
   t1 = TEXT(0.1, 0.09, 'time ='+strcompress(string((t-250)*2.24,FORMAT='(F7.1)')+'s'),$
             FONT_SIZE=12)

   ; Add frame to video
   timestamp = video.put(stream, s.copywindow())
   s.erase

   ; print the frame number
   print, t-250

endfor

; Save video to file
  device, /close
  set_plot, strlowcase(!version.os_family) eq 'windows' ? 'win' : 'x'
  video.cleanup
  print, 'File "' + video_file + '" written to current directory.'

end
