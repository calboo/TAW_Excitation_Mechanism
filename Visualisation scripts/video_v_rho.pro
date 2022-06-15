pro video_v_rho

;; DESCRIPTION

; IDL script used to produce supplemenary videos 1,2 and 3
; given with the paper,
; Torsional oscillations within a magnetic pore in the solar
; photosphere, M.Stangalini, R.Erdelyi, C.Boocock et al.

; This script produces a video that shows velocity streamlines
; from a Lare3d simulation plotted across a plane at a chosen height,
; z, at different times throughout the simulation.
; The streamlines are plotted over a coloured contour of the
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
  video_file = 'V_rho_500km.mp4'
  video = idlffvideowrite(video_file)
  framerate = 50
  framedims = [640,512]
  stream = video.addvideostream(framedims[0], framedims[1], framerate)
  set_plot, 'z', /copy
  device, set_resolution=framedims, Z_Buffer=0, set_pixel_depth=24, decomposed=0
  
; Select number of frames to be included in the video
  nframes = 1250

; Set contour levels
  clevels =  findgen(100)*0.3/100

; Begin loop for adding frames to video
for t=250,250+nframes do begin

   ; Load output at time t
   ds = getdata(t,/vx,/vy,/rho,/grid)

   ; Load positions and velocities, positions are scaled into Mm,
   ; velocities are normalised into kms-1 and then scaled for aesthetics.
   x = ds.x*0.2
   y = ds.y*0.2
   xdot = ds.vx[1:*,1:*,z]*8.9206/250
   ydot = ds.vy[1:*,1:*,z]*8.9206/250

   ; Plot density contour
   loadct, 1
   cn = CONTOUR(ds.rho[*,*,z]*0.1,x,y,/fill, $
                POSITION=[0.1,0.22,0.95,0.9], $ 
                XTITLE='X / Mm', YTITLE='Y / Mm', $         
                N_LEVELS=100,c_value=clevels,/buffer)

   ; Plot partial streamlines
   loadct, 3
   s = STREAMLINE(xdot, ydot, x, y, $
                  STREAMLINE_STEPSIZE=10.0, $
                  POSITION=[0.1,0.22,0.95,0.9], $
                  X_STREAMPARTICLES=20, Y_STREAMPARTICLES=20, $
                  XTITLE='X / Mm', YTITLE='Y / Mm', $
                  TITLE='Title of Plot',$
                  THICK=2,AUTO_COLOR = 1,AUTO_RANGE=[0.0,0.02],$
                  RGB_TABLE=3,/overplot,/buffer)

   ; Add a colorbar.
   c = COLORBAR(TARGET=s, MINOR=0, $
                POSITION=[0.35,0.1,0.7,0.125], $
                TITLE='Velocity / kms^-1',$
                TICKNAME=['0','0.1','0.2','0.3','0.4','0.5'])

   ; Add timer
   t1 = TEXT(0.1, 0.09, 'time ='+strcompress(string((t-250)*2.24,FORMAT='(F7.1)')+'s'),$
             FONT_SIZE=12)

   ; Add frame to video
   timestamp = video.put(stream, s.copywindow())
   s.erase

   ; Print the frame number
   print, t-250

endfor

  ; Save video to file
  device, /close
  set_plot, strlowcase(!version.os_family) eq 'windows' ? 'win' : 'x'
  video.cleanup
  print, 'File "' + video_file + '" written to current directory.'

end
