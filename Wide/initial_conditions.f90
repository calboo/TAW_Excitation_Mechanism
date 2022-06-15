MODULE initial_conditions

  USE shared_data
  USE neutral

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: set_initial_conditions

CONTAINS

  SUBROUTINE set_initial_conditions
    
    INTEGER :: ix, iy, iz
    REAL(num) :: beta, xpos, ypos, zpos, r, theta, H

  ! Gravity and plasma beta
    grav = 0.0_num
    beta = 1.0_num

  ! Velocities
  ! Static domain
    vx = 0.0_num
    vy = 0.0_num
    vz = 0.0_num

  ! Bx Field

    DO ix = -2, nx+2
       DO iy = -1, ny+2
          DO iz = -1, nz+2
             xpos = xb(ix)
             ypos = yc(iy)
             zpos = zc(iz) 
             r = sqrt(xpos**2+ypos**2)
             theta = atan2(ypos,xpos)
             bx(ix,iy,iz) = exp(-zpos/H)*bessel_j1(r/H)*cos(theta)
          END DO
       END DO
    END DO

  ! By Field

    DO ix = -1, nx+2
       DO iy = -2, ny+2
          DO iz = -1, nz+2
             xpos = xc(ix)
             ypos = yb(iy)
             zpos = zc(iz)
             r = sqrt(xpos**2+ypos**2)
             theta = atan2(ypos,xpos)
             by(ix,iy,iz) = exp(-zpos/H)*bessel_j1(r/H)*sin(theta)
          END DO
       END DO
    END DO

  ! Bz field

    DO ix = -1, nx+2
       DO iy = -1, ny+2
          DO iz = -2, nz+2
             xpos = xc(ix)
             ypos = yc(iy)
             zpos = zb(iz) 
             r = sqrt(xpos**2+ypos**2)
             bz(ix,iy,iz) = exp(-zpos/H)*bessel_j0(r/H)
          END DO
       END DO
END DO

  ! Density

    DO iy= -1,ny+2 
       DO iz = -1,nz+2 
          DO ix = -1,nx+2 
             rho(ix,iy,iz) = exp(-zc(iz)/H) *&
             (2.0_num - tanh(rampsteep*(sqrt(xc(ix)**2+yc(iy)**2)-vortexrad)))
          END DO
       END DO
    END DO

  ! Energy

    energy=0.5_num*(beta*1.0_num) / ((rho)*(gamma-1.0_num))

  ! Store Initial Values
  ! Initial values are stored in these arrays

    ALLOCATE(rho0(-1:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(bx0 (-2:nx+2, -1:ny+2, -1:nz+2))
    ALLOCATE(by0 (-1:nx+2, -2:ny+2, -1:nz+2))
    ALLOCATE(bz0 (-1:nx+2, -1:ny+2, -2:nz+2))
    ALLOCATE(energy0(-1:nx+2, -1:ny+2, -1:nz+2))

    bx0 = bx
    by0 = by
    bz0 = bz
    rho0 = rho
    energy0 = energy

  END SUBROUTINE set_initial_conditions

END MODULE initial_conditions
