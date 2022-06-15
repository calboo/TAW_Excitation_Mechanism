!******************************************************************************
! This module contains the boundary conditions for the entire code
! Any new boundary conditions should be added here
!******************************************************************************

MODULE boundary

  USE shared_data
  USE mpiboundary

  IMPLICIT NONE

CONTAINS

  !****************************************************************************
  ! Set up any necessary variables for the chosen boundary conditions
  !****************************************************************************

  SUBROUTINE set_boundary_conditions

    any_open = .FALSE.
    IF (xbc_min == BC_OPEN .OR. xbc_max == BC_OPEN &
        .OR. ybc_min == BC_OPEN .OR. ybc_max == BC_OPEN &
        .OR. zbc_min == BC_OPEN .OR. zbc_max == BC_OPEN) any_open = .TRUE.

  END SUBROUTINE set_boundary_conditions


  !****************************************************************************
  ! Call all of the boundaries needed by the core Lagrangian solver
  !****************************************************************************

  SUBROUTINE boundary_conditions

    CALL bfield_bcs
    CALL energy_bcs
    CALL density_bcs
    CALL velocity_bcs
    CALL damp_boundaries

  END SUBROUTINE boundary_conditions


  !****************************************************************************
  ! Boundary conditions for magnetic field through plane
  !****************************************************************************

  SUBROUTINE bfield_bcs

    CALL bfield_mpi

    ! We have used custom boundary conditions that force the ghost cells
    ! to vary relative to their initial values by the same amount as
    ! the mirrored cells within the domain.

    ! At the lower boundary we drive a shear Alfven wave in by with an amplitude a0,
    ! frequency components at multiples 1,3 an 5 times that of the frequency omega
    ! and a rampup period of t0.

   IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      bx(-1,:,:) = bx0(-1,:,:)+(bx(1,:,:)-bx0(1,:,:))
      bx(-2,:,:) = bx0(-2,:,:)+(bx(2,:,:)-bx0(2,:,:))
      by( 0,:,:) = by0( 0,:,:)+(by(1,:,:)-by0(1,:,:))
      by(-1,:,:) = by0(-1,:,:)+(by(2,:,:)-by0(2,:,:))
      bz( 0,:,:) = bz0( 0,:,:)+(bz(1,:,:)-bz0(1,:,:))
      bz(-1,:,:) = bz0(-1,:,:)+(bz(2,:,:)-bz0(2,:,:))
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      bx(nx+1,:,:) = bx0(nx+1,:,:)+(bx(nx-1,:,:)-bx0(nx-1,:,:))
      bx(nx+2,:,:) = bx0(nx+2,:,:)+(bx(nx-2,:,:)-bx0(nx-2,:,:))
      by(nx+1,:,:) = by0(nx+1,:,:)+(by(nx  ,:,:)-by0(nx  ,:,:))
      by(nx+2,:,:) = by0(nx+2,:,:)+(by(nx-1,:,:)-by0(nx-1,:,:))
      bz(nx+1,:,:) = bz0(nx+1,:,:)+(bz(nx  ,:,:)-bz0(nx  ,:,:))
      bz(nx+2,:,:) = bz0(nx+2,:,:)+(bz(nx-1,:,:)-bz0(nx-1,:,:))
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      bx(:, 0,:) = bx0(:, 0,:)+(bx(:,1,:)-bx0(:,1,:))
      bx(:,-1,:) = bx0(:,-1,:)+(bx(:,2,:)-bx0(:,2,:))
      by(:,-1,:) = by0(:,-1,:)+(by(:,1,:)-by0(:,1,:))
      by(:,-2,:) = by0(:,-2,:)+(by(:,2,:)-by0(:,2,:))
      bz(:, 0,:) = bz0(:, 0,:)+(bz(:,1,:)-bz0(:,1,:))
      bz(:,-1,:) = bz0(:,-1,:)+(bz(:,2,:)-bz0(:,2,:))
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      bx(:,ny+1,:) = bx0(:,ny+1,:)+(bx(:,ny  ,:)-bx0(:,ny  ,:))
      bx(:,ny+2,:) = bx0(:,ny+2,:)+(bx(:,ny-1,:)-bx0(:,ny-1,:))
      by(:,ny+1,:) = by0(:,ny+1,:)+(by(:,ny-1,:)-by0(:,ny-1,:))
      by(:,ny+2,:) = by0(:,ny+2,:)+(by(:,ny-2,:)-by0(:,ny-2,:))
      bz(:,ny+1,:) = bz0(:,ny+1,:)+(bz(:,ny  ,:)-bz0(:,ny  ,:))
      bz(:,ny+2,:) = bz0(:,ny+2,:)+(bz(:,ny-1,:)-bz0(:,ny-1,:))
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      bx(:,:,-1) = bx0(:,:,-1)+(bx(:,:,1)-bx0(:,:,1))
      bx(:,:,-2) = bx0(:,:,-2)+(bx(:,:,2)-bx0(:,:,2))
        DO ix = -1,nx+2
         DO iy = -2, ny+2
            by(ix,iy,-1:0) = by0(ix,iy,-1:0)- &
            a0*(sin(omega*time)+sin(3.0*omega*time)+sin(5.0*omega*time))* &
            (1-exp(-(time/t0)**3))*sqrt(rho0(ix,iy,-1:0))
         ENDDO
      ENDDO  
      bz(:,:,-1) = bz0(:,:,-1)+(bz(:,:,1)-bz0(:,:,1))
      bz(:,:,-2) = bz0(:,:,-2)+(bz(:,:,2)-bz0(:,:,2))
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      bx(:,:,nz+1) = bx0(:,:,nz+1)+(bx(:,:,nz  )-bx0(:,:,nz  ))
      bx(:,:,nz+2) = bx0(:,:,nz+2)+(bx(:,:,nz-1)-bx0(:,:,nz-1))
      by(:,:,nz+1) = by0(:,:,nz+1)+(by(:,:,nz  )-by0(:,:,nz  ))
      by(:,:,nz+2) = by0(:,:,nz+2)+(by(:,:,nz-1)-by0(:,:,nz-1))
      bz(:,:,nz+1) = bz0(:,:,nz+1)+(bz(:,:,nz-1)-bz0(:,:,nz-1))
      bz(:,:,nz+2) = bz0(:,:,nz+2)+(bz(:,:,nz-2)-bz0(:,:,nz-2))
    END IF

  END SUBROUTINE bfield_bcs


  !****************************************************************************
  ! Boundary conditions for specific internal energy
  !****************************************************************************

  SUBROUTINE energy_bcs

    CALL energy_mpi

    ! We used fixed boundary conditions for the internal energy.

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      energy( 0,:,:) = energy0( 0,:,:)
      energy(-1,:,:) = energy0(-1,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      energy(nx+1,:,:) = energy0(nx+1,:,:)
      energy(nx+2,:,:) = energy0(nx+2,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      energy(:, 0,:) = energy0(:, 0,:)
      energy(:,-1,:) = energy0(:,-1,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      energy(:,ny+1,:) = energy0(:,ny+1,:)
      energy(:,ny+2,:) = energy0(:,ny+2,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      energy(:,:, 0) = energy0(:,:, 0)
      energy(:,:,-1) = energy0(:,:,-1)
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      energy(:,:,nz+1) = energy0(:,:,nz+1)
      energy(:,:,nz+2) = energy0(:,:,nz+2)
    END IF

  END SUBROUTINE energy_bcs


  !****************************************************************************
  ! Boundary conditions for density
  !****************************************************************************

  SUBROUTINE density_bcs

    CALL density_mpi

    ! We used fixed boundary conditions for the density.

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      rho( 0,:,:) = rho0( 0,:,:)
      rho(-1,:,:) = rho0(-1,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      rho(nx+1,:,:) = rho0(nx+1,:,:)
      rho(nx+2,:,:) = rho0(nx+2,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      rho(:, 0,:) = rho0(:, 0,:)
      rho(:,-1,:) = rho0(:,-1,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      rho(:,ny+1,:) = rho0(:,ny+1,:)
      rho(:,ny+2,:) = rho0(:,ny+2,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      rho(:,:, 0) = rho0(:,:, 0)
      rho(:,:,-1) = rho0(:,:,-1)
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      rho(:,:,nz+1) = rho0(:,:,nz+1)
      rho(:,:,nz+2) = rho0(:,:,nz+2)
    END IF

  END SUBROUTINE density_bcs


  !****************************************************************************
  ! Full timestep velocity boundary conditions
  !****************************************************************************

  SUBROUTINE velocity_bcs

    CALL velocity_mpi

    ! We use reflective boundary conditions for the velocity.

    ! At the lower boundary we drive a shear Alfven wave in vy with an amplitude a0,
    ! frequency components at multiples 1,3 an 5 times that of the frequency omega
    ! and a rampup period of t0.

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      vx(-2,:,:) = vx(2,:,:)
      vx(-1,:,:) = vx(1,:,:)      
      vy(-2,:,:) = vy(2,:,:)
      vy(-1,:,:) = vy(1,:,:)      
      vz(-2,:,:) = vz(2,:,:)
      vz(-1,:,:) = vz(1,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      vx(nx+1,:,:) = vx(nx-1,:,:)
      vx(nx+2,:,:) = vx(nx-2,:,:)
      vy(nx+1,:,:) = vy(nx-1,:,:)
      vy(nx+2,:,:) = vy(nx-2,:,:)
      vz(nx+1,:,:) = vz(nx-1,:,:)
      vz(nx+2,:,:) = vz(nx-2,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      vx(:,-2,:) = vx(:,2,:)
      vx(:,-1,:) = vx(:,1,:)
      vy(:,-2,:) = vy(:,2,:)
      vy(:,-1,:) = vy(:,1,:)
      vz(:,-2,:) = vz(:,2,:)
      vz(:,-1,:) = vz(:,1,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      vx(:,ny+1,:) = vx(:,ny-1,:)
      vx(:,ny+2,:) = vx(:,ny-2,:)
      vy(:,ny+1,:) = vy(:,ny-1,:)
      vy(:,ny+2,:) = vy(:,ny-2,:)
      vz(:,ny+1,:) = vz(:,ny-1,:)
      vz(:,ny+2,:) = vz(:,ny-2,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      DO ix = -2,nx+2
         DO iy = -2, ny+2
            vy(ix,iy,-2:0) = a0*(sin(omega*time)+sin(3.0*omega*time)+sin(5.0*omega*time))* &
            (1-exp(-(time/t0)**3))
         ENDDO
      ENDDO
      vx(:,:,-2:0) = 0
      vz(:,:,-2:0) = 0
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      vx(:,:,nz+1) = vx(:,:,nz-1)
      vx(:,:,nz+2) = vx(:,:,nz-2)
      vy(:,:,nz+1) = vy(:,:,nz-1)
      vy(:,:,nz+2) = vy(:,:,nz-2)
      vz(:,:,nz+1) = vz(:,:,nz-1)
      vz(:,:,nz+2) = vz(:,:,nz-2)
    END IF

  END SUBROUTINE velocity_bcs


  !****************************************************************************
  ! Half timestep velocity boundary conditions
  !****************************************************************************

  SUBROUTINE remap_v_bcs

    CALL remap_v_mpi

    ! We use reflective boundary conditions for the half timestep velocity.

    ! At the lower boundary we drive a shear Alfven wave in vy1 with an amplitude a0,
    ! frequency components at multiples 1,3 an 5 times that of the frequency omega
    ! and a rampup period of t0.

    IF (proc_x_min == MPI_PROC_NULL .AND. xbc_min == BC_USER) THEN
      vx1(-2,:,:) = vx1(2,:,:)
      vx1(-1,:,:) = vx1(1,:,:)      
      vy1(-2,:,:) = vy1(2,:,:)
      vy1(-1,:,:) = vy1(1,:,:)      
      vz1(-2,:,:) = vz1(2,:,:)
      vz1(-1,:,:) = vz1(1,:,:)
    END IF

    IF (proc_x_max == MPI_PROC_NULL .AND. xbc_max == BC_USER) THEN
      vx1(nx+1,:,:) = vx1(nx-1,:,:)
      vx1(nx+2,:,:) = vx1(nx-2,:,:)
      vy1(nx+1,:,:) = vy1(nx-1,:,:)
      vy1(nx+2,:,:) = vy1(nx-2,:,:)
      vz1(nx+1,:,:) = vz1(nx-1,:,:)
      vz1(nx+2,:,:) = vz1(nx-2,:,:)
    END IF

    IF (proc_y_min == MPI_PROC_NULL .AND. ybc_min == BC_USER) THEN
      vx1(:,-2,:) = vx1(:,2,:)
      vx1(:,-1,:) = vx1(:,1,:)
      vy1(:,-2,:) = vy1(:,2,:)
      vy1(:,-1,:) = vy1(:,1,:)
      vz1(:,-2,:) = vz1(:,2,:)
      vz1(:,-1,:) = vz1(:,1,:)
    END IF

    IF (proc_y_max == MPI_PROC_NULL .AND. ybc_max == BC_USER) THEN
      vx1(:,ny+1,:) = vx1(:,ny-1,:)
      vx1(:,ny+2,:) = vx1(:,ny-2,:)
      vy1(:,ny+1,:) = vy1(:,ny-1,:)
      vy1(:,ny+2,:) = vy1(:,ny-2,:)
      vz1(:,ny+1,:) = vz1(:,ny-1,:)
      vz1(:,ny+2,:) = vz1(:,ny-2,:)
    END IF

    IF (proc_z_min == MPI_PROC_NULL .AND. zbc_min == BC_USER) THEN
      DO ix = -2,nx+2
         DO iy = -2, ny+2
            vy1(ix,iy,-2:0) = a0*(sin(omega*time)+sin(3.0*omega*time)+sin(5.0*omega*time))* &
            (1-exp(-(time/t0)**3))
         ENDDO
      ENDDO
      vx1(:,:,-2:0) = 0
      vz1(:,:,-2:0) = 0
    END IF

    IF (proc_z_max == MPI_PROC_NULL .AND. zbc_max == BC_USER) THEN
      vx1(:,:,nz+1) = vx1(:,:,nz-1)
      vx1(:,:,nz+2) = vx1(:,:,nz-2)
      vy1(:,:,nz+1) = vy1(:,:,nz-1)
      vy1(:,:,nz+2) = vy1(:,:,nz-2)
      vz1(:,:,nz+1) = vz1(:,:,nz-1)
      vz1(:,:,nz+2) = vz1(:,:,nz-2)
    END IF

  END SUBROUTINE remap_v_bcs


  !****************************************************************************
  ! Damped boundary conditions
  !****************************************************************************

  SUBROUTINE damp_boundaries

    ! We have applied exponential damping to the velocity in the upper domain,
    ! the damping begins at a height d.

    REAL(num) :: a, d

    IF (.NOT.damping) RETURN

      d = 10.0_num
      DO iz = -1, nz + 1
        DO iy = -1, ny + 1
          DO ix = -1, nx + 1
            IF (zb(iz) > d) THEN
              a = dt * (zb(iz) - d) / (z_max - d) +1.0_num
              vx(ix,iy,iz) = vx(ix,iy,iz) / a
              vy(ix,iy,iz) = vy(ix,iy,iz) / a
              vz(ix,iy,iz) = vz(ix,iy,iz) / a
            END IF
          END DO
        END DO
      END DO

  END SUBROUTINE damp_boundaries

END MODULE boundary
