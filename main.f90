module shock_utils
  implicit none
contains
  ! Function to compute Q value
  function Q_value(Nx, Ny, rho1, u1, v1, e1, gamma2) result(Q_val)
    implicit none
    integer, intent(in) :: Nx, Ny
    real, intent(in) :: rho1(Ny,Nx), u1(Ny,Nx), v1(Ny,Nx), e1(Ny,Nx)
    real, intent(in) :: gamma2
    real :: Q_val(Ny,Nx,4)

    Q_val(:,:,1) = rho1(:,:)
    Q_val(:,:,2) = rho1(:,:) * u1(:,:)
    Q_val(:,:,3) = rho1(:,:) * v1(:,:)
    Q_val(:,:,4) = rho1(:,:) * (e1(:,:) + 0.5 * (u1(:,:)**2 + v1(:,:)**2))
  end function Q_value

  function F_value(Nx,Ny, rho1, u1, v1, e1, p1, gamma2) result(F_val)
        implicit none
        integer, intent(in) :: Nx,Ny
        real, intent(in) :: rho1(Ny,Nx), u1(Ny,Nx),v1(Ny, Nx), e1(Ny,Nx), p1(Ny,Nx)
        real, intent(in) :: gamma2
        real ::F_val(Ny,Nx,4)

        F_val(:,:,1) = rho1(:,:) * u1(:,:)
        F_val(:,:,2) = rho1(:,:) * u1(:,:)**2 + p1(:,:)
        F_val(:,:,3) = rho1(:,:) * u1(:,:) * v1(:,:)
        F_val(:,:,4) = rho1(:,:) * u1(:,:) * (e1(:,:) + 0.5 * (u1(:,:)**2 + v1(:,:)**2)) + p1(:,:) * u1(:,:)
    end function F_value

    !Function to compute G value
   function G_value(Nx,Ny, rho1, u1, v1, e1, p1, gamma2) result(G_val)
        implicit none
        integer, intent(in) :: Nx,Ny
        real, intent(in) :: rho1(Ny,Nx), u1(Ny,Nx),v1(Ny, Nx), e1(Ny,Nx), p1(Ny,Nx)
        real, intent(in) :: gamma2
        real ::G_val(Ny,Nx,4)

        G_val(:,:,1) = rho1(:,:) * v1(:,:)
        G_val(:,:,2) = rho1(:,:) * u1(:,:) * v1(:,:)
        G_val(:,:,3) = rho1(:,:) * v1(:,:)**2 + p1(:,:)
        G_val(:,:,4) = rho1(:,:) * v1(:,:) * (e1(:,:) + 0.5 * (u1(:,:)**2 + v1(:,:)**2)) + p1(:,:) * v1(:,:)
    end function G_value

     !Function to compute A_roe
    function A_roe(ul,ur,vl,vr,rhol,rhor,hl,hr,al,ar) result(A_flux)
    implicit none

    real,intent(in)::rhor,ur,hr,rhol,ul,hl,vl,vr,al,ar
    real::u_av,Rho_av,alpha,E_av,v_av,h_av,a_av
    real,parameter::gamma1 = 1.4
    real::A_flux(4,4),P(4,4),Lambda(4,4),P_inv(4,4),A_r1(4,4)
    Integer::j,i
    	alpha = sqrt(rhol) / (sqrt(rhol) + sqrt(rhor))
    	rho_av = (1 - alpha) * rhol + alpha * rhor
    	u_av = alpha * ul + (1 - alpha) * ur
    	v_av = alpha * vl + (1- alpha) * vr
    	a_av = alpha * al + (1 - alpha) * ar
    	h_av = alpha * hl + (1 - alpha) * hr

        !Making a matrix P
        P(1,:) =(/1.0,1.0,1.0,0.0/)
        P(2,:) = (/u_av - a_av, u_av , u_av + a_av, 0.0 /)
        P(3,:) = (/v_av,v_av,v_av,-1.0/)
        P(4,:) = (/h_av - u_av * a_av,0.5 * (u_av**2 + a_av**2), h_av + u_av * a_av,-v_av/)

        !Making a matrix lambda
        Lambda(1,:) = (/u_av - a_av,0.0,0.0, 0.0/)
        Lambda(2,:) = (/0.0,u_av,0.0, 0.0/)
        Lambda(3,:) = (/0.0,0.0,u_av + a_av, 0.0/)
        Lambda(4,:) = (/0.0,0.0,0.0, u_av/)

        Lambda = abs(lambda)

        !Making a matrix P_inv
        P_inv(1,:) = (/((gamma1 - 1) * 0.5 * (u_av**2 + v_av**2) + a_av * u_av )/(2*a_av**2),&
                            ((1 - gamma1) * u_av - a_av )/ (2 * a_av**2),&
                            ((1 - gamma1) * v_av)/(2 * a_av**2),&
                            (gamma1 - 1)/(2 * a_av**2) /)
        P_inv(2,:) = (/ (a_av**2 - (gamma1 - 1) * 0.5 * (u_av**2 + v_av**2)) / (a_av**2), &
                            ((gamma1 - 1) * u_av) / a_av**2,&
                            ((gamma1 - 1) * v_av) / a_av**2, (1 - gamma1) / a_av**2/)
        P_inv(3,:) = (/((gamma1 - 1) * 0.5 * (u_av**2 + v_av**2) - a_av * u_av )/(2*a_av**2),&
                                ((1 - gamma1) * u_av + a_av) / (2*a_av**2),&
                            (1- gamma1) * v_av / (2 * a_av**2), (gamma1 - 1) / (2 * a_av**2)/)
        P_inv(4,:) = (/0.0,0.0,-1.0,0.0/)

    	!Making Array for A_roe
    	A_r1 = Matmul(P,lambda)
        A_flux = Matmul(A_r1,P_inv)

  end function A_roe

  !Function to compute B_roe
    function B_roe(ul,ur,vl,vr,rhol,rhor,hl,hr,al,ar) result(B_flux)
    implicit none

    real,intent(in)::rhor,ur,hr,rhol,ul,hl,vl,vr,al,ar
    real::u_av,Rho_av,R,E_av,v_av,alpha,h_av,a_av
    real,parameter::gamma1 = 1.4
    real::B_flux(4,4),P(4,4),Lambda(4,4),P_inv(4,4),B_r1(4,4)
    Integer::j,i


    	alpha = sqrt(rhol) / (sqrt(rhol) + sqrt(rhor))
    	rho_av = (1 - alpha) * rhol + alpha * rhor
    	u_av = alpha * ul + (1 - alpha) * ur
    	v_av = alpha * vl + (1- alpha) * vr
    	h_av = alpha * hl + (1 - alpha) * hr
    	a_av = alpha * al + (1 - alpha) * ar
    	E_av = (h_av + 0.5 * (gamma1 - 1) * (u_av**2 + v_av**2)) / gamma1


    	!Making Array for eigen value matrix
    	  !Making a matrix P
        P(1,:) =(/1.0,1.0,1.0,0.0/)
        P(2,:) = (/u_av, u_av , u_av, 1.0 /)
        P(3,:) = (/v_av - a_av,v_av,v_av + a_av, 0.0/)
        P(4,:) = (/h_av - v_av * a_av, 0.5 * (u_av**2 + a_av**2), h_av + v_av * a_av, u_av/)

 !Making a matrix lambda
        Lambda(1,:) = (/v_av - a_av,0.0,0.0, 0.0/)
        Lambda(2,:) = (/0.0,v_av,0.0, 0.0/)
        Lambda(3,:) = (/0.0,0.0,v_av + a_av, 0.0/)
        Lambda(4,:) = (/0.0,0.0,0.0, v_av/)

        Lambda = abs(Lambda)

 !Making a matrix P_inv
        P_inv(1,:) = (/(gamma1 - 1) * (0.5 * (u_av**2 + v_av**2) + a_av * v_av )/(2*a_av**2),&
                            ((1 - gamma1) * u_av )/ (2 * a_av**2),&
                            ((1 - gamma1) * v_av - a_av)/(2 * a_av**2),&
                            (gamma1 - 1) / (2.0 * a_av**2) /)
        P_inv(2,:) = (/ (a_av**2 - (gamma1 - 1) * 0.5 * (u_av**2 + v_av**2)) / (a_av**2),&
                            ((gamma1 - 1) * u_av) / a_av**2,&
                            ((gamma1 - 1) * v_av) / a_av**2,&
                            (1 - gamma1) / a_av**2/)
        P_inv(3,:) = (/((gamma1 - 1) * 0.5 * (u_av**2 + v_av**2) - a_av * v_av )/(2*a_av**2),&
                            (1- gamma1) * u_av / (2 * a_av**2),&
                            ((1 - gamma1) * v_av + a_av) / (2*a_av**2),&
                            (gamma1 - 1) / (2 * a_av**2)/)
        P_inv(4,:) = (/-u_av,1.0,0.0,0.0/)

    B_r1 = Matmul(P,lambda)
    B_flux = Matmul(B_r1,P_inv)

  end function B_roe

end module shock_utils

! Main Program starts here --------------------------------------------------------------------------------------------------
program ROE_Wedge_equation
  use shock_utils
  implicit none

  ! Parameters
  real, parameter :: gamma1 = 1.4, CFL = 0.5, L = 1.0
  real :: max_a, max_u
  integer, parameter :: Nx = 201, Ny = 201, Nt = 2300

  ! Variables
  real :: x_ini, x_fin, dt, y_ini, y_fin,conv,C_d1
  real :: dx, dy, R_const, u_tot, theta, Mach,P_sum,sum1
  integer :: i, j, k, m,z
  real, dimension(:), allocatable :: x, y
  real, dimension(:,:), allocatable :: rho, p, u, v, e, a, T, E_t,h,u_tot1,rho_old,Mach1,P0,T0,R0,ent,us,vs,div,C_p,C_d
  real, dimension(:,:,:), allocatable :: Q, F, G, Q_n1, F_jph,G_jph, Q1, Q2,Q_f
  real, dimension(4,4) :: A_r1, B_r1,A_r2,B_r2
  real, dimension(4) :: dq1, dq2, final_value1, final_value2,dq3,dq4,final_value3,final_value4

  ! Domain setup
  dx = L / (Nx - 1)
  dy = L / (Ny - 1)
  x_ini = 0.0
  x_fin = 1.0
  y_ini = x_ini
  y_fin = x_fin
  R_const = 287.0

  ! Allocate arrays
  allocate(x(Nx), y(Ny))
  allocate(rho(Ny,Nx), p(Ny,Nx), u(Ny,Nx), v(Ny,Nx), e(Ny,Nx), a(Ny,Nx), T(Ny,Nx), E_t(Ny,Nx),h(Ny,Nx))
  allocate(u_tot1(Ny,Nx),rho_old(Ny,Nx),Mach1(Ny,Nx),p0(Ny,Nx),T0(Ny,Nx),R0(Ny,Nx),ent(Ny,Nx),us(Ny,Nx),vs(Ny,Nx))
  allocate(div(Ny,Nx),C_p(Ny,Nx),C_d(Ny,Nx))
  allocate(Q(Ny,Nx,4), F(Ny,Nx,4),G(Ny,Nx,4), Q_n1(Ny,Nx,4), F_jph(Ny,Nx,4), G_jph(Ny,Nx,4))
  allocate(Q1(Ny,Nx,4), Q2(Ny,Nx,4),Q_f(Ny,Nx,4))

  ! Create x and y grid
  do i = 1, Nx
    x(i) = x_ini + dx * (i - 1)
  end do

  do j = 1, Ny
    y(j) = y_ini + dy * (j - 1)
  end do

  ! Initial conditions
  theta = 15.0 * 22.0 / (7.0 * 180.0)           ! converting degrees to radians
  Mach = 2.0
  p(:,:) = 101325.0
  T(:,:) = 298.0
  rho(:,:) = p(:,:) / (R_const * T(:,:))
  a(:,:) = sqrt(gamma1 * R_const * T(:,:))
  u_tot = Mach * maxval(a)
  u(:,:) = u_tot * cos(theta)
  v(:,:) = -u_tot * sin(theta)
  e(:,:) = p(:,:) / (rho(:,:) * (gamma1 - 1.0))
  E_t(:,:) = 0.5 * (u(:,:)**2 + v(:,:)**2) + e(:,:)
  h(:,:) = gamma1 * E_t(:,:) - 0.5 * (gamma1 - 1) * (u(:,:)**2 + v(:,:)**2)

  ! Initialize Q
  Q_f(:,:,1) = rho(:,:)
  Q_f(:,:,2) = rho(:,:) * u(:,:)
  Q_f(:,:,3) = rho(:,:) * v(:,:)
  Q_f(:,:,4) = rho(:,:) * E_t(:,:)


!open a txt file to write the output
open(unit = 11, file = "output_ROE.txt")
open(unit = 12, file="output_sgen.csv")
open(unit = 13, file = "shock_width.txt")
open(unit = 14, file ="divergence.txt")

do k = 1,Nt
  ! Update speed of sound and time step
  a(:,:) = sqrt(gamma1 * p(:,:) / rho(:,:))
  dt = CFL * dx / (maxval(abs(u(:,:)) + a(:,:)) + maxval(abs(v(:,:)) + a(:,:)))

  ! Calculate Q using the function
  Q = Q_value(Nx, Ny, rho, u, v, e, gamma1)
  F = F_value(Nx,Ny, rho, u,v, e, p, gamma1)
  G = G_value(Nx,Ny, rho, u,v, e, p, gamma1)

do i = 2,Ny-1
    do j = 2,Nx-1
        A_r1 = A_roe(u(i,j),u(i,j+1),v(i,j),v(i,j+1),rho(i,j),rho(i,j+1),h(i,j),h(i,j+1),a(i,j),a(i,j+1))
        A_r2 = A_roe(u(i,j-1),u(i,j),v(i,j-1),v(i,j),rho(i,j-1),rho(i,j),h(i,j-1),h(i,j),a(i,j-1),a(i,j))
       ! if (i==2 .and. j==2) then
           ! do m = 1,4
           ! print*,(A_r1(m,z),z=1,4)
          !  end do
        !end if
        B_r1 = B_roe(u(i,j),u(i+1,j),v(i,j),v(i+1,j),rho(i,j),rho(i+1,j),h(i,j),h(i+1,j),a(i,j),a(i+1,j))
        B_r2 = B_roe(u(i-1,j),u(i,j),v(i-1,j),v(i,j),rho(i-1,j),rho(i,j),h(i-1,j),h(i,j),a(i-1,j),a(i,j))
        dq1 =  Q(i,j,:) - Q(i,j+1,:)
        dq2 =  Q(i,j,:) - Q(i+1,j,:)
        dq3 =  Q(i,j-1,:) - Q(i,j,:)
        dq4 =  Q(i-1,j,:) - Q(i,j,:)
        Final_value1 = Matmul(A_r1,dq1)
        Final_value2 = Matmul(B_r1,dq2)
        Final_value3 = Matmul(A_r2,dq3)
        Final_value4 = Matmul(B_r2,dq4)
        F_jph (i,j,:) = 0.5 * (F(i,j,:) + F(i,j+1,:)) + 0.5 * Final_value1(:)
  	    F_jph (i,j-1,:) = 0.5 * (F(i,j-1,:) + F(i,j,:)) + 0.5 * Final_Value3(:)
  	    G_jph (i,j,:) = 0.5 * (G(i,j,:) + G(i+1,j,:)) + 0.5 * Final_value2(:)
  	    G_jph (i-1,j,:) = 0.5 * (G(i-1,j,:) + G(i,j,:)) + 0.5 * Final_Value4(:)

  	    Q_n1(i,j,:) = Q(i,j,:) - dt / dx * (F_jph(i,j,:) - F_jph(i,j-1,:)) - dt / dy * (G_jph(i,j,:) - G_jph(i-1,j,:))                                       !Main equation
    end do
end do

!Applying the Boundary conditions
    !-------------------------Top BC--------------------------
    Q_n1(Ny,:,1) = Q_f(Ny,:,1)
    Q_n1(Ny,:,2) = Q_f(Ny,:,2)
    Q_n1(Ny,:,3) = Q_f(Ny,:,3)
    Q_n1(Ny,:,4) = Q_f(Ny,:,4)

    !------------------------Inlet Boundary conditions
    Q_n1(2:Ny-1,:1,1) = Q_f(2:Ny-1,:1,1)
    Q_n1(2:Ny-1,:1,2) = Q_f(2:Ny-1,:1,2)
    Q_n1(2:Ny-1,:1,3) = Q_f(2:Ny-1,:1,3)
    Q_n1(2:Ny-1,:1,4) = Q_f(2:Ny-1,:1,4)

    !-----------------------Right Boundary conditions
    Q_n1(:,Nx,1) = Q_n1(:,Nx-1,1)
    Q_n1(:,Nx,2) = Q_n1(:,Nx-1,2)
    Q_n1(:,Nx,3) = Q_n1(:,Nx-1,3)
    Q_n1(:,Nx,4) = Q_n1(:,Nx-1,4)

    !----------------------Bottom Boundary conditions
    Q_n1(1,:,1) = Q_n1(2,:,1)
    Q_n1(1,:,2) = Q_n1(2,:,2)
    Q_n1(1,:,3) = 0
    Q_n1(1,:,4) = Q_n1(2,:,4)



    Q(:,:,:) = Q_n1(:,:,:)		!copying all values from n+1 time step to n




    rho(:,:) = Q(:,:,1)
	u(:,:) = Q(:,:,2) / rho(:,:)
	v(:,:) = Q(:,:,3) / rho(:,:)
	E_t(:,:) = Q(:,:,4) / rho(:,:)
	e(:,:) = E_t(:,:) - 0.5 * (u(:,:)**2 + v(:,:)**2)
	h(:,:) = gamma1 * E_t - 0.5 * (gamma1 - 1) * (u(:,:)**2 + v(:,:)**2)
	p(:,:) = rho(:,:) * e(:,:) * (gamma1 - 1)
	T(:,:) = p(:,:) / (rho(:,:) * R_const)
	a(:,:) = sqrt(gamma1 * R_const * T(:,:))
    u_tot1(:,:) = sqrt(u(:,:)**2 + v(:,:)**2)
    Mach1(:,:) = u_tot1(:,:) / a(:,:)
    P0(:,:) = p(:,:) * (1 + 0.5 * (gamma1 - 1) * Mach1(:,:)**2)**(gamma1 / (gamma1 - 1))
    T0(:,:) = T(:,:) * ( 1 + 0.5 * (gamma1 - 1) * Mach1**2)
    R0(:,:) = P0(:,:) / (R_const * T0(:,:))
    ent(:,:) = rho(:,:) *  (log(p(:,:)/101325.0) - gamma1 * log(rho(:,:) / 1.184727451)) / (gamma1 - 1)
    C_p(:,:) = (p(:,:) - 101325.0)/(0.5 * 1.184727451 * 2.0**2 *gamma1 * R_const * 298.0)
    C_d(:,:) = 4 * (p(:,:) / 101325.0 - 1.0) * tan(theta) / (gamma1 * 2.0**2)
    us = u(:,:) * ent(:,:)
    vs = v(:,:) * ent(:,:)




   do i =1,Ny
        do j =1,Nx
           if (us(i,j) < 0) then
                    us(i,j) = 0
                else
                    us(i,j) = us(i,j)
           end if
              if (vs(i,j) < 0) then
                    vs(i,j) = 0
                else
                    vs(i,j) = vs(i,j)
           end if
        end do

   end do

    If (k == Nt -1) then
        rho_old = rho
    End If

!writing data to a file
                if (k == Nt) then
                    do i=1,Ny
                        do j = 1,Nx
                           write(11,*)x(j),y(i),p(i,j)
                           end do
                        write(11,*)""
                        write(11,*)""
                    end do
                end if

                 if (k == Nt) then
                    do i=1,Ny
                        do j = 1,Nx
                           write(12,*)x(j),y(i),us(i,j),vs(i,j)
                           end do
                        write(12,*)""
                        write(12,*)""
                    end do
                end if

!Writing data file for Shock_width Calculation
 !if (k == Nt) then
                    !do i=1,Ny
                        !do j = 1,Nx
                           !write(11,*)x(j),y(i),p(i,j)
                           !end do
                        !write(11,*)""
                        !write(11,*)""
                    !end do
!end if

 if (k == Nt) then
                        do j = 1,Nx
                           write(13,*)x(j),p(31,j),Mach1(31,j)
                        end do
end if


end do

!Computation of divergence----------------------------------------------------------------------------

!For interior points
do i=2,Ny-1
    do j=2,Nx-1
        div(i,j) = (us(i,j+1) - us(i,j-1) ) / (2 * dx) + (vs(i+1,j) - vs(i-1,j) ) / (2 * dy)
    end do
end do

!For left boundary
do i =2,Ny-1
    div(i,1) = (us(i,2) - us(i,1) ) /  dx + (vs(i+1,1) - vs(i,1) ) /  dy
end do

!for right boundary
do i =2,Ny-1
    div(i,Nx) = (us(i,Nx) - us(i,Nx-1) ) /  dx + (vs(i+1,Nx) - vs(i,Nx) ) /  dy
end do

!for Upper boundary
do j =2,Nx-1
    div(Ny,j) = (us(Ny,j+1) - us(Ny,j) ) /  dx + (vs(Ny,j) - vs(Ny-1,j) ) /  dy
end do

!For Lower boundary
do j =2,Nx-1
    div(1,j) = (us(1,j+1) - us(1,j) ) /  dx + (vs(2,j) - vs(1,j) ) /  dy
end do

!For corner
div(1,1) = (us(1,2) - us(1,1)) / dx + (vs(2,1) - vs(1,1)) / dy
div(Ny,1) = (us(Ny,2) - us(Ny,1)) / dx + (vs(Ny,2) - vs(Ny,1)) / dy
div(1,Nx) =  (us(1,Nx) - us(1,Nx-1)) / dx + (vs(2,Nx) - vs(1,Nx)) / dy
div(Ny,Nx) = (us(Ny,Nx) - us(Ny,Nx-1)) / dx + (vs(Ny,Nx) - vs(Ny-1,Nx)) / dy

!divergence computation ends here---------------------------------------------------------------------------------------

do i=1,Ny
    do j=1,Nx
        write(14,*)x(j),y(i),div(i,j)
    end do
end do


 conv = sqrt(sum((rho_old - rho)**2 )/ size(rho))

Print*,conv
print*,Mach1(20,81)
do i=1,Nx
    P_sum = P_sum + P(1,i) * dx * sin(theta)
end do

C_d1 = P_sum / (0.5 * 1.184727451 * 2.0**2 * gamma1 * R_const * 298.0)
print*,"C_d=",C_d1

!call system('gnuplot -p graph.plt')
!call system('gnuplot -p graph_div.plt')
close(11)
close(12)
close(13)
close(14)
end program ROE_Wedge_equation


