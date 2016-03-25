PROGRAM supernova

! AUTHOR: Federico Esposito
!   DATE: 2016/03/08
!
!   grid xa -> m, v
!   grid xb -> d, p, e, q


IMPLICIT NONE

! precision parameter (use 4 or 8)
 INTEGER, PARAMETER :: b = 8
 REAL(KIND=b), PARAMETER :: pi = 3.141592653589793
 REAL(KIND=b), PARAMETER :: pc2cm = 3.086d18
 REAL(KIND=b), PARAMETER :: yr2sec = 3.1556926d7
 REAL(KIND=b), PARAMETER :: km2cm = 1.d5

 INTEGER :: j, n, coord, bc, order, counter

 REAL(KIND=b) :: gamma, dt, c2, K0, E0, dt_sf
 REAL(KIND=b) :: d0, v0, xmin, xmax, t_end, c_v
 REAL(KIND=b) :: delta_p, delta_q, delta_v
 REAL(KIND=b) :: alpha, v_mean, d_mean
 REAL(KIND=b) :: t, v2_bw, v2_fw, num, den
 REAL(KIND=b) :: first, second, dt_sf_i, dt_sf_f
 REAL(KIND=b), ALLOCATABLE :: xa(:), xb(:), dxa(:), dxb(:)
 REAL(KIND=b), ALLOCATABLE :: g2a(:), g2b(:), g3a(:)
 REAL(KIND=b), ALLOCATABLE :: g3b(:), dvla(:), dvlb(:)
 REAL(KIND=b), ALLOCATABLE :: dsa(:), dsb(:), sigma(:)
 REAL(KIND=b), ALLOCATABLE :: v(:), d(:), p(:), e(:), K(:)
 REAL(KIND=b), ALLOCATABLE :: dv(:), s(:), q(:), F(:)
 REAL(KIND=b), ALLOCATABLE :: delta_s(:), d_s(:)


! input file
 OPEN(11,file="sninput.dat")

! output file
 OPEN(12,file="sn_a.dat")
 OPEN(13,file="sn_b.dat")
 

! input file -> setting space & time
 READ(11,"(3(/), T53, I7)") n    !grid points for space
 READ(11,"(T53, F8.3)") xmin     !infimum of space
 READ(11,"(T53, F8.3)") xmax     !supremum of space
 READ(11,"(T53, F8.3)") t_end    !integration time
 
! input file -> initial conditions
 READ(11,"(T53, F8.3)") E0       !t=0 sn energy [erg]
 READ(11,"(T53, F8.3)") d0       !t=0 ism density [g/cm3]
 READ(11,"(T53, F8.3)") K0       !t=0 ism temperature
 READ(11,"(T53, F8.3)") v0       !t=0 ism velocity
 
! input file -> constants
 READ(11,"(T53, F8.3)") dt_sf_i  !t=0 safety factor for dt
 READ(11,"(T53, F8.3)") dt_sf_f  !t_end   "    "    for dt
 READ(11,"(T53, F8.3)") gamma    !adiabatic coefficient
 READ(11,"(T53, F8.3)") c_v      !specific heat [cgs]
 READ(11,"(T53, F8.3)") c2       !viscosity coefficient
 
! input file -> choices
 READ(11,"(T53, I5)") coord      !cartesian/spherical
 READ(11,"(T53, I5)") bc         !outflow/reflection
 READ(11,"(T53, I5)") order      !1st/2nd order for upwind
 
! arrange the measure units
 xmin = xmin*pc2cm
 xmax = xmax*pc2cm
 t_end = t_end*yr2sec
 
! allocating the arrays
 ALLOCATE (xa(n), xb(n), g2a(n), g2b(n), g3a(n), g3b(n))
 ALLOCATE (dxa(n-1), dxb(2:n), dvla(n-1), dvlb(2:n))
 ALLOCATE (dsa(n), dsb(n), delta_s(n), d_s(n), sigma(n))
 ALLOCATE (v(n), d(n), p(n), e(n), dv(n), s(n))
 ALLOCATE (q(n), F(n), K(n))
 
! defining the grid
 CALL grid (n, xmin, xmax, coord, xa, xb, dxa, dxb, &
          & g2a, g2b, g3a, g3b, dvla, dvlb, dsa, dsb)
 

! setting initial conditions
 DO j=1,n
   v(j) = v0
   d(j) = d0
   IF (j.GE.2) THEN
     s(j) = v(j)*0.5*(d(j-1) + d(j))
   END IF
 END DO
 s(1) = s(2)

 e = 0.
 e(2) = 3.*E0/(4.*pi*xb(4)**3)
 e(3) = e(2)
 
 DO j=1,n
   p(j) = (gamma - 1.)*e(j)
   K(j) = e(j)/(c_v*d(j))
 END DO
 
 
 
! ------------------ BIG TEMPORAL CYCLE ------------------
 t = 0.
 counter = 0
 dt_sf = dt_sf_i
 DO WHILE (t.LE.t_end)
   CALL delta_t (dt, n, dxa, dxb, p, d, v, gamma, c2)
   IF (dt_sf.LE.dt_sf_f) dt_sf = dt_sf + 0.1*dt_sf
   dt = dt*dt_sf
   t = t + dt
   counter = counter + 1
   
  ! ---------------- SOURCE STEP ----------------
  ! 1: update momentum for grad(pressure)
   DO j=2,n-1
     d_mean = 0.5*(d(j) + d(j-1))
     delta_p = p(j) - p(j-1)
     v(j) = v(j) - delta_p*dt/(dxb(j)*d_mean)
   END DO
   CALL bound_cond (n, v, 2)
   
  ! 2a: artificial viscosity q(j)
   DO j=2,n-1
     IF ((v(j+1) - v(j)).LT.0.) THEN
       q(j) = c2*d(j)*(v(j+1) - v(j))**2
     ELSE
       q(j) = 0.
     END IF
   END DO
   CALL bound_cond (n, q, bc)
   
  ! 2b: update velocity with viscosity
   DO j=2,n-1
     d_mean = 0.5*(d(j) + d(j-1))
     delta_q = q(j) - q(j-1)
     v(j) = v(j) - delta_q*dt/(dxb(j)*d_mean)
   END DO
   CALL bound_cond (n, v, 2)
   
  ! 2c: update energy with viscosity
   DO j=2,n-1
     delta_v = v(j+1) - v(j)
     e(j) = e(j) - q(j)*delta_v*dt/dxa(j)
   END DO
   CALL bound_cond (n, e, bc)
   
  ! 3a: velocity divergence dv(j)
   DO j=2,n-1
     v2_fw = g2a(j+1)*g3a(j+1)*v(j+1)
     v2_bw = g2a(j)*g3a(j)*v(j)
     dv(j) = (v2_fw - v2_bw)/dvla(j)
   END DO
   CALL bound_cond (n, dv, bc)
   
  ! 3b: update energy for compressional heating
   DO j=2,n-1
     alpha = 0.5*dt*(gamma - 1.)
     e(j) = e(j)*(1 - alpha*dv(j))/(1 + alpha*dv(j))
     p(j) = e(j)*(gamma - 1.)
   END DO
   CALL bound_cond (n, e, bc)
 
  ! -------------- TRANSPORT STEP ---------------
  ! 1: define momentum
   DO j=2,n-1
     s(j) = v(j)*0.5*(d(j-1) + d(j))
   END DO
   CALL bound_cond (n, s, 2)
   
  ! 2: apply upwind 1st or 2nd order
   SELECT CASE (order)
     CASE (1) ! upwind 1st order
      ! update d & e
       CALL upwind1 (d, dsa, dvla, v, dt, n, bc)
       CALL upwind1 (e, dsa, dvla, v, dt, n, bc)
       p = e*(gamma - 1.)
       
      ! momentum flux
       DO j=2,n-1
         v_mean = 0.5*(v(j) + v(j+1))
         IF (v_mean.GE.0.) THEN
           F(j) = s(j)*v_mean*dsb(j)
         ELSE
           F(j) = s(j+1)*v_mean*dsb(j)
         END IF
       END DO
       CALL bound_cond (n, F, bc)
       
      ! update momentum & velocity
       DO j=2,n-1
         s(j) = s(j) - (F(j) - F(j-1))*dt/dvlb(j)
         v(j) = s(j)*2./(d(j-1) + d(j))
       END DO
       CALL bound_cond (n, s, 2)
       CALL bound_cond (n, v, 2)
       
     CASE (2) ! upwind 2nd order (van leer method)
      ! update d & e
       CALL vanleer (d, dxa, v, dt, n, bc)
       CALL vanleer (e, dxa, v, dt, n, bc)
       
      ! update temperature K
       DO j=1,n
         K(j) = e(j)/(c_v*d(j))
       END DO
       
      ! define courant number sigma
       DO j=2,n-1
         v_mean = 0.5*(v(j+1) + v(j))
         sigma(j) = v_mean*dt/dxb(j)
       END DO
       CALL bound_cond (n, sigma, bc)
       
      ! define delta_s
       DO j=2,n-1
         delta_s(j) = s(j+1) - s(j)
       END DO
       CALL bound_cond (n, delta_s, bc)
   
      ! define d_s
       DO j=2,n-1
         num = delta_s(j-1)*delta_s(j)
         den = s(j+1) - s(j-1)
         IF (num.GT.0.) THEN
           d_s(j) = 2*num/den
         ELSE
           d_s(j) = 0.
         END IF
       END DO
       CALL bound_cond (n, d_s, bc)
       
       ! update s & v
       DO j=2,n-1
         v_mean = 0.5*(v(j) + v(j+1))
         IF (v_mean.GE.0.) THEN
           first = sigma(j)*(s(j+1) - s(j))
           second = (1. - sigma(j))*(d_s(j+1) - d_s(j))
         ELSE
           first = sigma(j)*(s(j) - s(j-1))
           second = (1. - sigma(j))*(d_s(j) - d_s(j-1))
         END IF
         s(j) = s(j) - first - (sigma(j)/2.)*second
         v(j) = s(j)*2./(d(j-1) + d(j))
       END DO
       CALL bound_cond (n, s, 2)
       CALL bound_cond (n, v, 2)
   END SELECT
 END DO
 
 11 CONTINUE
 120 FORMAT (X, 3(ES10.3, 2X))
 130 FORMAT (X, 4(ES10.3, 2X))
 DO j=1,n
   WRITE(12,120) xa(j)/pc2cm, v(j)/km2cm, s(j)
   WRITE(13,130) xb(j)/pc2cm, d(j), p(j), e(j)
 END DO
 
 
! closing files
 CLOSE(11)
 CLOSE(12)
 CLOSE(13)
 
 
 
 CONTAINS
 
SUBROUTINE grid (n, xmin, xmax, coord, xa, xb, dxa, dxb, &
               & g2a, g2b, g3a, g3b, dvla, dvlb, dsa, dsb)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, coord
   REAL(KIND=b), INTENT(IN) :: xmin, xmax
   REAL(KIND=b), INTENT(OUT) :: xa(n), xb(n), dxa(n-1)
   REAL(KIND=b), INTENT(OUT) :: dxb(2:n), g2a(n), g2b(n)
   REAL(KIND=b), INTENT(OUT) :: g3a(n), g3b(n)
   REAL(KIND=b), INTENT(OUT) :: dvla(n-1), dvlb(2:n)
   REAL(KIND=b), INTENT(OUT) :: dsa(n), dsb(n)
   REAL(KIND=b) :: xa_nplus1
   INTEGER :: j
   
  !I want xa to start on -dx
   DO j=1,n
     xa(j) = xmin + (j-2.)*(xmax - xmin)/(n-1.)
   END DO
  !I want xb to fall ahead xa
   DO j=1,n-1
     xb(j) = 0.5*(xa(j) + xa(j+1))
   END DO
   
  !I need a xa(n+1) point to define xb(n)
   xa_nplus1 = xmax
   xb(n) = 0.5*(xa(n) + xa_nplus1)
   
  !Intervals will be n-1 of course
   DO j=1,n-1
     dxa(j) = xa(j+1) - xa(j) !centered on xb(j)
   END DO
   DO j=2,n
     dxb(j) = xb(j) - xb(j-1) !centered on xa(j)
   END DO
 
  !Metric scale factors
   SELECT CASE (coord)
     CASE (1) !cartesian
       g2a = 1.
       g2b = 1.
       g3a = 1.
       g3b = 1.
       dvla = dxa
       dvlb = dxb
     CASE (2) !spherical
       g2a = xa
       g2b = xb
       g3a = xa
       g3b = xb
       DO j=1,n-1
         dvla(j) = (xa(j+1)**3 - xa(j)**3)/3.
       END DO
       DO j=2,n
         dvlb(j) = (xb(j)**3 - xb(j-1)**3)/3.
       END DO
   END SELECT
   
  !Surface elements
   DO j=1,n
     dsa(j) = g2a(j)*g3a(j)
     dsb(j) = g2b(j)*g3b(j)
   END DO
END SUBROUTINE grid


SUBROUTINE bound_cond (n, array, bc)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, bc
   REAL(KIND=b), INTENT(INOUT) :: array(n)
   
   SELECT CASE (bc)
     CASE (1) !outflow
       array(1) = array(2)
       array(n) = array(n-1)
     CASE (2) !reflection
       array(1) = -array(3)
       array(2) = 0.
       array(n) = array(n-2)
       array(n-1) = 0.       
   END SELECT
END SUBROUTINE bound_cond


SUBROUTINE delta_t (dt, n, dxa, dxb, p, d, v, gamma, c2)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n
   REAL(KIND=b), INTENT(OUT) :: dt
   REAL(KIND=b), INTENT(IN) :: dxa(n-1), dxb(2:n), c2
   REAL(KIND=b), INTENT(IN) :: p(n), d(n), v(n), gamma
   REAL(KIND=b) :: dx_min, c_s, dv(n-1)
   REAL(KIND=b) :: dt_1(2:n-1), dt_2(2:n-1)
   INTEGER :: j
   
   DO j=1,n-1
     dv(j) = ABS(v(j+1) - v(j))
   END DO
   
   DO j=2,n-1
     dx_min = MIN(dxa(j), dxb(j))
     c_s = SQRT(gamma*p(j)/d(j))
     dt_1(j) = dx_min/(c_s + ABS(v(j)))
     dt_2(j) = dx_min/(4*c2*dv(j))
   END DO
   dt = MIN(MINVAL(dt_1), MINVAL(dt_2))
END SUBROUTINE delta_t


SUBROUTINE upwind1 (q, ds, dvl, v, dt, n, bc)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, bc
   REAL(KIND=b), INTENT(INOUT) :: q(n)
   REAL(KIND=b), INTENT(IN) :: ds(n), v(n)
   REAL(KIND=b), INTENT(IN) :: dvl(n), dt
   REAL(KIND=b) :: F(n)
   INTEGER :: j
   
  ! define flux F
   DO j=2,n-1
     IF (v(j).GE.0.) THEN
       F(j) = q(j-1)*v(j)*ds(j)
     ELSE
       F(j) = q(j)*v(j)*ds(j)
     END IF
   END DO
   CALL bound_cond (n, F, bc)
   
  ! update q values
   DO j=2,n-1
     q(j) = q(j) - (F(j+1) - F(j))*dt/dvl(j)
   END DO
   CALL bound_cond (n, q, bc)
END SUBROUTINE upwind1


SUBROUTINE vanleer (q, dx, v, dt, n, bc)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n, bc
   REAL(KIND=b), INTENT(INOUT) :: q(n)
   REAL(KIND=b), INTENT(IN) :: dx(n), v(n)
   REAL(KIND=b), INTENT(IN) :: dt
   REAL(KIND=b) :: delta_q(n), d_q(n), sigma(n)
   REAL(KIND=b) :: num, den, first, second
   INTEGER :: j
   
  ! define courant number sigma
   DO j=2,n-1
     sigma(j) = v(j)*dt/dx(j)
   END DO
   CALL bound_cond (n, sigma, bc)
   
  ! define delta_quantities
   DO j=2,n-1
     delta_q(j) = q(j+1) - q(j)
   END DO
   CALL bound_cond (n, delta_q, bc)
   
  ! define d_quantities
   DO j=2,n-1
     num = delta_q(j-1)*delta_q(j)
     den = q(j+1) - q(j-1)
     IF (num.GT.0.) THEN
       d_q(j) = 2*num/den
     ELSE
       d_q(j) = 0.
     END IF
   END DO
   CALL bound_cond (n, d_q, bc)
   
   ! update q values
   DO j=2,n-1
     IF (v(j).GE.0.) THEN
       first = sigma(j)*(q(j) - q(j-1))
       second = (1. - sigma(j))*(d_q(j) - d_q(j-1))
     ELSE
       first = sigma(j)*(q(j+1) - q(j))
       second = (1. - sigma(j))*(d_q(j+1) - d_q(j))
     END IF
     q(j) = q(j) - first - (sigma(j)/2.)*second
   END DO
   CALL bound_cond (n, q, bc)
END SUBROUTINE vanleer


END PROGRAM
