program DebyeFunction
  use iso_fortran_env, only: real64
  implicit none

  !==================================================================
  ! Parámetros: define aquí los valores de T y ThetaD según necesites
  !==================================================================
  real(real64), parameter :: T      = 308.0_real64      ! Temperatura en K
  real(real64), parameter :: ThetaD = 315.0_real64      ! Temperatura de Debye en K
  real(real64), parameter :: R = 8.31446_real64         ! Constante de los gases ideales en J/molK
  real(real64), parameter :: S = 33.3_real64            ! Valor nominal de entropia del cobre en J/molK
  !==================================================================
  ! Parámetros de la cuadratura
  !==================================================================
  integer, parameter       :: N = 10000                ! Número de subintervalos
  real(real64)             :: xD, upper, dx, sum, x, xE
  real(real64)             :: result, SD, ThetaE, SE, ErrorE, ErrorD
  integer                  :: i

  !==================================================================
  ! Cálculo de la integral por regla del trapecio (punto medio)
  !==================================================================
  xD    = ThetaD / T
  upper = xD
  dx    = upper / real(N, real64)
  sum   = 0.0_real64

  do i = 1, N
    x   = (i - 0.5_real64) * dx
    sum = sum + x**3 / (exp(x) - 1.0_real64)
  end do

  result = 3.0_real64 * (T**3) / (ThetaD**3) * sum * dx ! Función de Debye
  
  SD=R*(4*result-3*log(1-exp(-xD))) !Entropía modelo de Debye
  
  ThetaE=sqrt(0.6)*ThetaD           !Temperatura de Einstein en K
  
  xE=ThetaE / T
  
  SE=3*R*((xE/(exp(xE)-1))-log(1-exp(-xE)))
  
  ErrorE=(abs(SE-S)/S)*100
  
  ErrorD=(abs(SD-S)/S)*100
  !==================================================================
  ! Salida
  !==================================================================
  print '(A, F12.6)', 'D(ThetaD/T) = ', result
  print '(A, F12.6,A)', 'SD = ', SD, 'J/molK '
  print '(A, F12.6,A)', 'ThetaE = ', ThetaE, 'K '
  print '(A, F12.6,A)', 'SE = ', SE, 'J/molK '
  print '(A, F12.6)', 'Error relativo porcentual modelo de Debye = ', ErrorD
  print '(A, F12.6)', 'Error relativo porcentual modelo de Einstein = ', ErrorE
end program DebyeFunction
