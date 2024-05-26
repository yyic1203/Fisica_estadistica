PROGRAM P3FE
	IMPLICIT NONE
	integer :: iter
	DOUBLE PRECISION, PARAMETER :: hbar=1.064D-34, NA = 6.022D+23, pi= 3.141592653589793d0, m= 9.10D-31, kb= 1.3806D-23
	CHARACTER(len=5) :: RB, ELE
	DOUBLE PRECISION :: PR, T
	DOUBLE PRECISION :: D, Ef, Vf, Mf, Mu, E, NAP, crit, Delta, mu_old
	
	WRITE(*,'(A)') "Este programa busca calcular las propiedades termodinámicas de redes cristalinas:  "
  WRITE(*,'(A)', ADVANCE= 'NO') "Densidad (p), Energía de Fermi (Ef), la velocidad de los electrones en la capa (Vf),  "
  WRITE(*,'(A)') "El potencial y la energía a cualquier temperatura T dada por el usuario "
  WRITE(*,'(A)') "Digite la nomenclatura del átomo: "
  READ *, ELE
  WRITE(*,'(A)') "Existen distintos tipos de redes de Bravais se digitan sus siglas de la siguiente manera: "
  WRITE(*,'(A)') "CS: Cúbica simple"
  WRITE(*,'(A)') "FCC: Cúbica centrada en las caras "
  WRITE(*,'(A)') "BCC: Cúbica centrada en el cuerpo "
  WRITE(*,'(A)',ADVANCE='NO') "Cuál es la red de Bravais de "
  WRITE(*,'(A)',ADVANCE='NO')  ELE
  WRITE(*,'(A)')  ":"
  READ *, RB
  
  IF (TRIM(	ADJUSTL(RB)) == "CS") THEN
  	NAP = 1.0D0
  ELSE IF (TRIM(	ADJUSTL(RB)) == "FCC") THEN
  	NAP = 4.0D0
  ELSE IF (TRIM(	ADJUSTL(RB)) == "BCC") THEN
  	NAP = 2.0D0
  ELSE
		  PRINT *, "Digitaste mal la red de Bravais"
	END IF
	
	WRITE(*,'(A)', ADVANCE='NO') "Introduzca el parámetro de red de "
	WRITE(*,'(A)', ADVANCE='NO') ELE
	WRITE(*,'(A)') ":"
  READ *, PR
  
  WRITE(*,'(A)') "¿A qué temperatura esta el cristal?:"
  READ *, T
  
  D = NAP/PR**3.0d0
  Mf=((hbar**2.0d0)/(2.0D0*m))*(3.0d0*pi**2.0d0*D)**0.666d0
  vf=(hbar/m)*(3.0d0*pi**2.0d0*D)**0.333d0
  Ef= (3.0d0/5.0d0)*NAP*Mf
  
  iter=0
  crit= 1.0D-25!!Criterio de estabilidad
  mu_old = 1 
  1 iter=iter+1
  Mu = Mf * (1.0d0 - (pi**2 / 12.0d0) * (kb**2 * T**2 / mu_old**2))
  Delta= ABS(Mu-mu_old)
  mu_old = Mu
  
  
  if (Delta .gt. crit) then
		go to 1 !!Retorna el ciclo
	end if
	E = EF * (1.0d0 - (5.0d0*pi**2 / 12.0d0) * (kb**2 * T**2 / mu_old**2))
	
  WRITE(*,'(A,E9.3,A)') "La densidad es: ", D, "Electones/m3"
  WRITE(*,'(A,E9.3,A)') "El potencial químico de Fermi o la energía de Fermi: ", Mf, "J"
  WRITE(*,'(A,E9.3,A)') "La energía del punto cero es: ", Ef, "J"
  WRITE(*,'(A,E9.3,A)') "La velocidad de Fermi es: ", vf, "m/s"
  WRITE(*,'(A,F8.3,A,E9.2,A)') "El potencia químico : ", T, "K es:", mu_old, "J"
  WRITE(*,'(A,F8.3,A,E9.2,A)') "La energía interna del cristal a: ", T, "K es:", E, "J"
END PROGRAM P3FE
