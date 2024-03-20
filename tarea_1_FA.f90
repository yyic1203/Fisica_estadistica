PROGRAM tarea_1_FA
  IMPLICIT NONE
  INTEGER :: numa, n
  REAL :: realNumber
  DOUBLE PRECISION, PARAMETER :: T=298.0D0, P=1.0D0, NA = 6.022D+23, R=0.083140D0 
  DOUBLE PRECISION, PARAMETER :: pi = 3.1416D0, h = 6.62615D-34, kb = 1.3806D-23
  DOUBLE PRECISION, PARAMETER :: RJ = 8.3145D0
  CHARACTER(len=20), ALLOCATABLE :: sima(:), confa(:)
  DOUBLE PRECISION, ALLOCATABLE :: masa(:), vm(:), qt(:), ut(:), gt(:), st(:), deg(:), ve(:), eto(:)
  DOUBLE PRECISION, ALLOCATABLE :: qe(:), ue(:), ge(:), se(:), uto(:), gto(:), sto(:), et(:)

  WRITE(*,'(A)', ADVANCE='NO') "Este programa busca calcular las propiedades termodinámicas:  "
  WRITE(*,'(A)', ADVANCE='NO') "Energía libre de Gibbs (G), energía interna (U) y entropia (S)  "
  WRITE(*,'(A)') "a condiciones estándar en átomos con la función de partición traslacional "
  WRITE(*,'(A)') "y electrónica, además de compararla con los valores reportados. "
  WRITE(*,'(A)', ADVANCE='NO') "Digite cuantos átomos piensa analizar: "
  READ *, numa

  ALLOCATE(sima(numa), confa(numa), masa(numa), vm(numa), qt(numa), ut(numa), gt(numa), sto(numa), et(numa), eto(numa))
  ALLOCATE(st(numa), deg(numa), qe(numa), ue(numa), ge(numa), se(numa), uto(numa), gto(numa), ve(numa))
  
  DO n = 1, numa
    WRITE(*,'(A,I3,A)', ADVANCE='NO') "Introduzca el símbolo del átomo numero ", n, ": "
  READ *, sima(n)
  WRITE(*,'(A,A)', ADVANCE='NO') "Introduzca la configuración electrónica de los electrones en el átomo ", sima(n)
  WRITE(*,'(A)', ADVANCE='NO') "el último subnivel del átomo (lx) l(minus)=subnivel, x=num.electrones: "
  READ *, confa(n)
  WRITE(*,'(A)', ADVANCE='NO') "Introduzca la masa atómica del átomo " // TRIM(sima(n)) // ": "
  READ *, masa(n)
  WRITE(*,'(A)', ADVANCE='NO') "Introduzca el valor experimental de entropía del átomo " // TRIM(sima(n)) // ": "
  READ *, ve(n)
  END DO
  
 !! Función de partición traslacional
 	DO n=1, numa
		vm(n)= R*T/(NA*P*1000)
		qt(n)= ((2.0D0 * pi * masa(n) * kb * T / (h**2))**(1.5D0))*vm(n)
		st(n)= RJ*dlog(qt(n)*exp(1.5D0))
		ut(n)= 1.5D0*RJ*T
		gt(n)= RJ*T*dlog(exp(1.0D0)/qt(n))
		et(n)= (DABS(ve(n)-st(n))/ve(n))*100.0D0
		WRITE(*,'(A,F9.3,A)') "La entropía traslacional de " // TRIM(sima(n)) // " es: ", st(n), " J/K mol"
		WRITE(*,'(A,F9.3,A)') "La energía interna traslacional de " // TRIM(sima(n)) // " es: ", ut(n), " J/mol"
		WRITE(*,'(A,E10.3,A)') "La energía libre de Gibbs traslacional de " // TRIM(sima(n)) // " es: ", gt(n), " J/mol"
		WRITE(*,'(A,F9.3,A)') "El error relativo porcentual traslacional de " // TRIM(sima(n)) // " es: ", et(n), "%"
		!! WRITE(*,'(A,A,A,E10.3)', ADVANCE='NO') "La función de partición traslacional de ", TRIM(sima(n)), " es: ", qt(n)
	End DO
	!!Calculo de la degenerancia
	DO n=1, numa
		IF (TRIM(ADJUSTL(confa(n))) == "s2" .OR. TRIM(ADJUSTL(confa(n))) == "p6") THEN
		  deg(n)=1.0D0
		ELSE IF (TRIM(ADJUSTL(confa(n))) == "d10" .OR. TRIM(ADJUSTL(confa(n))) == "p2" .OR. TRIM(ADJUSTL(confa(n))) == "d4") THEN
		  deg(n)=1.0D0
		ELSE IF (TRIM(ADJUSTL(confa(n))) == "p1") THEN
		  deg(n)=2.0D0
		ELSE IF (TRIM(ADJUSTL(confa(n))) == "p5" .OR. TRIM(ADJUSTL(confa(n))) == "d1") THEN
			deg(n)=4.0D0
		ELSE IF (TRIM(ADJUSTL(confa(n))) == "p4" .OR. TRIM(ADJUSTL(confa(n))) == "d2") THEN
			deg(n)=5.0D0
		ELSE IF (TRIM(ADJUSTL(confa(n))) == "p3" .OR. TRIM(ADJUSTL(confa(n))) == "d9") THEN
			deg(n)=6.0D0
		ELSE IF (TRIM(ADJUSTL(confa(n))) == "d8" .OR. TRIM(ADJUSTL(confa(n))) == "d6") THEN 
			deg(n)=9.0D0
		ELSE IF (TRIM(ADJUSTL(confa(n))) == "d7" .OR. TRIM(ADJUSTL(confa(n))) == "d5") THEN 
			deg(n)=10.0D0
		ELSE
		  PRINT *, "Digitaste mal la configuración electrónica"
		END IF
		!!WRITE(*,'(A,A,A,F9.3)', ADVANCE='NO') "La degenerancia del átomo ", TRIM(sima(n)), " es: ", deg(n)
	End DO
	!!Función de partición electrónica
	DO n=1, numa
		qe(n)=deg(n)
		se=RJ*dlog(qe(n))
		ue=0.0D0
		ge=RJ*T*dlog(exp(1.0D0)/qe(n))
		WRITE(*,'(A,F9.3,A)') "La entropía electrónica de " // TRIM(sima(n)) // " es: ", se(n), " J/K mol"
		WRITE(*,'(A,F9.3,A)') "La energía interna electrónica de " // TRIM(sima(n)) // " es: ", ue(n), " J/mol"
		WRITE(*,'(A,E10.3,A)') "La energía libre de Gibbs electrónica de " // TRIM(sima(n)) // " es: ", ge(n), " J/mol"
	End DO
	!!Usando las 2 funciones
	DO n=1, numa
		sto=2.5D0*RJ+RJ*dlog(qt(n)*qe(n))
		uto=1.5D0*RJ*T
		gto=RJ*T*dlog(qt(n)*qe(n)/NA)
		eto(n)= (DABS(ve(n)-sto(n))/ve(n))*100.0D0
		WRITE(*,'(A,F9.3,A)') "La entropía total de " // TRIM(sima(n)) // " es: ", sto(n), " J/K mol"
		WRITE(*,'(A,F9.3,A)') "La energía interna total de " // TRIM(sima(n)) // " es: ", uto(n), " J/mol"
		WRITE(*,'(A,E10.3,A)') "La energía libre de Gibbs total de " // TRIM(sima(n)) // " es: ", gto(n), " J/mol"
		WRITE(*,'(A,F9.3,A)') "El error relativo porcentual total de " // TRIM(sima(n)) // " es: ", eto(n), "%"
	End DO
		DEALLOCATE(sima, confa, masa, vm, qt, ut, gt, st, deg, qe,se,ue,ge,sto,gto,uto)
END PROGRAM tarea_1_FA



