program virial_coefficient
    implicit none

    ! Variables de entrada
    real(8) :: T        ! Temperatua en K
    real(8) :: sigma    ! Sigma en cm
    real(8) :: epk      ! Energía sobre la constate de boltzmann en K
    real(8) :: B_exp    ! Valor experimental del segundo coeficiente del virial en cm^3/mol

    ! Variables calculadas
    real(8) :: T_star       ! Temperatura reducida
    real(8) :: B2_star      ! Segundo coeficiente del virial reducido
    real(8) :: B2           ! Segundo coeficiente del virial
    real(8) :: error        ! Error relativo porcentual

    ! Valores de entrada 
    T=222
    sigma = 3.783d-8       
    epk = 148.9       
    B_exp = -97.7

  ! Calcular la temperatura reducida
  T_star =  T / epk

  ! Calcular el segundo coeficiente del virial reducido
  B2_star = (0.7071 / (T_star**(0.25))) * (-2.451 + (3.62256 / (T_star**(0.5))))

  ! Calcular el coeficiente del virial en unidades cm^3/mol
  B2 = B2_star * (2.0/3.0) * 3.14159265359 * (sigma**3) * 6.02214d23

  ! Calcular el error porcentual respecto al experimental
  error = abs((B2 - B_exp) / B_exp) * 100.0

  ! Imprimir resultados
  print *, 'Temperatura reducida (T*):', T_star
  print *, 'Coeficiente del virial reducido (B2*):', B2_star
  print *, 'Coeficiente del virial (B2 en cm^3/mol):', B2
  print *, 'Error porcentual respecto al experimental:', error, '%'

end program virial_coefficient