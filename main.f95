program calculate_cv
    implicit none

    ! Constantes
    real(8), parameter :: R = 8.31446261815324d0  ! Constante de los gases en J/(mol·K)
    real(8), parameter :: NA = 6.02214076d23      ! Número de Avogadro

    ! Variables de entrada
    real(8) :: epsilon_kB  ! epsilon / kB en Kelvin
    real(8) :: sigma       ! sigma en metros
    real(8) :: lambda      ! lambda (Adimensional)
    real(8) :: T           ! Temperatura en Kelvin
    real(8) :: V           ! Volumen en metros cúbicos por mol
    real(8) :: CV_experimental  ! Valor experimental de CV en J/(mol·K)

    ! Variables calculadas
    real(8) :: b0          ! b0
    real(8) :: CV          ! Capacidad calorífica a volumen constante
    real(8) :: error_relativo  ! Error relativo porcentual

    ! Valores de entrada 
    epsilon_kB = 284.0d0   
    sigma = 3.57d-10       
    lambda = 1.44d0        
    T = 300.0d0           
    V = 1.1d-4        
    CV_experimental = 20.8d0  

    ! Calcular b0
    b0 = (2.0d0 * 3.141592653589793d0 * NA * sigma**3) / 3.0d0

    ! Calcular CV
    CV = (R * b0 / V) * (lambda**3 - 1.0d0) * (epsilon_kB / T)**2 * exp(epsilon_kB / T)

    ! Calcular el error relativo porcentual
    error_relativo = abs(CV - CV_experimental) / CV_experimental * 100.0d0

    ! Imprimir resultados
    print *, "Valor calculado de CV: ", CV, " J/mol·K"
    print *, "Valor experimental de CV: ", CV_experimental, " J/mol·K"
    print *, "Error relativo porcentual: ", error_relativo, "%"

end program calculate_cv