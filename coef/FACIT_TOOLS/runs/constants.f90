module constants

  implicit none

  ! double precision floats
  integer, parameter :: rkind = kind(1.d0)

  ! mathematical constants
  real(rkind), parameter :: pi    = 4.d0*ATAN(1.d0)
  real(rkind), parameter :: pio2  = 2.d0*ATAN(1.d0)
  real(rkind), parameter :: twopi = 8.d0*ATAN(1.d0)
  real(rkind), parameter :: sqrt2 = SQRT(2.d0)

  ! physical constants
  real(rkind), parameter :: mp         = 1.67262192369e-27_rkind         ! proton mass [kg]
  real(rkind), parameter :: me         = 9.1093837015e-31_rkind          ! electron mass [kg]
  real(rkind), parameter :: q_e        = 1.602176634e-19_rkind           ! electron charge [C]
  real(rkind), parameter :: eps0       = 8.8541878128e-12_rkind          ! vacuum permittivity [F/m]
  real(rkind), parameter :: eps_pi_fac = 3.0*eps0**2*(2.0*pi/q_e)**1.5/q_e ! to be used in Coulomb logarithms [F^2/(m^2 C^(5/2))]

end module constants
