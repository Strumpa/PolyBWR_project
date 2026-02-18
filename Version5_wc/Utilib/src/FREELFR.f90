!
!-----------------------------------------------------------------------
!
!Purpose:
! Lead properties
!
!Author(s): Francesco Pepe
!
!Reference: 
! C. Fazio et Al. Handbook on Lead-bismuth Eutectic Alloy and Lead Properties, 
! Materials Compatibility, Thermal-hydraulics and Technologies. Technical 
! Report No. 7268, OECD/NEA, Paris, France, 2015.
!-----------------------------------------------------------------------

module FREELFR
    use iso_fortran_env, only: sp => real32, dp => real64, int32
    implicit none
    
    private
    public :: THMLT, THMLH, THMLMBT, CHECKLIM

    ! Lead melting and boiling temperatures at 1 atm in (K)
    real(sp), parameter :: T_MELTING = 600.6
    real(sp), parameter :: T_BOILING = 2021.0
    ! Maximum valid temperature for lead properties, corresponding to
    ! thermal conductivity
    real(sp), parameter :: T_MAX = 1300.0
    ! Coefficients for density, viscosity, conductivity, specific heat, enthalpy
    real(sp), parameter :: A_DENS = 11441.0, B_DENS = 1.2795
    real(sp), parameter :: A_VISC = 4.55e-4, B_VISC = 1069.0
    real(sp), parameter :: A_COND = 9.2, B_COND = 0.011
    real(sp), parameter :: A_CAP = 176.2, B_CAP = 4.923e-2, &
        C_CAP = 1.544e-5, D_CAP = 1.524e6 
    real(sp), parameter :: A_ENT = A_CAP, B_ENT = 2.4615e-2, &
        C_ENT = 5.147e-6, D_ENT = D_CAP
    ! Polynomial coefficients for T(h)
    real(dp), parameter :: COEFFS_H2T(0:12) = [ &
        -1.4386933d-62, 1.8610280d-56, -1.0864315d-50, 3.8078996d-45, & 
        -9.0392430d-40, 1.5620228d-34, -2.0613129d-29, 2.2198098d-24, & 
        -2.3788483d-19, 1.5494688d-14, 2.5647972d-09, 6.7577808d-03, & 
        6.0060000d+02] ! from degree 12 to degree 0
    ! Maximum number of iterations and tolerance for Tfromh
    integer(int32), parameter :: MAXITER = 50
    real(sp), parameter :: TOL = 1.0e-8

contains
    pure elemental real(sp) function dens(t) result(res)
        ! Computes density (kg/m3) as a function of temperature (K)
        real(sp), intent(in) :: t ! Temperature (K)
        res =  A_DENS - B_DENS * t
    end function

    pure elemental real(sp) function visc(t) result(res)
        ! Computes dynamic viscosity (Pa*s) as a function of temperature (K)
        real(sp), intent(in) :: t ! Temperature (K)
        res = A_VISC * exp(B_VISC / t)
    end function

    pure elemental real(sp) function cond(t) result(res)
        ! Computes thermal conductivity (W/m/K) as a function of temperature (K)
        real(sp), intent(in) :: t ! Temperature (K)
        res = A_COND + B_COND * t
    end function

    pure elemental real(sp) function cap(t) result(res)
        ! Computes heat capacity (J/kg/K) as a function of temperature (K)
        real(sp), intent(in) :: t ! Temperature (K)
        res = A_CAP - t * (B_CAP - C_CAP * t) - D_CAP / t / t
    end function

    pure elemental real(sp) function ent(t) result(res)
        ! Computes specific enthalpy (J/kg) as a function of temperature (K)
        real(sp), intent(in) :: t ! Temperature (K)
        res = t * (176.2 - t * (2.4615e-2 - 5.147e-6 * t)) &
            - 600.6 * (176.2 - 600.6 * (2.4615e-2 - 5.147e-6 * 600.6))&
            + 1.524e6 * (1 / t - 1 / 600.6)
        !res = t * (A_ENT - t * (B_ENT - C_ENT * t)) - T_MELTING * &
        !    (A_ENT - T_MELTING * (B_ENT - C_ENT * T_MELTING)) &
        !    + D_ENT * (1 / t - 1 / T_MELTING)
    end function

    pure elemental real(sp) function Tfromh(h) result(res)
        ! Computes temperature (K) as a function of specific enthalpy (J/kg)
        real(sp), intent(in) :: h ! Specific enthalpy (J/kg)
        real(dp) :: h_dp, res_dp
        integer(int32) :: i
        ! Cast h to double precision for polynomial evaluation
        h_dp = real(h, dp)
        ! Horner's scheme
        res_dp = COEFFS_H2T(0)
        do i = 1, 12
            res_dp = res_dp * h_dp + COEFFS_H2T(i)
        end do
        res = real(res_dp, sp)
    end function Tfromh

! -----------------------------------------------------------------------
! Public Subroutines
! -----------------------------------------------------------------------
    subroutine CHECKLIM(t) ! Check temperature limits for lead properties
        real(sp), intent(in) :: t ! Temperature (K)
        if (t .le. T_MELTING .or. t .ge. T_BOILING) then
            !call xabort('Temperature out of validity range for molten lead properties.')
            print*,'Temperature out of validity range for molten lead properties.'
            stop 
        end if
        if (t .gt. T_MAX) then
            print *, 'Warning: Temperature exceeds maximum valid temperature for lead properties correlations.'
        end if
    end subroutine CHECKLIM

    pure subroutine THMLMBT(tm, tb, tmax) ! THM Lead Melting, Boiling Temperature, maximum valid temperature
        ! return the lead melting and boiling temperatures at 1 atm, and the maximum valid temperature (K)
        real(sp), intent(out) :: tm, tb, tmax ! Melting temp (K), Boiling temp (K), Max valid temp (K)
    
        tm = T_MELTING
        tb = T_BOILING
        tmax = T_MAX
    end subroutine THMLMBT

    pure elemental subroutine THMLT(t, rho, h, k, mu, cp) ! THM Lead properties at temperature T
        ! return the density (kg/m3), enthalpy (J/kg), thermal conductivity (W/m/K),
        ! dynamic viscosity (kg/m/s) and specific heat capacity (J/kg/K)
        ! as a function of the temperature (K)
        real(sp), intent(in) :: t ! Temperature (K)
        real(sp), intent(out) :: rho, h, k, mu, cp ! Properties to be returned

        rho = dens(t)
        k = cond(t)
        mu = visc(t)
        cp = cap(t)
        h = ent(t)
    end subroutine THMLT

    pure elemental subroutine THMLH(h, rho, t) ! THM Lead properties at enthalpy H
        ! return the density (kg/m3) and temperature (K) as a function of 
        ! the enthalpy (J/kg)
        real(sp), intent(in) :: h ! Enthalpy (J/kg)
        real(sp), intent(out) :: rho, t ! Properties to be returned
        
        t = Tfromh(h)
        rho = dens(t)
    end subroutine THMLH
end module FREELFR
