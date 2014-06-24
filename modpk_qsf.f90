module modpk_qsf
  use modpkparams
  use modpk_numerics
  implicit none

  !Need one choice for every heavy direction (assumes one light direction)
  !Ref functions turning_function and derivs
  integer, dimension(:), allocatable :: turning_choice

  !For turning trajs
  integer :: effective_V_choice
  integer, dimension(:), allocatable :: number_knots_qsfrandom
  real(dp), dimension(:), allocatable :: stand_dev_qsfrandom
  real(dp), dimension(:,:,:), allocatable :: knot_positions
  real(dp), dimension(:), allocatable :: knot_range_min, knot_range_max
  logical :: custom_knot_range

  !For a given qsf trajectory defined by a parametric curve, gives
  !reference variables for calls to zero finding, etc
  type :: qsf_reference
    real(dp), dimension(:), allocatable :: phi
    real(dp) :: param
    integer :: hunt_guess
    logical :: traj_init = .false.
    real(dp), dimension(:,:), allocatable :: phi_light_vs_param
    contains
      procedure :: get_init_param => choose_initial_parameter
      procedure :: initialize_traj => integrate_through_trajectory
      procedure :: phi_light => get_phi_light
  end type

  type(qsf_reference) :: qsf_runref


contains

  !Function that parameterizes the turn for quasi--single-field trajectories
  !funct_i(phi_light)
  !Function parameters passed here via vparams(4+,:)
  real(dp) function turning_function(phi, heavy_field_index, turning_choice)
    real(dp), intent(in) :: phi
    integer, intent(in) :: turning_choice, heavy_field_index

    !Hyperbola params
    real(dp) :: offset_phi, asympt_angle, focal_length

    !hyperbolic tan params
    real(dp) :: turn_magnitude, turn_sharpness

    !random traj
    integer :: i
    real(dp) :: phi_knot1, phi_knot2
    real(dp) :: chi_knot1, chi_knot2
    real(dp) :: low_endpoint, high_endpoint
    real(dp) :: slope
    real(dp) :: range
    !heavy_field_index => which heavy field? for picking (light,heavy_i) direction
    !  This is the "i" in funct_i => need heavy_field_index>1

    !turning_choice => pick a function for turn in (light,heavy_i) direction

    !Parameters for turning_function are set in:
    !vparams(j+4,heavy_field_index) = j^th param for turn in heavy_field_index
    !   direction

    turning_function = 0e0_dp

    if (heavy_field_index <2) then
      print*, "Set heavy_field_index>1. heavy_field_index=", heavy_field_index
      stop
    end if

    select case(turning_choice)
    case(1)
      !No turn, effectively single-field
      turning_function = 0e0_dp
    case(2)
      !North-facing hyperbola, symm around y-axis, min at phi=0

      if (size(vparams,1) <6) then
        print*, "Not enough vparams to set turning_choice=", turning_choice
        stop
      end if

      offset_phi   = vparams(4,heavy_field_index) !min turning_function in phi
      asympt_angle = vparams(5,heavy_field_index) !turn angle
      focal_length = vparams(6,heavy_field_index) !dist focal length (sharpness)

      turning_function = sqrt( (focal_length**2)*(sin(asympt_angle)**2) + &
        ((phi - offset_phi)**2) * (tan(asympt_angle)**2) ) &
        -focal_length*sin(asympt_angle)
    case(3)
      !hyperbolic tangent

      if (size(vparams,1) <6) then
        print*, "Not enough vparams to set turning_choice=", turning_choice
        stop
      end if

      offset_phi     = vparams(4,heavy_field_index) !position of the turn in phi direction
      turn_magnitude = vparams(5,heavy_field_index) !the asymptotic displacement in the heavy direction
      turn_sharpness = vparams(6,heavy_field_index) !the sharpness of the turn

      turning_function = turn_magnitude*tanh(turn_sharpness*(phi - offset_phi))

    case(4)
      !Random turning trajectory

      call find_bounding_knots(phi, heavy_field_index, &
        phi_knot1, phi_knot2, &
        chi_knot1, chi_knot2)

      slope = (chi_knot2 - chi_knot1)/(phi_knot2 - phi_knot1)

      turning_function = slope*(phi - phi_knot1) + chi_knot1

    case(5)
      !Two-knot trajectory
       if (size(vparams,1) <6) then
        print*, "Not enough vparams to set turning_choice=", turning_choice
        stop
       end if

       !vparams(4,:) = phi position of turn midpoint
       !vparams(5,:) = angle from horizontal
       !vparams(6,:) = deviation in heavy field from horiz

       offset_phi = vparams(4,heavy_field_index)
       slope = tan(vparams(5,heavy_field_index))

       chi_knot2 = vparams(6,heavy_field_index)
       chi_knot1 = -chi_knot2

       range = 2.0e0_dp*vparams(6,heavy_field_index)/ &
         tan(vparams(5,heavy_field_index))

       phi_knot2 = offset_phi + range/2.0e0_dp
       phi_knot1 = offset_phi - range/2.0e0_dp

       if (phi > phi_knot2) then
         turning_function = chi_knot2
       else if (phi > phi_knot1) then
         turning_function = slope*(phi-offset_phi)
       else 
         turning_function = chi_knot1
       end if 

    case default
       write(*,*) 'QSF: Need to set turning_function in modpk_potential.f90 for turning_choice =',turning_choice
    end select

  end function turning_function

  !Function that parameterizes the turn for quasi--single-field trajectories
  real(dp) function dturndphi(phi, heavy_field_index, turning_choice)
    real(dp), intent(in) :: phi
    integer, intent(in) :: turning_choice, heavy_field_index

    !Hyperbola params
    real(dp) :: offset_phi, asympt_angle, focal_length

    !hyperbolic tan params
    real(dp) :: turn_magnitude, turn_sharpness

    !Random traj
    real(dp) :: phi_knot1, phi_knot2
    real(dp) :: chi_knot1, chi_knot2

    dturndphi = 0e0_dp


    select case(turning_choice)
    case(1)
      dturndphi = 0e0_dp
    case(2)
      !North-facing hyperbola, symm around y-axis, min at phi=0

      offset_phi   = vparams(4,heavy_field_index) !min turning_function in phi
      asympt_angle = vparams(5,heavy_field_index) !turn angle
      focal_length = vparams(6,heavy_field_index) !dist focal length (sharpness)

#define funct turning_function(phi,heavy_field_index,turning_choice)
      dturndphi = (tan(asympt_angle)**2)* &
        (phi - offset_phi)/(funct+focal_length*sin(asympt_angle))
#undef funct

    case(3)
      !hyperbolic tan, centred at chi_i = 0 and phi= offset_phi

      offset_phi     = vparams(4,heavy_field_index) !position of the turn in phi direction
      turn_magnitude = vparams(5,heavy_field_index) !the asymptotic displacement in the heavy direction
      turn_sharpness = vparams(6,heavy_field_index) !the sharpness of the turn

      dturndphi = turn_magnitude*turn_sharpness/(cosh(turn_sharpness*(phi - offset_phi))**2)

    case(4)
      !Random turning trajectories

      call find_bounding_knots(phi, heavy_field_index, &
        phi_knot1, phi_knot2, &
        chi_knot1, chi_knot2)

      dturndphi = (chi_knot2 - chi_knot1)/(phi_knot2 - phi_knot1)
   
   case(5)
      !Two-knot trajectory
      
       if (phi>phi_knot2) then
          dturndphi = 0e0_dp
        else if (phi>phi_knot1) then
           dturndphi = tan(vparams(5,heavy_field_index))
        else 
           dturndphi = 0e0_dp
       end if 

    case default
       write(*,*) 'QSF: Need to set turning_function in modpk_potential.f90 for turning_choice =',turning_choice
    end select


  end function dturndphi

  !Function that parameterizes the turn for quasi--single-field trajectories
  real(dp) function d2turndphi2(phi, heavy_field_index, turning_choice)
    real(dp), intent(in) :: phi
    integer, intent(in) :: turning_choice, heavy_field_index

    !Hyperbola params
    real(dp) :: offset_phi, asympt_angle, focal_length

    !hyperbolic tan params
    real(dp) :: turn_magnitude, turn_sharpness

    d2turndphi2=0e0_dp

    select case(turning_choice)
    case(1)
      d2turndphi2 = 0e0_dp

    case(2)
      !North-facing hyperbola, symm around y-axis, min at phi=0

      offset_phi   = vparams(4,heavy_field_index) !min turning_function in phi
      asympt_angle = vparams(5,heavy_field_index) !turn angle
      focal_length = vparams(6,heavy_field_index) !dist focal length (sharpness)

#define funct turning_function(phi,heavy_field_index,turning_choice)
#define dfunct dturndphi(phi,heavy_field_index,turning_choice)

      d2turndphi2 = ((tan(asympt_angle)**2)/&
        (funct+focal_length*sin(asympt_angle))) * &
        (1.0e0_dp  - (phi-offset_phi)*&
        (dfunct/ (funct &
        +focal_length*sin(asympt_angle))))

#undef funct
#undef dfunct

    case(3)
      !hyperbolic tan, centred at chi_i = 0 and phi= offset_phi

      offset_phi     = vparams(4,heavy_field_index) !position of the turn in phi direction
      turn_magnitude = vparams(5,heavy_field_index) !the asymptotic displacement in the heavy direction
      turn_sharpness = vparams(6,heavy_field_index) !the sharpness of the turn

      d2turndphi2 = -2.0e0_dp*turn_magnitude*(turn_sharpness**2)* &
                     sinh(turn_sharpness*(phi-offset_phi))/(cosh(turn_sharpness*(phi-offset_phi))**3)
    case(4)
      !Random turning trajectories
      d2turndphi2 = 0e0_dp

    case(5)
       !two-knot trajectory
      d2turndphi2 = 0e0_dp
    case default
       write(*,*) 'QSF: Need to set turning_function in modpk_potential.f90 for turning_choice =',turning_choice
    end select


  end function d2turndphi2


  !For qsf_random
  !Find which knots phi is between
  subroutine find_bounding_knots(phi, heavy_field_index, &
      phiknot_low, phiknot_high, &
      heavyknot_low, heavyknot_high)
    real(dp), intent(in) :: phi
    real(dp), intent(out) :: phiknot_low, phiknot_high
    real(dp), intent(out) :: heavyknot_low, heavyknot_high
    integer, intent(in) :: heavy_field_index
    integer :: i

#define HEAVY_INDX (heavy_field_index-1)

    phiknot_low = knot_positions(HEAVY_INDX,&
      number_knots_qsfrandom(HEAVY_INDX),&
      1)
    phiknot_high = phi_init0(1)

    heavyknot_low = knot_positions(HEAVY_INDX,&
      number_knots_qsfrandom(HEAVY_INDX),&
      2)
    heavyknot_high = 0e0_dp

    do i=1,number_knots_qsfrandom(HEAVY_INDX)
      if (phi < knot_positions(HEAVY_INDX,i,1)) then
        if (i==1) then
          phiknot_low = 0e0_dp  !ending condition
          heavyknot_low = 0e0_dp  !ending condition
        else
          phiknot_low = knot_positions(HEAVY_INDX,i-1,1)
          heavyknot_low = knot_positions(HEAVY_INDX,i-1,2)
        end if
        phiknot_high = knot_positions(HEAVY_INDX,i,1)
        heavyknot_high = knot_positions(HEAVY_INDX,i,2)
        exit

      !elseif (phi == knot_positions(HEAVY_INDX,i,1)) then
      !  print*, "ERROR: In turning_function, exactly at knot-point."
        !DEBUG
      !  stop
      end if
    end do

#undef HEAVY_INDX

  end subroutine find_bounding_knots

  !For numerical QSF trajectory potential
  !This is the parametric vector function in field space
  !that one defines, around which we construct a valley
  function turning_function_parametric(param) result(funct)
    implicit none
    real(dp), intent(in) :: param
    real(dp), dimension(num_inflaton) :: funct

    select case(turning_choice(1))
    case(0)
      !Line
      funct = param
    case(1)
      !Parabola
      funct(1) = param
      funct(2) = param**2
    case(2)
      !Helix
      funct(1) = sin(param)
      funct(2) = cos(param)
    case default
      print*, "ERROR: turning_function_parametric not implemented for"
      print*, "ERROR: turning_choice = ", turning_choice
      stop
    end select

  end function turning_function_parametric

  !Deriv wrt parameter
  function dturndparam(param) result(funct)
    implicit none
    real(dp), intent(in) :: param
    real(dp), dimension(num_inflaton) :: funct

    select case(turning_choice(1))
    case(0)
      !Line
      funct = 1e0_dp
    case(1)
      !Parabola
      funct(1) = 1e0_dp
      funct(2) = 2e0_dp*param
    case(2)
      !Helix
      funct(1) = cos(param)
      funct(2) = -sin(param)
    case default
      print*, "ERROR: turning_function_parametric not implemented for"
      print*, "ERROR: turning_choice = ", turning_choice
      stop
    end select

  end function dturndparam

  !Second deriv wrt parameter
  function d2turndparam2(param) result(funct)
    implicit none
    real(dp), intent(in) :: param
    real(dp), dimension(num_inflaton) :: funct

    select case(turning_choice(1))
    case(0)
      !Line
      funct = 0e0_dp
    case(1)
      !Parabola
      funct(1) = 0e0_dp
      funct(2) = 2e0_dp
    case(2)
      !Helix
      funct(1) = -sin(param)
      funct(2) = -cos(param)
    case default
      print*, "ERROR: turning_function_parametric not implemented for"
      print*, "ERROR: turning_choice = ", turning_choice
      stop
    end select

  end function d2turndparam2

  function distance(t)
    use modpkparams
    implicit none

    real(dp), intent(in) :: t
    real(dp) :: distance

    distance = sqrt(sum((qsf_runref%phi - turning_function_parametric(t))**2))

  end function distance

  !Modified deriv for minimizing
  function distance_deriv(t)
    use modpkparams
    implicit none

    real(dp), intent(in) :: t
    real(dp) :: distance_deriv

    distance_deriv = sum((qsf_runref%phi - turning_function_parametric(t))*&
      dturndparam(t))

  end function distance_deriv

  !Modified deriv for minimizing
  function distance_2deriv(t)
    use modpkparams
    implicit none

    real(dp), intent(in) :: t
    real(dp) :: distance_2deriv

    distance_2deriv = sum(&
          d2turndparam2(t)*(qsf_runref%phi-turning_function_parametric(t)) -&
          dturndparam(t)**2)

  end function distance_2deriv

  !Segment of parametric path to obtain integrated distance
  function path_segment(param) result(length)
    implicit none

    real(dp) :: length
    real(dp), intent(in) :: param

    length = sqrt(sum( (dturndparam(param))**2))

  end function path_segment

  !When initializing the integration around the parametrized curve,
  !this function will get the first guess for the parameter that puts
  !phi closest to the line
  subroutine choose_initial_parameter(this, phi_light0)
    implicit none

    class(qsf_reference) :: this
    integer :: ii
    real(dp), intent(in) :: phi_light0


    if (.not. allocated(this%phi_light_vs_param)) then
      print*, "QSF: Trying to set initial parameter guess,"
      print*, "QSF: but haven't integrated trajectory yet."
      stop
    end if

#define LIGHT (this%phi_light_vs_param(:,1))
    ii= locate(LIGHT, phi_light0)
    this%hunt_guess = ii
    this%param = this%phi_light_vs_param(ii,2)
#undef LIGHT


  end subroutine choose_initial_parameter

  !Since the valley is defined parametrically in field-space,
  !we need to find how far up the valley we are for any
  !given parameter.  This integrates the path until length>phi_start
  subroutine integrate_through_trajectory(this, phi_start, &
      param0, stepsize)
    implicit none

    class(qsf_reference) :: this
    real(dp), intent(in) :: param0, phi_start
    real(dp), intent(in), optional :: stepsize

    real(dp) :: dx
    real(dp), dimension(:,:), allocatable :: length
    integer(dp), parameter :: maxsteps=1e6

    integer :: ii, counter
    real(dp) :: x_a, x_b

    if (this%traj_init) return

    if (present(stepsize)) then
      dx = stepsize
    else
      dx = 1e-4_dp
    end if

    !Integrate along parametric curve and record values for interpolation
    allocate(length(maxsteps, 2))
    length(1,:) = (/ 0.0e0_dp, param0 /) !length=(/path_length,param/)

    counter=1
    x_a = param0
    do ii=1, maxsteps-1
      counter = counter + 1

      x_b = x_a + dx
      length(counter,1) = length(counter-1,1)&
        + integrate_1d(path_segment, x_a, x_b, 10)
      length(counter,2) = x_b
      x_a = x_b

      !Check if have gone far enough
      if (length(counter,1)>1.25*phi_start) exit

      if (ii==maxsteps-1) then
        print*, "QSF: Couldn't integrate until length of curve=", phi_start
        print*, "QSF: Either lower phi_start or make maxsteps or stepsize bigger."
        stop
      end if

    end do

    allocate(this%phi_light_vs_param(counter,2))
    this%phi_light_vs_param = length(1:counter,:)
    deallocate(length)

    this%traj_init = .true.

  end subroutine integrate_through_trajectory

  !After the trajectory has been integrated, this will give the
  !distance along the trajectory (from param0) to a reference param
  function get_phi_light(this, param) result(phi_light)
    implicit none

    class(qsf_reference) :: this
    real(dp), intent(in) :: param
    integer :: param_guess
    real(dp) :: phi_light
    real(dp) :: del_phi

    integer :: ii, jj

    if (.not. this%traj_init) then
      print*, "QSF: Trying to get phi_light before integrating traj."
      stop
    end if

    param_guess = this%hunt_guess

    call locator()
    call interpolator()

    this%hunt_guess = param_guess


    contains

      subroutine locator()
        integer :: low, high

#define INTLEN size(this%phi_light_vs_param,1)

        low = max(1, param_guess-100)
        high = min(INTLEN, param_guess+100)

!DEBUG
!#define LIGHT (this%phi_light_vs_param(:,1))
!#define P_ARR (this%phi_light_vs_param(:,2))
#define LIGHT (this%phi_light_vs_param(low:high,1))
#define P_ARR (this%phi_light_vs_param(low:high,2))

        !P_ARR is monotonic
        param_guess = param_guess - low
        call hunt(P_ARR, param, param_guess)
        param_guess = param_guess + low
        ii = param_guess
        !ii= locate(P_ARR, param)

        jj=min(max(ii-(4-1)/2,1),INTLEN+1-4)

#undef P_ARR
#undef LIGHT

      end subroutine

      subroutine interpolator()
#define LIGHT (this%phi_light_vs_param(jj:jj+4,1))
#define P_ARR (this%phi_light_vs_param(jj:jj+4,2))

        call polint(&
          P_ARR, &
          LIGHT,&
          param, phi_light, del_phi)

          !Check for interpolation errors
          if(del_phi > 0.1) then
            print*,'QSF: The interpolation in get_phi_light/locator has suspiciously large'
            print*,'QSF: QUITTING'
            print*,"QSF: del_phi", del_phi
            print*,"QSF: P_ARR", P_ARR, "param", param
            print*,"QSF: LIGHT", LIGHT
            stop
          else if (phi_light > this%phi_light_vs_param(INTLEN,1)) then
            print*, "QSF: phi_light > LIGHT(MAX)"
            print*, "QSF: phi_light =", phi_light
            print*, "QSF: LIGHT(MAX) =", this%phi_light_vs_param(INTLEN,1)
            stop
          end if

#undef LIGHT
#undef P_ARR
#undef INTLEN

      end subroutine

  end function get_phi_light

end module modpk_qsf
