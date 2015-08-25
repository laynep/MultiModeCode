module modpk_qsf
  use modpkparams
  use modpk_numerics
  use internals, only : pi
  use modpk_errorhandling, only : raise
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
      procedure, public :: get_param => choose_parameter
      procedure, public :: initialize_traj => integrate_through_trajectory
      procedure, public :: phi_light => get_phi_light
      procedure, public :: dphi_light_dparam => dphi_light_dparam
      procedure, public :: d2phi_light_dparam2 => d2phi_light_dparam2
      procedure, public :: min_dist => distance_minimizer
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
      print*, "MODECODE: heavy_field_index=", heavy_field_index
      call raise%fatal_code(&
              "Set heavy_field_index>1.", __FILE__, __LINE__)
    end if

    select case(turning_choice)
    case(1)
      !No turn, effectively single-field
      turning_function = 0e0_dp
    case(2)
      !North-facing hyperbola, symm around y-axis, min at phi=0

      if (size(vparams,1) <6) then
        print*, "MODECODE: turning_choice=", turning_choice
        call raise%fatal_code(&
                "Not enough vparams to set this turning choice.",&
                __FILE__, __LINE__)
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
        print*, "MODECODE: turning_choice=", turning_choice
        call raise%fatal_code(&
                "Not enough vparams to set this turning choice.",&
                __FILE__, __LINE__)
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
        print*, "MODECODE: turning_choice=", turning_choice
        call raise%fatal_code(&
                "Not enough vparams to set this turning choice.",&
                __FILE__, __LINE__)
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
      !  print*, "MODECODE: In turning_function, exactly at knot-point."
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

    real(dp) :: phi_turn, slope, f_decay

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
      f_decay = vparams(4,1)

      funct(1) = f_decay*sin(param)
      funct(2) =  f_decay*cos(param)
    case(3)
      !One sharp turn
      phi_turn = vparams(4,1)
      slope = vparams(4,2)

      if (param < phi_turn) then
        funct(1) = param
        funct(2) = 0e0_dp
      else
        funct(1) = slope*(param-phi_turn) + phi_turn
        funct(2) = (param-phi_turn)
      end if

    case default
      print*, "MODECODE: turning_choice = ", turning_choice
      call raise%fatal_code(&
              "MODECODE: turning_function_parametric not implemented.",&
              __FILE__, __LINE__)
    end select

  end function turning_function_parametric

  !Deriv wrt parameter
  function dturndparam(param) result(funct)
    implicit none

    real(dp), intent(in) :: param
    real(dp), dimension(num_inflaton) :: funct

    real(dp) :: phi_turn, slope, f_decay

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
      f_decay = vparams(4,1)

      funct(1) = f_decay*cos(param)
      funct(2) = -f_decay*sin(param)
    case(3)
      !One sharp turn
      phi_turn = vparams(4,1)
      slope = vparams(4,2)

      if (param < phi_turn) then
        funct(1) = 1e0_dp
        funct(2) = 0e0_dp
      else
        funct(1) = slope
        funct(2) = 1e0_dp
      end if
    case default
      print*, "MODECODE: turning_choice = ", turning_choice
      call raise%fatal_code(&
              "MODECODE: dturndparam not implemented.",&
              __FILE__, __LINE__)
    end select

  end function dturndparam

  !Second deriv wrt parameter
  function d2turndparam2(param) result(funct)
    implicit none
    real(dp), intent(in) :: param
    real(dp), dimension(num_inflaton) :: funct

    real(dp) :: phi_turn, f_decay

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
      f_decay = vparams(4,1)

      funct(1) = -f_decay*sin(param)
      funct(2) = -f_decay*cos(param)
    case(3)
      !One sharp turn
      phi_turn = vparams(4,1)

      if (param < phi_turn) then
        funct(1) = 0e0_dp
        funct(2) = 0e0_dp
      else
        funct(1) = 0e0_dp
        funct(2) = 0e0_dp
      end if
    case default
      print*, "MODECODE: turning_choice = ", turning_choice
      call raise%fatal_code(&
              "MODECODE: d2turndparam2 not implemented.",&
              __FILE__, __LINE__)
    end select

  end function d2turndparam2

  !Partial derivative of parameter of closest approach for
  !turning_function_parametric with respect to the field values.
  function dparam_closest_dphi(param_closest) result(funct)
    implicit none
    real(dp), intent(in) :: param_closest
    real(dp), dimension(num_inflaton) :: funct

    real(dp) :: phi_turn, slope

    select case(turning_choice(1))
    case(3)
      !One sharp turn
      phi_turn = vparams(4,1)
      slope = vparams(4,2)

      if (param_closest < phi_turn) then
        funct(1) = 1.0e0_dp
        funct(2) = 0.0e0_dp
      else if (param_closest > phi_turn) then
        funct(1) = ((1.0e0_dp+slope**2)**(-1))*slope
        funct(2) = ((1.0e0_dp+slope**2)**(-1))
      else
        funct(1) = 0e0_dp
        funct(2) = 0e0_dp
      end if

    case default

      print*, "MODECODE: turning_choice = ", turning_choice
      call raise%fatal_code(&
              "dparam_closest_dphi not implemented.",&
              __FILE__, __LINE__)

    end select


  end function dparam_closest_dphi

  !Second partial derivative of parameter of closest approach for
  !turning_function_parametric with respect to the field values.
  function d2param_closest_dphi2(param_closest) result(funct)
    implicit none
    real(dp), intent(in) :: param_closest
    real(dp), dimension(num_inflaton, num_inflaton) :: funct

    real(dp) :: phi_turn, slope

    select case(turning_choice(1))
    case(3)
      !One sharp turn
      phi_turn = vparams(4,1)
      slope = vparams(4,2)

      funct = 0e0_dp

    case default

      print*, "MODECODE: turning_choice = ", turning_choice
      call raise%fatal_code(&
              "MODECODE: d2param_closest_dphi2 not implemented.",&
              __FILE__, __LINE__)

    end select


  end function d2param_closest_dphi2

  !Euclidean distance function
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
  subroutine choose_parameter(self, phi_light, param)
    implicit none

    class(qsf_reference) :: self
    integer :: ii
    real(dp), intent(in), optional :: phi_light, param

    !Check for functionality
    if (.not. allocated(self%phi_light_vs_param)) then
      call raise%fatal_code(&
        "Trying to set initial parameter guess, &
        QSF: but haven't integrated trajectory yet.", &
        __FILE__, __LINE__)

    else if (present(phi_light) .and. present(param) ) then
      call raise%fatal_code("Can't have both phi_light and param.",&
        __FILE__, __LINE__)
    end if

    if (present(phi_light)) then

#define LIGHT (self%phi_light_vs_param(:,1))
      ii= locate(LIGHT, phi_light)
      self%hunt_guess = ii
      self%param = self%phi_light_vs_param(ii,2)
#undef LIGHT
    else if (present(param)) then
#define PK_ARR (self%phi_light_vs_param(:,2))
      ii= locate(PK_ARR, param)
      self%hunt_guess = ii
#undef PK_ARR

    else

      call raise%fatal_code(&
        "Give choose_parameter either phi_light or param.  &
        Neither is present.", &
        __FILE__, __LINE__)
    end if

  end subroutine choose_parameter

  !Since the valley is defined parametrically in field-space,
  !we need to find how far up the valley we are for any
  !given parameter.  This integrates the path until length>phi_start
  subroutine integrate_through_trajectory(self, phi_start, &
      param0, stepsize)
    implicit none

    class(qsf_reference) :: self
    real(dp), intent(in) :: param0, phi_start
    real(dp), intent(in), optional :: stepsize

    real(dp) :: dx
    real(dp), dimension(:,:), allocatable :: length
    integer(dp), parameter :: maxsteps=1e7

    integer :: ii, counter
    real(dp) :: x_a, x_b

    !If already initialized, then reset everything
    if (self%traj_init) then
      self%traj_init = .false.
      if (allocated(self%phi_light_vs_param)) then
        deallocate(self%phi_light_vs_param)
      end if
    end if

    if (present(stepsize)) then
      dx = stepsize
    else
      dx = 1e-4_dp
    end if

    !Integrate along parametric curve and record values for interpolation
    allocate(length(maxsteps, 2))
    length(1,:) = (/ 0.0e0_dp, param0 /) !length=(/path_length,param/)

    counter=1
    x_a = length(1,2)
    do ii=1, maxsteps-1
      counter = counter + 1

      x_b = x_a + dx
      length(counter,1) = length(counter-1,1)&
        + integrate_1d(path_segment, x_a, x_b, 10)
      length(counter,2) = x_b
      x_a = x_b

      !Check if have gone far enough
      if (length(counter,1)>2.0e0_dp*phi_start) exit

      if (ii==maxsteps-1) then
        print*, "MODECODE: length of curve=", 2.0e0_dp*phi_start

        call raise%fatal_cosmo(&
          "Couldn't integrate this far. &
          Either lower phi_start or make maxsteps or stepsize bigger.",&
          __FILE__, __LINE__)
      end if

    end do

    allocate(self%phi_light_vs_param(counter,2))
    self%phi_light_vs_param = length(1:counter,:)
    deallocate(length)

    self%traj_init = .true.

  end subroutine integrate_through_trajectory

  !After the trajectory has been integrated, this will give the
  !distance along the trajectory (from param0) to a reference param
  function get_phi_light(self, param, use_locate) result(phi_light)
    implicit none

    class(qsf_reference) :: self
    real(dp), intent(in) :: param
    integer :: param_guess
    real(dp) :: phi_light
    real(dp) :: del_phi

    logical, intent(in), optional :: use_locate

    integer :: ii, jj

    if (.not. self%traj_init) then
      call raise%fatal_code(&
        "QSF: Trying to get phi_light before integrating traj.", &
        __FILE__, __LINE__)
    end if

    param_guess = self%hunt_guess

    if (present(use_locate) .and. use_locate) then
      call locator()
    else
      call hunter()
    end if
    call interpolator()

#define INTLEN size(self%phi_light_vs_param,1)

    !Check for interpolation errors
    if(abs(del_phi) > 0.1 .or. &
      phi_light<self%phi_light_vs_param(1,1) .or. &
      phi_light > self%phi_light_vs_param(INTLEN,1)) then

      if (.not. present(use_locate) .or. .not. use_locate) then
        !Try locate, instead of hunt
        call locator()
        call interpolator()
      end if

    end if

    if(abs(del_phi) > 0.1) then
      print*,"MODECODE: del_phi", del_phi
      call raise%fatal_code(&
        'The interpolation in get_phi_light/locator has &
        suspiciously large errors',&
        __FILE__, __LINE__)

    else if (phi_light > self%phi_light_vs_param(INTLEN,1) .or. &
      phi_light < self%phi_light_vs_param(1,1)) then
      print*, "MODECODE: phi_light =", phi_light
      print*, "MODECODE: LIGHT(MAX) =", self%phi_light_vs_param(INTLEN,1)
      print*, "MODECODE: LIGHT(MIN) =", self%phi_light_vs_param(1,1)
      print*, "MODECODE: param_guess =", param_guess

      call raise%fatal_code(&
        "phi_light is out of bounds.  &
        phi_light > LIGHT(MAX) or < LIGHT(MIN)", &
        __FILE__, __LINE__)
    end if

    self%hunt_guess = param_guess


    contains

      subroutine locator()


#define LIGHT (self%phi_light_vs_param(:,1))
#define P_ARR (self%phi_light_vs_param(:,2))

        !P_ARR is monotonic
        ii= locate(P_ARR, param)

        jj=min(max(ii-(4-1)/2,1),INTLEN+1-4)

#undef P_ARR
#undef LIGHT

      end subroutine

      !Major speed-up
      subroutine hunter()
        integer :: low, high

        low = max(1, param_guess-50)
        high = min(INTLEN, param_guess+50)

#define LIGHT (self%phi_light_vs_param(low:high,1))
#define P_ARR (self%phi_light_vs_param(low:high,2))

        !P_ARR is monotonic
        param_guess = param_guess - low
        call hunt(P_ARR, param, param_guess)
        param_guess = param_guess + low
        ii = param_guess

        jj=min(max(ii-(4-1)/2,1),INTLEN+1-4)

#undef P_ARR
#undef LIGHT

      end subroutine hunter

      subroutine interpolator()
#define LIGHT (self%phi_light_vs_param(jj:jj+4,1))
#define P_ARR (self%phi_light_vs_param(jj:jj+4,2))

        call polint(&
          P_ARR, &
          LIGHT,&
          param, phi_light, del_phi)

      end subroutine
#undef LIGHT
#undef P_ARR
#undef INTLEN

  end function get_phi_light


  !Finds param_closest which minimizes the distance between
  !a field space point and the parametric curve
  !that is used to define the minimum of the QSF valley.  For some parametrized
  !curves this can be done analytically, else we use Newtonian optimization.
  function distance_minimizer(self) result(param_closest)

    class(qsf_reference) :: self
    real(dp) :: param_closest

    real(dp) :: phi_turn, slope
    real(dp) :: t_less, t_great
    real(dp) :: d_less, d_great, d_turn
    logical :: at_corner, fixing

    real(dp) :: radius

    select case(turning_choice(1))
    !case(1)
      !Parabola
    !case(2)
    !  !Helix
    !  if (abs(self%phi(1)) < 1e-5_dp) then
    !    !Protect against divide by zero
    !    if (self%phi(2) > 0e0_dp) then
    !      param_closest = pi/2.0e0_dp
    !    else
    !      param_closest = 3.0e0_dp*pi/2.0e0_dp
    !    end if
    !  else
    !    param_closest = atan(self%phi(2)/self%phi(1))
    !  end if

    !  !Use the param guess to guess the winding number
    !  param_closest = (3.0e0_dp/2.0e0_dp)*pi &
    !    -param_closest &
    !    + (2.0e0_dp*pi)*floor(self%param/2.0e0_dp/pi)

    !  if (abs(param_closest-self%param)>1e0_dp) then

    !    !Try fixing winding number
    !    param_closest = param_closest + pi

    !    if (abs(param_closest-self%param)>1e0_dp) then
    !      param_closest = param_closest - 2.0e0_dp*pi
    !    end if

    !  end if

    !case(2)
    !  !Helix

    !  if (num_inflaton > 2) then
    !    call raise%fatal_code(&
    !            "Analytic solution for closest parameter &
    !            not implemented with more than two fields.", &
    !            __FILE__, __LINE__)
    !  end if

    !  radius = sqrt(self%phi(1)**2 + self%phi(2)**2)
    !  if (radius < 0.0e0_dp) then
    !    call raise%fatal_code(&
    !      "The radius in field space is negative, which is weird.",&
    !      __FILE__, __LINE__)
    !  !else if (radius < 1e-9_dp) then
    !  !  !Divide by zero protection
    !  !  radius = 1e-9_dp
    !  end if
    !  param_closest = acos(self%phi(1)/radius)
    !  param_closest = (pi/2.0e0_dp) - param_closest !Since param=0 at f*(0,1)

    !  !Guess winding number
    !  param_closest = 2.0e0_dp*pi*ceiling(self%param/2.0e0_dp/pi) + param_closest

    !  !DEBUG
    !  !print*, "this is phi", self%phi, param_closest
    !  print*, self%phi, param_closest
    !  !print*, "this is param_closest", param_closest
    !  !print*, "this is recnostruct", radius*cos(param_closest)
    !  !print*, "this is self%param", self%param, radius*cos(self%param)
    !  !stop


    case(3)

      !One sharp turn
      phi_turn = vparams(4,1)
      slope = vparams(4,2)

      t_less = self%phi(1)
      t_great = (1.0e0_dp + slope**2)**(-1)*&
        (slope*self%phi(1) +&
        slope**2*phi_turn - &
        slope*phi_turn + &
        self%phi(2) + phi_turn)

      !Try different distances
      d_less = distance(t_less)
      d_great = distance(t_great)
      d_turn = distance(phi_turn)
      if (d_less <= d_great .and. d_less <=d_turn) then
        param_closest = t_less
      else if (d_great <= d_less .and. d_great <=d_turn) then
        param_closest = t_great
      else
        param_closest = phi_turn
      end if


    case default
      !Use Newtonian optimization
      param_closest = zero_finder(distance_deriv, &
        distance_2deriv, self%param)
    end select

  end function distance_minimizer

  !Derivative of the light field with respect to the parameter at which we
  !perform the interpolation of the pre-integrated trajectory.
  function dphi_light_dparam(self,param) result(dot_phi)
    implicit none

    class(qsf_reference) :: self
    real(dp) :: dot_phi
    real(dp), intent(in) :: param
    real(dp), dimension(num_inflaton) :: dfunct

    integer :: i

    dfunct = dturndparam(param)
    dot_phi = sqrt( sum( dfunct**2))

  end function dphi_light_dparam

  !Second derivative of the light field with respect to the parameter at which we
  !perform the interpolation of the pre-integrated trajectory.
  function d2phi_light_dparam2(self,param) result(ddot_phi)
    implicit none

    class(qsf_reference) :: self
    real(dp) :: ddot_phi
    real(dp), intent(in) :: param
    real(dp), dimension(num_inflaton) :: dfunct
    real(dp), dimension(num_inflaton) :: d2funct
    real(dp) :: dot_phi

    integer :: i

    dfunct = dturndparam(param)
    d2funct = d2turndparam2(param)
    dot_phi = self%dphi_light_dparam(param)

    ddot_phi = (1.0e0_dp/dot_phi)* sum( dfunct*d2funct)

  end function d2phi_light_dparam2


end module modpk_qsf
