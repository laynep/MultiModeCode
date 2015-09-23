!Macros for use in global environment.

!-------------------
!Indices for reference in main y-vector for ODES with
!dydN = f(y)
!For all
#define IND_FIELDS 1:num_inflaton
#define IND_VEL num_inflaton+1:2*num_inflaton

!For the background when including radn fluid
#define IND_RADN 2*num_inflaton+2:3*num_inflaton+1
#define IND_EFOLDS 2*num_inflaton+1

!For the modes
#define IND_MODES 2*num_inflaton+1:2*num_inflaton+num_inflaton**2
#define IND_MODES_VEL 2*num_inflaton+1+num_inflaton**2:2*num_inflaton+2*num_inflaton**2+1
#define IND_TENSOR 2*num_inflaton+2*num_inflaton**2+1
#define IND_TENSOR_VEL 2*num_inflaton+2*num_inflaton**2+2
#define IND_UZETA 2*num_inflaton+2*num_inflaton**2 + 3
#define IND_UZETA_VEL 2*num_inflaton+2*num_inflaton**2 + 4

!For looping
#define LOOP_MODES 2*num_inflaton+1,2*num_inflaton+num_inflaton**2

!For variables that have constraints, but not typically evolved
!in the ODES, e.g., 0<epsilon<3
!For use with dvode integrator only
#define IND_CONST_EPS_MODES 2*num_inflaton+2*num_inflaton**2 + 5
#define IND_CONST_EPS_BACK 2*num_inflaton+1
#define IND_CONST_EPS_RADN 3*num_inflaton+2
