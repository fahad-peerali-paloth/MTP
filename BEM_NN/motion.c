#include "udf.h"
#include "dynamesh_tools.h"
#include "storage.h"
DEFINE_CG_MOTION(object_motion, dt, vel, omega, time, dtime)
{
real a, b, w, pi;
real v = 0;
pi = 3.1415;

/* define motion variables */
a = -0;
b = -0;
w = 2*pi*0.4023;  /* 2Hz frequency */

/* define object movement law */

v = a*cos((w*time));
vel[0] = 0;
vel[1] = v;
vel[2] = 0;

omega[2] = b*cos((w*time));

}