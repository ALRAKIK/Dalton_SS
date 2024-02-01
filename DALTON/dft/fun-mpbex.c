/*


!
!  Dalton, a molecular electronic structure program
!  Copyright (C) The Dalton Authors (see AUTHORS file for details).
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!

!

*/
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-mpbex.c:

   Automatically generated code implementing MPBEX functional and
   its derivatives. It is generated by func-codegen.pl being a part of
   a "Automatic code generation framework for analytical functional
   derivative evaluation", Pawel Salek, 2005

    This functional is connected by making following changes:
    1. add "extern Functional mpbexFunctional;" to 'functionals.h'
    2. add "&mpbexFunctional," to 'functionals.c'
    3. add "fun-mpbex.c" to 'Makefile.am', 'Makefile.in' or 'Makefile'.

    This functional has been generated from following input:
    ------ cut here -------
xa:sqrt(grada*grada)/rhoa^(4/3);
xb:sqrt(gradb*gradb)/rhob^(4/3);
 
a:0.157;
C1:0.21951;
C2:-0.015;
R:0.804;
d:0.066725;
mu:d*%PI^2/3;
Sa:xa/(2*(6*%PI^2)^(1/3));
Sb:xb/(2*(6*%PI^2)^(1/3));
 
G(T):=T^2/(1+a*T^2);
F(S):=1+C1*G(S)+C2*G(S)^2;
Ea(n):=-3/(4*%PI)*(3*%PI^2)^(1/3)*n^(4/3)*F(Sa);
Eb(n):=-3/(4*%PI)*(3*%PI^2)^(1/3)*n^(4/3)*F(Sb);
 
K(rhoa,grada,rhob,gradb,gradab):=0.5*(Ea(2*rhoa)+Eb(2*rhob));


    ------ cut here -------
*/

 
/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas and pow() prototype
 * is not specified. */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif
#include <math.h>
#include <stddef.h>
#include "general.h"

#define __CVERSION__
 
#include "functionals.h"
 
/* INTERFACE PART */
static integer mpbex_isgga(void) { return 1; } /* FIXME: detect! */
static integer mpbex_read(const char *conf_line);
static real mpbex_energy(const FunDensProp* dp);
static void mpbex_first(FunFirstFuncDrv *ds,   real factor,
                         const FunDensProp* dp);
static void mpbex_second(FunSecondFuncDrv *ds, real factor,
                          const FunDensProp* dp);
static void mpbex_third(FunThirdFuncDrv *ds,   real factor,
                         const FunDensProp* dp);
static void mpbex_fourth(FunFourthFuncDrv *ds,   real factor,
                          const FunDensProp* dp);
 
Functional mPBExFunctional = {
  "mPBEx",       /* name */
  mpbex_isgga,   /* gga-corrected */
   1,
  mpbex_read,
  NULL,
  mpbex_energy,
  mpbex_first,
  mpbex_second,
  mpbex_third,
  mpbex_fourth
};
 
/* IMPLEMENTATION PART */
static integer
mpbex_read(const char *conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

static real
mpbex_energy(const FunDensProp *dp)
{
    real res;
    real rhoa = dp->rhoa, rhob = dp->rhob;
    real grada = dp->grada, gradb = dp->gradb, gradab = dp->gradab;

    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12, t13;

    t1 = pow(2.0,0.333333333333333);
    t2 = pow(3.0,0.333333333333333);
    t3 = 1/pow(3.141592653589793,0.333333333333333);
    t4 = 1/pow(6.0,0.333333333333333);
    t5 = 1/pow(3.141592653589793,2.666666666666667);
    t6 = 1/pow(6.0,0.666666666666667);
    t7 = 1/pow(3.141592653589793,1.333333333333333);
    t8 = pow(grada,2.0);
    t9 = 1/pow(rhoa,2.666666666666667);
    t10 = 0.03925*t6*t7*t8*t9+1.0;
    t11 = pow(gradb,2.0);
    t12 = 1/pow(rhob,2.666666666666667);
    t13 = 0.03925*t6*t7*t11*t12+1.0;

   /* code */
    res = 0.5*(-1.5*t1*t2*t3*pow(rhob,1.333333333333333)*
        (-1.5624999999999998E-4*t4*t5*pow(gradb,4.0)/(pow(t13,2.0)*
        pow(rhob,5.333333333333333))+0.0548775*t11*t12*t6*t7/t13+1.0)-
        1.5*t1*t2*t3*pow(rhoa,1.333333333333333)*(-1.5624999999999998E-4*
        t4*t5*pow(grada,4.0)/(pow(t10,2.0)*pow(rhoa,5.333333333333333))+
        0.0548775*t6*t7*t8*t9/t10+1.0));

    return res;
}

static void
mpbex_first_helper(real rhoa, real grada, real *res)
{    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12, t13, t14, t15, t16, t17;

    t1 = pow(2.0,0.333333333333333);
    t2 = pow(3.0,0.333333333333333);
    t3 = 1/pow(3.141592653589793,0.333333333333333);
    t4 = 1/pow(6.0,0.333333333333333);
    t5 = 1/pow(3.141592653589793,2.666666666666667);
    t6 = pow(grada,4.0);
    t7 = 1/pow(6.0,0.666666666666667);
    t8 = 1/pow(3.141592653589793,1.333333333333333);
    t9 = pow(grada,2.0);
    t10 = 1/pow(rhoa,2.666666666666667);
    t11 = 0.03925*t7*t8*t9*t10+1.0;
    t12 = 1/pow(t11,2.0);
    t13 = 1/pow(rhoa,5.333333333333333);
    t14 = 1/t11;
    t15 = 1/pow(3.141592653589793,4.0);
    t16 = 1/pow(t11,3.0);
    t17 = pow(rhoa,1.333333333333333);

   /* code */
    res[0] = 0.5*(-1.5*t1*t17*t2*t3*(-5.451388888888888E-6*
        t15*t16*pow(grada,6.0)/pow(rhoa,9.0)+0.001790640833333*t12*
        t4*t5*t6/pow(rhoa,6.333333333333333)-0.14634*t14*t7*t8*t9/
        pow(rhoa,3.666666666666667))-2.0*t1*(-1.5624999999999998E-4*
        t4*t5*t6*t12*t13+0.0548775*t7*t8*t9*t14*t10+1.0)*t2*t3*pow(rhoa,
        0.333333333333333));
    res[1] = -0.75*t1*t17*t2*t3*(4.088541666666666E-6*t15*
        t16*pow(grada,5.0)/pow(rhoa,8.0)-0.001342980625*t12*t13*t4*
        t5*pow(grada,3.0)+0.109755*t7*t8*grada*t14*t10);
}

static void
mpbex_first(FunFirstFuncDrv *ds, real factor, const FunDensProp *dp)
{
    real res[2];

    mpbex_first_helper(dp->rhoa, dp->grada, res);
   /* Final assignment */
    ds->df1000 += factor*res[0];
    ds->df0010 += factor*res[1];


    if(fabs(dp->rhoa-dp->rhob)>1e-13 ||
       fabs(dp->grada-dp->gradb)>1e-13)
        mpbex_first_helper(dp->rhob, dp->gradb, res);
    ds->df0100 += factor*res[0];
    ds->df0001 += factor*res[1];

}

static void
mpbex_second_helper(real rhoa, real grada, real *res)
{
    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12, t13, t14, t15, t16, t17, t18;
    real t19, t20, t21, t22, t23, t24, t25, t26;
    real t27, t28, t29, t30;

    t1 = pow(2.0,0.333333333333333);
    t2 = pow(3.0,0.333333333333333);
    t3 = 1/pow(3.141592653589793,0.333333333333333);
    t4 = 1/pow(6.0,0.333333333333333);
    t5 = 1/pow(3.141592653589793,2.666666666666667);
    t6 = pow(grada,4.0);
    t7 = 1/pow(6.0,0.666666666666667);
    t8 = 1/pow(3.141592653589793,1.333333333333333);
    t9 = pow(grada,2.0);
    t10 = 1/pow(rhoa,2.666666666666667);
    t11 = 0.03925*t7*t8*t9*t10+1.0;
    t12 = 1/pow(t11,2.0);
    t13 = 1/pow(rhoa,5.333333333333333);
    t14 = 1/t11;
    t15 = -1.5624999999999998E-4*t4*t5*t6*t12*t13+0.0548775*
        t7*t8*t9*t14*t10+1.0;
    t16 = pow(rhoa,0.333333333333333);
    t17 = 1/pow(3.141592653589793,4.0);
    t18 = pow(grada,6.0);
    t19 = 1/pow(t11,3.0);
    t20 = 1/pow(rhoa,9.0);
    t21 = 1/pow(rhoa,6.333333333333333);
    t22 = 1/pow(rhoa,3.666666666666667);
    t23 = -0.14634*t7*t8*t9*t14*t22+0.001790640833333*t4*
        t5*t6*t12*t21-5.451388888888888E-6*t17*t18*t19*t20;
    t24 = pow(rhoa,1.333333333333333);
    t25 = pow(grada,5.0);
    t26 = 1/pow(rhoa,8.0);
    t27 = pow(grada,3.0);
    t28 = 0.109755*t7*t8*grada*t14*t10-0.001342980625*t4*
        t5*t27*t12*t13+4.088541666666666E-6*t17*t25*t19*t26;
    t29 = 1/pow(3.141592653589793,5.333333333333333);
    t30 = 1/pow(t11,4.0);

   /* code */
    res[0] = 0.5*(-1.5*t1*t2*t23*t24*t3-2.0*t1*t15*t16*t2*
        t3);
    res[1] = -0.75*t1*t2*t3*t28*t24;
    res[2] = 0.5*(-1.5*t1*t2*t24*t3*(-1.7117361111111105E-6*
        t29*t30*t7*pow(grada,8.0)/pow(rhoa,12.66666666666667)+1.1153596907407406E-4*
        t17*t18*t19/pow(rhoa,10.0)-0.013893545277778*t12*t4*t5*t6/
        pow(rhoa,7.333333333333333)+0.53658*t14*t7*t8*t9/pow(rhoa,
        4.666666666666667))-0.666666666666667*t1*t15*t2*t3/pow(rhoa,
        0.666666666666667)-4.0*t1*t16*t2*t23*t3);
    res[3] = 0.5*(-1.5*t1*t2*t24*t3*(1.283802083333333E-6*
        t29*t30*t7*pow(grada,7.0)/pow(rhoa,11.66666666666667)-0.29268*
        t7*t8*grada*t14*t22+0.009077178333333*t4*t5*t27*t12*t21-7.956343513888887E-5*
        t17*t25*t19*t20)-2.0*t1*t16*t2*t28*t3);
    res[4] = -0.75*t1*t2*t24*t3*(-9.628515624999997E-7*t18*
        t29*t30*t7/pow(rhoa,10.66666666666667)+5.558403468749999E-5*
        t17*t6*t19*t26-0.005464903125*t4*t5*t9*t12*t13+0.109755*t7*
        t8*t14*t10);

}

static void
mpbex_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real res[5];
 
    mpbex_second_helper(dp->rhoa, dp->grada, res);

    ds->df1000 += factor*res[0];
    ds->df0010 += factor*res[1];

    ds->df2000 += factor*res[2];
    ds->df1010 += factor*res[3];
    ds->df0020 += factor*res[4];


    if(fabs(dp->rhoa-dp->rhob)>1e-13 ||
       fabs(dp->grada-dp->gradb)>1e-13)
        mpbex_second_helper(dp->rhob, dp->gradb, res);
    ds->df0100 += factor*res[0];
    ds->df0001 += factor*res[1];

    ds->df0200 += factor*res[2];
    ds->df0101 += factor*res[3];
    ds->df0002 += factor*res[4];

}

static void
mpbex_third_helper(real rhoa, real grada, real *res)
{
    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12, t13, t14, t15, t16, t17, t18;
    real t19, t20, t21, t22, t23, t24, t25, t26;
    real t27, t28, t29, t30, t31, t32, t33, t34;
    real t35, t36, t37, t38, t39, t40, t41, t42;
    real t43, t44;

    t1 = pow(2.0,0.333333333333333);
    t2 = pow(3.0,0.333333333333333);
    t3 = 1/pow(3.141592653589793,0.333333333333333);
    t4 = 1/pow(6.0,0.333333333333333);
    t5 = 1/pow(3.141592653589793,2.666666666666667);
    t6 = pow(grada,4.0);
    t7 = 1/pow(6.0,0.666666666666667);
    t8 = 1/pow(3.141592653589793,1.333333333333333);
    t9 = pow(grada,2.0);
    t10 = 1/pow(rhoa,2.666666666666667);
    t11 = 0.03925*t7*t8*t9*t10+1.0;
    t12 = 1/pow(t11,2.0);
    t13 = 1/pow(rhoa,5.333333333333333);
    t14 = 1/t11;
    t15 = -1.5624999999999998E-4*t4*t5*t6*t12*t13+0.0548775*
        t7*t8*t9*t14*t10+1.0;
    t16 = pow(rhoa,0.333333333333333);
    t17 = 1/pow(3.141592653589793,4.0);
    t18 = pow(grada,6.0);
    t19 = 1/pow(t11,3.0);
    t20 = 1/pow(rhoa,9.0);
    t21 = 1/pow(rhoa,6.333333333333333);
    t22 = 1/pow(rhoa,3.666666666666667);
    t23 = -0.14634*t7*t8*t9*t14*t22+0.001790640833333*t4*
        t5*t6*t12*t21-5.451388888888888E-6*t17*t18*t19*t20;
    t24 = pow(rhoa,1.333333333333333);
    t25 = pow(grada,5.0);
    t26 = 1/pow(rhoa,8.0);
    t27 = pow(grada,3.0);
    t28 = 0.109755*t7*t8*grada*t14*t10-0.001342980625*t4*
        t5*t27*t12*t13+4.088541666666666E-6*t17*t25*t19*t26;
    t29 = 1/pow(rhoa,0.666666666666667);
    t30 = 1/pow(3.141592653589793,5.333333333333333);
    t31 = pow(grada,8.0);
    t32 = 1/pow(t11,4.0);
    t33 = 1/pow(rhoa,12.66666666666667);
    t34 = 1/pow(rhoa,10.0);
    t35 = 1/pow(rhoa,7.333333333333333);
    t36 = 1/pow(rhoa,4.666666666666667);
    t37 = 0.53658*t7*t8*t9*t14*t36-0.013893545277778*t4*t5*
        t6*t12*t35+1.1153596907407406E-4*t17*t18*t19*t34-1.7117361111111105E-6*
        t7*t30*t31*t32*t33;
    t38 = pow(grada,7.0);
    t39 = 1/pow(rhoa,11.66666666666667);
    t40 = -0.29268*t7*t8*grada*t14*t22+0.009077178333333*
        t4*t5*t27*t12*t21-7.956343513888887E-5*t17*t25*t19*t20+1.283802083333333E-6*
        t7*t30*t38*t32*t39;
    t41 = 1/pow(rhoa,10.66666666666667);
    t42 = 0.109755*t7*t8*t14*t10-0.005464903125*t4*t5*t9*
        t12*t13+5.558403468749999E-5*t17*t6*t19*t26-9.628515624999997E-7*
        t7*t30*t18*t32*t41;
    t43 = 1/pow(3.141592653589793,6.666666666666667);
    t44 = 1/pow(t11,5.0);

   /* code */
    res[0] = 0.5*(-1.5*t1*t2*t23*t24*t3-2.0*t1*t15*t16*t2*
        t3);
    res[1] = -0.75*t1*t2*t3*t28*t24;
    res[2] = 0.5*(-1.5*t1*t2*t24*t3*t37-0.666666666666667*
        t1*t15*t2*t29*t3-4.0*t1*t16*t2*t23*t3);
    res[3] = 0.5*(-1.5*t1*t2*t24*t3*t40-2.0*t1*t16*t2*t28*
        t3);
    res[4] = -0.75*t1*t2*t3*t42*t24;
    res[5] = 0.5*(-1.5*t1*t2*t24*t3*(-1.194411419753086E-7*
        t4*t43*t44*pow(grada,10.0)/pow(rhoa,16.33333333333333)+5.670428502999997E-5*
        t30*t31*t32*t7/pow(rhoa,13.66666666666667)-0.00160009004821*
        t17*t18*t19/pow(rhoa,11.0)+0.111246338703704*t12*t4*t5*t6/
        pow(rhoa,8.333333333333334)-2.50404*t14*t7*t8*t9/pow(rhoa,
        5.666666666666667))+0.444444444444444*t1*t15*t2*t3/pow(rhoa,
        1.666666666666667)-6.0*t1*t16*t2*t3*t37-2.0*t1*t2*t23*t29*
        t3);
    res[6] = 0.5*(-1.5*t1*t2*t24*t3*(8.958085648148145E-8*
        t4*t43*t44*pow(grada,9.0)/pow(rhoa,15.33333333333333)+1.07316*
        t7*t8*grada*t14*t36-0.062594436111111*t4*t5*t27*t12*t35+0.001032763582546*
        t17*t25*t19*t34-3.996060960583333E-5*t7*t30*t38*t32*t33)-4.0*
        t1*t16*t2*t3*t40-0.666666666666667*t1*t2*t28*t29*t3);
    res[7] = 0.5*(-1.5*t1*t2*t24*t3*(-6.718564236111108E-8*
        t31*t4*t43*t44/pow(rhoa,14.33333333333333)+2.772380355854166E-5*
        t7*t30*t18*t32*t39-0.29268*t7*t8*t14*t22+0.031060765*t4*t5*
        t9*t12*t21-6.353366754166664E-4*t17*t6*t19*t20)-2.0*t1*t16*
        t2*t3*t42);
    res[8] = -0.75*t1*t2*t24*t3*(5.038923177083331E-8*t38*
        t4*t43*t44/pow(rhoa,13.33333333333333)-1.8867149543906244E-5*
        t7*t30*t25*t32*t41+3.653344371874999E-4*t17*t27*t19*t26-0.0123657675*
        t4*t5*grada*t12*t13);

}

static void
mpbex_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real res[9];
 
    mpbex_third_helper(dp->rhoa, dp->grada, res);

    ds->df1000 += factor*res[0];
    ds->df0010 += factor*res[1];

    ds->df2000 += factor*res[2];
    ds->df1010 += factor*res[3];
    ds->df0020 += factor*res[4];

    ds->df3000 += factor*res[5];
    ds->df2010 += factor*res[6];
    ds->df1020 += factor*res[7];
    ds->df0030 += factor*res[8];


    if(fabs(dp->rhoa-dp->rhob)>1e-13 ||
       fabs(dp->grada-dp->gradb)>1e-13)
        mpbex_third_helper(dp->rhob, dp->gradb, res);

    ds->df0100 += factor*res[0];
    ds->df0001 += factor*res[1];

    ds->df0200 += factor*res[2];
    ds->df0101 += factor*res[3];
    ds->df0002 += factor*res[4];

    ds->df0300 += factor*res[5];
    ds->df0201 += factor*res[6];
    ds->df0102 += factor*res[7];
    ds->df0003 += factor*res[8];

}

static void
mpbex_fourth_helper(real rhoa, real grada, real *res)
{
    real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    real t11, t12, t13, t14, t15, t16, t17, t18;
    real t19, t20, t21, t22, t23, t24, t25, t26;
    real t27, t28, t29, t30, t31, t32, t33, t34;
    real t35, t36, t37, t38, t39, t40, t41, t42;
    real t43, t44, t45, t46, t47, t48, t49, t50;
    real t51, t52, t53, t54, t55, t56, t57, t58;
    real t59, t60, t61;

    t1 = pow(2.0,0.333333333333333);
    t2 = pow(3.0,0.333333333333333);
    t3 = 1/pow(3.141592653589793,0.333333333333333);
    t4 = 1/pow(6.0,0.333333333333333);
    t5 = 1/pow(3.141592653589793,2.666666666666667);
    t6 = pow(grada,4.0);
    t7 = 1/pow(6.0,0.666666666666667);
    t8 = 1/pow(3.141592653589793,1.333333333333333);
    t9 = pow(grada,2.0);
    t10 = 1/pow(rhoa,2.666666666666667);
    t11 = 0.03925*t7*t8*t9*t10+1.0;
    t12 = 1/pow(t11,2.0);
    t13 = 1/pow(rhoa,5.333333333333333);
    t14 = 1/t11;
    t15 = -1.5624999999999998E-4*t4*t5*t6*t12*t13+0.0548775*
        t7*t8*t9*t14*t10+1.0;
    t16 = pow(rhoa,0.333333333333333);
    t17 = 1/pow(3.141592653589793,4.0);
    t18 = pow(grada,6.0);
    t19 = 1/pow(t11,3.0);
    t20 = 1/pow(rhoa,9.0);
    t21 = 1/pow(rhoa,6.333333333333333);
    t22 = 1/pow(rhoa,3.666666666666667);
    t23 = -0.14634*t7*t8*t9*t14*t22+0.001790640833333*t4*
        t5*t6*t12*t21-5.451388888888888E-6*t17*t18*t19*t20;
    t24 = pow(rhoa,1.333333333333333);
    t25 = pow(grada,5.0);
    t26 = 1/pow(rhoa,8.0);
    t27 = pow(grada,3.0);
    t28 = 0.109755*t7*t8*grada*t14*t10-0.001342980625*t4*
        t5*t27*t12*t13+4.088541666666666E-6*t17*t25*t19*t26;
    t29 = 1/pow(rhoa,0.666666666666667);
    t30 = 1/pow(3.141592653589793,5.333333333333333);
    t31 = pow(grada,8.0);
    t32 = 1/pow(t11,4.0);
    t33 = 1/pow(rhoa,12.66666666666667);
    t34 = 1/pow(rhoa,10.0);
    t35 = 1/pow(rhoa,7.333333333333333);
    t36 = 1/pow(rhoa,4.666666666666667);
    t37 = 0.53658*t7*t8*t9*t14*t36-0.013893545277778*t4*t5*
        t6*t12*t35+1.1153596907407406E-4*t17*t18*t19*t34-1.7117361111111105E-6*
        t7*t30*t31*t32*t33;
    t38 = pow(grada,7.0);
    t39 = 1/pow(rhoa,11.66666666666667);
    t40 = -0.29268*t7*t8*grada*t14*t22+0.009077178333333*
        t4*t5*t27*t12*t21-7.956343513888887E-5*t17*t25*t19*t20+1.283802083333333E-6*
        t7*t30*t38*t32*t39;
    t41 = 1/pow(rhoa,10.66666666666667);
    t42 = 0.109755*t7*t8*t14*t10-0.005464903125*t4*t5*t9*
        t12*t13+5.558403468749999E-5*t17*t6*t19*t26-9.628515624999997E-7*
        t7*t30*t18*t32*t41;
    t43 = 1/pow(rhoa,1.666666666666667);
    t44 = 1/pow(3.141592653589793,6.666666666666667);
    t45 = pow(grada,10.0);
    t46 = 1/pow(t11,5.0);
    t47 = 1/pow(rhoa,16.33333333333333);
    t48 = 1/pow(rhoa,13.66666666666667);
    t49 = 1/pow(rhoa,11.0);
    t50 = 1/pow(rhoa,8.333333333333334);
    t51 = 1/pow(rhoa,5.666666666666667);
    t52 = -2.50404*t7*t8*t9*t14*t51+0.111246338703704*t4*
        t5*t6*t12*t50-0.00160009004821*t17*t18*t19*t49+5.670428502999997E-5*
        t7*t30*t31*t32*t48-1.194411419753086E-7*t4*t44*t45*t46*t47;
    t53 = pow(grada,
        9.0);
    t54 = 1/pow(rhoa,15.33333333333333);
    t55 = 1.07316*t7*t8*grada*t14*t36-0.062594436111111*t4*
        t5*t27*t12*t35+0.001032763582546*t17*t25*t19*t34-3.996060960583333E-5*
        t7*t30*t38*t32*t33+8.958085648148145E-8*t4*t44*t53*t46*t54;
    t56 = 1/
        pow(rhoa,14.33333333333333);
    t57 = -0.29268*t7*t8*t14*t22+0.031060765*t4*t5*t9*t12*
        t21-6.353366754166664E-4*t17*t6*t19*t20+2.772380355854166E-5*
        t7*t30*t18*t32*t39-6.718564236111108E-8*t4*t44*t31*t46*t56;
    t58 = 1/
        pow(rhoa,13.33333333333333);
    t59 = -0.0123657675*t4*t5*grada*t12*t13+3.653344371874999E-4*
        t17*t27*t19*t26-1.8867149543906244E-5*t7*t30*t25*t32*t41+5.038923177083331E-8*
        t4*t44*t38*t46*t58;
    t60 = 1/pow(3.141592653589793,8.0);
    t61 = 1/pow(t11,6.0);

   /* code */
    res[0] = 0.5*(-1.5*t1*t2*t23*t24*t3-2.0*t1*t15*t16*t2*
        t3);
    res[1] = -0.75*t1*t2*t3*t28*t24;
    res[2] = 0.5*(-1.5*t1*t2*t24*t3*t37-0.666666666666667*
        t1*t15*t2*t29*t3-4.0*t1*t16*t2*t23*t3);
    res[3] = 0.5*(-1.5*t1*t2*t24*t3*t40-2.0*t1*t16*t2*t28*
        t3);
    res[4] = -0.75*t1*t2*t3*t42*t24;
    res[5] = 0.5*(-1.5*t1*t2*t24*t3*t52+0.444444444444444*
        t1*t15*t2*t3*t43-6.0*t1*t16*t2*t3*t37-2.0*t1*t2*t23*t29*t3);
    res[6] = 0.5*(-1.5*t1*t2*t24*t3*t55-4.0*t1*t16*t2*t3*
        t40-0.666666666666667*t1*t2*t28*t29*t3);
    res[7] = 0.5*(-1.5*t1*t2*t24*t3*t57-2.0*t1*t16*t2*t3*
        t42);
    res[8] = -0.75*t1*t2*t3*t59*t24;
    res[9] = 0.5*(-1.5*t1*t2*t24*t3*(-1.041792182784636E-8*
        t60*t61*pow(grada,12.0)/pow(rhoa,20.0)+5.907570985467816E-6*
        t4*t44*t45*t46/pow(rhoa,17.33333333333333)-0.001277386837215*
        t30*t31*t32*t7/pow(rhoa,14.66666666666667)+0.021482251680638*
        t17*t18*t19/pow(rhoa,12.0)-0.970734409197531*t12*t4*t5*t6/
        pow(rhoa,9.333333333333334)+14.18956*t14*t7*t8*t9/pow(rhoa,
        6.666666666666667))-8.0*t1*t16*t2*t3*t52+1.777777777777778*
        t1*t2*t23*t3*t43-4.0*t1*t2*t29*t3*t37-0.740740740740741*t1*
        t10*t15*t2*t3);
    res[10] = 0.5*(-1.5*t1*t2*t24*t3*(7.81344137088477E-9*
        t60*t61*pow(grada,11.0)/pow(rhoa,19.0)-5.00808*t7*t8*grada*
        t14*t51+0.477746544814815*t4*t5*t27*t12*t50-0.012511486152006*
        t17*t25*t19*t49+8.304554865934256E-4*t7*t30*t38*t32*t48-4.161935669656417E-6*
        t4*t44*t53*t46*t47)-6.0*t1*t16*t2*t3*t55+0.444444444444444*
        t1*t2*t28*t3*t43-2.0*t1*t2*t29*t3*t40);
    res[11] = 0.5*(-1.5*t1*t2*t24*t3*(-5.860081028163578E-9*
        t45*t60*t61/pow(rhoa,18.0)+2.89749961103861E-6*t4*t44*t31*
        t46*t54+1.07316*t7*t8*t14*t36-0.201823818333333*t4*t5*t9*t12*
        t35+0.006801705657639*t17*t6*t19*t34-5.22940090930486E-4*t7*
        t30*t18*t32*t33)-4.0*t1*t16*t2*t3*t57-0.666666666666667*t1*
        t2*t29*t3*t42);
    res[12] = 0.5*(-1.5*t1*t2*t24*t3*(4.395060771122683E-9*
        t53*t60*t61/pow(rhoa,17.0)-1.988364191785902E-6*t4*t44*t38*
        t46*t56+3.159646084118749E-4*t7*t30*t25*t32*t39+0.06595076*
        t4*t5*grada*t12*t21-0.003354103385833*t17*t27*t19*t20)-2.0*
        t1*t16*t2*t3*t59);
    res[13] = -0.75*t1*t2*t24*t3*(-3.296295578342012E-9*t31*
        t60*t61/pow(rhoa,16.0)+1.3401054485269264E-6*t4*t44*t18*t46*
        t58-1.8037200767718743E-4*t7*t30*t6*t32*t41+0.001419574227812*
        t17*t9*t19*t26-0.0123657675*t4*t5*t12*t13);

}

static void
mpbex_fourth(FunFourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real res[14];
 
    mpbex_fourth_helper(dp->rhoa, dp->grada, res);

    ds->df1000 += factor*res[0];
    ds->df0010 += factor*res[1];

    ds->df2000 += factor*res[2];
    ds->df1010 += factor*res[3];
    ds->df0020 += factor*res[4];

    ds->df3000 += factor*res[5];
    ds->df2010 += factor*res[6];
    ds->df1020 += factor*res[7];
    ds->df0030 += factor*res[8];

    ds->df4000 += factor*res[9];
    ds->df3010 += factor*res[10];
    ds->df2020 += factor*res[11];
    ds->df1030 += factor*res[12];
    ds->df0040 += factor*res[13];


    if(fabs(dp->rhoa-dp->rhob)>1e-13 ||
       fabs(dp->grada-dp->gradb)>1e-13)
        mpbex_fourth_helper(dp->rhob, dp->gradb, res);

    ds->df0100 += factor*res[0];
    ds->df0001 += factor*res[1];

    ds->df0200 += factor*res[2];
    ds->df0101 += factor*res[3];
    ds->df0002 += factor*res[4];

    ds->df0300 += factor*res[5];
    ds->df0201 += factor*res[6];
    ds->df0102 += factor*res[7];
    ds->df0003 += factor*res[8];

    ds->df0400 += factor*res[9];
    ds->df0301 += factor*res[10];
    ds->df0202 += factor*res[11];
    ds->df0103 += factor*res[12];
    ds->df0004 += factor*res[13];

}
