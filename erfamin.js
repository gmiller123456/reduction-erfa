'use strict';;
var ERFA = {};;
(function(era) {;
/* Pi */
var ERFA_DPI = 3.141592653589793238462643;

/* 2Pi */
var ERFA_D2PI = 6.283185307179586476925287;

/* Radians to degrees */
var ERFA_DR2D = 57.29577951308232087679815;

/* Degrees to radians */
var ERFA_DD2R = 1.745329251994329576923691e-2;

/* Radians to arcseconds */
var ERFA_DR2AS = 206264.8062470963551564734;

/* Arcseconds to radians */
var ERFA_DAS2R = 4.848136811095359935899141e-6;

/* Seconds of time to radians */
var ERFA_DS2R = 7.272205216643039903848712e-5;

/* Arcseconds in a full circle */
var ERFA_TURNAS = 1296000.0;

/* Milliarcseconds to radians */
var ERFA_DMAS2R = ERFA_DAS2R / 1e3;

/* Length of tropical year B1900 (days) */
var ERFA_DTY = 365.242198781;

/* Seconds per day. */
var ERFA_DAYSEC = 86400.0;

/* Days per Julian year */
var ERFA_DJY = 365.25;

/* Days per Julian century */
var ERFA_DJC = 36525.0;

/* Days per Julian millennium */
var ERFA_DJM = 365250.0;

/* Reference epoch = J2000.0;, Julian Date */
var ERFA_DJ00 = 2451545.0;

/* Julian Date of Modified Julian Date zero */
var ERFA_DJM0 = 2400000.5;

/* Reference epoch (J2000.0), Modified Julian Date */
var ERFA_DJM00 = 51544.5;

/* 1977 Jan 1.0 as MJD */
var ERFA_DJM77 = 43144.0;

/* TT minus TAI (s) */
var ERFA_TTMTAI = 32.184;

/* Astronomical unit (m) */
var ERFA_DAU = 149597870.7e3;

/* Speed of light (m/s) */
var ERFA_CMPS = 299792458.0;

/* Light time for 1 au (s) */
var ERFA_AULT = ERFA_DAU / ERFA_CMPS;

/* Speed of light (AU per day) */
var ERFA_DC = ERFA_DAYSEC / ERFA_AULT;

/* L_G = 1 - d(TT)/d(TCG) */
var ERFA_ELG = 6.969290134e-10;

/* L_B = 1 - d(TDB)/d(TCB), and TDB (s) at TAI 1977/1/1.0 */
var ERFA_ELB = 1.550519768e-8;
var ERFA_TDB0 = -6.55e-5;

/* Schwarzschild radius of the Sun (au) */
/* = 2 * 1.32712440041e20 / (2.99792458e8)^2 / 1.49597870700e11 */
var ERFA_SRS = 1.97412574336e-8;


/* ERFA_DINT(A) - truncate to nearest whole number towards zero (double) */
function ERFA_DINT(X) { return X < 0 ? Math.ceil(X) : Math.floor(X); }

/* ERFA_DNINT(A) - round to nearest whole number (double) */
function ERFA_DNINT(X) { return X < 0 ? Math.ceil(X-.5) : Math.floor(X+.5); }

/* ERFA_DSIGN(A,B) - magnitude of A with sign of B (double) */
function ERFA_DSIGN(A,B) { return ((B)<0.0?-Math.abs(A):Math.abs(A)); }

/* Reference ellipsoids */
var ERFA_WGS84 = 1;
var ERFA_GRS80 = 2;
var ERFA_WGS72 = 3;

/* Dates and Delta(AT)s */
var ERFA_LEAPSEC = [
    [ 1960,  1,  1.4178180 ],
    [ 1961,  1,  1.4228180 ],
    [ 1961,  8,  1.3728180 ],
    [ 1962,  1,  1.8458580 ],
    [ 1963, 11,  1.9458580 ],
    [ 1964,  1,  3.2401300 ],
    [ 1964,  4,  3.3401300 ],
    [ 1964,  9,  3.4401300 ],
    [ 1965,  1,  3.5401300 ],
    [ 1965,  3,  3.6401300 ],
    [ 1965,  7,  3.7401300 ],
    [ 1965,  9,  3.8401300 ],
    [ 1966,  1,  4.3131700 ],
    [ 1968,  2,  4.2131700 ],
    [ 1972,  1, 10.0       ],
    [ 1972,  7, 11.0       ],
    [ 1973,  1, 12.0       ],
    [ 1974,  1, 13.0       ],
    [ 1975,  1, 14.0       ],
    [ 1976,  1, 15.0       ],
    [ 1977,  1, 16.0       ],
    [ 1978,  1, 17.0       ],
    [ 1979,  1, 18.0       ],
    [ 1980,  1, 19.0       ],
    [ 1981,  7, 20.0       ],
    [ 1982,  7, 21.0       ],
    [ 1983,  7, 22.0       ],
    [ 1985,  7, 23.0       ],
    [ 1988,  1, 24.0       ],
    [ 1990,  1, 25.0       ],
    [ 1991,  1, 26.0       ],
    [ 1992,  7, 27.0       ],
    [ 1993,  7, 28.0       ],
    [ 1994,  7, 29.0       ],
    [ 1996,  1, 30.0       ],
    [ 1997,  7, 31.0       ],
    [ 1999,  1, 32.0       ],
    [ 2006,  1, 33.0       ],
    [ 2009,  1, 34.0       ],
    [ 2012,  7, 35.0       ],
    [ 2015,  7, 36.0       ],
    [ 2017,  1, 37.0       ]
];

/* Release year for eraDat */
var ERFA_IYV = 2019;

function eraAb(pnat, v, s, bm1)
/*
**  - - - - - -
**   e r a A b
**  - - - - - -
**
**  Apply aberration to transform natural direction into proper
**  direction.
**
**  Given:
**    pnat    double[3]   natural direction to the source (unit vector)
**    v       double[3]   observer barycentric velocity in units of c
**    s       double      distance between the Sun and the observer (au)
**    bm1     double      sqrt(1-|v|^2): reciprocal of Lorenz factor
**
**  Returned:
**    ppr     double[3]   proper direction to source (unit vector)
**
**  Notes:
**
**  1) The algorithm is based on Expr. (7.40) in the Explanatory
**     Supplement (Urban & Seidelmann 2013), but with the following
**     changes:
**
**     o  Rigorous rather than approximate normalization is applied.
**
**     o  The gravitational potential term from Expr. (7) in
**        Klioner (2003) is added, taking into account only the Sun's
**        contribution.  This has a maximum effect of about
**        0.4 microarcsecond.
**
**  2) In almost all cases, the maximum accuracy will be limited by the
**     supplied velocity.  For example, if the ERFA eraEpv00 function is
**     used, errors of up to 5 microarcseconds could occur.
**
**  References:
**
**     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
**     the Astronomical Almanac, 3rd ed., University Science Books
**     (2013).
**
**     Klioner, Sergei A., "A practical relativistic model for micro-
**     arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).
**
**  Called:
**     eraPdp       scalar product of two p-vectors
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var ppr = [0, 0, 0];;


   var i;
   var pdv, w1, w2, r2, w, p = [], r;


   pdv = eraPdp(pnat, v);
   w1 = 1.0 + pdv/(1.0 + bm1);
   w2 = ERFA_SRS/s;
   r2 = 0.0;
   for (i = 0; i < 3; i++) {
      w = pnat[i]*bm1 + w1*v[i] + w2*(v[i] - pdv*pnat[i]);
      p[i] = w;
      r2 = r2 + w*w;
   }
   r = Math.sqrt(r2);
   for (i = 0; i < 3; i++) {
      ppr[i] = p[i]/r;
   }

/* Finished. */

return ppr;
};
function eraAnp(a)
/*
**  - - - - - - -
**   e r a A n p
**  - - - - - - -
**
**  Normalize angle into the range 0 <= a < 2pi.
**
**  Given:
**     a        double     angle (radians)
**
**  Returned (function value):
**              double     angle in range 0-2pi
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var w;


   w = ((a) % (ERFA_D2PI));
   if (w < 0) w += ERFA_D2PI;

   return w;

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraCal2jd(iy, im, id)
/*
**  - - - - - - - - - -
**   e r a C a l 2 j d
**  - - - - - - - - - -
**
**  Gregorian Calendar to Julian Date.
**
**  Given:
**     iy,im,id  int     year, month, day in Gregorian calendar (Note 1)
**
**  Returned:
**     djm0      double  MJD zero-point: always 2400000.5
**     djm       double  Modified Julian Date for 0 hrs
**
**  Returned (function value):
**               int     status:
**                           0 = OK
**                          -1 = bad year   (Note 3: JD not computed)
**                          -2 = bad month  (JD not computed)
**                          -3 = bad day    (JD computed)
**
**  Notes:
**
**  1) The algorithm used is valid from -4800 March 1, but this
**     implementation rejects dates before -4799 January 1.
**
**  2) The Julian Date is returned in two pieces, in the usual ERFA
**     manner, which is designed to preserve time resolution.  The
**     Julian Date is available as a single number by adding djm0 and
**     djm.
**
**  3) In early eras the conversion is from the "Proleptic Gregorian
**     Calendar";  no account is taken of the date(s) of adoption of
**     the Gregorian Calendar, nor is the AD/BC numbering convention
**     observed.
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 12.92 (p604).
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var djm0 = 0.0;;
   var djm = 0.0;;


   var j, ly, my;
   var iypmy;

/* Earliest year allowed (4800BC) */
   var IYMIN = -4799;

/* Month lengths in days */
   var                 mtab = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];


/* Preset status. */
   j = 0;

/* Validate year and month. */
   if (iy < IYMIN) return [ -1, djm0, djm ];
   if (im < 1 || im > 12) return [ -2, djm0, djm ];

/* If February in a leap year, 1, otherwise 0. */
   ly = ~~(((im == 2) && !(iy%4) && (iy%100 || !(iy%400))));

/* Validate day, taking into account leap years. */
   if ( (id < 1) || (id > (mtab[im-1] + ly))) j = -3;

/* Return result. */
   my = ~~((im - 14) / 12);
   iypmy = ~~(~~(iy + my));
   djm0 = ERFA_DJM0;
   djm = (~~((1461 * (iypmy + 4800)) / 4)
                 + ~~((367 * (im - 2 - 12 * my)) / 12)
                 - ~~((3 * ~~(((iypmy + 4900)) / 100)) / 4)
                 + ~~id - 2432076);

/* Return status. */
   return [ j, djm0, djm ];

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraCp(p)
/*
**  - - - - - -
**   e r a C p
**  - - - - - -
**
**  Copy a p-vector.
**
**  Given:
**     p        double[3]     p-vector to be copied
**
**  Returned:
**     c        double[3]     copy
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var c = [0, 0, 0];;


   c[0] = p[0];
   c[1] = p[1];
   c[2] = p[2];

   return c;

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraCr(r)
/*
**  - - - - - -
**   e r a C r
**  - - - - - -
**
**  Copy an r-matrix.
**
**  Given:
**     r        double[3][3]    r-matrix to be copied
**
**  Returned:
**     c        double[3][3]    copy
**
**  Called:
**     eraCp        copy p-vector
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var c = [ [0,0,0], [0,0,0], [0,0,0] ];;
   if (typeof c == 'undefined') {
      c = [ [0,0,0], [0,0,0], [0,0,0] ];
   }

   c[0] = eraCp(r[0]);
   c[1] = eraCp(r[1]);
   c[2] = eraCp(r[2]);

   return c;

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraDat(iy, im, id, fd)
/*
**  - - - - - - -
**   e r a D a t
**  - - - - - - -
**
**  For a given UTC date, calculate Delta(AT) = TAI-UTC.
**
**     :------------------------------------------:
**     :                                          :
**     :                 IMPORTANT                :
**     :                                          :
**     :  A new version of this function must be  :
**     :  produced whenever a new leap second is  :
**     :  announced.  There are four items to     :
**     :  change on each such occasion:           :
**     :                                          :
**     :  1) A new line must be added to the set  :
**     :     of statements that initialize the    :
**     :     array "ERFA_LEAPSEC" in "erfam.js".  :
**     :                                          :
**     :  2) The constant ERFA_IYV must be set    :
**     :     to the current year in "erfam.js".   :
**     :                                          :
**     :  3) The "Latest leap second" comment     :
**     :     below must be set to the new leap    :
**     :     second date.                         :
**     :                                          :
**     :  4) The "This revision" comment, later,  :
**     :     must be set to the current date.     :
**     :                                          :
**     :  Change (2) must also be carried out     :
**     :  whenever the function is re-issued,     :
**     :  even if no leap seconds have been       :
**     :  added.                                  :
**     :                                          :
**     :  Latest leap second:  2016 December 31   :
**     :                                          :
**     :__________________________________________:
**
**  Given:
**     iy     int      UTC:  year (Notes 1 and 2)
**     im     int            month (Note 2)
**     id     int            day (Notes 2 and 3)
**     fd     double         fraction of day (Note 4)
**
**  Returned:
**     deltat double   TAI minus UTC, seconds
**
**  Returned (function value):
**            int      status (Note 5):
**                       1 = dubious year (Note 1)
**                       0 = OK
**                      -1 = bad year
**                      -2 = bad month
**                      -3 = bad day (Note 3)
**                      -4 = bad fraction (Note 4)
**                      -5 = internal error (Note 5)
**
**  Notes:
**
**  1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
**     to call the function with an earlier date.  If this is attempted,
**     zero is returned together with a warning status.
**
**     Because leap seconds cannot, in principle, be predicted in
**     advance, a reliable check for dates beyond the valid range is
**     impossible.  To guard against gross errors, a year five or more
**     after the release year of the present function (see the constant
**     ERFA_IYV) is considered dubious.  In this case a warning status
**     is returned but the result is computed in the normal way.
**
**     For both too-early and too-late years, the warning status is +1.
**     This is distinct from the error status -1, which signifies a year
**     so early that JD could not be computed.
**
**  2) If the specified date is for a day which ends with a leap second,
**     the TAI-UTC value returned is for the period leading up to the
**     leap second.  If the date is for a day which begins as a leap
**     second ends, the TAI-UTC returned is for the period following the
**     leap second.
**
**  3) The day number must be in the normal calendar range, for example
**     1 through 30 for April.  The "almanac" convention of allowing
**     such dates as January 0 and December 32 is not supported in this
**     function, in order to avoid confusion near leap seconds.
**
**  4) The fraction of day is used only for dates before the
**     introduction of leap seconds, the first of which occurred at the
**     end of 1971.  It is tested for validity (0 to 1 is the valid
**     range) even if not used;  if invalid, zero is used and status -4
**     is returned.  For many applications, setting fd to zero is
**     acceptable;  the resulting error is always less than 3 ms (and
**     occurs only pre-1972).
**
**  5) The status value returned in the case where there are multiple
**     errors refers to the first error detected.  For example, if the
**     month and day are 13 and 32 respectively, status -2 (bad month)
**     will be returned.  The "internal error" status refers to a
**     case that is impossible but causes some compilers to issue a
**     warning.
**
**  6) In cases where a valid result is not available, zero is returned.
**
**  References:
**
**  1) For dates from 1961 January 1 onwards, the expressions from the
**     file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.
**
**  2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
**     the 1992 Explanatory Supplement.
**
**  Called:
**     eraCal2jd    Gregorian calendar to JD
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var deltat = 0.0;
   var _rv2;

/* Reference dates (MJD) and drift rates (s/day), pre leap seconds */
   var                 drift = [
      [ 37300.0, 0.0012960 ],
      [ 37300.0, 0.0012960 ],
      [ 37300.0, 0.0012960 ],
      [ 37665.0, 0.0011232 ],
      [ 37665.0, 0.0011232 ],
      [ 38761.0, 0.0012960 ],
      [ 38761.0, 0.0012960 ],
      [ 38761.0, 0.0012960 ],
      [ 38761.0, 0.0012960 ],
      [ 38761.0, 0.0012960 ],
      [ 38761.0, 0.0012960 ],
      [ 38761.0, 0.0012960 ],
      [ 39126.0, 0.0025920 ],
      [ 39126.0, 0.0025920 ]
   ];

/* Number of Delta(AT) expressions before leap seconds were introduced */
   var NERA1 = drift.length;

/* Miscellaneous local variables */
   var j, i, m;
   var da, djm0, djm;


/* Initialize the result to zero. */
   deltat = da = 0.0;

/* If invalid fraction of a day, set error status and give up. */
   if (fd < 0.0 || fd > 1.0) return [ -4, deltat ];

/* Convert the date into an MJD. */
   j = ~~((_rv2 = eraCal2jd(iy, im, id))[0]);
   djm0 = _rv2[1];
   djm = _rv2[2];

/* If invalid year, month, or day, give up. */
   if (j < 0) return [ j, deltat ];

/* If pre-UTC year, set warning status and give up. */
   if (iy < ERFA_LEAPSEC[0][0]) return [ 1, deltat ];

/* If suspiciously late year, set warning status but proceed. */
   if (iy > ERFA_IYV + 5) j = 1;

/* Combine year and month to form a date-ordered integer... */
   m = ~~(12*iy + im);

/* ...and use it to find the preceding table entry. */
   for (i = ERFA_LEAPSEC.length-1; i >=0; i--) {
      if (m >= (12 * ERFA_LEAPSEC[i][0] + ERFA_LEAPSEC[i][1])) break;
   }

/* Prevent underflow warnings. */
   if (i < 0) return [ -5, deltat ];

/* Get the Delta(AT). */
   da = ERFA_LEAPSEC[i][2];

/* If pre-1972, adjust for drift. */
   if (i < NERA1) da += (djm + fd - drift[i][0]) * drift[i][1];

/* Return the Delta(AT) value. */
   deltat = da;

/* Return the status. */
   return [ j, deltat ];

}
/*
 *+----------------------------------------------------------------------
 *
 *  IAU SOFA functions converted to JS
 *  http:://www.github.com/mgreter/sofa.js
 *  2016 by Marcel Greter
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Please read notice below, as all rights go to the Standards
 *  Of Fundamental Astronomy (SOFA) Review Board of the International
 *  Astronomical Union, as far as applicable. There is no guarantee
 *  that the conversion is bug free and I give no warranty of
 *  usability or correctness whatsoever.
 *
 *  The agreement below (3c/d) says that functions should
 *  be renamed. From the preface I guess this only applies
 *  if the function behavior was changed in any way. Since
 *  this is a one-to-one conversion, it shouldn't apply?
 *
 *+----------------------------------------------------------------------
 * SOFA-Issue: 2016-05-03
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2016
 *  Standards Of Fundamental Astronomy Review Board
 *  of the International Astronomical Union.
 *
 *  =====================
 *  SOFA Software License
 *  =====================
 *
 *  NOTICE TO USER:
 *
 *  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
 *  WHICH APPLY TO ITS USE.
 *
 *  1. The Software is owned by the IAU SOFA Review Board ("the Board").
 *
 *  2. Permission is granted to anyone to use the SOFA software for any
 *     purpose, including commercial applications, free of charge and
 *     without payment of royalties, subject to the conditions and
 *     restrictions listed below.
 *
 *  3. You (the user) may copy and adapt the SOFA software and its
 *     algorithms for your own purposes and you may copy and distribute
 *     a resulting "derived work" to others on a world-wide, royalty-free
 *     basis, provided that the derived work complies with the following
 *     requirements:
 *
 *     a) Your work shall be marked or carry a statement that it (i) uses
 *        routines and computations derived by you from software provided
 *        by SOFA under license to you; and (ii) does not contain
 *        software provided by SOFA or software that has been distributed
 *        by or endorsed by SOFA.
 *
 *     b) The source code of your derived work must contain descriptions
 *        of how the derived work is based upon and/or differs from the
 *        original SOFA software.
 *
 *     c) The name(s) of all routine(s) that you distribute shall differ
 *        from the SOFA names, even when the SOFA content has not been
 *        otherwise changed.
 *
 *     d) The routine-naming prefix "iau" shall not be used.
 *
 *     e) The origin of the SOFA components of your derived work must not
 *        be misrepresented;  you must not claim that you wrote the
 *        original software, nor file a patent application for SOFA
 *        software or algorithms embedded in the SOFA software.
 *
 *     f) These requirements must be reproduced intact in any source
 *        distribution and shall apply to anyone to whom you have granted
 *        a further right to modify the source code of your derived work.
 *
 *  4. In any published work or commercial products which includes
 *     results achieved by using the SOFA software, you shall acknowledge
 *     that the SOFA software was used in obtaining those results.
 *
 *  5. You shall not cause the SOFA software to be brought into
 *     disrepute, either by misuse, or use for inappropriate tasks, or by
 *     inappropriate modification.
 *
 *  6. The SOFA software is provided "as is" and the Board makes no
 *     warranty as to its use or performance.   The Board does not and
 *     cannot warrant the performance or results which the user may obtain
 *     by using the SOFA software.  The Board makes no warranties, express
 *     or implied, as to non-infringement of third party rights,
 *     merchantability, or fitness for any particular purpose.  In no
 *     event will the Board be liable to the user for any consequential,
 *     incidental, or special damages, including any lost profits or lost
 *     savings, even if a Board representative has been advised of such
 *     damages, or for any claim by any third party.
 *
 *  7. The provision of any version of the SOFA software under the terms
 *     and conditions specified herein does not imply that future
 *     versions will also be made available under the same terms and
 *     conditions.

 *  Correspondence concerning SOFA software should be addressed as
 *  follows:
 *
 *     Internet email: sofa@rl.ac.uk
 *     Postal address: IAU SOFA Center
 *                     Rutherford Appleton Laboratory
 *                     Chilton, Didcot, Oxon OX11 0QX
 *                     United Kingdom
 *
 *-----------------------------------------------------------------------
*/
;
function eraEform( n)
/*
**  - - - - - - - - -
**   e r a E f o r m
**  - - - - - - - - -
**
**  Earth reference ellipsoids.
**
**  Given:
**     n    int         ellipsoid identifier (Note 1)
**
**  Returned:
**     a    double      equatorial radius (meters, Note 2)
**     f    double      flattening (Note 2)
**
**  Returned (function value):
**          int         status:  0 = OK
**                              -1 = illegal identifier (Note 3)
**
**  Notes:
**
**  1) The identifier n is a number that specifies the choice of
**     reference ellipsoid.  The following are supported:
**
**        n    ellipsoid
**
**        1     ERFA_WGS84
**        2     ERFA_GRS80
**        3     ERFA_WGS72
**
**     The n value has no significance outside the ERFA software.  For
**     convenience, symbols ERFA_WGS84 etc. are defined in erfam.h.
**
**  2) The ellipsoid parameters are returned in the form of equatorial
**     radius in meters (a) and flattening (f).  The latter is a number
**     around 0.00335, i.e. around 1/298.
**
**  3) For the case where an unsupported n value is supplied, zero a and
**     f are returned, as well as error status.
**
**  References:
**
**     Department of Defense World Geodetic System 1984, National
**     Imagery and Mapping Agency Technical Report 8350.2, Third
**     Edition, p3-2.
**
**     Moritz, H., Bull. Geodesique 66-2, 187 (1992).
**
**     The Department of Defense World Geodetic System 1972, World
**     Geodetic System Committee, May 1974.
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     p220.
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a = 0.0;;
   var f = 0.0;;


/* Look up a and f for the specified reference ellipsoid. */
   switch ( n ) {

   case ERFA_WGS84:
      a = 6378137.0;
      f = 1.0 / 298.257223563;
      break;

   case ERFA_GRS80:
      a = 6378137.0;
      f = 1.0 / 298.257222101;
      break;

   case ERFA_WGS72:
      a = 6378135.0;
      f = 1.0 / 298.26;
      break;

   default:

   /* Invalid identifier. */
      a = 0.0;
      f = 0.0;
      return [ -1, a, f  ];

   }

/* OK status. */
   return [ 0, a, f  ];

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraEra00(dj1, dj2)
/*
**  - - - - - - - - -
**   e r a E r a 0 0
**  - - - - - - - - -
**
**  Earth rotation angle (IAU 2000 model).
**
**  Given:
**     dj1,dj2   double    UT1 as a 2-part Julian Date (see note)
**
**  Returned (function value):
**               double    Earth rotation angle (radians), range 0-2pi
**
**  Notes:
**
**  1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
**     convenient way between the arguments dj1 and dj2.  For example,
**     JD(UT1)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**             dj1            dj2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  The date & time method is
**     best matched to the algorithm used:  maximum precision is
**     delivered when the dj1 argument is for 0hrs UT1 on the day in
**     question and the dj2 argument lies in the range 0 to 1, or vice
**     versa.
**
**  2) The algorithm is adapted from Expression 22 of Capitaine et al.
**     2000.  The time argument has been expressed in days directly,
**     and, to retain precision, integer contributions have been
**     eliminated.  The same formulation is given in IERS Conventions
**     (2003), Chap. 5, Eq. 14.
**
**  Called:
**     eraAnp       normalize angle into range 0 to 2pi
**
**  References:
**
**     Capitaine N., Guinot B. and McCarthy D.D, 2000, Astron.
**     Astrophys., 355, 398-405.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var d1, d2, t, f, theta;


/* Days since fundamental epoch. */
   if (dj1 < dj2) {
      d1 = dj1;
      d2 = dj2;
   } else {
      d1 = dj2;
      d2 = dj1;
   }
   t = d1 + (d2- ERFA_DJ00);

/* Fractional part of T (days). */
   f = ((d1) % (1.0)) + ((d2) % (1.0));

/* Earth rotation angle at this UT1. */
   theta = eraAnp(ERFA_D2PI * (f + 0.7790572732640
                            + 0.00273781191135448 * t));

   return theta;

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraFw2m(gamb, phib, psi, eps)
/*
**  - - - - - - - -
**   e r a F w 2 m
**  - - - - - - - -
**
**  Form rotation matrix given the Fukushima-Williams angles.
**
**  Given:
**     gamb     double         F-W angle gamma_bar (radians)
**     phib     double         F-W angle phi_bar (radians)
**     psi      double         F-W angle psi (radians)
**     eps      double         F-W angle epsilon (radians)
**
**  Returned:
**     r        double[3][3]   rotation matrix
**
**  Notes:
**
**  1) Naming the following points:
**
**           e = J2000.0 ecliptic pole,
**           p = GCRS pole,
**           E = ecliptic pole of date,
**     and   P = CIP,
**
**     the four Fukushima-Williams angles are as follows:
**
**        gamb = gamma = epE
**        phib = phi = pE
**        psi = psi = pEP
**        eps = epsilon = EP
**
**  2) The matrix representing the combined effects of frame bias,
**     precession and nutation is:
**
**        NxPxB = R_1(-eps).R_3(-psi).R_1(phib).R_3(gamb)
**
**  3) Three different matrices can be constructed, depending on the
**     supplied angles:
**
**     o  To obtain the nutation x precession x frame bias matrix,
**        generate the four precession angles, generate the nutation
**        components and add them to the psi_bar and epsilon_A angles,
**        and call the present function.
**
**     o  To obtain the precession x frame bias matrix, generate the
**        four precession angles and call the present function.
**
**     o  To obtain the frame bias matrix, generate the four precession
**        angles for date J2000.0 and call the present function.
**
**     The nutation-only and precession-only matrices can if necessary
**     be obtained by combining these three appropriately.
**
**  Called:
**     eraIr        initialize r-matrix to identity
**     eraRz        rotate around Z-axis
**     eraRx        rotate around X-axis
**
**  Reference:
**
**     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var r = [ [0,0,0], [0,0,0], [0,0,0] ];;
   if (typeof r == 'undefined') {
      r = [ [0,0,0], [0,0,0], [0,0,0] ];
   }

/* Construct the matrix. */
   r = eraIr();
   r = eraRz(gamb, r);
   r = eraRx(phib, r);
   r = eraRz(-psi, r);
   r = eraRx(-eps, r);

   return r;

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraGd2gc( n, elong, phi, height)
/*
**  - - - - - - - - -
**   e r a G d 2 g c
**  - - - - - - - - -
**
**  Transform geodetic coordinates to geocentric using the specified
**  reference ellipsoid.
**
**  Given:
**     n       int        ellipsoid identifier (Note 1)
**     elong   double     longitude (radians, east +ve)
**     phi     double     latitude (geodetic, radians, Note 3)
**     height  double     height above ellipsoid (geodetic, Notes 2,3)
**
**  Returned:
**     xyz     double[3]  geocentric vector (Note 2)
**
**  Returned (function value):
**             int        status:  0 = OK
**                                -1 = illegal identifier (Note 3)
**                                -2 = illegal case (Note 3)
**
**  Notes:
**
**  1) The identifier n is a number that specifies the choice of
**     reference ellipsoid.  The following are supported:
**
**        n    ellipsoid
**
**        1     ERFA_WGS84
**        2     ERFA_GRS80
**        3     ERFA_WGS72
**
**     The n value has no significance outside the ERFA software.  For
**     convenience, symbols ERFA_WGS84 etc. are defined in erfam.h.
**
**  2) The height (height, given) and the geocentric vector (xyz,
**     returned) are in meters.
**
**  3) No validation is performed on the arguments elong, phi and
**     height.  An error status -1 means that the identifier n is
**     illegal.  An error status -2 protects against cases that would
**     lead to arithmetic exceptions.  In all error cases, xyz is set
**     to zeros.
**
**  4) The inverse transformation is performed in the function eraGc2gd.
**
**  Called:
**     eraEform     Earth reference ellipsoids
**     eraGd2gce    geodetic to geocentric transformation, general
**     eraZp        zero p-vector
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var xyz = [0, 0, 0];;
   var _rv1, _rv2;

   var j;
   var a, f;


/* Obtain reference ellipsoid parameters. */
   j = ~~((_rv1 = eraEform(n))[0]);
   a = _rv1[1];
   f = _rv1[2];

/* If OK, transform longitude, geodetic latitude, height to x,y,z. */
   if ( j == 0 ) {
      j = ~~((_rv2 = eraGd2gce(a, f, elong, phi, height))[0]);
   xyz = _rv2[1];
      if ( j != 0 ) j = -2;
   }

/* Deal with any errors. */
   if ( j != 0 ) xyz = eraZp(xyz);

/* Return the status. */
   return [ j, xyz  ];

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraGd2gce( a, f, elong, phi, height)
/*
**  - - - - - - - - - -
**   e r a G d 2 g c e
**  - - - - - - - - - -
**
**  Transform geodetic coordinates to geocentric for a reference
**  ellipsoid of specified form.
**
**  Given:
**     a       double     equatorial radius (Notes 1,4)
**     f       double     flattening (Notes 2,4)
**     elong   double     longitude (radians, east +ve)
**     phi     double     latitude (geodetic, radians, Note 4)
**     height  double     height above ellipsoid (geodetic, Notes 3,4)
**
**  Returned:
**     xyz     double[3]  geocentric vector (Note 3)
**
**  Returned (function value):
**             int        status:  0 = OK
**                                -1 = illegal case (Note 4)
**  Notes:
**
**  1) The equatorial radius, a, can be in any units, but meters is
**     the conventional choice.
**
**  2) The flattening, f, is (for the Earth) a value around 0.00335,
**     i.e. around 1/298.
**
**  3) The equatorial radius, a, and the height, height, must be
**     given in the same units, and determine the units of the
**     returned geocentric vector, xyz.
**
**  4) No validation is performed on individual arguments.  The error
**     status -1 protects against (unrealistic) cases that would lead
**     to arithmetic exceptions.  If an error occurs, xyz is unchanged.
**
**  5) The inverse transformation is performed in the function
**     eraGc2gde.
**
**  6) The transformation for a standard ellipsoid (such as ERFA_WGS84) can
**     more conveniently be performed by calling eraGd2gc,  which uses a
**     numerical code to identify the required a and f values.
**
**  References:
**
**     Green, R.M., Spherical Astronomy, Cambridge University Press,
**     (1985) Section 4.5, p96.
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 4.22, p202.
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var xyz = [0, 0, 0];;


   var sp, cp, w, d, ac, as, r;


/* Functions of geodetic latitude. */
   sp = Math.sin(phi);
   cp = Math.cos(phi);
   w = 1.0 - f;
   w = w * w;
   d = cp*cp + w*sp*sp;
   if ( d <= 0.0 ) return [ -1, xyz  ];
   ac = a / Math.sqrt(d);
   as = w * ac;

/* Geocentric vector. */
   r = (ac + height) * cp;
   xyz[0] = r * Math.cos(elong);
   xyz[1] = r * Math.sin(elong);
   xyz[2] = (as + height) * sp;

/* Success. */
   return [ 0, xyz  ];

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraGmst00(uta, utb, tta, ttb)
/*
**  - - - - - - - - - -
**   e r a G m s t 0 0
**  - - - - - - - - - -
**
**  Greenwich mean sidereal time (model consistent with IAU 2000
**  resolutions).
**
**  Given:
**     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
**     tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)
**
**  Returned (function value):
**                double    Greenwich mean sidereal time (radians)
**
**  Notes:
**
**  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
**     Julian Dates, apportioned in any convenient way between the
**     argument pairs.  For example, JD=2450123.7 could be expressed in
**     any of these ways, among others:
**
**            Part A         Part B
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable (in the case of UT;  the TT is not at all critical
**     in this respect).  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  For UT, the date & time
**     method is best matched to the algorithm that is used by the Earth
**     Rotation Angle function, called internally:  maximum precision is
**     delivered when the uta argument is for 0hrs UT1 on the day in
**     question and the utb argument lies in the range 0 to 1, or vice
**     versa.
**
**  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
**     and TT to predict the effects of precession.  If UT1 is used for
**     both purposes, errors of order 100 microarcseconds result.
**
**  3) This GMST is compatible with the IAU 2000 resolutions and must be
**     used only in conjunction with other IAU 2000 compatible
**     components such as precession-nutation and equation of the
**     equinoxes.
**
**  4) The result is returned in the range 0 to 2pi.
**
**  5) The algorithm is from Capitaine et al. (2003) and IERS
**     Conventions 2003.
**
**  Called:
**     eraEra00     Earth rotation angle, IAU 2000
**     eraAnp       normalize angle into range 0 to 2pi
**
**  References:
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astronomy &
**     Astrophysics, 406, 1135-1149 (2003)
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var t, gmst;


/* TT Julian centuries since J2000.0. */
   t = ((tta - ERFA_DJ00) + ttb) / ERFA_DJC;

/* Greenwich Mean Sidereal Time, IAU 2000. */
   gmst = eraAnp(eraEra00(uta, utb) +
                   (     0.014506   +
                   (  4612.15739966 +
                   (     1.39667721 +
                   (    -0.00009344 +
                   (     0.00001882 )
          * t) * t) * t) * t) * ERFA_DAS2R);

   return gmst;

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraIr()
/*
**  - - - - - -
**   e r a I r
**  - - - - - -
**
**  Initialize an r-matrix to the identity matrix.
**
**  Returned:
**     r       double[3][3]    r-matrix
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var r = [ [0,0,0], [0,0,0], [0,0,0] ];;
   if (typeof r == 'undefined') {
      r = [ [0,0,0], [0,0,0], [0,0,0] ];
   }

   r[0][0] = 1.0;
   r[0][1] = 0.0;
   r[0][2] = 0.0;
   r[1][0] = 0.0;
   r[1][1] = 1.0;
   r[1][2] = 0.0;
   r[2][0] = 0.0;
   r[2][1] = 0.0;
   r[2][2] = 1.0;

   return r;

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraJd2cal(dj1, dj2)
/*
**  - - - - - - - - - -
**   e r a J d 2 c a l
**  - - - - - - - - - -
**
**  Julian Date to Gregorian year, month, day, and fraction of a day.
**
**  Given:
**     dj1,dj2   double   Julian Date (Notes 1, 2)
**
**  Returned (arguments):
**     iy        int      year
**     im        int      month
**     id        int      day
**     fd        double   fraction of day
**
**  Returned (function value):
**               int      status:
**                           0 = OK
**                          -1 = unacceptable date (Note 1)
**
**  Notes:
**
**  1) The earliest valid date is -68569.5 (-4900 March 1).  The
**     largest value accepted is 1e9.
**
**  2) The Julian Date is apportioned in any convenient way between
**     the arguments dj1 and dj2.  For example, JD=2450123.7 could
**     be expressed in any of these ways, among others:
**
**            dj1             dj2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**  3) In early eras the conversion is from the "proleptic Gregorian
**     calendar";  no account is taken of the date(s) of adoption of
**     the Gregorian calendar, nor is the AD/BC numbering convention
**     observed.
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 12.92 (p604).
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var iy = 0;;
   var im = 0;;
   var id = 0;;
   var fd = 0.0;;


/* Minimum and maximum allowed JD */
   var DJMIN = -68569.5;
   var DJMAX = 1e9;

   var jd, l, n, i, k;
   var dj, d1, d2, f1, f2, f, d;


/* Verify date is acceptable. */
   dj = dj1 + dj2;
   if (dj < DJMIN || dj > DJMAX) return [ -1, iy, im, id, fd ];

/* Copy the date, big then small, and re-align to midnight. */
   if (Math.abs(dj1) >= Math.abs(dj2)) {
      d1 = dj1;
      d2 = dj2;
   } else {
      d1 = dj2;
      d2 = dj1;
   }
   d2 -= 0.5;

/* Separate day and fraction. */
   f1 = ((d1) % (1.0));
   f2 = ((d2) % (1.0));
   f = ((f1 + f2) % (1.0));
   if (f < 0.0) f += 1.0;
   d = ERFA_DNINT(d1-f1) + ERFA_DNINT(d2-f2) + ERFA_DNINT(f1+f2-f);
   jd = ~~(~~ERFA_DNINT(d) + 1);

/* Express day in Gregorian calendar. */
   l = ~~(jd + 68569);
   n = ~~(~~((4 * l) / 146097));
   l -= ~~(~~((146097 * n + 3) / 4));
   i = ~~(~~((4000 * (l + 1)) / 1461001));
   l -= ~~(~~((1461 * i) / 4) - 31);
   k = ~~(~~((80 * l) / 2447));
   id = ~~(~~(l - ~~((2447 * k) / 80)));
   l = ~~(k / 11);
   im = ~~(~~(k + 2 - 12 * l));
   iy = ~~(~~(100 * (n - 49) + i + l));
   fd = f;

   return [ 0, iy, im, id, fd ];

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraNum00a(date1, date2)
/*
**  - - - - - - - - - -
**   e r a N u m 0 0 a
**  - - - - - - - - - -
**
**  Form the matrix of nutation for a given date, IAU 2000A model.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rmatn        double[3][3]    nutation matrix
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix operates in the sense V(true) = rmatn * V(mean), where
**     the p-vector V(true) is with respect to the true equatorial triad
**     of date and the p-vector V(mean) is with respect to the mean
**     equatorial triad of date.
**
**  3) A faster, but slightly less accurate result (about 1 mas), can be
**     obtained by using instead the eraNum00b function.
**
**  Called:
**     eraPn00a     bias/precession/nutation, IAU 2000A
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 3.222-3 (p114).
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var rmatn = [ [0,0,0], [0,0,0], [0,0,0] ];;
   if (typeof rmatn == 'undefined') {
      rmatn = [ [0,0,0], [0,0,0], [0,0,0] ];
   }   var _rv1;

   var dpsi, deps, epsa, rb = [[], [], []], rp = [[], [], []], rbp = [[], [], []], rbpn = [[], [], []];


/* Obtain the required matrix (discarding other results). */
   (_rv1 = eraPn00a(date1, date2))[0];
   dpsi = _rv1[0];
   deps = _rv1[1];
   epsa = _rv1[2];
   rb = _rv1[3];
   rp = _rv1[4];
   rbp = _rv1[5];
   rmatn = _rv1[6];
   rbpn = _rv1[7];

   return rmatn;

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraObl06(date1, date2)
/*
**  - - - - - - - - -
**   e r a O b l 0 6
**  - - - - - - - - -
**
**  Mean obliquity of the ecliptic, IAU 2006 precession model.
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double   obliquity of the ecliptic (radians, Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The result is the angle between the ecliptic and mean equator of
**     date date1+date2.
**
**  Reference:
**
**     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var t, eps0;


/* Interval between fundamental date J2000.0 and given date (JC). */
   t = ((date1 - ERFA_DJ00) + date2) / ERFA_DJC;

/* Mean obliquity. */
   eps0 = (84381.406     +
          (-46.836769    +
          ( -0.0001831   +
          (  0.00200340  +
          ( -0.000000576 +
          ( -0.0000000434) * t) * t) * t) * t) * t) * ERFA_DAS2R;

   return eps0;

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraPdp(a, b)
/*
**  - - - - - - -
**   e r a P d p
**  - - - - - - -
**
**  p-vector inner (=scalar=dot) product.
**
**  Given:
**     a      double[3]     first p-vector
**     b      double[3]     second p-vector
**
**  Returned (function value):
**            double        a . b
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var w;


   w  = a[0] * b[0]
      + a[1] * b[1]
      + a[2] * b[2];

   return w;

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraPfw06(date1, date2)
/*
**  - - - - - - - - -
**   e r a P f w 0 6
**  - - - - - - - - -
**
**  Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     gamb         double   F-W angle gamma_bar (radians)
**     phib         double   F-W angle phi_bar (radians)
**     psib         double   F-W angle psi_bar (radians)
**     epsa         double   F-W angle epsilon_A (radians)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) Naming the following points:
**
**           e = J2000.0 ecliptic pole,
**           p = GCRS pole,
**           E = mean ecliptic pole of date,
**     and   P = mean pole of date,
**
**     the four Fukushima-Williams angles are as follows:
**
**        gamb = gamma_bar = epE
**        phib = phi_bar = pE
**        psib = psi_bar = pEP
**        epsa = epsilon_A = EP
**
**  3) The matrix representing the combined effects of frame bias and
**     precession is:
**
**        PxB = R_1(-epsa).R_3(-psib).R_1(phib).R_3(gamb)
**
**  4) The matrix representing the combined effects of frame bias,
**     precession and nutation is simply:
**
**        NxPxB = R_1(-epsa-dE).R_3(-psib-dP).R_1(phib).R_3(gamb)
**
**     where dP and dE are the nutation components with respect to the
**     ecliptic of date.
**
**  Reference:
**
**     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
**
**  Called:
**     eraObl06     mean obliquity, IAU 2006
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var gamb = 0.0;;
   var phib = 0.0;;
   var psib = 0.0;;
   var epsa = 0.0;;


   var t;


/* Interval between fundamental date J2000.0 and given date (JC). */
   t = ((date1 - ERFA_DJ00) + date2) / ERFA_DJC;

/* P03 bias+precession angles. */
   gamb = (    -0.052928     +
           (    10.556378     +
           (     0.4932044    +
           (    -0.00031238   +
           (    -0.000002788  +
           (     0.0000000260 )
           * t) * t) * t) * t) * t) * ERFA_DAS2R;
   phib = ( 84381.412819     +
           (   -46.811016     +
           (     0.0511268    +
           (     0.00053289   +
           (    -0.000000440  +
           (    -0.0000000176 )
           * t) * t) * t) * t) * t) * ERFA_DAS2R;
   psib = (    -0.041775     +
           (  5038.481484     +
           (     1.5584175    +
           (    -0.00018522   +
           (    -0.000026452  +
           (    -0.0000000148 )
           * t) * t) * t) * t) * t) * ERFA_DAS2R;
   epsa = eraObl06(date1, date2);

   return [gamb, phib, psib, epsa];

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraPnm06a(date1, date2)
/*
**  - - - - - - - - - -
**   e r a P n m 0 6 a
**  - - - - - - - - - -
**
**  Form the matrix of precession-nutation for a given date (including
**  frame bias), IAU 2006 precession and IAU 2000A nutation models.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rnpb        double[3][3] bias-precession-nutation matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix operates in the sense V(date) = rnpb * V(GCRS), where
**     the p-vector V(date) is with respect to the true equatorial triad
**     of date date1+date2 and the p-vector V(GCRS) is with respect to
**     the Geocentric Celestial Reference System (IAU, 2000).
**
**  Called:
**     eraPfw06     bias-precession F-W angles, IAU 2006
**     eraNut06a    nutation, IAU 2006/2000A
**     eraFw2m      F-W angles to r-matrix
**
**  Reference:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855.
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var rnpb = [ [0,0,0], [0,0,0], [0,0,0] ];;
   if (typeof rnpb == 'undefined') {
      rnpb = [ [0,0,0], [0,0,0], [0,0,0] ];
   }   var _rv1, _rv2;

   var gamb, phib, psib, epsa, dp, de;


/* Fukushima-Williams angles for frame bias and precession. */
   (_rv1 = eraPfw06(date1, date2))[0];
   gamb = _rv1[0];
   phib = _rv1[1];
   psib = _rv1[2];
   epsa = _rv1[3];

/* Nutation components. */
   (_rv2 = eraNut06a(date1, date2))[0];
   dp = _rv2[0];
   de = _rv2[1];

/* Equinox based nutation x precession x bias matrix. */
   rnpb = eraFw2m(gamb, phib, psib + dp, epsa + de);

   return rnpb;

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraPom00(xp, yp, sp)
/*
**  - - - - - - - - - -
**   e r a P o m 0 0
**  - - - - - - - - - -
**
**  Form the matrix of polar motion for a given date, IAU 2000.
**
**  Given:
**     xp,yp    double    coordinates of the pole (radians, Note 1)
**     sp       double    the TIO locator s' (radians, Note 2)
**
**  Returned:
**     rpom     double[3][3]   polar-motion matrix (Note 3)
**
**  Notes:
**
**  1) The arguments xp and yp are the coordinates (in radians) of the
**     Celestial Intermediate Pole with respect to the International
**     Terrestrial Reference System (see IERS Conventions 2003),
**     measured along the meridians to 0 and 90 deg west respectively.
**
**  2) The argument sp is the TIO locator s', in radians, which
**     positions the Terrestrial Intermediate Origin on the equator.  It
**     is obtained from polar motion observations by numerical
**     integration, and so is in essence unpredictable.  However, it is
**     dominated by a secular drift of about 47 microarcseconds per
**     century, and so can be taken into account by using s' = -47*t,
**     where t is centuries since J2000.0.  The function eraSp00
**     implements this approximation.
**
**  3) The matrix operates in the sense V(TRS) = rpom * V(CIP), meaning
**     that it is the final rotation when computing the pointing
**     direction to a celestial source.
**
**  Called:
**     eraIr        initialize r-matrix to identity
**     eraRz        rotate around Z-axis
**     eraRy        rotate around Y-axis
**     eraRx        rotate around X-axis
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var rpom = [ [0,0,0], [0,0,0], [0,0,0] ];;
   if (typeof rpom == 'undefined') {
      rpom = [ [0,0,0], [0,0,0], [0,0,0] ];
   }

/* Construct the matrix. */
   rpom = eraIr();
   rpom = eraRz(sp, rpom);
   rpom = eraRy(-xp, rpom);
   rpom = eraRx(-yp, rpom);

   return rpom;

}
/*
 *+----------------------------------------------------------------------
 *
 *  ERFA/SOFA functions converted to JS
 *  Copyright (C) 2020 by Marcel Greter
 *  http:://www.github.com/mgreter/sofa.js
 *
 *  The conversion is done by a custom hacked perl script.
 *  Automatically generates QUnit tests for all functions.
 *
 *  Conversion is made from liberfa sources:
 *  https://github.com/liberfa/erfa
 *
 *+----------------------------------------------------------------------
 *  THIS WORK IS RELEASED UNDER THE SAME TERMS AS ERFA:
 *+----------------------------------------------------------------------
 *
 *  Copyright (C) 2013-2014, NumFOCUS Foundation.
 *  All rights reserved.
 *  
 *  This library is derived, with permission, from the International
 *  Astronomical Union's "Standards of Fundamental Astronomy" library,
 *  available from http://www.iausofa.org.
 *  
 *  The ERFA version is intended to retain identical
 *  functionality to the SOFA library, but made distinct through
 *  different function and file names, as set out in the SOFA license
 *  conditions. The SOFA original has a role as a reference standard
 *  for the IAU and IERS, and consequently redistribution is permitted only
 *  in its unaltered state. The ERFA version is not subject to this
 *  restriction and therefore can be included in distributions which do not
 *  support the concept of "read only" software.
 *  
 *  Although the intent is to replicate the SOFA API (other than replacement of
 *  prefix names) and results (with the exception of bugs; any that are
 *  discovered will be fixed), SOFA is not responsible for any errors found
 *  in this version of the library.
 *  
 *  If you wish to acknowledge the SOFA heritage, please acknowledge that
 *  you are using a library derived from SOFA, rather than SOFA itself.
 *  
 *  
 *  TERMS AND CONDITIONS
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  
 *  1 Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *  
 *  2 Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *  
 *  3 Neither the name of the Standards Of Fundamental Astronomy Board, the
 *     International Astronomical Union nor the names of its contributors
 *     may be used to endorse or promote products derived from this software
 *     without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 *  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 *  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 *  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *-----------------------------------------------------------------------
*/
;
function eraPvtob(elong, phi, hm, xp, yp, sp, theta)
/*
**  - - - - - - - - -
**   e r a P v t o b
**  - - - - - - - - -
**
**  Position and velocity of a terrestrial observing station.
**
**  Given:
**     elong   double       longitude (radians, east +ve, Note 1)
**     phi     double       latitude (geodetic, radians, Note 1)
**     hm      double       height above ref. ellipsoid (geodetic, m)
**     xp,yp   double       coordinates of the pole (radians, Note 2)
**     sp      double       the TIO locator s' (radians, Note 2)
**     theta   double       Earth rotation angle (radians, Note 3)
**
**  Returned:
**     pv      double[2][3] position/velocity vector (m, m/s, CIRS)
**
**  Notes:
**
**  1) The terrestrial coordinates are with respect to the ERFA_WGS84
**     reference ellipsoid.
**
**  2) xp and yp are the coordinates (in radians) of the Celestial
**     Intermediate Pole with respect to the International Terrestrial
**     Reference System (see IERS Conventions), measured along the
**     meridians 0 and 90 deg west respectively.  sp is the TIO locator
**     s', in radians, which positions the Terrestrial Intermediate
**     Origin on the equator.  For many applications, xp, yp and
**     (especially) sp can be set to zero.
**
**  3) If theta is Greenwich apparent sidereal time instead of Earth
**     rotation angle, the result is with respect to the true equator
**     and equinox of date, i.e. with the x-axis at the equinox rather
**     than the celestial intermediate origin.
**
**  4) The velocity units are meters per UT1 second, not per SI second.
**     This is unlikely to have any practical consequences in the modern
**     era.
**
**  5) No validation is performed on the arguments.  Error cases that
**     could lead to arithmetic exceptions are trapped by the eraGd2gc
**     function, and the result set to zeros.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
**     the Astronomical Almanac, 3rd ed., University Science Books
**     (2013), Section 7.4.3.3.
**
**  Called:
**     eraGd2gc     geodetic to geocentric transformation
**     eraPom00     polar motion matrix
**     eraTrxp      product of transpose of r-matrix and p-vector
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var pv = [ [0,0,0], [0,0,0] ];;
   var _rv1;

/* Earth rotation rate in radians per UT1 second */
   var OM = 1.00273781191135448 * ERFA_D2PI / ERFA_DAYSEC;

   var xyzm = [], rpm = [[], [], []], xyz = [], x, y, z, s, c;


/* Geodetic to geocentric transformation (ERFA_WGS84). */
   (_rv1 = eraGd2gc(1, elong, phi, hm))[0];
   xyzm = _rv1[1];

/* Polar motion and TIO position. */
   rpm = eraPom00(xp, yp, sp);
   xyz = eraTrxp(rpm, xyzm);
   x = xyz[0];
   y = xyz[1];
   z = xyz[2];

/* Functions of ERA. */
   s = Math.sin(theta);
   c = Math.cos(theta);

/* Position. */
   pv[0][0] = c*x - s*y;
   pv[0][1] = s*x + c*y;
   pv[0][2] = z;

/* Velocity. */
   pv[1][0] = OM * ( -s*x - c*y );
   pv[1][1] = OM * (  c*x - s*y );
   pv[1][2] = 0.0;

/* Finished. */

return pv;
};
function eraRx(phi, r)
/*
**  - - - - - -
**   e r a R x
**  - - - - - -
**
**  Rotate an r-matrix about the x-axis.
**
**  Given:
**     phi    double          angle (radians)
**
**  Given and returned:
**     r      double[3][3]    r-matrix, rotated
**
**  Notes:
**
**  1) Calling this function with positive phi incorporates in the
**     supplied r-matrix r an additional rotation, about the x-axis,
**     anticlockwise as seen looking towards the origin from positive x.
**
**  2) The additional rotation can be represented by this matrix:
**
**         (  1        0            0      )
**         (                               )
**         (  0   + cos(phi)   + sin(phi)  )
**         (                               )
**         (  0   - sin(phi)   + cos(phi)  )
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   if (typeof r == 'undefined') {
      r = [ [0,0,0], [0,0,0], [0,0,0] ];
   }

   var s, c, a10, a11, a12, a20, a21, a22;


   s = Math.sin(phi);
   c = Math.cos(phi);

   a10 =   c*r[1][0] + s*r[2][0];
   a11 =   c*r[1][1] + s*r[2][1];
   a12 =   c*r[1][2] + s*r[2][2];
   a20 = - s*r[1][0] + c*r[2][0];
   a21 = - s*r[1][1] + c*r[2][1];
   a22 = - s*r[1][2] + c*r[2][2];

   r[1][0] = a10;
   r[1][1] = a11;
   r[1][2] = a12;
   r[2][0] = a20;
   r[2][1] = a21;
   r[2][2] = a22;

   return r;

};
function eraRy(theta, r)
/*
**  - - - - - -
**   e r a R y
**  - - - - - -
**
**  Rotate an r-matrix about the y-axis.
**
**  Given:
**     theta  double          angle (radians)
**
**  Given and returned:
**     r      double[3][3]    r-matrix, rotated
**
**  Notes:
**
**  1) Calling this function with positive theta incorporates in the
**     supplied r-matrix r an additional rotation, about the y-axis,
**     anticlockwise as seen looking towards the origin from positive y.
**
**  2) The additional rotation can be represented by this matrix:
**
**         (  + cos(theta)     0      - sin(theta)  )
**         (                                        )
**         (       0           1           0        )
**         (                                        )
**         (  + sin(theta)     0      + cos(theta)  )
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   if (typeof r == 'undefined') {
      r = [ [0,0,0], [0,0,0], [0,0,0] ];
   }

   var s, c, a00, a01, a02, a20, a21, a22;


   s = Math.sin(theta);
   c = Math.cos(theta);

   a00 = c*r[0][0] - s*r[2][0];
   a01 = c*r[0][1] - s*r[2][1];
   a02 = c*r[0][2] - s*r[2][2];
   a20 = s*r[0][0] + c*r[2][0];
   a21 = s*r[0][1] + c*r[2][1];
   a22 = s*r[0][2] + c*r[2][2];

   r[0][0] = a00;
   r[0][1] = a01;
   r[0][2] = a02;
   r[2][0] = a20;
   r[2][1] = a21;
   r[2][2] = a22;

   return r;

};
function eraRz(psi, r)
/*
**  - - - - - -
**   e r a R z
**  - - - - - -
**
**  Rotate an r-matrix about the z-axis.
**
**  Given:
**     psi    double          angle (radians)
**
**  Given and returned:
**     r      double[3][3]    r-matrix, rotated
**
**  Notes:
**
**  1) Calling this function with positive psi incorporates in the
**     supplied r-matrix r an additional rotation, about the z-axis,
**     anticlockwise as seen looking towards the origin from positive z.
**
**  2) The additional rotation can be represented by this matrix:
**
**         (  + cos(psi)   + sin(psi)     0  )
**         (                                 )
**         (  - sin(psi)   + cos(psi)     0  )
**         (                                 )
**         (       0            0         1  )
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   if (typeof r == 'undefined') {
      r = [ [0,0,0], [0,0,0], [0,0,0] ];
   }

   var s, c, a00, a01, a02, a10, a11, a12;


   s = Math.sin(psi);
   c = Math.cos(psi);

   a00 =   c*r[0][0] + s*r[1][0];
   a01 =   c*r[0][1] + s*r[1][1];
   a02 =   c*r[0][2] + s*r[1][2];
   a10 = - s*r[0][0] + c*r[1][0];
   a11 = - s*r[0][1] + c*r[1][1];
   a12 = - s*r[0][2] + c*r[1][2];

   r[0][0] = a00;
   r[0][1] = a01;
   r[0][2] = a02;
   r[1][0] = a10;
   r[1][1] = a11;
   r[1][2] = a12;

   return r;

};
function eraTaitt(tai1, tai2)
/*
**  - - - - - - - - -
**   e r a T a i t t
**  - - - - - - - - -
**
**  Time scale transformation:  International Atomic Time, TAI, to
**  Terrestrial Time, TT.
**
**  Given:
**     tai1,tai2  double    TAI as a 2-part Julian Date
**
**  Returned:
**     tt1,tt2    double    TT as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Note:
**
**     tai1+tai2 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where tai1 is the Julian
**     Day Number and tai2 is the fraction of a day.  The returned
**     tt1,tt2 follow suit.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var tt1 = 0.0;;
   var tt2 = 0.0;;


/* TT minus TAI (days). */
   var dtat = ERFA_TTMTAI/ERFA_DAYSEC;


/* Result, safeguarding precision. */
   if ( Math.abs(tai1) > Math.abs(tai2) ) {
      tt1 = tai1;
      tt2 = tai2 + dtat;
   } else {
      tt1 = tai1 + dtat;
      tt2 = tai2;
   }

/* Status (always OK). */
   return [ 0, tt1, tt2 ];

};
function eraTaiut1(tai1, tai2, dta)
/*
**  - - - - - - - - - -
**   e r a T a i u t 1
**  - - - - - - - - - -
**
**  Time scale transformation:  International Atomic Time, TAI, to
**  Universal Time, UT1.
**
**  Given:
**     tai1,tai2  double    TAI as a 2-part Julian Date
**     dta        double    UT1-TAI in seconds
**
**  Returned:
**     ut11,ut12  double    UT1 as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Notes:
**
**  1) tai1+tai2 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where tai1 is the Julian
**     Day Number and tai2 is the fraction of a day.  The returned
**     UT11,UT12 follow suit.
**
**  2) The argument dta, i.e. UT1-TAI, is an observed quantity, and is
**     available from IERS tabulations.
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var ut11 = 0.0;;
   var ut12 = 0.0;;


   var dtad;


/* Result, safeguarding precision. */
   dtad = dta / ERFA_DAYSEC;
   if ( Math.abs(tai1) > Math.abs(tai2) ) {
      ut11 = tai1;
      ut12 = tai2 + dtad;
   } else {
      ut11 = tai1 + dtad;
      ut12 = tai2;
   }

/* Status (always OK). */
   return [ 0, ut11, ut12 ];

};
function eraTr(r)
/*
**  - - - - - -
**   e r a T r
**  - - - - - -
**
**  Transpose an r-matrix.
**
**  Given:
**     r        double[3][3]    r-matrix
**
**  Returned:
**     rt       double[3][3]    transpose
**
**  Note:
**     It is permissible for r and rt to be the same array.
**
**  Called:
**     eraCr        copy r-matrix
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var rt = [ [0,0,0], [0,0,0], [0,0,0] ];;
   if (typeof rt == 'undefined') {
      rt = [ [0,0,0], [0,0,0], [0,0,0] ];
   }

   var wm = [[], [], []];
   var i, j;


   for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
         wm[i][j] = r[j][i];
      }
   }
   rt = eraCr(wm);

   return rt;

};
function eraTrxp(r, p)
/*
**  - - - - - - - -
**   e r a T r x p
**  - - - - - - - -
**
**  Multiply a p-vector by the transpose of an r-matrix.
**
**  Given:
**     r        double[3][3]   r-matrix
**     p        double[3]      p-vector
**
**  Returned:
**     trp      double[3]      r * p
**
**  Note:
**     It is permissible for p and trp to be the same array.
**
**  Called:
**     eraTr        transpose r-matrix
**     eraRxp       product of r-matrix and p-vector
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var trp = [0, 0, 0];;


   var tr = [[], [], []];


/* Transpose of matrix r. */
   tr = eraTr(r);

/* Matrix tr * vector p -> vector trp. */
   trp = eraRxp(tr, p);

   return trp;

}
;
function eraUtctai(utc1, utc2)
/*
**  - - - - - - - - - -
**   e r a U t c t a i
**  - - - - - - - - - -
**
**  Time scale transformation:  Coordinated Universal Time, UTC, to
**  International Atomic Time, TAI.
**
**  Given:
**     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-4)
**
**  Returned:
**     tai1,tai2  double   TAI as a 2-part Julian Date (Note 5)
**
**  Returned (function value):
**                int      status: +1 = dubious year (Note 3)
**                                  0 = OK
**                                 -1 = unacceptable date
**
**  Notes:
**
**  1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
**     convenient way between the two arguments, for example where utc1
**     is the Julian Day Number and utc2 is the fraction of a day.
**
**  2) JD cannot unambiguously represent UTC during a leap second unless
**     special measures are taken.  The convention in the present
**     function is that the JD day represents UTC days whether the
**     length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era
**     there were smaller jumps (in either direction) each time the
**     linear UTC(TAI) expression was changed, and these "mini-leaps"
**     are also included in the ERFA convention.
**
**  3) The warning status "dubious year" flags UTCs that predate the
**     introduction of the time scale or that are too far in the future
**     to be trusted.  See eraDat for further details.
**
**  4) The function eraDtf2d converts from calendar date and time of day
**     into 2-part Julian Date, and in the case of UTC implements the
**     leap-second-ambiguity convention described above.
**
**  5) The returned TAI1,TAI2 are such that their sum is the TAI Julian
**     Date.
**
**  Called:
**     eraJd2cal    JD to Gregorian calendar
**     eraDat       delta(AT) = TAI-UTC
**     eraCal2jd    Gregorian calendar to JD
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var tai1 = 0.0;;
   var tai2 = 0.0;;
   var _rv1, _rv2, _rv3, _rv4, _rv5, _rv6;

   var big1;
   var iy, im, id, j, iyt, imt, idt;
   var u1, u2, fd, dat0, dat12, w, dat24, dlod, dleap, z1, z2, a2;


/* Put the two parts of the UTC into big-first order. */
   big1 = ~~(( Math.abs(utc1) >= Math.abs(utc2) ));
   if ( big1 ) {
      u1 = utc1;
      u2 = utc2;
   } else {
      u1 = utc2;
      u2 = utc1;
   }

/* Get TAI-UTC at 0h today. */
   j = ~~((_rv1 = eraJd2cal(u1, u2))[0]);
   iy = ~~(_rv1[1]);
   im = ~~(_rv1[2]);
   id = ~~(_rv1[3]);
   fd = _rv1[4];
   if ( j ) return [ j, tai1, tai2 ];
   j = ~~((_rv2 = eraDat(iy, im, id, 0.0))[0]);
   dat0 = _rv2[1];
   if ( j < 0 ) return [ j, tai1, tai2 ];

/* Get TAI-UTC at 12h today (to detect drift). */
   j = ~~((_rv3 = eraDat(iy, im, id, 0.5))[0]);
   dat12 = _rv3[1];
   if ( j < 0 ) return [ j, tai1, tai2 ];

/* Get TAI-UTC at 0h tomorrow (to detect jumps). */
   j = ~~((_rv4 = eraJd2cal(u1+1.5, u2-fd))[0]);
   iyt = ~~(_rv4[1]);
   imt = ~~(_rv4[2]);
   idt = ~~(_rv4[3]);
   w = _rv4[4];
   if ( j ) return [ j, tai1, tai2 ];
   j = ~~((_rv5 = eraDat(iyt, imt, idt, 0.0))[0]);
   dat24 = _rv5[1];
   if ( j < 0 ) return [ j, tai1, tai2 ];

/* Separate TAI-UTC change into per-day (DLOD) and any jump (DLEAP). */
   dlod = 2.0 * (dat12 - dat0);
   dleap = dat24 - (dat0 + dlod);

/* Remove any scaling applied to spread leap into preceding day. */
   fd *= (ERFA_DAYSEC+dleap)/ERFA_DAYSEC;

/* Scale from (pre-1972) UTC seconds to SI seconds. */
   fd *= (ERFA_DAYSEC+dlod)/ERFA_DAYSEC;

/* Today's calendar date to 2-part JD. */
   var __rv__ = (_rv6 = eraCal2jd(iy, im, id))[0];
   z1 = _rv6[1];
   z2 = _rv6[2];
   if(__rv__) return [ -1, tai1, tai2 ];

/* Assemble the TAI result, preserving the UTC split and order. */
   a2 = z1 - u1;
   a2 += z2;
   a2 += fd + dat0/ERFA_DAYSEC;
   if ( big1 ) {
      tai1 = u1;
      tai2 = a2;
   } else {
      tai1 = a2;
      tai2 = u1;
   }

/* Status. */
   return [ j, tai1, tai2 ];

}
;
function eraUtcut1(utc1, utc2, dut1)
/*
**  - - - - - - - - - -
**   e r a U t c u t 1
**  - - - - - - - - - -
**
**  Time scale transformation:  Coordinated Universal Time, UTC, to
**  Universal Time, UT1.
**
**  Given:
**     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-4)
**     dut1       double   Delta UT1 = UT1-UTC in seconds (Note 5)
**
**  Returned:
**     ut11,ut12  double   UT1 as a 2-part Julian Date (Note 6)
**
**  Returned (function value):
**                int      status: +1 = dubious year (Note 3)
**                                  0 = OK
**                                 -1 = unacceptable date
**
**  Notes:
**
**  1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
**     convenient way between the two arguments, for example where utc1
**     is the Julian Day Number and utc2 is the fraction of a day.
**
**  2) JD cannot unambiguously represent UTC during a leap second unless
**     special measures are taken.  The convention in the present
**     function is that the JD day represents UTC days whether the
**     length is 86399, 86400 or 86401 SI seconds.
**
**  3) The warning status "dubious year" flags UTCs that predate the
**     introduction of the time scale or that are too far in the future
**     to be trusted.  See eraDat for further details.
**
**  4) The function eraDtf2d converts from calendar date and time of
**     day into 2-part Julian Date, and in the case of UTC implements
**     the leap-second-ambiguity convention described above.
**
**  5) Delta UT1 can be obtained from tabulations provided by the
**     International Earth Rotation and Reference Systems Service.
**     It is the caller's responsibility to supply a dut1 argument
**     containing the UT1-UTC value that matches the given UTC.
**
**  6) The returned ut11,ut12 are such that their sum is the UT1 Julian
**     Date.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Called:
**     eraJd2cal    JD to Gregorian calendar
**     eraDat       delta(AT) = TAI-UTC
**     eraUtctai    UTC to TAI
**     eraTaiut1    TAI to UT1
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var ut11 = 0.0;;
   var ut12 = 0.0;;
   var _rv1, _rv2, _rv3, _rv4;

   var iy, im, id, js, jw;
   var w, dat, dta, tai1, tai2;


/* Look up TAI-UTC. */
   var __rv__ = (_rv1 = eraJd2cal(utc1, utc2))[0];
   iy = ~~(_rv1[1]);
   im = ~~(_rv1[2]);
   id = ~~(_rv1[3]);
   w = _rv1[4];
   if(__rv__) return [ -1, ut11, ut12 ];
   js = ~~((_rv2 = eraDat(iy, im, id, 0.0))[0]);
   dat = _rv2[1];
   if ( js < 0 ) return [ -1, ut11, ut12 ];

/* Form UT1-TAI. */
   dta = dut1 - dat;

/* UTC to TAI to UT1. */
   jw = ~~((_rv3 = eraUtctai(utc1, utc2))[0]);
   tai1 = _rv3[1];
   tai2 = _rv3[2];
   if ( jw < 0 ) {
      return [ -1, ut11, ut12 ];
   } else if ( jw > 0 ) {
      js = ~~(jw);
   }
   var __rv__ = (_rv4 = eraTaiut1(tai1, tai2, dta))[0];
   ut11 = _rv4[1];
   ut12 = _rv4[2];
   if(__rv__) return [ -1, ut11, ut12 ];

/* Status. */
   return [ js, ut11, ut12 ];

}
;
function eraZp(p)
/*
**  - - - - - -
**   e r a Z p
**  - - - - - -
**
**  Zero a p-vector.
**
**  Returned:
**     p        double[3]      p-vector
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   if (typeof p == 'undefined') {
      p = [0, 0, 0];
   }

   p[0] = 0.0;
   p[1] = 0.0;
   p[2] = 0.0;

   return p;

}
;
function eraGst00a(uta, utb, tta, ttb)
/*
**  - - - - - - - - - -
**   e r a G s t 0 0 a
**  - - - - - - - - - -
**
**  Greenwich apparent sidereal time (consistent with IAU 2000
**  resolutions).
**
**  Given:
**     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
**     tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)
**
**  Returned (function value):
**                double    Greenwich apparent sidereal time (radians)
**
**  Notes:
**
**  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
**     Julian Dates, apportioned in any convenient way between the
**     argument pairs.  For example, JD=2450123.7 could be expressed in
**     any of these ways, among others:
**
**            Part A        Part B
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable (in the case of UT;  the TT is not at all critical
**     in this respect).  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  For UT, the date & time
**     method is best matched to the algorithm that is used by the Earth
**     Rotation Angle function, called internally:  maximum precision is
**     delivered when the uta argument is for 0hrs UT1 on the day in
**     question and the utb argument lies in the range 0 to 1, or vice
**     versa.
**
**  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
**     and TT to predict the effects of precession-nutation.  If UT1 is
**     used for both purposes, errors of order 100 microarcseconds
**     result.
**
**  3) This GAST is compatible with the IAU 2000 resolutions and must be
**     used only in conjunction with other IAU 2000 compatible
**     components such as precession-nutation.
**
**  4) The result is returned in the range 0 to 2pi.
**
**  5) The algorithm is from Capitaine et al. (2003) and IERS
**     Conventions 2003.
**
**  Called:
**     eraGmst00    Greenwich mean sidereal time, IAU 2000
**     eraEe00a     equation of the equinoxes, IAU 2000A
**     eraAnp       normalize angle into range 0 to 2pi
**
**  References:
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astronomy &
**     Astrophysics, 406, 1135-1149 (2003)
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var gmst00, ee00a, gst;


   gmst00 = eraGmst00(uta, utb, tta, ttb);
   ee00a = eraEe00a(tta, ttb);
   gst = eraAnp(gmst00 + ee00a);

   return gst;

}
;

function eraEe00a(date1, date2)
/*
**  - - - - - - - - -
**   e r a E e 0 0 a
**  - - - - - - - - -
**
**  Equation of the equinoxes, compatible with IAU 2000 resolutions.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    equation of the equinoxes (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The result, which is in radians, operates in the following sense:
**
**        Greenwich apparent ST = GMST + equation of the equinoxes
**
**  3) The result is compatible with the IAU 2000 resolutions.  For
**     further details, see IERS Conventions 2003 and Capitaine et al.
**     (2002).
**
**  Called:
**     eraPr00      IAU 2000 precession adjustments
**     eraObl80     mean obliquity, IAU 1980
**     eraNut00a    nutation, IAU 2000A
**     eraEe00      equation of the equinoxes, IAU 2000
**
**  References:
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astronomy &
**     Astrophysics, 406, 1135-1149 (2003).
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004).
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var _rv1, _rv3;

   var dpsipr, depspr, epsa, dpsi, deps, ee;


/* IAU 2000 precession-rate adjustments. */
   (_rv1 = eraPr00(date1, date2))[0];
   dpsipr = _rv1[0];
   depspr = _rv1[1];

/* Mean obliquity, consistent with IAU 2000 precession-nutation. */
   epsa = eraObl80(date1, date2) + depspr;

/* Nutation in longitude. */
   (_rv3 = eraNut00a(date1, date2))[0];
   dpsi = _rv3[0];
   deps = _rv3[1];

/* Equation of the equinoxes. */
   ee = eraEe00(date1, date2, epsa, dpsi);

   return ee;

}
;
function eraPr00(date1, date2)
/*
**  - - - - - - - -
**   e r a P r 0 0
**  - - - - - - - -
**
**  Precession-rate part of the IAU 2000 precession-nutation models
**  (part of MHB2000).
**
**  Given:
**     date1,date2    double  TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsipr,depspr  double  precession corrections (Notes 2,3)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The precession adjustments are expressed as "nutation
**     components", corrections in longitude and obliquity with respect
**     to the J2000.0 equinox and ecliptic.
**
**  3) Although the precession adjustments are stated to be with respect
**     to Lieske et al. (1977), the MHB2000 model does not specify which
**     set of Euler angles are to be used and how the adjustments are to
**     be applied.  The most literal and straightforward procedure is to
**     adopt the 4-rotation epsilon_0, psi_A, omega_A, xi_A option, and
**     to add dpsipr to psi_A and depspr to both omega_A and eps_A.
**
**  4) This is an implementation of one aspect of the IAU 2000A nutation
**     model, formally adopted by the IAU General Assembly in 2000,
**     namely MHB2000 (Mathews et al. 2002).
**
**  References:
**
**     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B., "Expressions
**     for the precession quantities based upon the IAU (1976) System of
**     Astronomical Constants", Astron.Astrophys., 58, 1-16 (1977)
**
**     Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
**     and precession   New nutation series for nonrigid Earth and
**     insights into the Earth's interior", J.Geophys.Res., 107, B4,
**     2002.  The MHB2000 code itself was obtained on 9th September 2002
**     from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
**
**     Wallace, P.T., "Software for Implementing the IAU 2000
**     Resolutions", in IERS Workshop 5.1 (2002).
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var dpsipr = 0.0;;
   var depspr = 0.0;;


   var t;

/* Precession and obliquity corrections (radians per century) */
   var PRECOR = -0.29965 * ERFA_DAS2R,
                       OBLCOR = -0.02524 * ERFA_DAS2R;


/* Interval between fundamental epoch J2000.0 and given date (JC). */
   t = ((date1 - ERFA_DJ00) + date2) / ERFA_DJC;

/* Precession rate contributions with respect to IAU 1976/80. */
   dpsipr = PRECOR * t;
   depspr = OBLCOR * t;

   return [dpsipr, depspr];

}
;
function eraObl80(date1, date2)
/*
**  - - - - - - - - -
**   e r a O b l 8 0
**  - - - - - - - - -
**
**  Mean obliquity of the ecliptic, IAU 1980 model.
**
**  Given:
**     date1,date2   double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                   double    obliquity of the ecliptic (radians, Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The result is the angle between the ecliptic and mean equator of
**     date date1+date2.
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Expression 3.222-1 (p114).
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var t, eps0;


/* Interval between fundamental epoch J2000.0 and given date (JC). */
   t = ((date1 - ERFA_DJ00) + date2) / ERFA_DJC;

/* Mean obliquity of date. */
   eps0 = ERFA_DAS2R * (84381.448  +
                  (-46.8150   +
                  (-0.00059   +
                  ( 0.001813) * t) * t) * t);

   return eps0;

}
;
function eraNut00a(date1, date2)
/*
**  - - - - - - - - - -
**   e r a N u t 0 0 a
**  - - - - - - - - - -
**
**  Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation
**  with free core nutation omitted).
**
**  Given:
**     date1,date2   double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi,deps     double   nutation, luni-solar + planetary (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The nutation components in longitude and obliquity are in radians
**     and with respect to the equinox and ecliptic of date.  The
**     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
**     value of 84381.448 arcsec.
**
**     Both the luni-solar and planetary nutations are included.  The
**     latter are due to direct planetary nutations and the
**     perturbations of the lunar and terrestrial orbits.
**
**  3) The function computes the MHB2000 nutation series with the
**     associated corrections for planetary nutations.  It is an
**     implementation of the nutation part of the IAU 2000A precession-
**     nutation model, formally adopted by the IAU General Assembly in
**     2000, namely MHB2000 (Mathews et al. 2002), but with the free
**     core nutation (FCN - see Note 4) omitted.
**
**  4) The full MHB2000 model also contains contributions to the
**     nutations in longitude and obliquity due to the free-excitation
**     of the free-core-nutation during the period 1979-2000.  These FCN
**     terms, which are time-dependent and unpredictable, are NOT
**     included in the present function and, if required, must be
**     independently computed.  With the FCN corrections included, the
**     present function delivers a pole which is at current epochs
**     accurate to a few hundred microarcseconds.  The omission of FCN
**     introduces further errors of about that size.
**
**  5) The present function provides classical nutation.  The MHB2000
**     algorithm, from which it is adapted, deals also with (i) the
**     offsets between the GCRS and mean poles and (ii) the adjustments
**     in longitude and obliquity due to the changed precession rates.
**     These additional functions, namely frame bias and precession
**     adjustments, are supported by the ERFA functions eraBi00  and
**     eraPr00.
**
**  6) The MHB2000 algorithm also provides "total" nutations, comprising
**     the arithmetic sum of the frame bias, precession adjustments,
**     luni-solar nutation and planetary nutation.  These total
**     nutations can be used in combination with an existing IAU 1976
**     precession implementation, such as eraPmat76,  to deliver GCRS-
**     to-true predictions of sub-mas accuracy at current dates.
**     However, there are three shortcomings in the MHB2000 model that
**     must be taken into account if more accurate or definitive results
**     are required (see Wallace 2002):
**
**       (i) The MHB2000 total nutations are simply arithmetic sums,
**           yet in reality the various components are successive Euler
**           rotations.  This slight lack of rigor leads to cross terms
**           that exceed 1 mas after a century.  The rigorous procedure
**           is to form the GCRS-to-true rotation matrix by applying the
**           bias, precession and nutation in that order.
**
**      (ii) Although the precession adjustments are stated to be with
**           respect to Lieske et al. (1977), the MHB2000 model does
**           not specify which set of Euler angles are to be used and
**           how the adjustments are to be applied.  The most literal
**           and straightforward procedure is to adopt the 4-rotation
**           epsilon_0, psi_A, omega_A, xi_A option, and to add DPSIPR
**           to psi_A and DEPSPR to both omega_A and eps_A.
**
**     (iii) The MHB2000 model predates the determination by Chapront
**           et al. (2002) of a 14.6 mas displacement between the
**           J2000.0 mean equinox and the origin of the ICRS frame.  It
**           should, however, be noted that neglecting this displacement
**           when calculating star coordinates does not lead to a
**           14.6 mas change in right ascension, only a small second-
**           order distortion in the pattern of the precession-nutation
**           effect.
**
**     For these reasons, the ERFA functions do not generate the "total
**     nutations" directly, though they can of course easily be
**     generated by calling eraBi00, eraPr00 and the present function
**     and adding the results.
**
**  7) The MHB2000 model contains 41 instances where the same frequency
**     appears multiple times, of which 38 are duplicates and three are
**     triplicates.  To keep the present code close to the original MHB
**     algorithm, this small inefficiency has not been corrected.
**
**  Called:
**     eraFal03     mean anomaly of the Moon
**     eraFaf03     mean argument of the latitude of the Moon
**     eraFaom03    mean longitude of the Moon's ascending node
**     eraFame03    mean longitude of Mercury
**     eraFave03    mean longitude of Venus
**     eraFae03     mean longitude of Earth
**     eraFama03    mean longitude of Mars
**     eraFaju03    mean longitude of Jupiter
**     eraFasa03    mean longitude of Saturn
**     eraFaur03    mean longitude of Uranus
**     eraFapa03    general accumulated precession in longitude
**
**  References:
**
**     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
**     Astron.Astrophys. 387, 700
**
**     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
**     Astron.Astrophys. 58, 1-16
**
**     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
**     107, B4.  The MHB_2000 code itself was obtained on 9th September
**     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**     Wallace, P.T., "Software for Implementing the IAU 2000
**     Resolutions", in IERS Workshop 5.1 (2002)
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var dpsi = 0.0;;
   var deps = 0.0;;


   var i;
   var t, el, elp, f, d, om, arg, dp, de, sarg, carg, al, af, ad, aom, alme, alve, alea, alma, alju, alsa, alur, alne, apa, dpsils, depsls, dpsipl, depspl;

/* Units of 0.1 microarcsecond to radians */
   var U2R = ERFA_DAS2R / 1e7;

/* ------------------------- */
/* Luni-Solar nutation model */
/* ------------------------- */

/* The units for the sine and cosine coefficients are */
/* 0.1 microarcsecond and the same per Julian century */

   var xls = [

   /* 1- 10 */
      [ 0, 0, 0, 0, 1,
         -172064161.0, -174666.0, 33386.0, 92052331.0, 9086.0, 15377.0],
      [ 0, 0, 2,-2, 2,
           -13170906.0, -1675.0, -13696.0, 5730336.0, -3015.0, -4587.0],
      [ 0, 0, 2, 0, 2,-2276413.0,-234.0,2796.0,978459.0,-485.0, 1374.0],
      [ 0, 0, 0, 0, 2,2074554.0, 207.0, -698.0,-897492.0,470.0, -291.0],
      [ 0, 1, 0, 0, 0,1475877.0,-3633.0,11817.0,73871.0,-184.0,-1924.0],
      [ 0, 1, 2,-2, 2,-516821.0,1226.0, -524.0,224386.0,-677.0, -174.0],
      [ 1, 0, 0, 0, 0, 711159.0,  73.0, -872.0,  -6750.0,  0.0,  358.0],
      [ 0, 0, 2, 0, 1,-387298.0,-367.0,  380.0, 200728.0, 18.0,  318.0],
      [ 1, 0, 2, 0, 2,-301461.0, -36.0,  816.0, 129025.0,-63.0,  367.0],
      [ 0,-1, 2,-2, 2, 215829.0,-494.0,  111.0, -95929.0,299.0,  132.0],

   /* 11-20 */
      [ 0, 0, 2,-2, 1, 128227.0, 137.0,  181.0, -68982.0, -9.0,   39.0],
      [-1, 0, 2, 0, 2, 123457.0,  11.0,   19.0, -53311.0, 32.0,   -4.0],
      [-1, 0, 0, 2, 0, 156994.0,  10.0, -168.0,  -1235.0,  0.0,   82.0],
      [ 1, 0, 0, 0, 1,  63110.0,  63.0,   27.0, -33228.0,  0.0,   -9.0],
      [-1, 0, 0, 0, 1, -57976.0, -63.0, -189.0,  31429.0,  0.0,  -75.0],
      [-1, 0, 2, 2, 2, -59641.0, -11.0,  149.0,  25543.0,-11.0,   66.0],
      [ 1, 0, 2, 0, 1, -51613.0, -42.0,  129.0,  26366.0,  0.0,   78.0],
      [-2, 0, 2, 0, 1,  45893.0,  50.0,   31.0, -24236.0,-10.0,   20.0],
      [ 0, 0, 0, 2, 0,  63384.0,  11.0, -150.0,  -1220.0,  0.0,   29.0],
      [ 0, 0, 2, 2, 2, -38571.0,  -1.0,  158.0,  16452.0,-11.0,   68.0],

   /* 21-30 */
      [ 0,-2, 2,-2, 2,  32481.0,   0.0,    0.0, -13870.0,  0.0,    0.0],
      [-2, 0, 0, 2, 0, -47722.0,   0.0,  -18.0,    477.0,  0.0,  -25.0],
      [ 2, 0, 2, 0, 2, -31046.0,  -1.0,  131.0,  13238.0,-11.0,   59.0],
      [ 1, 0, 2,-2, 2,  28593.0,   0.0,   -1.0, -12338.0, 10.0,   -3.0],
      [-1, 0, 2, 0, 1,  20441.0,  21.0,   10.0, -10758.0,  0.0,   -3.0],
      [ 2, 0, 0, 0, 0,  29243.0,   0.0,  -74.0,   -609.0,  0.0,   13.0],
      [ 0, 0, 2, 0, 0,  25887.0,   0.0,  -66.0,   -550.0,  0.0,   11.0],
      [ 0, 1, 0, 0, 1, -14053.0, -25.0,   79.0,   8551.0, -2.0,  -45.0],
      [-1, 0, 0, 2, 1,  15164.0,  10.0,   11.0,  -8001.0,  0.0,   -1.0],
      [ 0, 2, 2,-2, 2, -15794.0,  72.0,  -16.0,   6850.0,-42.0,   -5.0],

   /* 31-40 */
      [ 0, 0,-2, 2, 0,  21783.0,   0.0,   13.0,   -167.0,  0.0,   13.0],
      [ 1, 0, 0,-2, 1, -12873.0, -10.0,  -37.0,   6953.0,  0.0,  -14.0],
      [ 0,-1, 0, 0, 1, -12654.0,  11.0,   63.0,   6415.0,  0.0,   26.0],
      [-1, 0, 2, 2, 1, -10204.0,   0.0,   25.0,   5222.0,  0.0,   15.0],
      [ 0, 2, 0, 0, 0,  16707.0, -85.0,  -10.0,    168.0, -1.0,   10.0],
      [ 1, 0, 2, 2, 2,  -7691.0,   0.0,   44.0,   3268.0,  0.0,   19.0],
      [-2, 0, 2, 0, 0, -11024.0,   0.0,  -14.0,    104.0,  0.0,    2.0],
      [ 0, 1, 2, 0, 2,   7566.0, -21.0,  -11.0,  -3250.0,  0.0,   -5.0],
      [ 0, 0, 2, 2, 1,  -6637.0, -11.0,   25.0,   3353.0,  0.0,   14.0],
      [ 0,-1, 2, 0, 2,  -7141.0,  21.0,    8.0,   3070.0,  0.0,    4.0],

   /* 41-50 */
      [ 0, 0, 0, 2, 1,  -6302.0, -11.0,    2.0,   3272.0,  0.0,    4.0],
      [ 1, 0, 2,-2, 1,   5800.0,  10.0,    2.0,  -3045.0,  0.0,   -1.0],
      [ 2, 0, 2,-2, 2,   6443.0,   0.0,   -7.0,  -2768.0,  0.0,   -4.0],
      [-2, 0, 0, 2, 1,  -5774.0, -11.0,  -15.0,   3041.0,  0.0,   -5.0],
      [ 2, 0, 2, 0, 1,  -5350.0,   0.0,   21.0,   2695.0,  0.0,   12.0],
      [ 0,-1, 2,-2, 1,  -4752.0, -11.0,   -3.0,   2719.0,  0.0,   -3.0],
      [ 0, 0, 0,-2, 1,  -4940.0, -11.0,  -21.0,   2720.0,  0.0,   -9.0],
      [-1,-1, 0, 2, 0,   7350.0,   0.0,   -8.0,    -51.0,  0.0,    4.0],
      [ 2, 0, 0,-2, 1,   4065.0,   0.0,    6.0,  -2206.0,  0.0,    1.0],
      [ 1, 0, 0, 2, 0,   6579.0,   0.0,  -24.0,   -199.0,  0.0,    2.0],

   /* 51-60 */
      [ 0, 1, 2,-2, 1,   3579.0,   0.0,    5.0,  -1900.0,  0.0,    1.0],
      [ 1,-1, 0, 0, 0,   4725.0,   0.0,   -6.0,    -41.0,  0.0,    3.0],
      [-2, 0, 2, 0, 2,  -3075.0,   0.0,   -2.0,   1313.0,  0.0,   -1.0],
      [ 3, 0, 2, 0, 2,  -2904.0,   0.0,   15.0,   1233.0,  0.0,    7.0],
      [ 0,-1, 0, 2, 0,   4348.0,   0.0,  -10.0,    -81.0,  0.0,    2.0],
      [ 1,-1, 2, 0, 2,  -2878.0,   0.0,    8.0,   1232.0,  0.0,    4.0],
      [ 0, 0, 0, 1, 0,  -4230.0,   0.0,    5.0,    -20.0,  0.0,   -2.0],
      [-1,-1, 2, 2, 2,  -2819.0,   0.0,    7.0,   1207.0,  0.0,    3.0],
      [-1, 0, 2, 0, 0,  -4056.0,   0.0,    5.0,     40.0,  0.0,   -2.0],
      [ 0,-1, 2, 2, 2,  -2647.0,   0.0,   11.0,   1129.0,  0.0,    5.0],

   /* 61-70 */
      [-2, 0, 0, 0, 1,  -2294.0,   0.0,  -10.0,   1266.0,  0.0,   -4.0],
      [ 1, 1, 2, 0, 2,   2481.0,   0.0,   -7.0,  -1062.0,  0.0,   -3.0],
      [ 2, 0, 0, 0, 1,   2179.0,   0.0,   -2.0,  -1129.0,  0.0,   -2.0],
      [-1, 1, 0, 1, 0,   3276.0,   0.0,    1.0,     -9.0,  0.0,    0.0],
      [ 1, 1, 0, 0, 0,  -3389.0,   0.0,    5.0,     35.0,  0.0,   -2.0],
      [ 1, 0, 2, 0, 0,   3339.0,   0.0,  -13.0,   -107.0,  0.0,    1.0],
      [-1, 0, 2,-2, 1,  -1987.0,   0.0,   -6.0,   1073.0,  0.0,   -2.0],
      [ 1, 0, 0, 0, 2,  -1981.0,   0.0,    0.0,    854.0,  0.0,    0.0],
      [-1, 0, 0, 1, 0,   4026.0,   0.0, -353.0,   -553.0,  0.0, -139.0],
      [ 0, 0, 2, 1, 2,   1660.0,   0.0,   -5.0,   -710.0,  0.0,   -2.0],

   /* 71-80 */
      [-1, 0, 2, 4, 2,  -1521.0,   0.0,    9.0,    647.0,  0.0,    4.0],
      [-1, 1, 0, 1, 1,   1314.0,   0.0,    0.0,   -700.0,  0.0,    0.0],
      [ 0,-2, 2,-2, 1,  -1283.0,   0.0,    0.0,    672.0,  0.0,    0.0],
      [ 1, 0, 2, 2, 1,  -1331.0,   0.0,    8.0,    663.0,  0.0,    4.0],
      [-2, 0, 2, 2, 2,   1383.0,   0.0,   -2.0,   -594.0,  0.0,   -2.0],
      [-1, 0, 0, 0, 2,   1405.0,   0.0,    4.0,   -610.0,  0.0,    2.0],
      [ 1, 1, 2,-2, 2,   1290.0,   0.0,    0.0,   -556.0,  0.0,    0.0],
      [-2, 0, 2, 4, 2,  -1214.0,   0.0,    5.0,    518.0,  0.0,    2.0],
      [-1, 0, 4, 0, 2,   1146.0,   0.0,   -3.0,   -490.0,  0.0,   -1.0],
      [ 2, 0, 2,-2, 1,   1019.0,   0.0,   -1.0,   -527.0,  0.0,   -1.0],

   /* 81-90 */
      [ 2, 0, 2, 2, 2,  -1100.0,   0.0,    9.0,    465.0,  0.0,    4.0],
      [ 1, 0, 0, 2, 1,   -970.0,   0.0,    2.0,    496.0,  0.0,    1.0],
      [ 3, 0, 0, 0, 0,   1575.0,   0.0,   -6.0,    -50.0,  0.0,    0.0],
      [ 3, 0, 2,-2, 2,    934.0,   0.0,   -3.0,   -399.0,  0.0,   -1.0],
      [ 0, 0, 4,-2, 2,    922.0,   0.0,   -1.0,   -395.0,  0.0,   -1.0],
      [ 0, 1, 2, 0, 1,    815.0,   0.0,   -1.0,   -422.0,  0.0,   -1.0],
      [ 0, 0,-2, 2, 1,    834.0,   0.0,    2.0,   -440.0,  0.0,    1.0],
      [ 0, 0, 2,-2, 3,   1248.0,   0.0,    0.0,   -170.0,  0.0,    1.0],
      [-1, 0, 0, 4, 0,   1338.0,   0.0,   -5.0,    -39.0,  0.0,    0.0],
      [ 2, 0,-2, 0, 1,    716.0,   0.0,   -2.0,   -389.0,  0.0,   -1.0],

   /* 91-100 */
      [-2, 0, 0, 4, 0,   1282.0,   0.0,   -3.0,    -23.0,  0.0,    1.0],
      [-1,-1, 0, 2, 1,    742.0,   0.0,    1.0,   -391.0,  0.0,    0.0],
      [-1, 0, 0, 1, 1,   1020.0,   0.0,  -25.0,   -495.0,  0.0,  -10.0],
      [ 0, 1, 0, 0, 2,    715.0,   0.0,   -4.0,   -326.0,  0.0,    2.0],
      [ 0, 0,-2, 0, 1,   -666.0,   0.0,   -3.0,    369.0,  0.0,   -1.0],
      [ 0,-1, 2, 0, 1,   -667.0,   0.0,    1.0,    346.0,  0.0,    1.0],
      [ 0, 0, 2,-1, 2,   -704.0,   0.0,    0.0,    304.0,  0.0,    0.0],
      [ 0, 0, 2, 4, 2,   -694.0,   0.0,    5.0,    294.0,  0.0,    2.0],
      [-2,-1, 0, 2, 0,  -1014.0,   0.0,   -1.0,      4.0,  0.0,   -1.0],
      [ 1, 1, 0,-2, 1,   -585.0,   0.0,   -2.0,    316.0,  0.0,   -1.0],

   /* 101-110 */
      [-1, 1, 0, 2, 0,   -949.0,   0.0,    1.0,      8.0,  0.0,   -1.0],
      [-1, 1, 0, 1, 2,   -595.0,   0.0,    0.0,    258.0,  0.0,    0.0],
      [ 1,-1, 0, 0, 1,    528.0,   0.0,    0.0,   -279.0,  0.0,    0.0],
      [ 1,-1, 2, 2, 2,   -590.0,   0.0,    4.0,    252.0,  0.0,    2.0],
      [-1, 1, 2, 2, 2,    570.0,   0.0,   -2.0,   -244.0,  0.0,   -1.0],
      [ 3, 0, 2, 0, 1,   -502.0,   0.0,    3.0,    250.0,  0.0,    2.0],
      [ 0, 1,-2, 2, 0,   -875.0,   0.0,    1.0,     29.0,  0.0,    0.0],
      [-1, 0, 0,-2, 1,   -492.0,   0.0,   -3.0,    275.0,  0.0,   -1.0],
      [ 0, 1, 2, 2, 2,    535.0,   0.0,   -2.0,   -228.0,  0.0,   -1.0],
      [-1,-1, 2, 2, 1,   -467.0,   0.0,    1.0,    240.0,  0.0,    1.0],

   /* 111-120 */
      [ 0,-1, 0, 0, 2,    591.0,   0.0,    0.0,   -253.0,  0.0,    0.0],
      [ 1, 0, 2,-4, 1,   -453.0,   0.0,   -1.0,    244.0,  0.0,   -1.0],
      [-1, 0,-2, 2, 0,    766.0,   0.0,    1.0,      9.0,  0.0,    0.0],
      [ 0,-1, 2, 2, 1,   -446.0,   0.0,    2.0,    225.0,  0.0,    1.0],
      [ 2,-1, 2, 0, 2,   -488.0,   0.0,    2.0,    207.0,  0.0,    1.0],
      [ 0, 0, 0, 2, 2,   -468.0,   0.0,    0.0,    201.0,  0.0,    0.0],
      [ 1,-1, 2, 0, 1,   -421.0,   0.0,    1.0,    216.0,  0.0,    1.0],
      [-1, 1, 2, 0, 2,    463.0,   0.0,    0.0,   -200.0,  0.0,    0.0],
      [ 0, 1, 0, 2, 0,   -673.0,   0.0,    2.0,     14.0,  0.0,    0.0],
      [ 0,-1,-2, 2, 0,    658.0,   0.0,    0.0,     -2.0,  0.0,    0.0],

   /* 121-130 */
      [ 0, 3, 2,-2, 2,   -438.0,   0.0,    0.0,    188.0,  0.0,    0.0],
      [ 0, 0, 0, 1, 1,   -390.0,   0.0,    0.0,    205.0,  0.0,    0.0],
      [-1, 0, 2, 2, 0,    639.0, -11.0,   -2.0,    -19.0,  0.0,    0.0],
      [ 2, 1, 2, 0, 2,    412.0,   0.0,   -2.0,   -176.0,  0.0,   -1.0],
      [ 1, 1, 0, 0, 1,   -361.0,   0.0,    0.0,    189.0,  0.0,    0.0],
      [ 1, 1, 2, 0, 1,    360.0,   0.0,   -1.0,   -185.0,  0.0,   -1.0],
      [ 2, 0, 0, 2, 0,    588.0,   0.0,   -3.0,    -24.0,  0.0,    0.0],
      [ 1, 0,-2, 2, 0,   -578.0,   0.0,    1.0,      5.0,  0.0,    0.0],
      [-1, 0, 0, 2, 2,   -396.0,   0.0,    0.0,    171.0,  0.0,    0.0],
      [ 0, 1, 0, 1, 0,    565.0,   0.0,   -1.0,     -6.0,  0.0,    0.0],

   /* 131-140 */
      [ 0, 1, 0,-2, 1,   -335.0,   0.0,   -1.0,    184.0,  0.0,   -1.0],
      [-1, 0, 2,-2, 2,    357.0,   0.0,    1.0,   -154.0,  0.0,    0.0],
      [ 0, 0, 0,-1, 1,    321.0,   0.0,    1.0,   -174.0,  0.0,    0.0],
      [-1, 1, 0, 0, 1,   -301.0,   0.0,   -1.0,    162.0,  0.0,    0.0],
      [ 1, 0, 2,-1, 2,   -334.0,   0.0,    0.0,    144.0,  0.0,    0.0],
      [ 1,-1, 0, 2, 0,    493.0,   0.0,   -2.0,    -15.0,  0.0,    0.0],
      [ 0, 0, 0, 4, 0,    494.0,   0.0,   -2.0,    -19.0,  0.0,    0.0],
      [ 1, 0, 2, 1, 2,    337.0,   0.0,   -1.0,   -143.0,  0.0,   -1.0],
      [ 0, 0, 2, 1, 1,    280.0,   0.0,   -1.0,   -144.0,  0.0,    0.0],
      [ 1, 0, 0,-2, 2,    309.0,   0.0,    1.0,   -134.0,  0.0,    0.0],

   /* 141-150 */
      [-1, 0, 2, 4, 1,   -263.0,   0.0,    2.0,    131.0,  0.0,    1.0],
      [ 1, 0,-2, 0, 1,    253.0,   0.0,    1.0,   -138.0,  0.0,    0.0],
      [ 1, 1, 2,-2, 1,    245.0,   0.0,    0.0,   -128.0,  0.0,    0.0],
      [ 0, 0, 2, 2, 0,    416.0,   0.0,   -2.0,    -17.0,  0.0,    0.0],
      [-1, 0, 2,-1, 1,   -229.0,   0.0,    0.0,    128.0,  0.0,    0.0],
      [-2, 0, 2, 2, 1,    231.0,   0.0,    0.0,   -120.0,  0.0,    0.0],
      [ 4, 0, 2, 0, 2,   -259.0,   0.0,    2.0,    109.0,  0.0,    1.0],
      [ 2,-1, 0, 0, 0,    375.0,   0.0,   -1.0,     -8.0,  0.0,    0.0],
      [ 2, 1, 2,-2, 2,    252.0,   0.0,    0.0,   -108.0,  0.0,    0.0],
      [ 0, 1, 2, 1, 2,   -245.0,   0.0,    1.0,    104.0,  0.0,    0.0],

   /* 151-160 */
      [ 1, 0, 4,-2, 2,    243.0,   0.0,   -1.0,   -104.0,  0.0,    0.0],
      [-1,-1, 0, 0, 1,    208.0,   0.0,    1.0,   -112.0,  0.0,    0.0],
      [ 0, 1, 0, 2, 1,    199.0,   0.0,    0.0,   -102.0,  0.0,    0.0],
      [-2, 0, 2, 4, 1,   -208.0,   0.0,    1.0,    105.0,  0.0,    0.0],
      [ 2, 0, 2, 0, 0,    335.0,   0.0,   -2.0,    -14.0,  0.0,    0.0],
      [ 1, 0, 0, 1, 0,   -325.0,   0.0,    1.0,      7.0,  0.0,    0.0],
      [-1, 0, 0, 4, 1,   -187.0,   0.0,    0.0,     96.0,  0.0,    0.0],
      [-1, 0, 4, 0, 1,    197.0,   0.0,   -1.0,   -100.0,  0.0,    0.0],
      [ 2, 0, 2, 2, 1,   -192.0,   0.0,    2.0,     94.0,  0.0,    1.0],
      [ 0, 0, 2,-3, 2,   -188.0,   0.0,    0.0,     83.0,  0.0,    0.0],

   /* 161-170 */
      [-1,-2, 0, 2, 0,    276.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 2, 1, 0, 0, 0,   -286.0,   0.0,    1.0,      6.0,  0.0,    0.0],
      [ 0, 0, 4, 0, 2,    186.0,   0.0,   -1.0,    -79.0,  0.0,    0.0],
      [ 0, 0, 0, 0, 3,   -219.0,   0.0,    0.0,     43.0,  0.0,    0.0],
      [ 0, 3, 0, 0, 0,    276.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 0, 0, 2,-4, 1,   -153.0,   0.0,   -1.0,     84.0,  0.0,    0.0],
      [ 0,-1, 0, 2, 1,   -156.0,   0.0,    0.0,     81.0,  0.0,    0.0],
      [ 0, 0, 0, 4, 1,   -154.0,   0.0,    1.0,     78.0,  0.0,    0.0],
      [-1,-1, 2, 4, 2,   -174.0,   0.0,    1.0,     75.0,  0.0,    0.0],
      [ 1, 0, 2, 4, 2,   -163.0,   0.0,    2.0,     69.0,  0.0,    1.0],

   /* 171-180 */
      [-2, 2, 0, 2, 0,   -228.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [-2,-1, 2, 0, 1,     91.0,   0.0,   -4.0,    -54.0,  0.0,   -2.0],
      [-2, 0, 0, 2, 2,    175.0,   0.0,    0.0,    -75.0,  0.0,    0.0],
      [-1,-1, 2, 0, 2,   -159.0,   0.0,    0.0,     69.0,  0.0,    0.0],
      [ 0, 0, 4,-2, 1,    141.0,   0.0,    0.0,    -72.0,  0.0,    0.0],
      [ 3, 0, 2,-2, 1,    147.0,   0.0,    0.0,    -75.0,  0.0,    0.0],
      [-2,-1, 0, 2, 1,   -132.0,   0.0,    0.0,     69.0,  0.0,    0.0],
      [ 1, 0, 0,-1, 1,    159.0,   0.0,  -28.0,    -54.0,  0.0,   11.0],
      [ 0,-2, 0, 2, 0,    213.0,   0.0,    0.0,     -4.0,  0.0,    0.0],
      [-2, 0, 0, 4, 1,    123.0,   0.0,    0.0,    -64.0,  0.0,    0.0],

   /* 181-190 */
      [-3, 0, 0, 0, 1,   -118.0,   0.0,   -1.0,     66.0,  0.0,    0.0],
      [ 1, 1, 2, 2, 2,    144.0,   0.0,   -1.0,    -61.0,  0.0,    0.0],
      [ 0, 0, 2, 4, 1,   -121.0,   0.0,    1.0,     60.0,  0.0,    0.0],
      [ 3, 0, 2, 2, 2,   -134.0,   0.0,    1.0,     56.0,  0.0,    1.0],
      [-1, 1, 2,-2, 1,   -105.0,   0.0,    0.0,     57.0,  0.0,    0.0],
      [ 2, 0, 0,-4, 1,   -102.0,   0.0,    0.0,     56.0,  0.0,    0.0],
      [ 0, 0, 0,-2, 2,    120.0,   0.0,    0.0,    -52.0,  0.0,    0.0],
      [ 2, 0, 2,-4, 1,    101.0,   0.0,    0.0,    -54.0,  0.0,    0.0],
      [-1, 1, 0, 2, 1,   -113.0,   0.0,    0.0,     59.0,  0.0,    0.0],
      [ 0, 0, 2,-1, 1,   -106.0,   0.0,    0.0,     61.0,  0.0,    0.0],

   /* 191-200 */
      [ 0,-2, 2, 2, 2,   -129.0,   0.0,    1.0,     55.0,  0.0,    0.0],
      [ 2, 0, 0, 2, 1,   -114.0,   0.0,    0.0,     57.0,  0.0,    0.0],
      [ 4, 0, 2,-2, 2,    113.0,   0.0,   -1.0,    -49.0,  0.0,    0.0],
      [ 2, 0, 0,-2, 2,   -102.0,   0.0,    0.0,     44.0,  0.0,    0.0],
      [ 0, 2, 0, 0, 1,    -94.0,   0.0,    0.0,     51.0,  0.0,    0.0],
      [ 1, 0, 0,-4, 1,   -100.0,   0.0,   -1.0,     56.0,  0.0,    0.0],
      [ 0, 2, 2,-2, 1,     87.0,   0.0,    0.0,    -47.0,  0.0,    0.0],
      [-3, 0, 0, 4, 0,    161.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [-1, 1, 2, 0, 1,     96.0,   0.0,    0.0,    -50.0,  0.0,    0.0],
      [-1,-1, 0, 4, 0,    151.0,   0.0,   -1.0,     -5.0,  0.0,    0.0],

   /* 201-210 */
      [-1,-2, 2, 2, 2,   -104.0,   0.0,    0.0,     44.0,  0.0,    0.0],
      [-2,-1, 2, 4, 2,   -110.0,   0.0,    0.0,     48.0,  0.0,    0.0],
      [ 1,-1, 2, 2, 1,   -100.0,   0.0,    1.0,     50.0,  0.0,    0.0],
      [-2, 1, 0, 2, 0,     92.0,   0.0,   -5.0,     12.0,  0.0,   -2.0],
      [-2, 1, 2, 0, 1,     82.0,   0.0,    0.0,    -45.0,  0.0,    0.0],
      [ 2, 1, 0,-2, 1,     82.0,   0.0,    0.0,    -45.0,  0.0,    0.0],
      [-3, 0, 2, 0, 1,    -78.0,   0.0,    0.0,     41.0,  0.0,    0.0],
      [-2, 0, 2,-2, 1,    -77.0,   0.0,    0.0,     43.0,  0.0,    0.0],
      [-1, 1, 0, 2, 2,      2.0,   0.0,    0.0,     54.0,  0.0,    0.0],
      [ 0,-1, 2,-1, 2,     94.0,   0.0,    0.0,    -40.0,  0.0,    0.0],

   /* 211-220 */
      [-1, 0, 4,-2, 2,    -93.0,   0.0,    0.0,     40.0,  0.0,    0.0],
      [ 0,-2, 2, 0, 2,    -83.0,   0.0,   10.0,     40.0,  0.0,   -2.0],
      [-1, 0, 2, 1, 2,     83.0,   0.0,    0.0,    -36.0,  0.0,    0.0],
      [ 2, 0, 0, 0, 2,    -91.0,   0.0,    0.0,     39.0,  0.0,    0.0],
      [ 0, 0, 2, 0, 3,    128.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [-2, 0, 4, 0, 2,    -79.0,   0.0,    0.0,     34.0,  0.0,    0.0],
      [-1, 0,-2, 0, 1,    -83.0,   0.0,    0.0,     47.0,  0.0,    0.0],
      [-1, 1, 2, 2, 1,     84.0,   0.0,    0.0,    -44.0,  0.0,    0.0],
      [ 3, 0, 0, 0, 1,     83.0,   0.0,    0.0,    -43.0,  0.0,    0.0],
      [-1, 0, 2, 3, 2,     91.0,   0.0,    0.0,    -39.0,  0.0,    0.0],

   /* 221-230 */
      [ 2,-1, 2, 0, 1,    -77.0,   0.0,    0.0,     39.0,  0.0,    0.0],
      [ 0, 1, 2, 2, 1,     84.0,   0.0,    0.0,    -43.0,  0.0,    0.0],
      [ 0,-1, 2, 4, 2,    -92.0,   0.0,    1.0,     39.0,  0.0,    0.0],
      [ 2,-1, 2, 2, 2,    -92.0,   0.0,    1.0,     39.0,  0.0,    0.0],
      [ 0, 2,-2, 2, 0,    -94.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1,-1, 2,-1, 1,     68.0,   0.0,    0.0,    -36.0,  0.0,    0.0],
      [ 0,-2, 0, 0, 1,    -61.0,   0.0,    0.0,     32.0,  0.0,    0.0],
      [ 1, 0, 2,-4, 2,     71.0,   0.0,    0.0,    -31.0,  0.0,    0.0],
      [ 1,-1, 0,-2, 1,     62.0,   0.0,    0.0,    -34.0,  0.0,    0.0],
      [-1,-1, 2, 0, 1,    -63.0,   0.0,    0.0,     33.0,  0.0,    0.0],

   /* 231-240 */
      [ 1,-1, 2,-2, 2,    -73.0,   0.0,    0.0,     32.0,  0.0,    0.0],
      [-2,-1, 0, 4, 0,    115.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [-1, 0, 0, 3, 0,   -103.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-2,-1, 2, 2, 2,     63.0,   0.0,    0.0,    -28.0,  0.0,    0.0],
      [ 0, 2, 2, 0, 2,     74.0,   0.0,    0.0,    -32.0,  0.0,    0.0],
      [ 1, 1, 0, 2, 0,   -103.0,   0.0,   -3.0,      3.0,  0.0,   -1.0],
      [ 2, 0, 2,-1, 2,    -69.0,   0.0,    0.0,     30.0,  0.0,    0.0],
      [ 1, 0, 2, 1, 1,     57.0,   0.0,    0.0,    -29.0,  0.0,    0.0],
      [ 4, 0, 0, 0, 0,     94.0,   0.0,    0.0,     -4.0,  0.0,    0.0],
      [ 2, 1, 2, 0, 1,     64.0,   0.0,    0.0,    -33.0,  0.0,    0.0],

   /* 241-250 */
      [ 3,-1, 2, 0, 2,    -63.0,   0.0,    0.0,     26.0,  0.0,    0.0],
      [-2, 2, 0, 2, 1,    -38.0,   0.0,    0.0,     20.0,  0.0,    0.0],
      [ 1, 0, 2,-3, 1,    -43.0,   0.0,    0.0,     24.0,  0.0,    0.0],
      [ 1, 1, 2,-4, 1,    -45.0,   0.0,    0.0,     23.0,  0.0,    0.0],
      [-1,-1, 2,-2, 1,     47.0,   0.0,    0.0,    -24.0,  0.0,    0.0],
      [ 0,-1, 0,-1, 1,    -48.0,   0.0,    0.0,     25.0,  0.0,    0.0],
      [ 0,-1, 0,-2, 1,     45.0,   0.0,    0.0,    -26.0,  0.0,    0.0],
      [-2, 0, 0, 0, 2,     56.0,   0.0,    0.0,    -25.0,  0.0,    0.0],
      [-2, 0,-2, 2, 0,     88.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-1, 0,-2, 4, 0,    -75.0,   0.0,    0.0,      0.0,  0.0,    0.0],

   /* 251-260 */
      [ 1,-2, 0, 0, 0,     85.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 1, 0, 1, 1,     49.0,   0.0,    0.0,    -26.0,  0.0,    0.0],
      [-1, 2, 0, 2, 0,    -74.0,   0.0,   -3.0,     -1.0,  0.0,   -1.0],
      [ 1,-1, 2,-2, 1,    -39.0,   0.0,    0.0,     21.0,  0.0,    0.0],
      [ 1, 2, 2,-2, 2,     45.0,   0.0,    0.0,    -20.0,  0.0,    0.0],
      [ 2,-1, 2,-2, 2,     51.0,   0.0,    0.0,    -22.0,  0.0,    0.0],
      [ 1, 0, 2,-1, 1,    -40.0,   0.0,    0.0,     21.0,  0.0,    0.0],
      [ 2, 1, 2,-2, 1,     41.0,   0.0,    0.0,    -21.0,  0.0,    0.0],
      [-2, 0, 0,-2, 1,    -42.0,   0.0,    0.0,     24.0,  0.0,    0.0],
      [ 1,-2, 2, 0, 2,    -51.0,   0.0,    0.0,     22.0,  0.0,    0.0],

   /* 261-270 */
      [ 0, 1, 2, 1, 1,    -42.0,   0.0,    0.0,     22.0,  0.0,    0.0],
      [ 1, 0, 4,-2, 1,     39.0,   0.0,    0.0,    -21.0,  0.0,    0.0],
      [-2, 0, 4, 2, 2,     46.0,   0.0,    0.0,    -18.0,  0.0,    0.0],
      [ 1, 1, 2, 1, 2,    -53.0,   0.0,    0.0,     22.0,  0.0,    0.0],
      [ 1, 0, 0, 4, 0,     82.0,   0.0,    0.0,     -4.0,  0.0,    0.0],
      [ 1, 0, 2, 2, 0,     81.0,   0.0,   -1.0,     -4.0,  0.0,    0.0],
      [ 2, 0, 2, 1, 2,     47.0,   0.0,    0.0,    -19.0,  0.0,    0.0],
      [ 3, 1, 2, 0, 2,     53.0,   0.0,    0.0,    -23.0,  0.0,    0.0],
      [ 4, 0, 2, 0, 1,    -45.0,   0.0,    0.0,     22.0,  0.0,    0.0],
      [-2,-1, 2, 0, 0,    -44.0,   0.0,    0.0,     -2.0,  0.0,    0.0],

   /* 271-280 */
      [ 0, 1,-2, 2, 1,    -33.0,   0.0,    0.0,     16.0,  0.0,    0.0],
      [ 1, 0,-2, 1, 0,    -61.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [ 0,-1,-2, 2, 1,     28.0,   0.0,    0.0,    -15.0,  0.0,    0.0],
      [ 2,-1, 0,-2, 1,    -38.0,   0.0,    0.0,     19.0,  0.0,    0.0],
      [-1, 0, 2,-1, 2,    -33.0,   0.0,    0.0,     21.0,  0.0,    0.0],
      [ 1, 0, 2,-3, 2,    -60.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 1, 2,-2, 3,     48.0,   0.0,    0.0,    -10.0,  0.0,    0.0],
      [ 0, 0, 2,-3, 1,     27.0,   0.0,    0.0,    -14.0,  0.0,    0.0],
      [-1, 0,-2, 2, 1,     38.0,   0.0,    0.0,    -20.0,  0.0,    0.0],
      [ 0, 0, 2,-4, 2,     31.0,   0.0,    0.0,    -13.0,  0.0,    0.0],

   /* 281-290 */
      [-2, 1, 0, 0, 1,    -29.0,   0.0,    0.0,     15.0,  0.0,    0.0],
      [-1, 0, 0,-1, 1,     28.0,   0.0,    0.0,    -15.0,  0.0,    0.0],
      [ 2, 0, 2,-4, 2,    -32.0,   0.0,    0.0,     15.0,  0.0,    0.0],
      [ 0, 0, 4,-4, 4,     45.0,   0.0,    0.0,     -8.0,  0.0,    0.0],
      [ 0, 0, 4,-4, 2,    -44.0,   0.0,    0.0,     19.0,  0.0,    0.0],
      [-1,-2, 0, 2, 1,     28.0,   0.0,    0.0,    -15.0,  0.0,    0.0],
      [-2, 0, 0, 3, 0,    -51.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1, 0,-2, 2, 1,    -36.0,   0.0,    0.0,     20.0,  0.0,    0.0],
      [-3, 0, 2, 2, 2,     44.0,   0.0,    0.0,    -19.0,  0.0,    0.0],
      [-3, 0, 2, 2, 1,     26.0,   0.0,    0.0,    -14.0,  0.0,    0.0],

   /* 291-300 */
      [-2, 0, 2, 2, 0,    -60.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 2,-1, 0, 0, 1,     35.0,   0.0,    0.0,    -18.0,  0.0,    0.0],
      [-2, 1, 2, 2, 2,    -27.0,   0.0,    0.0,     11.0,  0.0,    0.0],
      [ 1, 1, 0, 1, 0,     47.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [ 0, 1, 4,-2, 2,     36.0,   0.0,    0.0,    -15.0,  0.0,    0.0],
      [-1, 1, 0,-2, 1,    -36.0,   0.0,    0.0,     20.0,  0.0,    0.0],
      [ 0, 0, 0,-4, 1,    -35.0,   0.0,    0.0,     19.0,  0.0,    0.0],
      [ 1,-1, 0, 2, 1,    -37.0,   0.0,    0.0,     19.0,  0.0,    0.0],
      [ 1, 1, 0, 2, 1,     32.0,   0.0,    0.0,    -16.0,  0.0,    0.0],
      [-1, 2, 2, 2, 2,     35.0,   0.0,    0.0,    -14.0,  0.0,    0.0],

   /* 301-310 */
      [ 3, 1, 2,-2, 2,     32.0,   0.0,    0.0,    -13.0,  0.0,    0.0],
      [ 0,-1, 0, 4, 0,     65.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 2,-1, 0, 2, 0,     47.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [ 0, 0, 4, 0, 1,     32.0,   0.0,    0.0,    -16.0,  0.0,    0.0],
      [ 2, 0, 4,-2, 2,     37.0,   0.0,    0.0,    -16.0,  0.0,    0.0],
      [-1,-1, 2, 4, 1,    -30.0,   0.0,    0.0,     15.0,  0.0,    0.0],
      [ 1, 0, 0, 4, 1,    -32.0,   0.0,    0.0,     16.0,  0.0,    0.0],
      [ 1,-2, 2, 2, 2,    -31.0,   0.0,    0.0,     13.0,  0.0,    0.0],
      [ 0, 0, 2, 3, 2,     37.0,   0.0,    0.0,    -16.0,  0.0,    0.0],
      [-1, 1, 2, 4, 2,     31.0,   0.0,    0.0,    -13.0,  0.0,    0.0],

   /* 311-320 */
      [ 3, 0, 0, 2, 0,     49.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [-1, 0, 4, 2, 2,     32.0,   0.0,    0.0,    -13.0,  0.0,    0.0],
      [ 1, 1, 2, 2, 1,     23.0,   0.0,    0.0,    -12.0,  0.0,    0.0],
      [-2, 0, 2, 6, 2,    -43.0,   0.0,    0.0,     18.0,  0.0,    0.0],
      [ 2, 1, 2, 2, 2,     26.0,   0.0,    0.0,    -11.0,  0.0,    0.0],
      [-1, 0, 2, 6, 2,    -32.0,   0.0,    0.0,     14.0,  0.0,    0.0],
      [ 1, 0, 2, 4, 1,    -29.0,   0.0,    0.0,     14.0,  0.0,    0.0],
      [ 2, 0, 2, 4, 2,    -27.0,   0.0,    0.0,     12.0,  0.0,    0.0],
      [ 1, 1,-2, 1, 0,     30.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-3, 1, 2, 1, 2,    -11.0,   0.0,    0.0,      5.0,  0.0,    0.0],

   /* 321-330 */
      [ 2, 0,-2, 0, 2,    -21.0,   0.0,    0.0,     10.0,  0.0,    0.0],
      [-1, 0, 0, 1, 2,    -34.0,   0.0,    0.0,     15.0,  0.0,    0.0],
      [-4, 0, 2, 2, 1,    -10.0,   0.0,    0.0,      6.0,  0.0,    0.0],
      [-1,-1, 0, 1, 0,    -36.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 0,-2, 2, 2,     -9.0,   0.0,    0.0,      4.0,  0.0,    0.0],
      [ 1, 0, 0,-1, 2,    -12.0,   0.0,    0.0,      5.0,  0.0,    0.0],
      [ 0,-1, 2,-2, 3,    -21.0,   0.0,    0.0,      5.0,  0.0,    0.0],
      [-2, 1, 2, 0, 0,    -29.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [ 0, 0, 2,-2, 4,    -15.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [-2,-2, 0, 2, 0,    -20.0,   0.0,    0.0,      0.0,  0.0,    0.0],

   /* 331-340 */
      [-2, 0,-2, 4, 0,     28.0,   0.0,    0.0,      0.0,  0.0,   -2.0],
      [ 0,-2,-2, 2, 0,     17.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1, 2, 0,-2, 1,    -22.0,   0.0,    0.0,     12.0,  0.0,    0.0],
      [ 3, 0, 0,-4, 1,    -14.0,   0.0,    0.0,      7.0,  0.0,    0.0],
      [-1, 1, 2,-2, 2,     24.0,   0.0,    0.0,    -11.0,  0.0,    0.0],
      [ 1,-1, 2,-4, 1,     11.0,   0.0,    0.0,     -6.0,  0.0,    0.0],
      [ 1, 1, 0,-2, 2,     14.0,   0.0,    0.0,     -6.0,  0.0,    0.0],
      [-3, 0, 2, 0, 0,     24.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-3, 0, 2, 0, 2,     18.0,   0.0,    0.0,     -8.0,  0.0,    0.0],
      [-2, 0, 0, 1, 0,    -38.0,   0.0,    0.0,      0.0,  0.0,    0.0],

   /* 341-350 */
      [ 0, 0,-2, 1, 0,    -31.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-3, 0, 0, 2, 1,    -16.0,   0.0,    0.0,      8.0,  0.0,    0.0],
      [-1,-1,-2, 2, 0,     29.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 1, 2,-4, 1,    -18.0,   0.0,    0.0,     10.0,  0.0,    0.0],
      [ 2, 1, 0,-4, 1,    -10.0,   0.0,    0.0,      5.0,  0.0,    0.0],
      [ 0, 2, 0,-2, 1,    -17.0,   0.0,    0.0,     10.0,  0.0,    0.0],
      [ 1, 0, 0,-3, 1,      9.0,   0.0,    0.0,     -4.0,  0.0,    0.0],
      [-2, 0, 2,-2, 2,     16.0,   0.0,    0.0,     -6.0,  0.0,    0.0],
      [-2,-1, 0, 0, 1,     22.0,   0.0,    0.0,    -12.0,  0.0,    0.0],
      [-4, 0, 0, 2, 0,     20.0,   0.0,    0.0,      0.0,  0.0,    0.0],

   /* 351-360 */
      [ 1, 1, 0,-4, 1,    -13.0,   0.0,    0.0,      6.0,  0.0,    0.0],
      [-1, 0, 2,-4, 1,    -17.0,   0.0,    0.0,      9.0,  0.0,    0.0],
      [ 0, 0, 4,-4, 1,    -14.0,   0.0,    0.0,      8.0,  0.0,    0.0],
      [ 0, 3, 2,-2, 2,      0.0,   0.0,    0.0,     -7.0,  0.0,    0.0],
      [-3,-1, 0, 4, 0,     14.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-3, 0, 0, 4, 1,     19.0,   0.0,    0.0,    -10.0,  0.0,    0.0],
      [ 1,-1,-2, 2, 0,    -34.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1,-1, 0, 2, 2,    -20.0,   0.0,    0.0,      8.0,  0.0,    0.0],
      [ 1,-2, 0, 0, 1,      9.0,   0.0,    0.0,     -5.0,  0.0,    0.0],
      [ 1,-1, 0, 0, 2,    -18.0,   0.0,    0.0,      7.0,  0.0,    0.0],

   /* 361-370 */
      [ 0, 0, 0, 1, 2,     13.0,   0.0,    0.0,     -6.0,  0.0,    0.0],
      [-1,-1, 2, 0, 0,     17.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1,-2, 2,-2, 2,    -12.0,   0.0,    0.0,      5.0,  0.0,    0.0],
      [ 0,-1, 2,-1, 1,     15.0,   0.0,    0.0,     -8.0,  0.0,    0.0],
      [-1, 0, 2, 0, 3,    -11.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [ 1, 1, 0, 0, 2,     13.0,   0.0,    0.0,     -5.0,  0.0,    0.0],
      [-1, 1, 2, 0, 0,    -18.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1, 2, 0, 0, 0,    -35.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1, 2, 2, 0, 2,      9.0,   0.0,    0.0,     -4.0,  0.0,    0.0],
      [-1, 0, 4,-2, 1,    -19.0,   0.0,    0.0,     10.0,  0.0,    0.0],

   /* 371-380 */
      [ 3, 0, 2,-4, 2,    -26.0,   0.0,    0.0,     11.0,  0.0,    0.0],
      [ 1, 2, 2,-2, 1,      8.0,   0.0,    0.0,     -4.0,  0.0,    0.0],
      [ 1, 0, 4,-4, 2,    -10.0,   0.0,    0.0,      4.0,  0.0,    0.0],
      [-2,-1, 0, 4, 1,     10.0,   0.0,    0.0,     -6.0,  0.0,    0.0],
      [ 0,-1, 0, 2, 2,    -21.0,   0.0,    0.0,      9.0,  0.0,    0.0],
      [-2, 1, 0, 4, 0,    -15.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-2,-1, 2, 2, 1,      9.0,   0.0,    0.0,     -5.0,  0.0,    0.0],
      [ 2, 0,-2, 2, 0,    -29.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1, 0, 0, 1, 1,    -19.0,   0.0,    0.0,     10.0,  0.0,    0.0],
      [ 0, 1, 0, 2, 2,     12.0,   0.0,    0.0,     -5.0,  0.0,    0.0],

   /* 381-390 */
      [ 1,-1, 2,-1, 2,     22.0,   0.0,    0.0,     -9.0,  0.0,    0.0],
      [-2, 0, 4, 0, 1,    -10.0,   0.0,    0.0,      5.0,  0.0,    0.0],
      [ 2, 1, 0, 0, 1,    -20.0,   0.0,    0.0,     11.0,  0.0,    0.0],
      [ 0, 1, 2, 0, 0,    -20.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0,-1, 4,-2, 2,    -17.0,   0.0,    0.0,      7.0,  0.0,    0.0],
      [ 0, 0, 4,-2, 4,     15.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [ 0, 2, 2, 0, 1,      8.0,   0.0,    0.0,     -4.0,  0.0,    0.0],
      [-3, 0, 0, 6, 0,     14.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1,-1, 0, 4, 1,    -12.0,   0.0,    0.0,      6.0,  0.0,    0.0],
      [ 1,-2, 0, 2, 0,     25.0,   0.0,    0.0,      0.0,  0.0,    0.0],

   /* 391-400 */
      [-1, 0, 0, 4, 2,    -13.0,   0.0,    0.0,      6.0,  0.0,    0.0],
      [-1,-2, 2, 2, 1,    -14.0,   0.0,    0.0,      8.0,  0.0,    0.0],
      [-1, 0, 0,-2, 2,     13.0,   0.0,    0.0,     -5.0,  0.0,    0.0],
      [ 1, 0,-2,-2, 1,    -17.0,   0.0,    0.0,      9.0,  0.0,    0.0],
      [ 0, 0,-2,-2, 1,    -12.0,   0.0,    0.0,      6.0,  0.0,    0.0],
      [-2, 0,-2, 0, 1,    -10.0,   0.0,    0.0,      5.0,  0.0,    0.0],
      [ 0, 0, 0, 3, 1,     10.0,   0.0,    0.0,     -6.0,  0.0,    0.0],
      [ 0, 0, 0, 3, 0,    -15.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1, 1, 0, 4, 0,    -22.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1,-1, 2, 2, 0,     28.0,   0.0,    0.0,     -1.0,  0.0,    0.0],

   /* 401-410 */
      [-2, 0, 2, 3, 2,     15.0,   0.0,    0.0,     -7.0,  0.0,    0.0],
      [ 1, 0, 0, 2, 2,     23.0,   0.0,    0.0,    -10.0,  0.0,    0.0],
      [ 0,-1, 2, 1, 2,     12.0,   0.0,    0.0,     -5.0,  0.0,    0.0],
      [ 3,-1, 0, 0, 0,     29.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [ 2, 0, 0, 1, 0,    -25.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [ 1,-1, 2, 0, 0,     22.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 0, 2, 1, 0,    -18.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1, 0, 2, 0, 3,     15.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [ 3, 1, 0, 0, 0,    -23.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 3,-1, 2,-2, 2,     12.0,   0.0,    0.0,     -5.0,  0.0,    0.0],

   /* 411-420 */
      [ 2, 0, 2,-1, 1,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0],
      [ 1, 1, 2, 0, 0,    -19.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 0, 4,-1, 2,    -10.0,   0.0,    0.0,      4.0,  0.0,    0.0],
      [ 1, 2, 2, 0, 2,     21.0,   0.0,    0.0,     -9.0,  0.0,    0.0],
      [-2, 0, 0, 6, 0,     23.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [ 0,-1, 0, 4, 1,    -16.0,   0.0,    0.0,      8.0,  0.0,    0.0],
      [-2,-1, 2, 4, 1,    -19.0,   0.0,    0.0,      9.0,  0.0,    0.0],
      [ 0,-2, 2, 2, 1,    -22.0,   0.0,    0.0,     10.0,  0.0,    0.0],
      [ 0,-1, 2, 2, 0,     27.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [-1, 0, 2, 3, 1,     16.0,   0.0,    0.0,     -8.0,  0.0,    0.0],

   /* 421-430 */
      [-2, 1, 2, 4, 2,     19.0,   0.0,    0.0,     -8.0,  0.0,    0.0],
      [ 2, 0, 0, 2, 2,      9.0,   0.0,    0.0,     -4.0,  0.0,    0.0],
      [ 2,-2, 2, 0, 2,     -9.0,   0.0,    0.0,      4.0,  0.0,    0.0],
      [-1, 1, 2, 3, 2,     -9.0,   0.0,    0.0,      4.0,  0.0,    0.0],
      [ 3, 0, 2,-1, 2,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0],
      [ 4, 0, 2,-2, 1,     18.0,   0.0,    0.0,     -9.0,  0.0,    0.0],
      [-1, 0, 0, 6, 0,     16.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [-1,-2, 2, 4, 2,    -10.0,   0.0,    0.0,      4.0,  0.0,    0.0],
      [-3, 0, 2, 6, 2,    -23.0,   0.0,    0.0,      9.0,  0.0,    0.0],
      [-1, 0, 2, 4, 0,     16.0,   0.0,    0.0,     -1.0,  0.0,    0.0],

   /* 431-440 */
      [ 3, 0, 0, 2, 1,    -12.0,   0.0,    0.0,      6.0,  0.0,    0.0],
      [ 3,-1, 2, 0, 1,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0],
      [ 3, 0, 2, 0, 0,     30.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 1, 0, 4, 0, 2,     24.0,   0.0,    0.0,    -10.0,  0.0,    0.0],
      [ 5, 0, 2,-2, 2,     10.0,   0.0,    0.0,     -4.0,  0.0,    0.0],
      [ 0,-1, 2, 4, 1,    -16.0,   0.0,    0.0,      7.0,  0.0,    0.0],
      [ 2,-1, 2, 2, 1,    -16.0,   0.0,    0.0,      7.0,  0.0,    0.0],
      [ 0, 1, 2, 4, 2,     17.0,   0.0,    0.0,     -7.0,  0.0,    0.0],
      [ 1,-1, 2, 4, 2,    -24.0,   0.0,    0.0,     10.0,  0.0,    0.0],
      [ 3,-1, 2, 2, 2,    -12.0,   0.0,    0.0,      5.0,  0.0,    0.0],

   /* 441-450 */
      [ 3, 0, 2, 2, 1,    -24.0,   0.0,    0.0,     11.0,  0.0,    0.0],
      [ 5, 0, 2, 0, 2,    -23.0,   0.0,    0.0,      9.0,  0.0,    0.0],
      [ 0, 0, 2, 6, 2,    -13.0,   0.0,    0.0,      5.0,  0.0,    0.0],
      [ 4, 0, 2, 2, 2,    -15.0,   0.0,    0.0,      7.0,  0.0,    0.0],
      [ 0,-1, 1,-1, 1,      0.0,   0.0,-1988.0,      0.0,  0.0,-1679.0],
      [-1, 0, 1, 0, 3,      0.0,   0.0,  -63.0,      0.0,  0.0,  -27.0],
      [ 0,-2, 2,-2, 3,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1, 0,-1, 0, 1,      0.0,   0.0,    5.0,      0.0,  0.0,    4.0],
      [ 2,-2, 0,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [-1, 0, 1, 0, 2,      0.0,   0.0,  364.0,      0.0,  0.0,  176.0],

   /* 451-460 */
      [-1, 0, 1, 0, 1,      0.0,   0.0,-1044.0,      0.0,  0.0, -891.0],
      [-1,-1, 2,-1, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [-2, 2, 0, 2, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [-1, 0, 1, 0, 0,      0.0,   0.0,  330.0,      0.0,  0.0,    0.0],
      [-4, 1, 2, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [-3, 0, 2, 1, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [-2,-1, 2, 0, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [ 1, 0,-2, 1, 1,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 2,-1,-2, 0, 1,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [-4, 0, 2, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0],

   /* 461-470 */
      [-3, 1, 0, 3, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1, 0,-1, 2, 0,      0.0,   0.0,    5.0,      0.0,  0.0,    0.0],
      [ 0,-2, 0, 0, 2,      0.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [ 0,-2, 0, 0, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [-3, 0, 0, 3, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-2,-1, 0, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [-1, 0,-2, 3, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-4, 0, 0, 4, 0,    -12.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 2, 1,-2, 0, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [ 2,-1, 0,-2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],

   /* 471-480 */
      [ 0, 0, 1,-1, 0,     -5.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1, 2, 0, 1, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-2, 1, 2, 0, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [ 1, 1, 0,-1, 1,      7.0,   0.0,    0.0,     -4.0,  0.0,    0.0],
      [ 1, 0, 1,-2, 1,      0.0,   0.0,  -12.0,      0.0,  0.0,  -10.0],
      [ 0, 2, 0, 0, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 1,-1, 2,-3, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [-1, 1, 2,-1, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-2, 0, 4,-2, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [-2, 0, 4,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],

   /* 481-490 */
      [-2,-2, 0, 2, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [-2, 0,-2, 4, 0,      0.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1, 2, 2,-4, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [ 1, 1, 2,-4, 2,      7.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [-1, 2, 2,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 2, 0, 0,-3, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [-1, 2, 0, 0, 1,     -5.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [ 0, 0, 0,-2, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1,-1, 2,-2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-1, 1, 0, 0, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0],

   /* 491-500 */
      [ 0, 0, 0,-1, 2,     -8.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [-2, 1, 0, 1, 0,      9.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1,-2, 0,-2, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [ 1, 0,-2, 0, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-3, 1, 0, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1, 1,-2, 2, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1,-1, 0, 0, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [-3, 0, 0, 2, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-3,-1, 0, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 2, 0, 2,-6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0],

   /* 501-510 */
      [ 0, 1, 2,-4, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 2, 0, 0,-4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [-2, 1, 2,-2, 1,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 0,-1, 2,-4, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 0, 1, 0,-2, 2,      9.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [-1, 0, 0,-2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 2, 0,-2,-2, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [-4, 0, 2, 0, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-1,-1, 0,-1, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 0, 0,-2, 0, 2,      9.0,   0.0,    0.0,     -3.0,  0.0,    0.0],

   /* 511-520 */
      [-3, 0, 0, 1, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1, 0,-2, 1, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-2, 0,-2, 2, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 0, 0,-4, 2, 0,      8.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-2,-1,-2, 2, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1, 0, 2,-6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-1, 0, 2,-4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [ 1, 0, 0,-4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [ 2, 1, 2,-4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [ 2, 1, 2,-4, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0],

   /* 521-530 */
      [ 0, 1, 4,-4, 4,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 1, 4,-4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [-1,-1,-2, 4, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1,-3, 0, 2, 0,      9.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1, 0,-2, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-2,-1, 0, 3, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 0,-2, 3, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-2, 0, 0, 3, 1,     -5.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [ 0,-1, 0, 1, 0,    -13.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-3, 0, 2, 2, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0],

   /* 531-540 */
      [ 1, 1,-2, 2, 0,     10.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1, 1, 0, 2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [ 1,-2, 2,-2, 1,     10.0,   0.0,   13.0,      6.0,  0.0,   -5.0],
      [ 0, 0, 1, 0, 2,      0.0,   0.0,   30.0,      0.0,  0.0,   14.0],
      [ 0, 0, 1, 0, 1,      0.0,   0.0, -162.0,      0.0,  0.0, -138.0],
      [ 0, 0, 1, 0, 0,      0.0,   0.0,   75.0,      0.0,  0.0,    0.0],
      [-1, 2, 0, 2, 1,     -7.0,   0.0,    0.0,      4.0,  0.0,    0.0],
      [ 0, 0, 2, 0, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-2, 0, 2, 0, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 2, 0, 0,-1, 1,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0],

   /* 541-550 */
      [ 3, 0, 0,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [ 1, 0, 2,-2, 3,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1, 2, 0, 0, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 2, 0, 2,-3, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-1, 1, 4,-2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-2,-2, 0, 4, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0,-3, 0, 2, 0,      9.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 0,-2, 4, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1,-1, 0, 3, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-2, 0, 0, 4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],

   /* 551-560 */
      [-1, 0, 0, 3, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 2,-2, 0, 0, 0,      7.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1,-1, 0, 1, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1, 0, 0, 2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0,-2, 2, 0, 1,     -6.0,   0.0,   -3.0,      3.0,  0.0,    1.0],
      [-1, 0, 1, 2, 1,      0.0,   0.0,   -3.0,      0.0,  0.0,   -2.0],
      [-1, 1, 0, 3, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1,-1, 2, 1, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [ 0,-1, 2, 0, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-2, 1, 2, 2, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0],

   /* 561-570 */
      [ 2,-2, 2,-2, 2,     -1.0,   0.0,    3.0,      3.0,  0.0,   -1.0],
      [ 1, 1, 0, 1, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 1, 0, 1, 0, 1,      0.0,   0.0,  -13.0,      0.0,  0.0,  -11.0],
      [ 1, 0, 1, 0, 0,      3.0,   0.0,    6.0,      0.0,  0.0,    0.0],
      [ 0, 2, 0, 2, 0,     -7.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 2,-1, 2,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [ 0,-1, 4,-2, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [ 0, 0, 4,-2, 3,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 1, 4,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [ 4, 0, 2,-4, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0],

   /* 571-580 */
      [ 2, 2, 2,-2, 2,      8.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [ 2, 0, 4,-4, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-1,-2, 0, 4, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1,-3, 2, 2, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [-3, 0, 2, 4, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [-3, 0, 2,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-1,-1, 0,-2, 1,      8.0,   0.0,    0.0,     -4.0,  0.0,    0.0],
      [-3, 0, 0, 0, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [-3, 0,-2, 2, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 1, 0,-4, 1,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0],

   /* 581-590 */
      [-2, 1, 0,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-4, 0, 0, 0, 1,     -8.0,   0.0,    0.0,      4.0,  0.0,    0.0],
      [-1, 0, 0,-4, 1,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [-3, 0, 0,-2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 0, 0, 0, 3, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [-1, 1, 0, 4, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [ 1,-2, 2, 0, 1,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [ 0, 1, 0, 3, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-1, 0, 2, 2, 3,      6.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [ 0, 0, 2, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0],

   /* 591-600 */
      [-2, 0, 2, 2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-1, 1, 2, 2, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 3, 0, 0, 0, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 2, 1, 0, 1, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 2,-1, 2,-1, 2,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [ 0, 0, 2, 0, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 0, 0, 3, 0, 3,      0.0,   0.0,  -26.0,      0.0,  0.0,  -11.0],
      [ 0, 0, 3, 0, 2,      0.0,   0.0,  -10.0,      0.0,  0.0,   -5.0],
      [-1, 2, 2, 2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [-1, 0, 4, 0, 0,    -13.0,   0.0,    0.0,      0.0,  0.0,    0.0],

   /* 601-610 */
      [ 1, 2, 2, 0, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 3, 1, 2,-2, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 1, 1, 4,-2, 2,      7.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [-2,-1, 0, 6, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0,-2, 0, 4, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-2, 0, 0, 6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-2,-2, 2, 4, 2,     -6.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 0,-3, 2, 2, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 0, 0, 0, 4, 2,     -7.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [-1,-1, 2, 3, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0],

   /* 611-620 */
      [-2, 0, 2, 4, 0,     13.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 2,-1, 0, 2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 1, 0, 0, 3, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 1, 0, 4, 1,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 0, 1, 0, 4, 0,    -11.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1,-1, 2, 1, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 0, 0, 2, 2, 3,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1, 0, 2, 2, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [-1, 0, 2, 2, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-2, 0, 4, 2, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0],

   /* 621-630 */
      [ 2, 1, 0, 2, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 2, 1, 0, 2, 0,    -12.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 2,-1, 2, 0, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1, 0, 2, 1, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 1, 2, 2, 0,     -4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 2, 0, 2, 0, 3,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 3, 0, 2, 0, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [ 1, 0, 2, 0, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [ 1, 0, 3, 0, 3,      0.0,   0.0,   -5.0,      0.0,  0.0,   -2.0],
      [ 1, 1, 2, 1, 1,     -7.0,   0.0,    0.0,      4.0,  0.0,    0.0],

   /* 631-640 */
      [ 0, 2, 2, 2, 2,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [ 2, 1, 2, 0, 0,     -3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 2, 0, 4,-2, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [ 4, 1, 2,-2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [-1,-1, 0, 6, 0,      3.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [-3,-1, 2, 6, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [-1, 0, 0, 6, 1,     -5.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [-3, 0, 2, 6, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 1,-1, 0, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 1,-1, 0, 4, 0,     12.0,   0.0,    0.0,      0.0,  0.0,    0.0],

   /* 641-650 */
      [-2, 0, 2, 5, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [ 1,-2, 2, 2, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 3,-1, 0, 2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1,-1, 2, 2, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 0, 2, 3, 1,      5.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [-1, 1, 2, 4, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 0, 1, 2, 3, 2,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [-1, 0, 4, 2, 1,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 2, 0, 2, 1, 1,      6.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [ 5, 0, 0, 0, 0,      6.0,   0.0,    0.0,      0.0,  0.0,    0.0],

   /* 651-660 */
      [ 2, 1, 2, 1, 2,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [ 1, 0, 4, 0, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 3, 1, 2, 0, 1,      7.0,   0.0,    0.0,     -4.0,  0.0,    0.0],
      [ 3, 0, 4,-2, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [-2,-1, 2, 6, 2,     -5.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 0, 0, 0, 6, 0,      5.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0,-2, 2, 4, 2,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [-2, 0, 2, 6, 1,     -6.0,   0.0,    0.0,      3.0,  0.0,    0.0],
      [ 2, 0, 0, 4, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 2, 0, 0, 4, 0,     10.0,   0.0,    0.0,      0.0,  0.0,    0.0],

   /* 661-670 */
      [ 2,-2, 2, 2, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 0, 0, 2, 4, 0,      7.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 1, 0, 2, 3, 2,      7.0,   0.0,    0.0,     -3.0,  0.0,    0.0],
      [ 4, 0, 0, 2, 0,      4.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 2, 0, 2, 2, 0,     11.0,   0.0,    0.0,      0.0,  0.0,    0.0],
      [ 0, 0, 4, 2, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 4,-1, 2, 0, 2,     -6.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 3, 0, 2, 1, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 2, 1, 2, 2, 1,      3.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 4, 1, 2, 0, 2,      5.0,   0.0,    0.0,     -2.0,  0.0,    0.0],

   /* 671-678 */
      [-1,-1, 2, 6, 2,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [-1, 0, 2, 6, 1,     -4.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 1,-1, 2, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0],
      [ 1, 1, 2, 4, 2,      4.0,   0.0,    0.0,     -2.0,  0.0,    0.0],
      [ 3, 1, 2, 2, 2,      3.0,   0.0,    0.0,     -1.0,  0.0,    0.0],
      [ 5, 0, 2, 0, 1,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [ 2,-1, 2, 4, 2,     -3.0,   0.0,    0.0,      1.0,  0.0,    0.0],
      [ 2, 0, 2, 4, 1,     -3.0,   0.0,    0.0,      2.0,  0.0,    0.0]
   ];

/* Number of terms in the luni-solar nutation model */
   var NLS = xls.length;

/* ------------------------ */
/* Planetary nutation model */
/* ------------------------ */

/* The units for the sine and cosine coefficients are */
/* 0.1 microarcsecond                                 */

   var xpl = [

   /* 1-10 */
      [ 0, 0, 0, 0, 0,  0,  8,-16, 4, 5, 0, 0, 0, 1440,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0, -8, 16,-4,-5, 0, 0, 2,   56,-117,  -42, -40],
      [ 0, 0, 0, 0, 0,  0,  8,-16, 4, 5, 0, 0, 2,  125, -43,    0, -54],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 0,-1, 2, 2,    0,   5,    0,   0],
      [ 0, 0, 0, 0, 0,  0, -4,  8,-1,-5, 0, 0, 2,    3,  -7,   -3,   0],
      [ 0, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 1,    3,   0,    0,  -2],
      [ 0, 1,-1, 1, 0,  0,  3, -8, 3, 0, 0, 0, 0, -114,   0,    0,  61],
      [-1, 0, 0, 0, 0, 10, -3,  0, 0, 0, 0, 0, 0, -219,  89,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  0,  0,-2, 6,-3, 0, 2,   -3,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0, -462,1604,    0,   0],

   /* 11-20 */
      [ 0, 1,-1, 1, 0,  0, -5,  8,-3, 0, 0, 0, 0,   99,   0,    0, -53],
      [ 0, 0, 0, 0, 0,  0, -4,  8,-3, 0, 0, 0, 1,   -3,   0,    0,   2],
      [ 0, 0, 0, 0, 0,  0,  4, -8, 1, 5, 0, 0, 2,    0,   6,    2,   0],
      [ 0, 0, 0, 0, 0, -5,  6,  4, 0, 0, 0, 0, 2,    3,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 2,-5, 0, 0, 2,  -12,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 2,-5, 0, 0, 1,   14,-218,  117,   8],
      [ 0, 1,-1, 1, 0,  0, -1,  0, 2,-5, 0, 0, 0,   31,-481, -257, -17],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 2,-5, 0, 0, 0, -491, 128,    0,   0],
      [ 0, 1,-1, 1, 0,  0, -1,  0,-2, 5, 0, 0, 0,-3084,5123, 2735,1647],
      [ 0, 0, 0, 0, 0,  0,  0,  0,-2, 5, 0, 0, 1,-1444,2409,-1286,-771],

   /* 21-30 */
      [ 0, 0, 0, 0, 0,  0,  0,  0,-2, 5, 0, 0, 2,   11, -24,  -11,  -9],
      [ 2,-1,-1, 0, 0,  0,  3, -7, 0, 0, 0, 0, 0,   26,  -9,    0,   0],
      [ 1, 0,-2, 0, 0, 19,-21,  3, 0, 0, 0, 0, 0,  103, -60,    0,   0],
      [ 0, 1,-1, 1, 0,  2, -4,  0,-3, 0, 0, 0, 0,    0, -13,   -7,   0],
      [ 1, 0,-1, 1, 0,  0, -1,  0, 2, 0, 0, 0, 0,  -26, -29,  -16,  14],
      [ 0, 1,-1, 1, 0,  0, -1,  0,-4,10, 0, 0, 0,    9, -27,  -14,  -5],
      [-2, 0, 2, 1, 0,  0,  2,  0, 0,-5, 0, 0, 0,   12,   0,    0,  -6],
      [ 0, 0, 0, 0, 0,  3, -7,  4, 0, 0, 0, 0, 0,   -7,   0,    0,   0],
      [ 0,-1, 1, 0, 0,  0,  1,  0, 1,-1, 0, 0, 0,    0,  24,    0,   0],
      [-2, 0, 2, 1, 0,  0,  2,  0,-2, 0, 0, 0, 0,  284,   0,    0,-151],

   /* 31-40 */
      [-1, 0, 0, 0, 0, 18,-16,  0, 0, 0, 0, 0, 0,  226, 101,    0,   0],
      [-2, 1, 1, 2, 0,  0,  1,  0,-2, 0, 0, 0, 0,    0,  -8,   -2,   0],
      [-1, 1,-1, 1, 0, 18,-17,  0, 0, 0, 0, 0, 0,    0,  -6,   -3,   0],
      [-1, 0, 1, 1, 0,  0,  2, -2, 0, 0, 0, 0, 0,    5,   0,    0,  -3],
      [ 0, 0, 0, 0, 0, -8, 13,  0, 0, 0, 0, 0, 2,  -41, 175,   76,  17],
      [ 0, 2,-2, 2, 0, -8, 11,  0, 0, 0, 0, 0, 0,    0,  15,    6,   0],
      [ 0, 0, 0, 0, 0, -8, 13,  0, 0, 0, 0, 0, 1,  425, 212, -133, 269],
      [ 0, 1,-1, 1, 0, -8, 12,  0, 0, 0, 0, 0, 0, 1200, 598,  319,-641],
      [ 0, 0, 0, 0, 0,  8,-13,  0, 0, 0, 0, 0, 0,  235, 334,    0,   0],
      [ 0, 1,-1, 1, 0,  8,-14,  0, 0, 0, 0, 0, 0,   11, -12,   -7,  -6],

   /* 41-50 */
      [ 0, 0, 0, 0, 0,  8,-13,  0, 0, 0, 0, 0, 1,    5,  -6,    3,   3],
      [-2, 0, 2, 1, 0,  0,  2,  0,-4, 5, 0, 0, 0,   -5,   0,    0,   3],
      [-2, 0, 2, 2, 0,  3, -3,  0, 0, 0, 0, 0, 0,    6,   0,    0,  -3],
      [-2, 0, 2, 0, 0,  0,  2,  0,-3, 1, 0, 0, 0,   15,   0,    0,   0],
      [ 0, 0, 0, 1, 0,  3, -5,  0, 2, 0, 0, 0, 0,   13,   0,    0,  -7],
      [-2, 0, 2, 0, 0,  0,  2,  0,-4, 3, 0, 0, 0,   -6,  -9,    0,   0],
      [ 0,-1, 1, 0, 0,  0,  0,  2, 0, 0, 0, 0, 0,  266, -78,    0,   0],
      [ 0, 0, 0, 1, 0,  0, -1,  2, 0, 0, 0, 0, 0, -460,-435, -232, 246],
      [ 0, 1,-1, 2, 0,  0, -2,  2, 0, 0, 0, 0, 0,    0,  15,    7,   0],
      [-1, 1, 0, 1, 0,  3, -5,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   2],

   /* 51-60 */
      [-1, 0, 1, 0, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0, 131,    0,   0],
      [-2, 0, 2, 0, 0,  0,  2,  0,-2,-2, 0, 0, 0,    4,   0,    0,   0],
      [-2, 2, 0, 2, 0,  0, -5,  9, 0, 0, 0, 0, 0,    0,   3,    0,   0],
      [ 0, 1,-1, 1, 0,  0, -1,  0, 0, 0,-1, 0, 0,    0,   4,    2,   0],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 1, 0, 0,    0,   3,    0,   0],
      [ 0, 1,-1, 1, 0,  0, -1,  0, 0, 0, 0, 2, 0,  -17, -19,  -10,   9],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 2, 1,   -9, -11,    6,  -5],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 2, 2,   -6,   0,    0,   3],
      [-1, 0, 1, 0, 0,  0,  3, -4, 0, 0, 0, 0, 0,  -16,   8,    0,   0],
      [ 0,-1, 1, 0, 0,  0,  1,  0, 0, 2, 0, 0, 0,    0,   3,    0,   0],

   /* 61-70 */
      [ 0, 1,-1, 2, 0,  0, -1,  0, 0, 2, 0, 0, 0,   11,  24,   11,  -5],
      [ 0, 0, 0, 1, 0,  0, -9, 17, 0, 0, 0, 0, 0,   -3,  -4,   -2,   1],
      [ 0, 0, 0, 2, 0, -3,  5,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1],
      [ 0, 1,-1, 1, 0,  0, -1,  0,-1, 2, 0, 0, 0,    0,  -8,   -4,   0],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 1,-2, 0, 0, 0,    0,   3,    0,   0],
      [ 1, 0,-2, 0, 0, 17,-16,  0,-2, 0, 0, 0, 0,    0,   5,    0,   0],
      [ 0, 1,-1, 1, 0,  0, -1,  0, 1,-3, 0, 0, 0,    0,   3,    2,   0],
      [-2, 0, 2, 1, 0,  0,  5, -6, 0, 0, 0, 0, 0,   -6,   4,    2,   3],
      [ 0,-2, 2, 0, 0,  0,  9,-13, 0, 0, 0, 0, 0,   -3,  -5,    0,   0],
      [ 0, 1,-1, 2, 0,  0, -1,  0, 0, 1, 0, 0, 0,   -5,   0,    0,   2],

   /* 71-80 */
      [ 0, 0, 0, 1, 0,  0,  0,  0, 0, 1, 0, 0, 0,    4,  24,   13,  -2],
      [ 0,-1, 1, 0, 0,  0,  1,  0, 0, 1, 0, 0, 0,  -42,  20,    0,   0],
      [ 0,-2, 2, 0, 0,  5, -6,  0, 0, 0, 0, 0, 0,  -10, 233,    0,   0],
      [ 0,-1, 1, 1, 0,  5, -7,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1],
      [-2, 0, 2, 0, 0,  6, -8,  0, 0, 0, 0, 0, 0,   78, -18,    0,   0],
      [ 2, 1,-3, 1, 0, -6,  7,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0],
      [ 0, 0, 0, 2, 0,  0,  0,  0, 1, 0, 0, 0, 0,    0,  -3,   -1,   0],
      [ 0,-1, 1, 1, 0,  0,  1,  0, 1, 0, 0, 0, 0,    0,  -4,   -2,   1],
      [ 0, 1,-1, 1, 0,  0, -1,  0, 0, 0, 2, 0, 0,    0,  -8,   -4,  -1],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 2, 0, 1,    0,  -5,    3,   0],

   /* 81-90 */
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 2, 0, 2,   -7,   0,    0,   3],
      [ 0, 0, 0, 0, 0,  0, -8, 15, 0, 0, 0, 0, 2,  -14,   8,    3,   6],
      [ 0, 0, 0, 0, 0,  0, -8, 15, 0, 0, 0, 0, 1,    0,   8,   -4,   0],
      [ 0, 1,-1, 1, 0,  0, -9, 15, 0, 0, 0, 0, 0,    0,  19,   10,   0],
      [ 0, 0, 0, 0, 0,  0,  8,-15, 0, 0, 0, 0, 0,   45, -22,    0,   0],
      [ 1,-1,-1, 0, 0,  0,  8,-15, 0, 0, 0, 0, 0,   -3,   0,    0,   0],
      [ 2, 0,-2, 0, 0,  2, -5,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0],
      [-2, 0, 2, 0, 0,  0,  2,  0,-5, 5, 0, 0, 0,    0,   3,    0,   0],
      [ 2, 0,-2, 1, 0,  0, -6,  8, 0, 0, 0, 0, 0,    3,   5,    3,  -2],
      [ 2, 0,-2, 1, 0,  0, -2,  0, 3, 0, 0, 0, 0,   89, -16,   -9, -48],

   /* 91-100 */
      [-2, 1, 1, 0, 0,  0,  1,  0,-3, 0, 0, 0, 0,    0,   3,    0,   0],
      [-2, 1, 1, 1, 0,  0,  1,  0,-3, 0, 0, 0, 0,   -3,   7,    4,   2],
      [-2, 0, 2, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0, -349, -62,    0,   0],
      [-2, 0, 2, 0, 0,  0,  6, -8, 0, 0, 0, 0, 0,  -15,  22,    0,   0],
      [-2, 0, 2, 0, 0,  0,  2,  0,-1,-5, 0, 0, 0,   -3,   0,    0,   0],
      [-1, 0, 1, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,  -53,   0,    0,   0],
      [-1, 1, 1, 1, 0,-20, 20,  0, 0, 0, 0, 0, 0,    5,   0,    0,  -3],
      [ 1, 0,-2, 0, 0, 20,-21,  0, 0, 0, 0, 0, 0,    0,  -8,    0,   0],
      [ 0, 0, 0, 1, 0,  0,  8,-15, 0, 0, 0, 0, 0,   15,  -7,   -4,  -8],
      [ 0, 2,-2, 1, 0,  0,-10, 15, 0, 0, 0, 0, 0,   -3,   0,    0,   1],

   /* 101-110 */
      [ 0,-1, 1, 0, 0,  0,  1,  0, 1, 0, 0, 0, 0,  -21, -78,    0,   0],
      [ 0, 0, 0, 1, 0,  0,  0,  0, 1, 0, 0, 0, 0,   20, -70,  -37, -11],
      [ 0, 1,-1, 2, 0,  0, -1,  0, 1, 0, 0, 0, 0,    0,   6,    3,   0],
      [ 0, 1,-1, 1, 0,  0, -1,  0,-2, 4, 0, 0, 0,    5,   3,    2,  -2],
      [ 2, 0,-2, 1, 0, -6,  8,  0, 0, 0, 0, 0, 0,  -17,  -4,   -2,   9],
      [ 0,-2, 2, 1, 0,  5, -6,  0, 0, 0, 0, 0, 0,    0,   6,    3,   0],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0,-1, 0, 0, 1,   32,  15,   -8,  17],
      [ 0, 1,-1, 1, 0,  0, -1,  0, 0,-1, 0, 0, 0,  174,  84,   45, -93],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 1, 0, 0, 0,   11,  56,    0,   0],
      [ 0, 1,-1, 1, 0,  0, -1,  0, 0, 1, 0, 0, 0,  -66, -12,   -6,  35],

   /* 111-120 */
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 1, 0, 0, 1,   47,   8,    4, -25],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 1, 0, 0, 2,    0,   8,    4,   0],
      [ 0, 2,-2, 1, 0,  0, -9, 13, 0, 0, 0, 0, 0,   10, -22,  -12,  -5],
      [ 0, 0, 0, 1, 0,  0,  7,-13, 0, 0, 0, 0, 0,   -3,   0,    0,   2],
      [-2, 0, 2, 0, 0,  0,  5, -6, 0, 0, 0, 0, 0,  -24,  12,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  9,-17, 0, 0, 0, 0, 0,    5,  -6,    0,   0],
      [ 0, 0, 0, 0, 0,  0, -9, 17, 0, 0, 0, 0, 2,    3,   0,    0,  -2],
      [ 1, 0,-1, 1, 0,  0, -3,  4, 0, 0, 0, 0, 0,    4,   3,    1,  -2],
      [ 1, 0,-1, 1, 0, -3,  4,  0, 0, 0, 0, 0, 0,    0,  29,   15,   0],
      [ 0, 0, 0, 2, 0,  0, -1,  2, 0, 0, 0, 0, 0,   -5,  -4,   -2,   2],

   /* 121-130 */
      [ 0,-1, 1, 1, 0,  0,  0,  2, 0, 0, 0, 0, 0,    8,  -3,   -1,  -5],
      [ 0,-2, 2, 0, 1,  0, -2,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0],
      [ 0, 0, 0, 0, 0,  3, -5,  0, 2, 0, 0, 0, 0,   10,   0,    0,   0],
      [-2, 0, 2, 1, 0,  0,  2,  0,-3, 1, 0, 0, 0,    3,   0,    0,  -2],
      [-2, 0, 2, 1, 0,  3, -3,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   3],
      [ 0, 0, 0, 1, 0,  8,-13,  0, 0, 0, 0, 0, 0,   46,  66,   35, -25],
      [ 0,-1, 1, 0, 0,  8,-12,  0, 0, 0, 0, 0, 0,  -14,   7,    0,   0],
      [ 0, 2,-2, 1, 0, -8, 11,  0, 0, 0, 0, 0, 0,    0,   3,    2,   0],
      [-1, 0, 1, 0, 0,  0,  2, -2, 0, 0, 0, 0, 0,   -5,   0,    0,   0],
      [-1, 0, 0, 1, 0, 18,-16,  0, 0, 0, 0, 0, 0,  -68, -34,  -18,  36],

   /* 131-140 */
      [ 0, 1,-1, 1, 0,  0, -1,  0,-1, 1, 0, 0, 0,    0,  14,    7,   0],
      [ 0, 0, 0, 1, 0,  3, -7,  4, 0, 0, 0, 0, 0,   10,  -6,   -3,  -5],
      [-2, 1, 1, 1, 0,  0, -3,  7, 0, 0, 0, 0, 0,   -5,  -4,   -2,   3],
      [ 0, 1,-1, 2, 0,  0, -1,  0,-2, 5, 0, 0, 0,   -3,   5,    2,   1],
      [ 0, 0, 0, 1, 0,  0,  0,  0,-2, 5, 0, 0, 0,   76,  17,    9, -41],
      [ 0, 0, 0, 1, 0,  0, -4,  8,-3, 0, 0, 0, 0,   84, 298,  159, -45],
      [ 1, 0, 0, 1, 0,-10,  3,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1],
      [ 0, 2,-2, 1, 0,  0, -2,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   2],
      [-1, 0, 0, 1, 0, 10, -3,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1],
      [ 0, 0, 0, 1, 0,  0,  4, -8, 3, 0, 0, 0, 0,  -82, 292,  156,  44],

   /* 141-150 */
      [ 0, 0, 0, 1, 0,  0,  0,  0, 2,-5, 0, 0, 0,  -73,  17,    9,  39],
      [ 0,-1, 1, 0, 0,  0,  1,  0, 2,-5, 0, 0, 0,   -9, -16,    0,   0],
      [ 2,-1,-1, 1, 0,  0,  3, -7, 0, 0, 0, 0, 0,    3,   0,   -1,  -2],
      [-2, 0, 2, 0, 0,  0,  2,  0, 0,-5, 0, 0, 0,   -3,   0,    0,   0],
      [ 0, 0, 0, 1, 0, -3,  7, -4, 0, 0, 0, 0, 0,   -9,  -5,   -3,   5],
      [-2, 0, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0, -439,   0,    0,   0],
      [ 1, 0, 0, 1, 0,-18, 16,  0, 0, 0, 0, 0, 0,   57, -28,  -15, -30],
      [-2, 1, 1, 1, 0,  0,  1,  0,-2, 0, 0, 0, 0,    0,  -6,   -3,   0],
      [ 0, 1,-1, 2, 0, -8, 12,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   2],
      [ 0, 0, 0, 1, 0, -8, 13,  0, 0, 0, 0, 0, 0,  -40,  57,   30,  21],

   /* 151-160 */
      [ 0, 0, 0, 0, 0,  0,  1, -2, 0, 0, 0, 0, 1,   23,   7,    3, -13],
      [ 0, 1,-1, 1, 0,  0,  0, -2, 0, 0, 0, 0, 0,  273,  80,   43,-146],
      [ 0, 0, 0, 0, 0,  0,  1, -2, 0, 0, 0, 0, 0, -449, 430,    0,   0],
      [ 0, 1,-1, 1, 0,  0, -2,  2, 0, 0, 0, 0, 0,   -8, -47,  -25,   4],
      [ 0, 0, 0, 0, 0,  0, -1,  2, 0, 0, 0, 0, 1,    6,  47,   25,  -3],
      [-1, 0, 1, 1, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0,  23,   13,   0],
      [-1, 0, 1, 1, 0,  0,  3, -4, 0, 0, 0, 0, 0,   -3,   0,    0,   2],
      [ 0, 1,-1, 1, 0,  0, -1,  0, 0,-2, 0, 0, 0,    3,  -4,   -2,  -2],
      [ 0, 1,-1, 1, 0,  0, -1,  0, 0, 2, 0, 0, 0,  -48,-110,  -59,  26],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 2, 0, 0, 1,   51, 114,   61, -27],

   /* 161-170 */
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 2, 0, 0, 2, -133,   0,    0,  57],
      [ 0, 1,-1, 0, 0,  3, -6,  0, 0, 0, 0, 0, 0,    0,   4,    0,   0],
      [ 0, 0, 0, 1, 0, -3,  5,  0, 0, 0, 0, 0, 0,  -21,  -6,   -3,  11],
      [ 0, 1,-1, 2, 0, -3,  4,  0, 0, 0, 0, 0, 0,    0,  -3,   -1,   0],
      [ 0, 0, 0, 1, 0,  0, -2,  4, 0, 0, 0, 0, 0,  -11, -21,  -11,   6],
      [ 0, 2,-2, 1, 0, -5,  6,  0, 0, 0, 0, 0, 0,  -18,-436, -233,   9],
      [ 0,-1, 1, 0, 0,  5, -7,  0, 0, 0, 0, 0, 0,   35,  -7,    0,   0],
      [ 0, 0, 0, 1, 0,  5, -8,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0],
      [-2, 0, 2, 1, 0,  6, -8,  0, 0, 0, 0, 0, 0,   11,  -3,   -1,  -6],
      [ 0, 0, 0, 1, 0,  0, -8, 15, 0, 0, 0, 0, 0,   -5,  -3,   -1,   3],

   /* 171-180 */
      [-2, 0, 2, 1, 0,  0,  2,  0,-3, 0, 0, 0, 0,  -53,  -9,   -5,  28],
      [-2, 0, 2, 1, 0,  0,  6, -8, 0, 0, 0, 0, 0,    0,   3,    2,   1],
      [ 1, 0,-1, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,    4,   0,    0,  -2],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 3,-5, 0, 0, 0,    0,  -4,    0,   0],
      [ 0, 1,-1, 1, 0,  0, -1,  0,-1, 0, 0, 0, 0,  -50, 194,  103,  27],
      [ 0, 0, 0, 0, 0,  0,  0,  0,-1, 0, 0, 0, 1,  -13,  52,   28,   7],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 0,  -91, 248,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 1,    6,  49,   26,  -3],
      [ 0, 1,-1, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,   -6, -47,  -25,   3],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 1,    0,   5,    3,   0],

   /* 181-190 */
      [ 0, 0, 0, 0, 0,  0,  0,  0, 1, 0, 0, 0, 2,   52,  23,   10, -23],
      [ 0, 1,-1, 2, 0,  0, -1,  0, 0,-1, 0, 0, 0,   -3,   0,    0,   1],
      [ 0, 0, 0, 1, 0,  0,  0,  0, 0,-1, 0, 0, 0,    0,   5,    3,   0],
      [ 0,-1, 1, 0, 0,  0,  1,  0, 0,-1, 0, 0, 0,   -4,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0, -7, 13, 0, 0, 0, 0, 2,   -4,   8,    3,   2],
      [ 0, 0, 0, 0, 0,  0,  7,-13, 0, 0, 0, 0, 0,   10,   0,    0,   0],
      [ 2, 0,-2, 1, 0,  0, -5,  6, 0, 0, 0, 0, 0,    3,   0,    0,  -2],
      [ 0, 2,-2, 1, 0,  0, -8, 11, 0, 0, 0, 0, 0,    0,   8,    4,   0],
      [ 0, 2,-2, 1,-1,  0,  2,  0, 0, 0, 0, 0, 0,    0,   8,    4,   1],
      [-2, 0, 2, 0, 0,  0,  4, -4, 0, 0, 0, 0, 0,   -4,   0,    0,   0],

   /* 191-200 */
      [ 0, 0, 0, 0, 0,  0,  0,  0, 2,-2, 0, 0, 0,   -4,   0,    0,   0],
      [ 0, 1,-1, 1, 0,  0, -1,  0, 0, 3, 0, 0, 0,   -8,   4,    2,   4],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 3, 0, 0, 1,    8,  -4,   -2,  -4],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 3, 0, 0, 2,    0,  15,    7,   0],
      [-2, 0, 2, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0, -138,   0,    0,   0],
      [ 0, 0, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,  -7,   -3,   0],
      [ 0, 0, 0, 2, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -7,   -3,   0],
      [ 2, 0,-2, 1, 0,  0, -2,  0, 2, 0, 0, 0, 0,   54,   0,    0, -29],
      [ 0, 1,-1, 2, 0,  0, -1,  0, 2, 0, 0, 0, 0,    0,  10,    4,   0],
      [ 0, 1,-1, 2, 0,  0,  0, -2, 0, 0, 0, 0, 0,   -7,   0,    0,   3],

   /* 201-210 */
      [ 0, 0, 0, 1, 0,  0,  1, -2, 0, 0, 0, 0, 0,  -37,  35,   19,  20],
      [ 0,-1, 1, 0, 0,  0,  2, -2, 0, 0, 0, 0, 0,    0,   4,    0,   0],
      [ 0,-1, 1, 0, 0,  0,  1,  0, 0,-2, 0, 0, 0,   -4,   9,    0,   0],
      [ 0, 2,-2, 1, 0,  0, -2,  0, 0, 2, 0, 0, 0,    8,   0,    0,  -4],
      [ 0, 1,-1, 1, 0,  3, -6,  0, 0, 0, 0, 0, 0,   -9, -14,   -8,   5],
      [ 0, 0, 0, 0, 0,  3, -5,  0, 0, 0, 0, 0, 1,   -3,  -9,   -5,   3],
      [ 0, 0, 0, 0, 0,  3, -5,  0, 0, 0, 0, 0, 0, -145,  47,    0,   0],
      [ 0, 1,-1, 1, 0, -3,  4,  0, 0, 0, 0, 0, 0,  -10,  40,   21,   5],
      [ 0, 0, 0, 0, 0, -3,  5,  0, 0, 0, 0, 0, 1,   11, -49,  -26,  -7],
      [ 0, 0, 0, 0, 0, -3,  5,  0, 0, 0, 0, 0, 2,-2150,   0,    0, 932],

   /* 211-220 */
      [ 0, 2,-2, 2, 0, -3,  3,  0, 0, 0, 0, 0, 0,  -12,   0,    0,   5],
      [ 0, 0, 0, 0, 0, -3,  5,  0, 0, 0, 0, 0, 2,   85,   0,    0, -37],
      [ 0, 0, 0, 0, 0,  0,  2, -4, 0, 0, 0, 0, 1,    4,   0,    0,  -2],
      [ 0, 1,-1, 1, 0,  0,  1, -4, 0, 0, 0, 0, 0,    3,   0,    0,  -2],
      [ 0, 0, 0, 0, 0,  0,  2, -4, 0, 0, 0, 0, 0,  -86, 153,    0,   0],
      [ 0, 0, 0, 0, 0,  0, -2,  4, 0, 0, 0, 0, 1,   -6,   9,    5,   3],
      [ 0, 1,-1, 1, 0,  0, -3,  4, 0, 0, 0, 0, 0,    9, -13,   -7,  -5],
      [ 0, 0, 0, 0, 0,  0, -2,  4, 0, 0, 0, 0, 1,   -8,  12,    6,   4],
      [ 0, 0, 0, 0, 0,  0, -2,  4, 0, 0, 0, 0, 2,  -51,   0,    0,  22],
      [ 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 2,  -11,-268, -116,   5],

   /* 221-230 */
      [ 0, 2,-2, 2, 0, -5,  6,  0, 0, 0, 0, 0, 0,    0,  12,    5,   0],
      [ 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 2,    0,   7,    3,   0],
      [ 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 1,   31,   6,    3, -17],
      [ 0, 1,-1, 1, 0, -5,  7,  0, 0, 0, 0, 0, 0,  140,  27,   14, -75],
      [ 0, 0, 0, 0, 0, -5,  8,  0, 0, 0, 0, 0, 1,   57,  11,    6, -30],
      [ 0, 0, 0, 0, 0,  5, -8,  0, 0, 0, 0, 0, 0,  -14, -39,    0,   0],
      [ 0, 1,-1, 2, 0,  0, -1,  0,-1, 0, 0, 0, 0,    0,  -6,   -2,   0],
      [ 0, 0, 0, 1, 0,  0,  0,  0,-1, 0, 0, 0, 0,    4,  15,    8,  -2],
      [ 0,-1, 1, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    0,   4,    0,   0],
      [ 0, 2,-2, 1, 0,  0, -2,  0, 1, 0, 0, 0, 0,   -3,   0,    0,   1],

   /* 231-240 */
      [ 0, 0, 0, 0, 0,  0, -6, 11, 0, 0, 0, 0, 2,    0,  11,    5,   0],
      [ 0, 0, 0, 0, 0,  0,  6,-11, 0, 0, 0, 0, 0,    9,   6,    0,   0],
      [ 0, 0, 0, 0,-1,  0,  4,  0, 0, 0, 0, 0, 2,   -4,  10,    4,   2],
      [ 0, 0, 0, 0, 1,  0, -4,  0, 0, 0, 0, 0, 0,    5,   3,    0,   0],
      [ 2, 0,-2, 1, 0, -3,  3,  0, 0, 0, 0, 0, 0,   16,   0,    0,  -9],
      [-2, 0, 2, 0, 0,  0,  2,  0, 0,-2, 0, 0, 0,   -3,   0,    0,   0],
      [ 0, 2,-2, 1, 0,  0, -7,  9, 0, 0, 0, 0, 0,    0,   3,    2,  -1],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 4,-5, 0, 0, 2,    7,   0,    0,  -3],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 0,  -25,  22,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 1,   42, 223,  119, -22],

   /* 241-250 */
      [ 0, 1,-1, 1, 0,  0, -1,  0, 2, 0, 0, 0, 0,  -27,-143,  -77,  14],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 1,    9,  49,   26,  -5],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 2, 0, 0, 0, 2,-1166,   0,    0, 505],
      [ 0, 2,-2, 2, 0,  0, -2,  0, 2, 0, 0, 0, 0,   -5,   0,    0,   2],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 5, 0, 0, 2,   -6,   0,    0,   3],
      [ 0, 0, 0, 1, 0,  3, -5,  0, 0, 0, 0, 0, 0,   -8,   0,    1,   4],
      [ 0,-1, 1, 0, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0],
      [ 0, 2,-2, 1, 0, -3,  3,  0, 0, 0, 0, 0, 0,  117,   0,    0, -63],
      [ 0, 0, 0, 1, 0,  0,  2, -4, 0, 0, 0, 0, 0,   -4,   8,    4,   2],
      [ 0, 2,-2, 1, 0,  0, -4,  4, 0, 0, 0, 0, 0,    3,   0,    0,  -2],

   /* 251-260 */
      [ 0, 1,-1, 2, 0, -5,  7,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   2],
      [ 0, 0, 0, 0, 0,  0,  3, -6, 0, 0, 0, 0, 0,    0,  31,    0,   0],
      [ 0, 0, 0, 0, 0,  0, -3,  6, 0, 0, 0, 0, 1,   -5,   0,    1,   3],
      [ 0, 1,-1, 1, 0,  0, -4,  6, 0, 0, 0, 0, 0,    4,   0,    0,  -2],
      [ 0, 0, 0, 0, 0,  0, -3,  6, 0, 0, 0, 0, 1,   -4,   0,    0,   2],
      [ 0, 0, 0, 0, 0,  0, -3,  6, 0, 0, 0, 0, 2,  -24, -13,   -6,  10],
      [ 0,-1, 1, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    3,   0,    0,   0],
      [ 0, 0, 0, 1, 0,  2, -3,  0, 0, 0, 0, 0, 0,    0, -32,  -17,   0],
      [ 0, 0, 0, 0, 0,  0, -5,  9, 0, 0, 0, 0, 2,    8,  12,    5,  -3],
      [ 0, 0, 0, 0, 0,  0, -5,  9, 0, 0, 0, 0, 1,    3,   0,    0,  -1],

   /* 261-270 */
      [ 0, 0, 0, 0, 0,  0,  5, -9, 0, 0, 0, 0, 0,    7,  13,    0,   0],
      [ 0,-1, 1, 0, 0,  0,  1,  0,-2, 0, 0, 0, 0,   -3,  16,    0,   0],
      [ 0, 2,-2, 1, 0,  0, -2,  0, 2, 0, 0, 0, 0,   50,   0,    0, -27],
      [-2, 1, 1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -5,   -3,   0],
      [ 0,-2, 2, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,   13,   0,    0,   0],
      [ 0, 0, 0, 0, 0, -6, 10,  0, 0, 0, 0, 0, 1,    0,   5,    3,   1],
      [ 0, 0, 0, 0, 0, -6, 10,  0, 0, 0, 0, 0, 2,   24,   5,    2, -11],
      [ 0, 0, 0, 0, 0, -2,  3,  0, 0, 0, 0, 0, 2,    5, -11,   -5,  -2],
      [ 0, 0, 0, 0, 0, -2,  3,  0, 0, 0, 0, 0, 1,   30,  -3,   -2, -16],
      [ 0, 1,-1, 1, 0, -2,  2,  0, 0, 0, 0, 0, 0,   18,   0,    0,  -9],

   /* 271-280 */
      [ 0, 0, 0, 0, 0,  2, -3,  0, 0, 0, 0, 0, 0,    8, 614,    0,   0],
      [ 0, 0, 0, 0, 0,  2, -3,  0, 0, 0, 0, 0, 1,    3,  -3,   -1,  -2],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 3, 0, 0, 0, 1,    6,  17,    9,  -3],
      [ 0, 1,-1, 1, 0,  0, -1,  0, 3, 0, 0, 0, 0,   -3,  -9,   -5,   2],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 3, 0, 0, 0, 1,    0,   6,    3,  -1],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 3, 0, 0, 0, 2, -127,  21,    9,  55],
      [ 0, 0, 0, 0, 0,  0,  4, -8, 0, 0, 0, 0, 0,    3,   5,    0,   0],
      [ 0, 0, 0, 0, 0,  0, -4,  8, 0, 0, 0, 0, 2,   -6, -10,   -4,   3],
      [ 0,-2, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,    5,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0, -4,  7, 0, 0, 0, 0, 2,   16,   9,    4,  -7],

   /* 281-290 */
      [ 0, 0, 0, 0, 0,  0, -4,  7, 0, 0, 0, 0, 1,    3,   0,    0,  -2],
      [ 0, 0, 0, 0, 0,  0,  4, -7, 0, 0, 0, 0, 0,    0,  22,    0,   0],
      [ 0, 0, 0, 1, 0, -2,  3,  0, 0, 0, 0, 0, 0,    0,  19,   10,   0],
      [ 0, 2,-2, 1, 0,  0, -2,  0, 3, 0, 0, 0, 0,    7,   0,    0,  -4],
      [ 0, 0, 0, 0, 0,  0, -5, 10, 0, 0, 0, 0, 2,    0,  -5,   -2,   0],
      [ 0, 0, 0, 1, 0, -1,  2,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0],
      [ 0, 0, 0, 0, 0,  0,  0,  0, 4, 0, 0, 0, 2,   -9,   3,    1,   4],
      [ 0, 0, 0, 0, 0,  0, -3,  5, 0, 0, 0, 0, 2,   17,   0,    0,  -7],
      [ 0, 0, 0, 0, 0,  0, -3,  5, 0, 0, 0, 0, 1,    0,  -3,   -2,  -1],
      [ 0, 0, 0, 0, 0,  0,  3, -5, 0, 0, 0, 0, 0,  -20,  34,    0,   0],

   /* 291-300 */
      [ 0, 0, 0, 0, 0,  1, -2,  0, 0, 0, 0, 0, 1,  -10,   0,    1,   5],
      [ 0, 1,-1, 1, 0,  1, -3,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   2],
      [ 0, 0, 0, 0, 0,  1, -2,  0, 0, 0, 0, 0, 0,   22, -87,    0,   0],
      [ 0, 0, 0, 0, 0, -1,  2,  0, 0, 0, 0, 0, 1,   -4,   0,    0,   2],
      [ 0, 0, 0, 0, 0, -1,  2,  0, 0, 0, 0, 0, 2,   -3,  -6,   -2,   1],
      [ 0, 0, 0, 0, 0, -7, 11,  0, 0, 0, 0, 0, 2,  -16,  -3,   -1,   7],
      [ 0, 0, 0, 0, 0, -7, 11,  0, 0, 0, 0, 0, 1,    0,  -3,   -2,   0],
      [ 0,-2, 2, 0, 0,  4, -4,  0, 0, 0, 0, 0, 0,    4,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  2, -3, 0, 0, 0, 0, 0,  -68,  39,    0,   0],
      [ 0, 2,-2, 1, 0, -4,  4,  0, 0, 0, 0, 0, 0,   27,   0,    0, -14],

   /* 301-310 */
      [ 0,-1, 1, 0, 0,  4, -5,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  1, -1, 0, 0, 0, 0, 0,  -25,   0,    0,   0],
      [ 0, 0, 0, 0, 0, -4,  7,  0, 0, 0, 0, 0, 1,  -12,  -3,   -2,   6],
      [ 0, 1,-1, 1, 0, -4,  6,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1],
      [ 0, 0, 0, 0, 0, -4,  7,  0, 0, 0, 0, 0, 2,    3,  66,   29,  -1],
      [ 0, 0, 0, 0, 0, -4,  6,  0, 0, 0, 0, 0, 2,  490,   0,    0,-213],
      [ 0, 0, 0, 0, 0, -4,  6,  0, 0, 0, 0, 0, 1,  -22,  93,   49,  12],
      [ 0, 1,-1, 1, 0, -4,  5,  0, 0, 0, 0, 0, 0,   -7,  28,   15,   4],
      [ 0, 0, 0, 0, 0, -4,  6,  0, 0, 0, 0, 0, 1,   -3,  13,    7,   2],
      [ 0, 0, 0, 0, 0,  4, -6,  0, 0, 0, 0, 0, 0,  -46,  14,    0,   0],

   /* 311-320 */
      [-2, 0, 2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  0,  1, 0, 0, 0, 0, 0,    2,   1,    0,   0],
      [ 0,-1, 1, 0, 0,  1,  0,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0],
      [ 0, 0, 0, 1, 0,  1, -1,  0, 0, 0, 0, 0, 0,  -28,   0,    0,  15],
      [ 0, 0, 0, 0, 0,  0, -1,  0, 5, 0, 0, 0, 2,    5,   0,    0,  -2],
      [ 0, 0, 0, 0, 0,  0,  1, -3, 0, 0, 0, 0, 0,    0,   3,    0,   0],
      [ 0, 0, 0, 0, 0,  0, -1,  3, 0, 0, 0, 0, 2,  -11,   0,    0,   5],
      [ 0, 0, 0, 0, 0,  0, -7, 12, 0, 0, 0, 0, 2,    0,   3,    1,   0],
      [ 0, 0, 0, 0, 0, -1,  1,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1],
      [ 0, 0, 0, 0, 0, -1,  1,  0, 0, 0, 0, 0, 1,   25, 106,   57, -13],

   /* 321-330 */
      [ 0, 1,-1, 1, 0, -1,  0,  0, 0, 0, 0, 0, 0,    5,  21,   11,  -3],
      [ 0, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0, 1485,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 1,   -7, -32,  -17,   4],
      [ 0, 1,-1, 1, 0,  1, -2,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0],
      [ 0, 0, 0, 0, 0,  0, -2,  5, 0, 0, 0, 0, 2,   -6,  -3,   -2,   3],
      [ 0, 0, 0, 0, 0,  0, -1,  0, 4, 0, 0, 0, 2,   30,  -6,   -2, -13],
      [ 0, 0, 0, 0, 0,  0,  1,  0,-4, 0, 0, 0, 0,   -4,   4,    0,   0],
      [ 0, 0, 0, 1, 0, -1,  1,  0, 0, 0, 0, 0, 0,  -19,   0,    0,  10],
      [ 0, 0, 0, 0, 0,  0, -6, 10, 0, 0, 0, 0, 2,    0,   4,    2,  -1],
      [ 0, 0, 0, 0, 0,  0, -6, 10, 0, 0, 0, 0, 0,    0,   3,    0,   0],

   /* 331-340 */
      [ 0, 2,-2, 1, 0,  0, -3,  0, 3, 0, 0, 0, 0,    4,   0,    0,  -2],
      [ 0, 0, 0, 0, 0,  0, -3,  7, 0, 0, 0, 0, 2,    0,  -3,   -1,   0],
      [-2, 0, 2, 0, 0,  4, -4,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0, -5,  8, 0, 0, 0, 0, 2,    5,   3,    1,  -2],
      [ 0, 0, 0, 0, 0,  0,  5, -8, 0, 0, 0, 0, 0,    0,  11,    0,   0],
      [ 0, 0, 0, 0, 0,  0, -1,  0, 3, 0, 0, 0, 2,  118,   0,    0, -52],
      [ 0, 0, 0, 0, 0,  0, -1,  0, 3, 0, 0, 0, 1,    0,  -5,   -3,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  0,-3, 0, 0, 0, 0,  -28,  36,    0,   0],
      [ 0, 0, 0, 0, 0,  2, -4,  0, 0, 0, 0, 0, 0,    5,  -5,    0,   0],
      [ 0, 0, 0, 0, 0, -2,  4,  0, 0, 0, 0, 0, 1,   14, -59,  -31,  -8],

   /* 341-350 */
      [ 0, 1,-1, 1, 0, -2,  3,  0, 0, 0, 0, 0, 0,    0,   9,    5,   1],
      [ 0, 0, 0, 0, 0, -2,  4,  0, 0, 0, 0, 0, 2, -458,   0,    0, 198],
      [ 0, 0, 0, 0, 0, -6,  9,  0, 0, 0, 0, 0, 2,    0, -45,  -20,   0],
      [ 0, 0, 0, 0, 0, -6,  9,  0, 0, 0, 0, 0, 1,    9,   0,    0,  -5],
      [ 0, 0, 0, 0, 0,  6, -9,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0],
      [ 0, 0, 0, 1, 0,  0,  1,  0,-2, 0, 0, 0, 0,    0,  -4,   -2,  -1],
      [ 0, 2,-2, 1, 0, -2,  2,  0, 0, 0, 0, 0, 0,   11,   0,    0,  -6],
      [ 0, 0, 0, 0, 0,  0, -4,  6, 0, 0, 0, 0, 2,    6,   0,    0,  -2],
      [ 0, 0, 0, 0, 0,  0,  4, -6, 0, 0, 0, 0, 0,  -16,  23,    0,   0],
      [ 0, 0, 0, 1, 0,  3, -4,  0, 0, 0, 0, 0, 0,    0,  -4,   -2,   0],

   /* 351-360 */
      [ 0, 0, 0, 0, 0,  0, -1,  0, 2, 0, 0, 0, 2,   -5,   0,    0,   2],
      [ 0, 0, 0, 0, 0,  0,  1,  0,-2, 0, 0, 0, 0, -166, 269,    0,   0],
      [ 0, 0, 0, 1, 0,  0,  1,  0,-1, 0, 0, 0, 0,   15,   0,    0,  -8],
      [ 0, 0, 0, 0, 0, -5,  9,  0, 0, 0, 0, 0, 2,   10,   0,    0,  -4],
      [ 0, 0, 0, 0, 0,  0,  3, -4, 0, 0, 0, 0, 0,  -78,  45,    0,   0],
      [ 0, 0, 0, 0, 0, -3,  4,  0, 0, 0, 0, 0, 2,    0,  -5,   -2,   0],
      [ 0, 0, 0, 0, 0, -3,  4,  0, 0, 0, 0, 0, 1,    7,   0,    0,  -4],
      [ 0, 0, 0, 0, 0,  3, -4,  0, 0, 0, 0, 0, 0,   -5, 328,    0,   0],
      [ 0, 0, 0, 0, 0,  3, -4,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -2],
      [ 0, 0, 0, 1, 0,  0,  2, -2, 0, 0, 0, 0, 0,    5,   0,    0,  -2],

   /* 361-370 */
      [ 0, 0, 0, 1, 0,  0, -1,  0, 2, 0, 0, 0, 0,    0,   3,    1,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  0, 0,-3, 0, 0, 0,   -3,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  0, 1,-5, 0, 0, 0,   -3,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0, -1,  0, 1, 0, 0, 0, 1,    0,  -4,   -2,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,-1223, -26,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 1,    0,   7,    3,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  0,-3, 5, 0, 0, 0,    3,   0,    0,   0],
      [ 0, 0, 0, 1, 0, -3,  4,  0, 0, 0, 0, 0, 0,    0,   3,    2,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  0, 0,-2, 0, 0, 0,   -6,  20,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  2, -2, 0, 0, 0, 0, 0, -368,   0,    0,   0],

   /* 371-380 */
      [ 0, 0, 0, 0, 0,  0,  1,  0, 0,-1, 0, 0, 0,  -75,   0,    0,   0],
      [ 0, 0, 0, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,   11,   0,    0,  -6],
      [ 0, 0, 0, 1, 0,  0, -2,  2, 0, 0, 0, 0, 0,    3,   0,    0,  -2],
      [ 0, 0, 0, 0, 0, -8, 14,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1],
      [ 0, 0, 0, 0, 0,  0,  1,  0, 2,-5, 0, 0, 0,  -13, -30,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  5, -8, 3, 0, 0, 0, 0,   21,   3,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  5, -8, 3, 0, 0, 0, 2,   -3,   0,    0,   1],
      [ 0, 0, 0, 0, 0,  0, -1,  0, 0, 0, 0, 0, 1,   -4,   0,    0,   2],
      [ 0, 0, 0, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    8, -27,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  3, -8, 3, 0, 0, 0, 0,  -19, -11,    0,   0],

   /* 381-390 */
      [ 0, 0, 0, 0, 0,  0, -3,  8,-3, 0, 0, 0, 2,   -4,   0,    0,   2],
      [ 0, 0, 0, 0, 0,  0,  1,  0,-2, 5, 0, 0, 2,    0,   5,    2,   0],
      [ 0, 0, 0, 0, 0, -8, 12,  0, 0, 0, 0, 0, 2,   -6,   0,    0,   2],
      [ 0, 0, 0, 0, 0, -8, 12,  0, 0, 0, 0, 0, 0,   -8,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  0, 1,-2, 0, 0, 0,   -1,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  0, 0, 1, 0, 0, 2,  -14,   0,    0,   6],
      [ 0, 0, 0, 0, 0,  0,  0,  2, 0, 0, 0, 0, 0,    6,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  0,  2, 0, 0, 0, 0, 2,  -74,   0,    0,  32],
      [ 0, 0, 0, 0, 0,  0,  1,  0, 0, 2, 0, 0, 2,    0,  -3,   -1,   0],
      [ 0, 2,-2, 1, 0, -5,  5,  0, 0, 0, 0, 0, 0,    4,   0,    0,  -2],

   /* 391-400 */
      [ 0, 0, 0, 0, 0,  0,  1,  0, 1, 0, 0, 0, 0,    8,  11,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  0, 1, 0, 0, 0, 1,    0,   3,    2,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  0, 1, 0, 0, 0, 2, -262,   0,    0, 114],
      [ 0, 0, 0, 0, 0,  3, -6,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0],
      [ 0, 0, 0, 0, 0, -3,  6,  0, 0, 0, 0, 0, 1,   -7,   0,    0,   4],
      [ 0, 0, 0, 0, 0, -3,  6,  0, 0, 0, 0, 0, 2,    0, -27,  -12,   0],
      [ 0, 0, 0, 0, 0,  0, -1,  4, 0, 0, 0, 0, 2,  -19,  -8,   -4,   8],
      [ 0, 0, 0, 0, 0, -5,  7,  0, 0, 0, 0, 0, 2,  202,   0,    0, -87],
      [ 0, 0, 0, 0, 0, -5,  7,  0, 0, 0, 0, 0, 1,   -8,  35,   19,   5],
      [ 0, 1,-1, 1, 0, -5,  6,  0, 0, 0, 0, 0, 0,    0,   4,    2,   0],

   /* 401-410 */
      [ 0, 0, 0, 0, 0,  5, -7,  0, 0, 0, 0, 0, 0,   16,  -5,    0,   0],
      [ 0, 2,-2, 1, 0,  0, -1,  0, 1, 0, 0, 0, 0,    5,   0,    0,  -3],
      [ 0, 0, 0, 0, 0,  0, -1,  0, 1, 0, 0, 0, 0,    0,  -3,    0,   0],
      [ 0, 0, 0, 0,-1,  0,  3,  0, 0, 0, 0, 0, 2,    1,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  0, 2, 0, 0, 0, 2,  -35, -48,  -21,  15],
      [ 0, 0, 0, 0, 0,  0, -2,  6, 0, 0, 0, 0, 2,   -3,  -5,   -2,   1],
      [ 0, 0, 0, 1, 0,  2, -2,  0, 0, 0, 0, 0, 0,    6,   0,    0,  -3],
      [ 0, 0, 0, 0, 0,  0, -6,  9, 0, 0, 0, 0, 2,    3,   0,    0,  -1],
      [ 0, 0, 0, 0, 0,  0,  6, -9, 0, 0, 0, 0, 0,    0,  -5,    0,   0],
      [ 0, 0, 0, 0, 0, -2,  2,  0, 0, 0, 0, 0, 1,   12,  55,   29,  -6],

   /* 411-420 */
      [ 0, 1,-1, 1, 0, -2,  1,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0],
      [ 0, 0, 0, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0, -598,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  2, -2,  0, 0, 0, 0, 0, 1,   -3, -13,   -7,   1],
      [ 0, 0, 0, 0, 0,  0,  1,  0, 3, 0, 0, 0, 2,   -5,  -7,   -3,   2],
      [ 0, 0, 0, 0, 0,  0, -5,  7, 0, 0, 0, 0, 2,    3,   0,    0,  -1],
      [ 0, 0, 0, 0, 0,  0,  5, -7, 0, 0, 0, 0, 0,    5,  -7,    0,   0],
      [ 0, 0, 0, 1, 0, -2,  2,  0, 0, 0, 0, 0, 0,    4,   0,    0,  -2],
      [ 0, 0, 0, 0, 0,  0,  4, -5, 0, 0, 0, 0, 0,   16,  -6,    0,   0],
      [ 0, 0, 0, 0, 0,  1, -3,  0, 0, 0, 0, 0, 0,    8,  -3,    0,   0],
      [ 0, 0, 0, 0, 0, -1,  3,  0, 0, 0, 0, 0, 1,    8, -31,  -16,  -4],

   /* 421-430 */
      [ 0, 1,-1, 1, 0, -1,  2,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0],
      [ 0, 0, 0, 0, 0, -1,  3,  0, 0, 0, 0, 0, 2,  113,   0,    0, -49],
      [ 0, 0, 0, 0, 0, -7, 10,  0, 0, 0, 0, 0, 2,    0, -24,  -10,   0],
      [ 0, 0, 0, 0, 0, -7, 10,  0, 0, 0, 0, 0, 1,    4,   0,    0,  -2],
      [ 0, 0, 0, 0, 0,  0,  3, -3, 0, 0, 0, 0, 0,   27,   0,    0,   0],
      [ 0, 0, 0, 0, 0, -4,  8,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1],
      [ 0, 0, 0, 0, 0, -4,  5,  0, 0, 0, 0, 0, 2,    0,  -4,   -2,   0],
      [ 0, 0, 0, 0, 0, -4,  5,  0, 0, 0, 0, 0, 1,    5,   0,    0,  -2],
      [ 0, 0, 0, 0, 0,  4, -5,  0, 0, 0, 0, 0, 0,    0,  -3,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  1, 0, 0, 0, 0, 2,  -13,   0,    0,   6],

   /* 431-440 */
      [ 0, 0, 0, 0, 0,  0, -2,  0, 5, 0, 0, 0, 2,    5,   0,    0,  -2],
      [ 0, 0, 0, 0, 0,  0,  0,  3, 0, 0, 0, 0, 2,  -18, -10,   -4,   8],
      [ 0, 0, 0, 0, 0,  1,  0,  0, 0, 0, 0, 0, 0,   -4, -28,    0,   0],
      [ 0, 0, 0, 0, 0,  1,  0,  0, 0, 0, 0, 0, 2,   -5,   6,    3,   2],
      [ 0, 0, 0, 0, 0, -9, 13,  0, 0, 0, 0, 0, 2,   -3,   0,    0,   1],
      [ 0, 0, 0, 0, 0,  0, -1,  5, 0, 0, 0, 0, 2,   -5,  -9,   -4,   2],
      [ 0, 0, 0, 0, 0,  0, -2,  0, 4, 0, 0, 0, 2,   17,   0,    0,  -7],
      [ 0, 0, 0, 0, 0,  0,  2,  0,-4, 0, 0, 0, 0,   11,   4,    0,   0],
      [ 0, 0, 0, 0, 0,  0, -2,  7, 0, 0, 0, 0, 2,    0,  -6,   -2,   0],
      [ 0, 0, 0, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   83,  15,    0,   0],

   /* 441-450 */
      [ 0, 0, 0, 0, 0, -2,  5,  0, 0, 0, 0, 0, 1,   -4,   0,    0,   2],
      [ 0, 0, 0, 0, 0, -2,  5,  0, 0, 0, 0, 0, 2,    0,-114,  -49,   0],
      [ 0, 0, 0, 0, 0, -6,  8,  0, 0, 0, 0, 0, 2,  117,   0,    0, -51],
      [ 0, 0, 0, 0, 0, -6,  8,  0, 0, 0, 0, 0, 1,   -5,  19,   10,   2],
      [ 0, 0, 0, 0, 0,  6, -8,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0],
      [ 0, 0, 0, 1, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   2],
      [ 0, 0, 0, 0, 0,  0, -3,  9, 0, 0, 0, 0, 2,    0,  -3,   -1,   0],
      [ 0, 0, 0, 0, 0,  0,  5, -6, 0, 0, 0, 0, 0,    3,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  5, -6, 0, 0, 0, 0, 2,    0,  -6,   -2,   0],
      [ 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,  393,   3,    0,   0],

   /* 451-460 */
      [ 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 1,   -4,  21,   11,   2],
      [ 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 2,   -6,   0,   -1,   3],
      [ 0, 0, 0, 0, 0, -5, 10,  0, 0, 0, 0, 0, 2,   -3,   8,    4,   1],
      [ 0, 0, 0, 0, 0,  0,  4, -4, 0, 0, 0, 0, 0,    8,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  4, -4, 0, 0, 0, 0, 2,   18, -29,  -13,  -8],
      [ 0, 0, 0, 0, 0, -3,  3,  0, 0, 0, 0, 0, 1,    8,  34,   18,  -4],
      [ 0, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,   89,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 1,    3,  12,    6,  -1],
      [ 0, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 2,   54, -15,   -7, -24],
      [ 0, 0, 0, 0, 0,  0,  2,  0, 0,-3, 0, 0, 0,    0,   3,    0,   0],

   /* 461-470 */
      [ 0, 0, 0, 0, 0,  0, -5, 13, 0, 0, 0, 0, 2,    3,   0,    0,  -1],
      [ 0, 0, 0, 0, 0,  0,  2,  0,-1, 0, 0, 0, 0,    0,  35,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  2,  0,-1, 0, 0, 0, 2, -154, -30,  -13,  67],
      [ 0, 0, 0, 0, 0,  0,  2,  0, 0,-2, 0, 0, 0,   15,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  2,  0, 0,-2, 0, 0, 1,    0,   4,    2,   0],
      [ 0, 0, 0, 0, 0,  0,  3, -2, 0, 0, 0, 0, 0,    0,   9,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  3, -2, 0, 0, 0, 0, 2,   80, -71,  -31, -35],
      [ 0, 0, 0, 0, 0,  0,  2,  0, 0,-1, 0, 0, 2,    0, -20,   -9,   0],
      [ 0, 0, 0, 0, 0,  0, -6, 15, 0, 0, 0, 0, 2,   11,   5,    2,  -5],
      [ 0, 0, 0, 0, 0, -8, 15,  0, 0, 0, 0, 0, 2,   61, -96,  -42, -27],

   /* 471-480 */
      [ 0, 0, 0, 0, 0, -3,  9, -4, 0, 0, 0, 0, 2,   14,   9,    4,  -6],
      [ 0, 0, 0, 0, 0,  0,  2,  0, 2,-5, 0, 0, 2,  -11,  -6,   -3,   5],
      [ 0, 0, 0, 0, 0,  0, -2,  8,-1,-5, 0, 0, 2,    0,  -3,   -1,   0],
      [ 0, 0, 0, 0, 0,  0,  6, -8, 3, 0, 0, 0, 2,  123,-415, -180, -53],
      [ 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 0,    0,   0,    0, -35],
      [ 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 1,    7, -32,  -17,  -4],
      [ 0, 1,-1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -9,   -5,   0],
      [ 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 1,    0,  -4,    2,   0],
      [ 0, 0, 0, 0, 0,  0,  2,  0, 0, 0, 0, 0, 2,  -89,   0,    0,  38],

   /* 481-490 */
      [ 0, 0, 0, 0, 0,  0, -6, 16,-4,-5, 0, 0, 2,    0, -86,  -19,  -6],
      [ 0, 0, 0, 0, 0,  0, -2,  8,-3, 0, 0, 0, 2,    0,   0,  -19,   6],
      [ 0, 0, 0, 0, 0,  0, -2,  8,-3, 0, 0, 0, 2, -123,-416, -180,  53],
      [ 0, 0, 0, 0, 0,  0,  6, -8, 1, 5, 0, 0, 2,    0,  -3,   -1,   0],
      [ 0, 0, 0, 0, 0,  0,  2,  0,-2, 5, 0, 0, 2,   12,  -6,   -3,  -5],
      [ 0, 0, 0, 0, 0,  3, -5,  4, 0, 0, 0, 0, 2,  -13,   9,    4,   6],
      [ 0, 0, 0, 0, 0, -8, 11,  0, 0, 0, 0, 0, 2,    0, -15,   -7,   0],
      [ 0, 0, 0, 0, 0, -8, 11,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -1],
      [ 0, 0, 0, 0, 0, -8, 11,  0, 0, 0, 0, 0, 2,  -62, -97,  -42,  27],
      [ 0, 0, 0, 0, 0,  0, 11,  0, 0, 0, 0, 0, 2,  -11,   5,    2,   5],

   /* 491-500 */
      [ 0, 0, 0, 0, 0,  0,  2,  0, 0, 1, 0, 0, 2,    0, -19,   -8,   0],
      [ 0, 0, 0, 0, 0,  3, -3,  0, 2, 0, 0, 0, 2,   -3,   0,    0,   1],
      [ 0, 2,-2, 1, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,   4,    2,   0],
      [ 0, 1,-1, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,   3,    0,   0],
      [ 0, 2,-2, 1, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,   4,    2,   0],
      [ 0, 0, 0, 0, 0,  0,  1,  2, 0, 0, 0, 0, 2,  -85, -70,  -31,  37],
      [ 0, 0, 0, 0, 0,  0,  2,  0, 1, 0, 0, 0, 2,  163, -12,   -5, -72],
      [ 0, 0, 0, 0, 0, -3,  7,  0, 0, 0, 0, 0, 2,  -63, -16,   -7,  28],
      [ 0, 0, 0, 0, 0,  0,  0,  4, 0, 0, 0, 0, 2,  -21, -32,  -14,   9],
      [ 0, 0, 0, 0, 0, -5,  6,  0, 0, 0, 0, 0, 2,    0,  -3,   -1,   0],

   /* 501-510 */
      [ 0, 0, 0, 0, 0, -5,  6,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -2],
      [ 0, 0, 0, 0, 0,  5, -6,  0, 0, 0, 0, 0, 0,    0,   8,    0,   0],
      [ 0, 0, 0, 0, 0,  5, -6,  0, 0, 0, 0, 0, 2,    3,  10,    4,  -1],
      [ 0, 0, 0, 0, 0,  0,  2,  0, 2, 0, 0, 0, 2,    3,   0,    0,  -1],
      [ 0, 0, 0, 0, 0,  0, -1,  6, 0, 0, 0, 0, 2,    0,  -7,   -3,   0],
      [ 0, 0, 0, 0, 0,  0,  7, -9, 0, 0, 0, 0, 2,    0,  -4,   -2,   0],
      [ 0, 0, 0, 0, 0,  2, -1,  0, 0, 0, 0, 0, 0,    6,  19,    0,   0],
      [ 0, 0, 0, 0, 0,  2, -1,  0, 0, 0, 0, 0, 2,    5,-173,  -75,  -2],
      [ 0, 0, 0, 0, 0,  0,  6, -7, 0, 0, 0, 0, 2,    0,  -7,   -3,   0],
      [ 0, 0, 0, 0, 0,  0,  5, -5, 0, 0, 0, 0, 2,    7, -12,   -5,  -3],

   /* 511-520 */
      [ 0, 0, 0, 0, 0, -1,  4,  0, 0, 0, 0, 0, 1,   -3,   0,    0,   2],
      [ 0, 0, 0, 0, 0, -1,  4,  0, 0, 0, 0, 0, 2,    3,  -4,   -2,  -1],
      [ 0, 0, 0, 0, 0, -7,  9,  0, 0, 0, 0, 0, 2,   74,   0,    0, -32],
      [ 0, 0, 0, 0, 0, -7,  9,  0, 0, 0, 0, 0, 1,   -3,  12,    6,   2],
      [ 0, 0, 0, 0, 0,  0,  4, -3, 0, 0, 0, 0, 2,   26, -14,   -6, -11],
      [ 0, 0, 0, 0, 0,  0,  3, -1, 0, 0, 0, 0, 2,   19,   0,    0,  -8],
      [ 0, 0, 0, 0, 0, -4,  4,  0, 0, 0, 0, 0, 1,    6,  24,   13,  -3],
      [ 0, 0, 0, 0, 0,  4, -4,  0, 0, 0, 0, 0, 0,   83,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  4, -4,  0, 0, 0, 0, 0, 1,    0, -10,   -5,   0],
      [ 0, 0, 0, 0, 0,  4, -4,  0, 0, 0, 0, 0, 2,   11,  -3,   -1,  -5],

   /* 521-530 */
      [ 0, 0, 0, 0, 0,  0,  2,  1, 0, 0, 0, 0, 2,    3,   0,    1,  -1],
      [ 0, 0, 0, 0, 0,  0, -3,  0, 5, 0, 0, 0, 2,    3,   0,    0,  -1],
      [ 0, 0, 0, 0, 0,  1,  1,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  1,  1,  0, 0, 0, 0, 0, 1,    5, -23,  -12,  -3],
      [ 0, 0, 0, 0, 0,  1,  1,  0, 0, 0, 0, 0, 2, -339,   0,    0, 147],
      [ 0, 0, 0, 0, 0, -9, 12,  0, 0, 0, 0, 0, 2,    0, -10,   -5,   0],
      [ 0, 0, 0, 0, 0,  0,  3,  0,-4, 0, 0, 0, 0,    5,   0,    0,   0],
      [ 0, 2,-2, 1, 0,  1, -1,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1],
      [ 0, 0, 0, 0, 0,  0,  7, -8, 0, 0, 0, 0, 2,    0,  -4,   -2,   0],
      [ 0, 0, 0, 0, 0,  0,  3,  0,-3, 0, 0, 0, 0,   18,  -3,    0,   0],

   /* 531-540 */
      [ 0, 0, 0, 0, 0,  0,  3,  0,-3, 0, 0, 0, 2,    9, -11,   -5,  -4],
      [ 0, 0, 0, 0, 0, -2,  6,  0, 0, 0, 0, 0, 2,   -8,   0,    0,   4],
      [ 0, 0, 0, 0, 0, -6,  7,  0, 0, 0, 0, 0, 1,    3,   0,    0,  -1],
      [ 0, 0, 0, 0, 0,  6, -7,  0, 0, 0, 0, 0, 0,    0,   9,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  6, -6, 0, 0, 0, 0, 2,    6,  -9,   -4,  -2],
      [ 0, 0, 0, 0, 0,  0,  3,  0,-2, 0, 0, 0, 0,   -4, -12,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  3,  0,-2, 0, 0, 0, 2,   67, -91,  -39, -29],
      [ 0, 0, 0, 0, 0,  0,  5, -4, 0, 0, 0, 0, 2,   30, -18,   -8, -13],
      [ 0, 0, 0, 0, 0,  3, -2,  0, 0, 0, 0, 0, 0,    0,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  3, -2,  0, 0, 0, 0, 0, 2,    0,-114,  -50,   0],

   /* 541-550 */
      [ 0, 0, 0, 0, 0,  0,  3,  0,-1, 0, 0, 0, 2,    0,   0,    0,  23],
      [ 0, 0, 0, 0, 0,  0,  3,  0,-1, 0, 0, 0, 2,  517,  16,    7,-224],
      [ 0, 0, 0, 0, 0,  0,  3,  0, 0,-2, 0, 0, 2,    0,  -7,   -3,   0],
      [ 0, 0, 0, 0, 0,  0,  4, -2, 0, 0, 0, 0, 2,  143,  -3,   -1, -62],
      [ 0, 0, 0, 0, 0,  0,  3,  0, 0,-1, 0, 0, 2,   29,   0,    0, -13],
      [ 0, 2,-2, 1, 0,  0,  1,  0,-1, 0, 0, 0, 0,   -4,   0,    0,   2],
      [ 0, 0, 0, 0, 0, -8, 16,  0, 0, 0, 0, 0, 2,   -6,   0,    0,   3],
      [ 0, 0, 0, 0, 0,  0,  3,  0, 2,-5, 0, 0, 2,    5,  12,    5,  -2],
      [ 0, 0, 0, 0, 0,  0,  7, -8, 3, 0, 0, 0, 2,  -25,   0,    0,  11],
      [ 0, 0, 0, 0, 0,  0, -5, 16,-4,-5, 0, 0, 2,   -3,   0,    0,   1],

   /* 551-560 */
      [ 0, 0, 0, 0, 0,  0,  3,  0, 0, 0, 0, 0, 2,    0,   4,    2,   0],
      [ 0, 0, 0, 0, 0,  0, -1,  8,-3, 0, 0, 0, 2,  -22,  12,    5,  10],
      [ 0, 0, 0, 0, 0, -8, 10,  0, 0, 0, 0, 0, 2,   50,   0,    0, -22],
      [ 0, 0, 0, 0, 0, -8, 10,  0, 0, 0, 0, 0, 1,    0,   7,    4,   0],
      [ 0, 0, 0, 0, 0, -8, 10,  0, 0, 0, 0, 0, 2,    0,   3,    1,   0],
      [ 0, 0, 0, 0, 0,  0,  2,  2, 0, 0, 0, 0, 2,   -4,   4,    2,   2],
      [ 0, 0, 0, 0, 0,  0,  3,  0, 1, 0, 0, 0, 2,   -5, -11,   -5,   2],
      [ 0, 0, 0, 0, 0, -3,  8,  0, 0, 0, 0, 0, 2,    0,   4,    2,   0],
      [ 0, 0, 0, 0, 0, -5,  5,  0, 0, 0, 0, 0, 1,    4,  17,    9,  -2],
      [ 0, 0, 0, 0, 0,  5, -5,  0, 0, 0, 0, 0, 0,   59,   0,    0,   0],

   /* 561-570 */
      [ 0, 0, 0, 0, 0,  5, -5,  0, 0, 0, 0, 0, 1,    0,  -4,   -2,   0],
      [ 0, 0, 0, 0, 0,  5, -5,  0, 0, 0, 0, 0, 2,   -8,   0,    0,   4],
      [ 0, 0, 0, 0, 0,  2,  0,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  2,  0,  0, 0, 0, 0, 0, 1,    4, -15,   -8,  -2],
      [ 0, 0, 0, 0, 0,  2,  0,  0, 0, 0, 0, 0, 2,  370,  -8,    0,-160],
      [ 0, 0, 0, 0, 0,  0,  7, -7, 0, 0, 0, 0, 2,    0,   0,   -3,   0],
      [ 0, 0, 0, 0, 0,  0,  7, -7, 0, 0, 0, 0, 2,    0,   3,    1,   0],
      [ 0, 0, 0, 0, 0,  0,  6, -5, 0, 0, 0, 0, 2,   -6,   3,    1,   3],
      [ 0, 0, 0, 0, 0,  7, -8,  0, 0, 0, 0, 0, 0,    0,   6,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  5, -3, 0, 0, 0, 0, 2,  -10,   0,    0,   4],

   /* 571-580 */
      [ 0, 0, 0, 0, 0,  4, -3,  0, 0, 0, 0, 0, 2,    0,   9,    4,   0],
      [ 0, 0, 0, 0, 0,  1,  2,  0, 0, 0, 0, 0, 2,    4,  17,    7,  -2],
      [ 0, 0, 0, 0, 0, -9, 11,  0, 0, 0, 0, 0, 2,   34,   0,    0, -15],
      [ 0, 0, 0, 0, 0, -9, 11,  0, 0, 0, 0, 0, 1,    0,   5,    3,   0],
      [ 0, 0, 0, 0, 0,  0,  4,  0,-4, 0, 0, 0, 2,   -5,   0,    0,   2],
      [ 0, 0, 0, 0, 0,  0,  4,  0,-3, 0, 0, 0, 2,  -37,  -7,   -3,  16],
      [ 0, 0, 0, 0, 0, -6,  6,  0, 0, 0, 0, 0, 1,    3,  13,    7,  -2],
      [ 0, 0, 0, 0, 0,  6, -6,  0, 0, 0, 0, 0, 0,   40,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  6, -6,  0, 0, 0, 0, 0, 1,    0,  -3,   -2,   0],
      [ 0, 0, 0, 0, 0,  0,  4,  0,-2, 0, 0, 0, 2, -184,  -3,   -1,  80],

   /* 581-590 */
      [ 0, 0, 0, 0, 0,  0,  6, -4, 0, 0, 0, 0, 2,   -3,   0,    0,   1],
      [ 0, 0, 0, 0, 0,  3, -1,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  3, -1,  0, 0, 0, 0, 0, 1,    0, -10,   -6,  -1],
      [ 0, 0, 0, 0, 0,  3, -1,  0, 0, 0, 0, 0, 2,   31,  -6,    0, -13],
      [ 0, 0, 0, 0, 0,  0,  4,  0,-1, 0, 0, 0, 2,   -3, -32,  -14,   1],
      [ 0, 0, 0, 0, 0,  0,  4,  0, 0,-2, 0, 0, 2,   -7,   0,    0,   3],
      [ 0, 0, 0, 0, 0,  0,  5, -2, 0, 0, 0, 0, 2,    0,  -8,   -4,   0],
      [ 0, 0, 0, 0, 0,  0,  4,  0, 0, 0, 0, 0, 0,    3,  -4,    0,   0],
      [ 0, 0, 0, 0, 0,  8, -9,  0, 0, 0, 0, 0, 0,    0,   4,    0,   0],
      [ 0, 0, 0, 0, 0,  5, -4,  0, 0, 0, 0, 0, 2,    0,   3,    1,   0],

   /* 591-600 */
      [ 0, 0, 0, 0, 0,  2,  1,  0, 0, 0, 0, 0, 2,   19, -23,  -10,   2],
      [ 0, 0, 0, 0, 0,  2,  1,  0, 0, 0, 0, 0, 1,    0,   0,    0, -10],
      [ 0, 0, 0, 0, 0,  2,  1,  0, 0, 0, 0, 0, 1,    0,   3,    2,   0],
      [ 0, 0, 0, 0, 0, -7,  7,  0, 0, 0, 0, 0, 1,    0,   9,    5,  -1],
      [ 0, 0, 0, 0, 0,  7, -7,  0, 0, 0, 0, 0, 0,   28,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 1,    0,  -7,   -4,   0],
      [ 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 2,    8,  -4,    0,  -4],
      [ 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 0,    0,   0,   -2,   0],
      [ 0, 0, 0, 0, 0,  4, -2,  0, 0, 0, 0, 0, 0,    0,   3,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  5,  0,-4, 0, 0, 0, 2,   -3,   0,    0,   1],

   /* 601-610 */
      [ 0, 0, 0, 0, 0,  0,  5,  0,-3, 0, 0, 0, 2,   -9,   0,    1,   4],
      [ 0, 0, 0, 0, 0,  0,  5,  0,-2, 0, 0, 0, 2,    3,  12,    5,  -1],
      [ 0, 0, 0, 0, 0,  3,  0,  0, 0, 0, 0, 0, 2,   17,  -3,   -1,   0],
      [ 0, 0, 0, 0, 0, -8,  8,  0, 0, 0, 0, 0, 1,    0,   7,    4,   0],
      [ 0, 0, 0, 0, 0,  8, -8,  0, 0, 0, 0, 0, 0,   19,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  5, -3,  0, 0, 0, 0, 0, 1,    0,  -5,   -3,   0],
      [ 0, 0, 0, 0, 0,  5, -3,  0, 0, 0, 0, 0, 2,   14,  -3,    0,  -1],
      [ 0, 0, 0, 0, 0, -9,  9,  0, 0, 0, 0, 0, 1,    0,   0,   -1,   0],
      [ 0, 0, 0, 0, 0, -9,  9,  0, 0, 0, 0, 0, 1,    0,   0,    0,  -5],
      [ 0, 0, 0, 0, 0, -9,  9,  0, 0, 0, 0, 0, 1,    0,   5,    3,   0],

   /* 611-620 */
      [ 0, 0, 0, 0, 0,  9, -9,  0, 0, 0, 0, 0, 0,   13,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  6, -4,  0, 0, 0, 0, 0, 1,    0,  -3,   -2,   0],
      [ 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 2,    2,   9,    4,   3],
      [ 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 0,    0,   0,    0,  -4],
      [ 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 0,    8,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 1,    0,   4,    2,   0],
      [ 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 2,    6,   0,    0,  -3],
      [ 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 0,    6,   0,    0,   0],
      [ 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 1,    0,   3,    1,   0],
      [ 0, 0, 0, 0, 0,  0,  6,  0, 0, 0, 0, 0, 2,    5,   0,    0,  -2],

   /* 621-630 */
      [ 0, 0, 0, 0, 0,  0,  0,  0, 0, 0, 0, 0, 2,    3,   0,    0,  -1],
      [ 1, 0,-2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   0],
      [ 1, 0,-2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    6,   0,    0,   0],
      [ 1, 0,-2, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    7,   0,    0,   0],
      [ 1, 0,-2, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   0],
      [-1, 0, 0, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,    4,   0,    0,   0],
      [-1, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,    6,   0,    0,   0],
      [-1, 0, 2, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -4,    0,   0],
      [ 1, 0,-2, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -4,    0,   0],
      [-2, 0, 2, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    5,   0,    0,   0],

   /* 631-640 */
      [-1, 0, 0, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   -3,   0,    0,   0],
      [-1, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    4,   0,    0,   0],
      [-1, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0,   -5,   0,    0,   0],
      [-1, 0, 2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    4,   0,    0,   0],
      [ 1,-1, 1, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,   3,    0,   0],
      [-1, 0, 2, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   13,   0,    0,   0],
      [-2, 0, 0, 0, 0,  0,  2,  0,-3, 0, 0, 0, 0,   21,  11,    0,   0],
      [ 1, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -5,    0,   0],
      [-1, 1,-1, 1, 0,  0, -1,  0, 0, 0, 0, 0, 0,    0,  -5,   -2,   0],
      [ 1, 1,-1, 1, 0,  0, -1,  0, 0, 0, 0, 0, 0,    0,   5,    3,   0],

   /* 641-650 */
      [-1, 0, 0, 0, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,  -5,    0,   0],
      [-1, 0, 2, 1, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   2],
      [ 0, 0, 0, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,   20,  10,    0,   0],
      [-1, 0, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,  -34,   0,    0,   0],
      [-1, 0, 2, 0, 0,  3, -3,  0, 0, 0, 0, 0, 0,  -19,   0,    0,   0],
      [ 1, 0,-2, 1, 0,  0, -2,  0, 2, 0, 0, 0, 0,    3,   0,    0,  -2],
      [ 1, 2,-2, 2, 0, -3,  3,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1],
      [ 1, 2,-2, 2, 0,  0, -2,  0, 2, 0, 0, 0, 0,   -6,   0,    0,   3],
      [ 1, 0, 0, 0, 0,  1, -1,  0, 0, 0, 0, 0, 0,   -4,   0,    0,   0],
      [ 1, 0, 0, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    3,   0,    0,   0],

   /* 651-660 */
      [ 0, 0,-2, 0, 0,  2, -2,  0, 0, 0, 0, 0, 0,    3,   0,    0,   0],
      [ 0, 0,-2, 0, 0,  0,  1,  0,-1, 0, 0, 0, 0,    4,   0,    0,   0],
      [ 0, 2, 0, 2, 0, -2,  2,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1],
      [ 0, 2, 0, 2, 0,  0, -1,  0, 1, 0, 0, 0, 0,    6,   0,    0,  -3],
      [ 0, 2, 0, 2, 0, -1,  1,  0, 0, 0, 0, 0, 0,   -8,   0,    0,   3],
      [ 0, 2, 0, 2, 0, -2,  3,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0],
      [ 0, 0, 2, 0, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   0],
      [ 0, 1, 1, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -3,   -2,   0],
      [ 1, 2, 0, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,  126, -63,  -27, -55],
      [-1, 2, 0, 2, 0, 10, -3,  0, 0, 0, 0, 0, 0,   -5,   0,    1,   2],

   /* 661-670 */
      [ 0, 1, 1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,   -3,  28,   15,   2],
      [ 1, 2, 0, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,    5,   0,    1,  -2],
      [ 0, 2, 0, 2, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,   9,    4,   1],
      [ 0, 2, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,   9,    4,  -1],
      [-1, 2, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0, -126, -63,  -27,  55],
      [ 2, 2,-2, 2, 0,  0, -2,  0, 3, 0, 0, 0, 0,    3,   0,    0,  -1],
      [ 1, 2, 0, 1, 0,  0, -2,  0, 3, 0, 0, 0, 0,   21, -11,   -6, -11],
      [ 0, 1, 1, 0, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,  -4,    0,   0],
      [-1, 2, 0, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,  -21, -11,   -6,  11],
      [-2, 2, 2, 2, 0,  0,  2,  0,-2, 0, 0, 0, 0,   -3,   0,    0,   1],

   /* 671-680 */
      [ 0, 2, 0, 2, 0,  2, -3,  0, 0, 0, 0, 0, 0,    0,   3,    1,   0],
      [ 0, 2, 0, 2, 0,  1, -1,  0, 0, 0, 0, 0, 0,    8,   0,    0,  -4],
      [ 0, 2, 0, 2, 0,  0,  1,  0,-1, 0, 0, 0, 0,   -6,   0,    0,   3],
      [ 0, 2, 0, 2, 0,  2, -2,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1],
      [-1, 2, 2, 2, 0,  0, -1,  0, 1, 0, 0, 0, 0,    3,   0,    0,  -1],
      [ 1, 2, 0, 2, 0, -1,  1,  0, 0, 0, 0, 0, 0,   -3,   0,    0,   1],
      [-1, 2, 2, 2, 0,  0,  2,  0,-3, 0, 0, 0, 0,   -5,   0,    0,   2],
      [ 2, 2, 0, 2, 0,  0,  2,  0,-3, 0, 0, 0, 0,   24, -12,   -5, -11],
      [ 1, 2, 0, 2, 0,  0, -4,  8,-3, 0, 0, 0, 0,    0,   3,    1,   0],
      [ 1, 2, 0, 2, 0,  0,  4, -8, 3, 0, 0, 0, 0,    0,   3,    1,   0],

   /* 681-687 */
      [ 1, 1, 1, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    0,   3,    2,   0],
      [ 0, 2, 0, 2, 0,  0,  1,  0, 0, 0, 0, 0, 0,  -24, -12,   -5,  10],
      [ 2, 2, 0, 1, 0,  0,  1,  0, 0, 0, 0, 0, 0,    4,   0,   -1,  -2],
      [-1, 2, 2, 2, 0,  0,  2,  0,-2, 0, 0, 0, 0,   13,   0,    0,  -6],
      [-1, 2, 2, 2, 0,  3, -3,  0, 0, 0, 0, 0, 0,    7,   0,    0,  -3],
      [ 1, 2, 0, 2, 0,  1, -1,  0, 0, 0, 0, 0, 0,    3,   0,    0,  -1],
      [ 0, 2, 2, 2, 0,  0,  2,  0,-2, 0, 0, 0, 0,    3,   0,    0,  -1]
   ];

/* Number of terms in the planetary nutation model */
   var NPL = xpl.length;

/* ------------------------------------------------------------------ */

/* Interval between fundamental date J2000.0 and given date (JC). */
   t = ((date1 - ERFA_DJ00) + date2) / ERFA_DJC;

/* ------------------- */
/* LUNI-SOLAR NUTATION */
/* ------------------- */

/* Fundamental (Delaunay) arguments */

/* Mean anomaly of the Moon (IERS 2003). */
   el = eraFal03(t);

/* Mean anomaly of the Sun (MHB2000). */
   elp = ((1287104.79305  +
            t * (129596581.0481  +
            t * (-0.5532  +
            t * (0.000136  +
            t * (-0.00001149))))) % (ERFA_TURNAS)) * ERFA_DAS2R;

/* Mean longitude of the Moon minus that of the ascending node */
/* (IERS 2003. */
   f = eraFaf03(t);

/* Mean elongation of the Moon from the Sun (MHB2000). */
   d = ((1072260.70369  +
          t * (1602961601.2090  +
          t * (-6.3706  +
          t * (0.006593  +
          t * (-0.00003169))))) % (ERFA_TURNAS)) * ERFA_DAS2R;

/* Mean longitude of the ascending node of the Moon (IERS 2003). */
   om = eraFaom03(t);

/* Initialize the nutation values. */
   dp = 0.0;
   de = 0.0;

/* Summation of luni-solar nutation series (in reverse order). */
   for (i = ~~(NLS-1); i >= 0; i--) {

   /* Argument and functions. */
      arg = ((xls[i][0]  * el +
                 xls[i][1] * elp +
                 xls[i][2]  * f +
                 xls[i][3]  * d +
                 xls[i][4] * om) % (ERFA_D2PI));
      sarg = Math.sin(arg);
      carg = Math.cos(arg);

   /* Term. */
      dp += (xls[i][5] + xls[i][6] * t) * sarg + xls[i][7] * carg;
      de += (xls[i][8] + xls[i][9] * t) * carg + xls[i][10] * sarg;
   }

/* Convert from 0.1 microarcsec units to radians. */
   dpsils = dp * U2R;
   depsls = de * U2R;

/* ------------------ */
/* PLANETARY NUTATION */
/* ------------------ */

/* n.b.  The MHB2000 code computes the luni-solar and planetary nutation */
/* in different functions, using slightly different Delaunay */
/* arguments in the two cases.  This behaviour is faithfully */
/* reproduced here.  Use of the IERS 2003 expressions for both */
/* cases leads to negligible changes, well below */
/* 0.1 microarcsecond. */

/* Mean anomaly of the Moon (MHB2000). */
   al = ((2.35555598 + 8328.6914269554 * t) % (ERFA_D2PI));

/* Mean longitude of the Moon minus that of the ascending node */
/*(MHB2000). */
   af = ((1.627905234 + 8433.466158131 * t) % (ERFA_D2PI));

/* Mean elongation of the Moon from the Sun (MHB2000). */
   ad = ((5.198466741 + 7771.3771468121 * t) % (ERFA_D2PI));

/* Mean longitude of the ascending node of the Moon (MHB2000). */
   aom = ((2.18243920 - 33.757045 * t) % (ERFA_D2PI));

/* General accumulated precession in longitude (IERS 2003). */
   apa = eraFapa03(t);

/* Planetary longitudes, Mercury through Uranus (IERS 2003). */
   alme = eraFame03(t);
   alve = eraFave03(t);
   alea = eraFae03(t);
   alma = eraFama03(t);
   alju = eraFaju03(t);
   alsa = eraFasa03(t);
   alur = eraFaur03(t);

/* Neptune longitude (MHB2000). */
   alne = ((5.321159000 + 3.8127774000 * t) % (ERFA_D2PI));

/* Initialize the nutation values. */
   dp = 0.0;
   de = 0.0;

/* Summation of planetary nutation series (in reverse order). */
   for (i = ~~(NPL-1); i >= 0; i--) {

   /* Argument and functions. */
      arg = ((xpl[i][0]  * al   +
                 xpl[i][1]  * af   +
                 xpl[i][2]  * ad   +
                 xpl[i][3] * aom  +
                 xpl[i][4] * alme +
                 xpl[i][5] * alve +
                 xpl[i][6] * alea +
                 xpl[i][7] * alma +
                 xpl[i][8] * alju +
                 xpl[i][9] * alsa +
                 xpl[i][10] * alur +
                 xpl[i][11] * alne +
                 xpl[i][12] * apa) % (ERFA_D2PI));
      sarg = Math.sin(arg);
      carg = Math.cos(arg);

   /* Term. */
      dp += xpl[i][13] * sarg + xpl[i][14] * carg;
      de += xpl[i][15] * sarg + xpl[i][16] * carg;

   }

/* Convert from 0.1 microarcsec units to radians. */
   dpsipl = dp * U2R;
   depspl = de * U2R;

/* ------- */
/* RESULTS */
/* ------- */

/* Add luni-solar and planetary components. */
   dpsi = dpsils + dpsipl;
   deps = depsls + depspl;

   return [dpsi, deps];

}
;
function eraFal03(t)
/*
**  - - - - - - - - -
**   e r a F a l 0 3
**  - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean anomaly of the Moon.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    l, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a;


/* Mean anomaly of the Moon (IERS Conventions 2003). */
   a = ((485868.249036  +
             t * ( 1717915923.2178 +
             t * (         31.8792 +
             t * (          0.051635 +
             t * (        - 0.00024470 ) ) ) )) % (ERFA_TURNAS)) * ERFA_DAS2R;

   return a;

}
;
function eraFaf03(t)
/*
**  - - - - - - - - -
**   e r a F a f 0 3
**  - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of the Moon minus mean longitude of the ascending
**  node.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    F, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a;


/* Mean longitude of the Moon minus that of the ascending node */
/* (IERS Conventions 2003).                                    */
   a = ((335779.526232 +
             t * ( 1739527262.8478 +
             t * (       - 12.7512 +
             t * (        - 0.001037 +
             t * (          0.00000417 ) ) ) )) % (ERFA_TURNAS)) * ERFA_DAS2R;

   return a;


}
;
function eraFaom03(t)
/*
**  - - - - - - - - - -
**   e r a F a o m 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of the Moon's ascending node.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    Omega, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a;


/* Mean longitude of the Moon's ascending node */
/* (IERS Conventions 2003).                    */
   a = ((450160.398036 +
             t * ( - 6962890.5431 +
             t * (         7.4722 +
             t * (         0.007702 +
             t * (       - 0.00005939 ) ) ) )) % (ERFA_TURNAS)) * ERFA_DAS2R;

   return a;

}
;
function eraFapa03(t)
/*
**  - - - - - - - - - -
**   e r a F a p a 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  general accumulated precession in longitude.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    general precession in longitude, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003).  It
**     is taken from Kinoshita & Souchay (1990) and comes originally
**     from Lieske et al. (1977).
**
**  References:
**
**     Kinoshita, H. and Souchay J. 1990, Celest.Mech. and Dyn.Astron.
**     48, 187
**
**     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
**     Astron.Astrophys. 58, 1-16
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a;


/* General accumulated precession in longitude. */
   a = (0.024381750 + 0.00000538691 * t) * t;

   return a;

}
;
function eraFame03(t)
/*
**  - - - - - - - - - -
**   e r a F a m e 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Mercury.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Mercury, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a;


/* Mean longitude of Mercury (IERS Conventions 2003). */
   a = ((4.402608842 + 2608.7903141574 * t) % (ERFA_D2PI));

   return a;

}
;
function eraFave03(t)
/*
**  - - - - - - - - - -
**   e r a F a v e 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Venus.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Venus, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a;


/* Mean longitude of Venus (IERS Conventions 2003). */
   a = ((3.176146697 + 1021.3285546211 * t) % (ERFA_D2PI));

   return a;

}
;
function eraFae03(t)
/*
**  - - - - - - - - -
**   e r a F a e 0 3
**  - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Earth.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Earth, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a;


/* Mean longitude of Earth (IERS Conventions 2003). */
   a = ((1.753470314 + 628.3075849991 * t) % (ERFA_D2PI));

   return a;

}
;
function eraFama03(t)
/*
**  - - - - - - - - - -
**   e r a F a m a 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Mars.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Mars, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a;


/* Mean longitude of Mars (IERS Conventions 2003). */
   a = ((6.203480913 + 334.0612426700 * t) % (ERFA_D2PI));

   return a;

}
;
function eraFaju03(t)
/*
**  - - - - - - - - - -
**   e r a F a j u 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Jupiter.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Jupiter, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a;


/* Mean longitude of Jupiter (IERS Conventions 2003). */
   a = ((0.599546497 + 52.9690962641 * t) % (ERFA_D2PI));

   return a;

}
;
function eraFasa03(t)
/*
**  - - - - - - - - - -
**   e r a F a s a 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Saturn.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Saturn, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a;


/* Mean longitude of Saturn (IERS Conventions 2003). */
   a = ((0.874016757 + 21.3299104960 * t) % (ERFA_D2PI));

   return a;

}
;
function eraFaur03(t)
/*
**  - - - - - - - - - -
**   e r a F a u r 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Uranus.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned  (function value):
**           double    mean longitude of Uranus, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is adapted from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a;


/* Mean longitude of Uranus (IERS Conventions 2003). */
   a = ((5.481293872 + 7.4781598567 * t) % (ERFA_D2PI));

   return a;

}
;
function eraEe00(date1, date2, epsa, dpsi)
/*
**  - - - - - - - -
**   e r a E e 0 0
**  - - - - - - - -
**
**  The equation of the equinoxes, compatible with IAU 2000 resolutions,
**  given the nutation in longitude and the mean obliquity.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**     epsa         double    mean obliquity (Note 2)
**     dpsi         double    nutation in longitude (Note 3)
**
**  Returned (function value):
**                  double    equation of the equinoxes (Note 4)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The obliquity, in radians, is mean of date.
**
**  3) The result, which is in radians, operates in the following sense:
**
**        Greenwich apparent ST = GMST + equation of the equinoxes
**
**  4) The result is compatible with the IAU 2000 resolutions.  For
**     further details, see IERS Conventions 2003 and Capitaine et al.
**     (2002).
**
**  Called:
**     eraEect00    equation of the equinoxes complementary terms
**
**  References:
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astronomy &
**     Astrophysics, 406, 1135-1149 (2003)
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var ee;


/* Equation of the equinoxes. */
   ee = dpsi * Math.cos(epsa) + eraEect00(date1, date2);

   return ee;

}
;
function eraEect00(date1, date2)
/*
**  - - - - - - - - - -
**   e r a E e c t 0 0
**  - - - - - - - - - -
**
**  Equation of the equinoxes complementary terms, consistent with
**  IAU 2000 resolutions.
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double   complementary terms (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The "complementary terms" are part of the equation of the
**     equinoxes (EE), classically the difference between apparent and
**     mean Sidereal Time:
**
**        GAST = GMST + EE
**
**     with:
**
**        EE = dpsi * cos(eps)
**
**     where dpsi is the nutation in longitude and eps is the obliquity
**     of date.  However, if the rotation of the Earth were constant in
**     an inertial frame the classical formulation would lead to
**     apparent irregularities in the UT1 timescale traceable to side-
**     effects of precession-nutation.  In order to eliminate these
**     effects from UT1, "complementary terms" were introduced in 1994
**     (IAU, 1994) and took effect from 1997 (Capitaine and Gontier,
**     1993):
**
**        GAST = GMST + CT + EE
**
**     By convention, the complementary terms are included as part of
**     the equation of the equinoxes rather than as part of the mean
**     Sidereal Time.  This slightly compromises the "geometrical"
**     interpretation of mean sidereal time but is otherwise
**     inconsequential.
**
**     The present function computes CT in the above expression,
**     compatible with IAU 2000 resolutions (Capitaine et al., 2002, and
**     IERS Conventions 2003).
**
**  Called:
**     eraFal03     mean anomaly of the Moon
**     eraFalp03    mean anomaly of the Sun
**     eraFaf03     mean argument of the latitude of the Moon
**     eraFad03     mean elongation of the Moon from the Sun
**     eraFaom03    mean longitude of the Moon's ascending node
**     eraFave03    mean longitude of Venus
**     eraFae03     mean longitude of Earth
**     eraFapa03    general accumulated precession in longitude
**
**  References:
**
**     Capitaine, N. & Gontier, A.-M., Astron.Astrophys., 275,
**     645-650 (1993)
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astron.Astrophys., 406,
**     1135-1149 (2003)
**
**     IAU Resolution C7, Recommendation 3 (1994)
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
/* Time since J2000.0, in Julian centuries */
   var t;

/* Miscellaneous */
   var i, j;
   var a, s0, s1;

/* Fundamental arguments */
   var fa = [];

/* Returned value. */
   var eect;

/* ----------------------------------------- */
/* The series for the EE complementary terms */
/* ----------------------------------------- */

   /* Terms of order t^0 */
   var e0 = [

   /* 1-10 */
      [ 0,  0,  0,  0,  1,  0,  0,  0,  2640.96e-6, -0.39e-6 ],
      [ 0,  0,  0,  0,  2,  0,  0,  0,    63.52e-6, -0.02e-6 ],
      [ 0,  0,  2, -2,  3,  0,  0,  0,    11.75e-6,  0.01e-6 ],
      [ 0,  0,  2, -2,  1,  0,  0,  0,    11.21e-6,  0.01e-6 ],
      [ 0,  0,  2, -2,  2,  0,  0,  0,    -4.55e-6,  0.00e-6 ],
      [ 0,  0,  2,  0,  3,  0,  0,  0,     2.02e-6,  0.00e-6 ],
      [ 0,  0,  2,  0,  1,  0,  0,  0,     1.98e-6,  0.00e-6 ],
      [ 0,  0,  0,  0,  3,  0,  0,  0,    -1.72e-6,  0.00e-6 ],
      [ 0,  1,  0,  0,  1,  0,  0,  0,    -1.41e-6, -0.01e-6 ],
      [ 0,  1,  0,  0, -1,  0,  0,  0,    -1.26e-6, -0.01e-6 ],

   /* 11-20 */
      [ 1,  0,  0,  0, -1,  0,  0,  0,    -0.63e-6,  0.00e-6 ],
      [ 1,  0,  0,  0,  1,  0,  0,  0,    -0.63e-6,  0.00e-6 ],
      [ 0,  1,  2, -2,  3,  0,  0,  0,     0.46e-6,  0.00e-6 ],
      [ 0,  1,  2, -2,  1,  0,  0,  0,     0.45e-6,  0.00e-6 ],
      [ 0,  0,  4, -4,  4,  0,  0,  0,     0.36e-6,  0.00e-6 ],
      [ 0,  0,  1, -1,  1, -8, 12,  0,    -0.24e-6, -0.12e-6 ],
      [ 0,  0,  2,  0,  0,  0,  0,  0,     0.32e-6,  0.00e-6 ],
      [ 0,  0,  2,  0,  2,  0,  0,  0,     0.28e-6,  0.00e-6 ],
      [ 1,  0,  2,  0,  3,  0,  0,  0,     0.27e-6,  0.00e-6 ],
      [ 1,  0,  2,  0,  1,  0,  0,  0,     0.26e-6,  0.00e-6 ],

   /* 21-30 */
      [ 0,  0,  2, -2,  0,  0,  0,  0,    -0.21e-6,  0.00e-6 ],
      [ 0,  1, -2,  2, -3,  0,  0,  0,     0.19e-6,  0.00e-6 ],
      [ 0,  1, -2,  2, -1,  0,  0,  0,     0.18e-6,  0.00e-6 ],
      [ 0,  0,  0,  0,  0,  8,-13, -1,    -0.10e-6,  0.05e-6 ],
      [ 0,  0,  0,  2,  0,  0,  0,  0,     0.15e-6,  0.00e-6 ],
      [ 2,  0, -2,  0, -1,  0,  0,  0,    -0.14e-6,  0.00e-6 ],
      [ 1,  0,  0, -2,  1,  0,  0,  0,     0.14e-6,  0.00e-6 ],
      [ 0,  1,  2, -2,  2,  0,  0,  0,    -0.14e-6,  0.00e-6 ],
      [ 1,  0,  0, -2, -1,  0,  0,  0,     0.14e-6,  0.00e-6 ],
      [ 0,  0,  4, -2,  4,  0,  0,  0,     0.13e-6,  0.00e-6 ],

   /* 31-33 */
      [ 0,  0,  2, -2,  4,  0,  0,  0,    -0.11e-6,  0.00e-6 ],
      [ 1,  0, -2,  0, -3,  0,  0,  0,     0.11e-6,  0.00e-6 ],
      [ 1,  0, -2,  0, -1,  0,  0,  0,     0.11e-6,  0.00e-6 ]
   ];


/* Terms of order t^1 */
   var e1 = [
      [ 0,  0,  0,  0,  1,  0,  0,  0,     -0.87e-6,  0.00e-6 ]
   ];


/* Number of terms in the series */
   var NE0 = e0.length;
   var NE1 = e1.length;

/* ------------------------------------------------------------------ */

/* Interval between fundamental epoch J2000.0 and current date (JC). */
   t = ((date1 - ERFA_DJ00) + date2) / ERFA_DJC;

/* Fundamental Arguments (from IERS Conventions 2003) */

/* Mean anomaly of the Moon. */
   fa[0] = eraFal03(t);

/* Mean anomaly of the Sun. */
   fa[1] = eraFalp03(t);

/* Mean longitude of the Moon minus that of the ascending node. */
   fa[2] = eraFaf03(t);

/* Mean elongation of the Moon from the Sun. */
   fa[3] = eraFad03(t);

/* Mean longitude of the ascending node of the Moon. */
   fa[4] = eraFaom03(t);

/* Mean longitude of Venus. */
   fa[5] = eraFave03(t);

/* Mean longitude of Earth. */
   fa[6] = eraFae03(t);

/* General precession in longitude. */
   fa[7] = eraFapa03(t);

/* Evaluate the EE complementary terms. */
   s0 = 0.0;
   s1 = 0.0;

   for (i = ~~(NE0-1); i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < 8; j++) {
         a += (e0[i][j]) * fa[j];
      }
      s0 += e0[i][8] * Math.sin(a) + e0[i][9] * Math.cos(a);
   }

   for (i = ~~(NE1-1); i >= 0; i--) {
      a = 0.0;
      for (j = 0; j < 8; j++) {
         a += (e1[i][j]) * fa[j];
      }
      s1 += e1[i][8] * Math.sin(a) + e1[i][9] * Math.cos(a);
   }

   eect = (s0 + s1 * t ) * ERFA_DAS2R;

   return eect;

}
;
function eraFalp03(t)
/*
**  - - - - - - - - - -
**   e r a F a l p 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean anomaly of the Sun.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    l', radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a;


/* Mean anomaly of the Sun (IERS Conventions 2003). */
   a = ((1287104.793048 +
             t * ( 129596581.0481 +
             t * (       - 0.5532 +
             t * (         0.000136 +
             t * (       - 0.00001149 ) ) ) )) % (ERFA_TURNAS)) * ERFA_DAS2R;

   return a;

}
;
function eraFad03(t)
/*
**  - - - - - - - - -
**   e r a F a d 0 3
**  - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean elongation of the Moon from the Sun.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    D, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var a;


/* Mean elongation of the Moon from the Sun (IERS Conventions 2003). */
   a = ((1072260.703692 +
             t * ( 1602961601.2090 +
             t * (        - 6.3706 +
             t * (          0.006593 +
             t * (        - 0.00003169 ) ) ) )) % (ERFA_TURNAS)) * ERFA_DAS2R;

   return a;

}
;
function eraRxp(r, p)
/*
**  - - - - - - -
**   e r a R x p
**  - - - - - - -
**
**  Multiply a p-vector by an r-matrix.
**
**  Given:
**     r        double[3][3]    r-matrix
**     p        double[3]       p-vector
**
**  Returned:
**     rp       double[3]       r * p
**
**  Note:
**     It is permissible for p and rp to be the same array.
**
**  Called:
**     eraCp        copy p-vector
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var rp = [0, 0, 0];;


   var w, wrp = [];
   var i, j;


/* Matrix r * vector p. */
   for (j = 0; j < 3; j++) {
       w = 0.0;
       for (i = 0; i < 3; i++) {
           w += r[j][i] * p[i];
       }
       wrp[j] = w;
   }

/* Return the result. */
   rp = eraCp(wrp);

   return rp;

}
;
function eraNut06a(date1, date2)
/*
**  - - - - - - - - - -
**   e r a N u t 0 6 a
**  - - - - - - - - - -
**
**  IAU 2000A nutation with adjustments to match the IAU 2006
**  precession.
**
**  Given:
**     date1,date2   double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi,deps     double   nutation, luni-solar + planetary (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The nutation components in longitude and obliquity are in radians
**     and with respect to the mean equinox and ecliptic of date,
**     IAU 2006 precession model (Hilton et al. 2006, Capitaine et al.
**     2005).
**
**  3) The function first computes the IAU 2000A nutation, then applies
**     adjustments for (i) the consequences of the change in obliquity
**     from the IAU 1980 ecliptic to the IAU 2006 ecliptic and (ii) the
**     secular variation in the Earth's dynamical form factor J2.
**
**  4) The present function provides classical nutation, complementing
**     the IAU 2000 frame bias and IAU 2006 precession.  It delivers a
**     pole which is at current epochs accurate to a few tens of
**     microarcseconds, apart from the free core nutation.
**
**  Called:
**     eraNut00a    nutation, IAU 2000A
**
**  References:
**
**     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
**     Astron.Astrophys. 387, 700
**
**     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
**     Astron.Astrophys. 58, 1-16
**
**     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
**     107, B4.  The MHB_2000 code itself was obtained on 9th September
**     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**     Wallace, P.T., "Software for Implementing the IAU 2000
**     Resolutions", in IERS Workshop 5.1 (2002)
**
**  Copyright (C) 2013-2019, NumFOCUS Foundation.
**  Derived, with permission, from the SOFA library.  See notes at end of file.
*/
{
   var dpsi = 0.0;;
   var deps = 0.0;;
   var _rv1;

   var t, fj2, dp, de;


/* Interval between fundamental date J2000.0 and given date (JC). */
   t = ((date1 - ERFA_DJ00) + date2) / ERFA_DJC;

/* Factor correcting for secular variation of J2. */
   fj2 = -2.7774e-6 * t;

/* Obtain IAU 2000A nutation. */
   (_rv1 = eraNut00a(date1, date2))[0];
   dp = _rv1[0];
   de = _rv1[1];

/* Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5). */
   dpsi = dp + dp * (0.4697e-6 + fj2);
   deps = de + de * fj2;

   return [dpsi, deps];

}
;

// export functions
//era.a2af    = eraA2af;
//era.a2tf    = eraA2tf;
era.ab      = eraAb;
//era.af2a    = eraAf2a;
era.anp     = eraAnp;
//era.anpm    = eraAnpm;
//era.apcg    = eraApcg;
//era.apcg13  = eraApcg13;
//era.apci    = eraApci;
//era.apci13  = eraApci13;
//era.apco    = eraApco;
//era.apco13  = eraApco13;
//era.apcs    = eraApcs;
//era.apcs13  = eraApcs13;
//era.aper    = eraAper;
//era.aper13  = eraAper13;
//era.apio    = eraApio;
//era.apio13  = eraApio13;
//era.atci13  = eraAtci13;
//era.atciq   = eraAtciq;
//era.atciqn  = eraAtciqn;
//era.atciqz  = eraAtciqz;
//era.atco13  = eraAtco13;
//era.atic13  = eraAtic13;
//era.aticq   = eraAticq;
//era.aticqn  = eraAticqn;
//era.atio13  = eraAtio13;
//era.atioq   = eraAtioq;
//era.atoc13  = eraAtoc13;
//era.atoi13  = eraAtoi13;
//era.atoiq   = eraAtoiq;
//era.bi00    = eraBi00;
//era.bp00    = eraBp00;
//era.bp06    = eraBp06;
//era.bpn2xy  = eraBpn2xy;
//era.c2i00a  = eraC2i00a;
//era.c2i00b  = eraC2i00b;
//era.c2i06a  = eraC2i06a;
//era.c2ibpn  = eraC2ibpn;
//era.c2ixy   = eraC2ixy;
//era.c2ixys  = eraC2ixys;
//era.c2s     = eraC2s;
//era.c2t00a  = eraC2t00a;
//era.c2t00b  = eraC2t00b;
//era.c2t06a  = eraC2t06a;
//era.c2tcio  = eraC2tcio;
//era.c2teqx  = eraC2teqx;
//era.c2tpe   = eraC2tpe;
//era.c2txy   = eraC2txy;
era.cal2jd  = eraCal2jd;
era.cp      = eraCp;
//era.cpv     = eraCpv;
era.cr      = eraCr;
//era.d2dtf   = eraD2dtf;
//era.d2tf    = eraD2tf;
era.dat     = eraDat;
//era.dtdb    = eraDtdb;
//era.dtf2d   = eraDtf2d;
//era.ee00    = eraEe00;
//era.ee00a   = eraEe00a;
//era.ee00b   = eraEe00b;
//era.ee06a   = eraEe06a;
//era.eect00  = eraEect00;
era.eform   = eraEform;
//era.eo06a   = eraEo06a;
//era.eors    = eraEors;
//era.epb     = eraEpb;
//era.epb2jd  = eraEpb2jd;
//era.epj     = eraEpj;
//era.epj2jd  = eraEpj2jd;
//era.epv00   = eraEpv00;
//era.eqeq94  = eraEqeq94;
era.era00   = eraEra00;
//era.fad03   = eraFad03;
//era.fae03   = eraFae03;
//era.faf03   = eraFaf03;
//era.faju03  = eraFaju03;
//era.fal03   = eraFal03;
//era.falp03  = eraFalp03;
//era.fama03  = eraFama03;
//era.fame03  = eraFame03;
//era.fane03  = eraFane03;
//era.faom03  = eraFaom03;
//era.fapa03  = eraFapa03;
//era.fasa03  = eraFasa03;
//era.faur03  = eraFaur03;
//era.fave03  = eraFave03;
//era.fk52h   = eraFk52h;
//era.fk5hip  = eraFk5hip;
//era.fk5hz   = eraFk5hz;
era.fw2m    = eraFw2m;
//era.fw2xy   = eraFw2xy;
//era.gc2gd   = eraGc2gd;
//era.gc2gde  = eraGc2gde;
era.gd2gc   = eraGd2gc;
era.gd2gce  = eraGd2gce;
era.gmst00  = eraGmst00;
//era.gmst06  = eraGmst06;
//era.gmst82  = eraGmst82;
era.gst00a  = eraGst00a;
//era.gst00b  = eraGst00b;
//era.gst06   = eraGst06;
//era.gst06a  = eraGst06a;
//era.gst94   = eraGst94;
//era.h2fk5   = eraH2fk5;
//era.hfk5z   = eraHfk5z;
era.ir      = eraIr;
era.jd2cal  = eraJd2cal;
//era.jdcalf  = eraJdcalf;
//era.ld      = eraLd;
//era.ldn     = eraLdn;
//era.ldsun   = eraLdsun;
era.num00a  = eraNum00a;
//era.num00b  = eraNum00b;
//era.num06a  = eraNum06a;
//era.numat   = eraNumat;
//era.nut00a  = eraNut00a;
//era.nut00b  = eraNut00b;
//era.nut06a  = eraNut06a;
//era.nut80   = eraNut80;
//era.nutm80  = eraNutm80;
era.obl06   = eraObl06;
//era.obl80   = eraObl80;
//era.p06e    = eraP06e;
//era.p2pv    = eraP2pv;
//era.p2s     = eraP2s;
//era.pap     = eraPap;
//era.pas     = eraPas;
//era.pb06    = eraPb06;
era.pdp     = eraPdp;
era.pfw06   = eraPfw06;
//era.plan94  = eraPlan94;
//era.pm      = eraPm;
//era.pmat00  = eraPmat00;
//era.pmat06  = eraPmat06;
//era.pmat76  = eraPmat76;
//era.pmp     = eraPmp;
//era.pmpx    = eraPmpx;
//era.pmsafe  = eraPmsafe;
//era.pn      = eraPn;
//era.pn00    = eraPn00;
//era.pn00a   = eraPn00a;
//era.pn00b   = eraPn00b;
//era.pn06    = eraPn06;
//era.pn06a   = eraPn06a;
//era.pnm00a  = eraPnm00a;
//era.pnm00b  = eraPnm00b;
era.pnm06a  = eraPnm06a;
//era.pnm80   = eraPnm80;
era.pom00   = eraPom00;
//era.ppp     = eraPpp;
//era.ppsp    = eraPpsp;
//era.pr00    = eraPr00;
//era.prec76  = eraPrec76;
//era.pv2p    = eraPv2p;
//era.pv2s    = eraPv2s;
//era.pvdpv   = eraPvdpv;
//era.pvm     = eraPvm;
//era.pvmpv   = eraPvmpv;
//era.pvppv   = eraPvppv;
//era.pvstar  = eraPvstar;
era.pvtob   = eraPvtob;
//era.pvu     = eraPvu;
//era.pvup    = eraPvup;
//era.pvxpv   = eraPvxpv;
//era.pxp     = eraPxp;
//era.refco   = eraRefco;
//era.rm2v    = eraRm2v;
//era.rv2m    = eraRv2m;
era.rx      = eraRx;
era.rxp     = eraRxp;
//era.rxpv    = eraRxpv;
//era.rxr     = eraRxr;
era.ry      = eraRy;
era.rz      = eraRz;
//era.s00     = eraS00;
//era.s00a    = eraS00a;
//era.s00b    = eraS00b;
//era.s06     = eraS06;
//era.s06a    = eraS06a;
//era.s2c     = eraS2c;
//era.s2p     = eraS2p;
//era.s2pv    = eraS2pv;
//era.s2xpv   = eraS2xpv;
//era.sepp    = eraSepp;
//era.seps    = eraSeps;
//era.sp00    = eraSp00;
//era.starpm  = eraStarpm;
//era.starpv  = eraStarpv;
//era.sxp     = eraSxp;
//era.sxpv    = eraSxpv;
era.taitt   = eraTaitt;
era.taiut1  = eraTaiut1;
//era.taiutc  = eraTaiutc;
//era.tcbtdb  = eraTcbtdb;
//era.tcgtt   = eraTcgtt;
//era.tdbtcb  = eraTdbtcb;
//era.tdbtt   = eraTdbtt;
//era.tf2a    = eraTf2a;
//era.tf2d    = eraTf2d;
era.tr      = eraTr;
era.trxp    = eraTrxp;
//era.trxpv   = eraTrxpv;
//era.tttai   = eraTttai;
//era.tttcg   = eraTttcg;
//era.tttdb   = eraTttdb;
//era.ttut1   = eraTtut1;
//era.ut1tai  = eraUt1tai;
//era.ut1tt   = eraUt1tt;
//era.ut1utc  = eraUt1utc;
era.utctai  = eraUtctai;
era.utcut1  = eraUtcut1;
//era.xy06    = eraXy06;
//era.xys00a  = eraXys00a;
//era.xys00b  = eraXys00b;
//era.xys06a  = eraXys06a;
era.zp      = eraZp;
//era.zpv     = eraZpv;
//era.zr      = eraZr;
//era.ltpequ  = eraLtpequ;
//era.ltpecl  = eraLtpecl;
//era.ltpb    = eraLtpb;
//era.ltp     = eraLtp;
//era.lteqec  = eraLteqec;
//era.ltecm   = eraLtecm;
//era.lteceq  = eraLteceq;
//era.icrs2g  = eraIcrs2g;
//era.g2icrs  = eraG2icrs;
//era.eqec06  = eraEqec06;
//era.ecm06   = eraEcm06;
//era.eceq06  = eraEceq06;
//// added with release 14
//era.ae2hd   = eraAe2hd;
//era.hd2ae   = eraHd2ae;
//era.hd2pa   = eraHd2pa;
//era.tpors   = eraTpors;
//era.tporv   = eraTporv;
//era.tpsts   = eraTpsts;
//era.tpstv   = eraTpstv;
//era.tpxes   = eraTpxes;
//era.tpxev   = eraTpxev;
//// added with release 15
//era.fk425   = eraFk425;
//era.fk45z   = eraFk45z;
//era.fk524   = eraFk524;
//era.fk54z   = eraFk54z;
			;
})(ERFA)
/* crc: 8271EBE5E498EC6B6FD45EEC44E0A4DA */
