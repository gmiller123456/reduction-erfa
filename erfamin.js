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
}

function eraAnp(a)
{
   var w;


   w = ((a) % (ERFA_D2PI));
   if (w < 0) w += ERFA_D2PI;

   return w;

}

function eraCal2jd(iy, im, id)
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

function eraCp(p)
{
	return [p[0],p[1],p[2]];
}

function eraCr(r)
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

function eraDat(iy, im, id, fd)
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

function eraEform( n)
{
   return [ 0, 6378137.0, 1.0 / 298.257223563  ];
}

function eraEra00(dj1, dj2)
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
;
function eraFw2m(gamb, phib, psi, eps)
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
;
function eraGd2gc( n, elong, phi, height)
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
;
function eraGd2gce( a, f, elong, phi, height)
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
;
function eraGmst00(uta, utb, tta, ttb)
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

function eraIr()
{
   return [ [1,0,0], [0,1,0], [0,0,1] ];;
}

function eraJd2cal(dj1, dj2)
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
;
function eraNum00a(date1, date2)
{
	return eraPn00a(date1, date2)[6];
}

function eraObl06(date1, date2)
{

/* Interval between fundamental date J2000.0 and given date (JC). */
   const t = ((date1 - ERFA_DJ00) + date2) / ERFA_DJC;

/* Mean obliquity. */
   const eps0 = (84381.406     +
          (-46.836769    +
          ( -0.0001831   +
          (  0.00200340  +
          ( -0.000000576 +
          ( -0.0000000434) * t) * t) * t) * t) * t) * ERFA_DAS2R;

   return eps0;

}
;
function eraPdp(a, b)
{
   var w;


   w  = a[0] * b[0]
      + a[1] * b[1]
      + a[2] * b[2];

   return w;

}
;
function eraPfw06(date1, date2)
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
;
function eraPnm06a(date1, date2)
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

function eraPom00(xp, yp, sp)
{
   var rpom = [ [0,0,0], [0,0,0], [0,0,0] ];;

/* Construct the matrix. */
   rpom = eraIr();
   rpom = eraRz(sp, rpom);
   rpom = eraRy(-xp, rpom);
   rpom = eraRx(-yp, rpom);

   return rpom;

}

function eraPvtob(elong, phi, hm, xp, yp, sp, theta)
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
}

function eraRx(phi, r)
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

}

function eraRy(theta, r)
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

}

function eraRz(psi, r)
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

}

function eraTaitt(tai1, tai2)
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

}
function eraTaiut1(tai1, tai2, dta)
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

}
function eraTr(r)
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
{
	return [0,0,0];
}

function eraGst00a(uta, utb, tta, ttb)
{
   var gmst00, ee00a, gst;


   gmst00 = eraGmst00(uta, utb, tta, ttb);
   ee00a = eraEe00a(tta, ttb);
   gst = eraAnp(gmst00 + ee00a);

   return gst;

}
;

function eraEe00a(date1, date2)
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

function eraFal03(t)
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
{
   var a;


/* General accumulated precession in longitude. */
   a = (0.024381750 + 0.00000538691 * t) * t;

   return a;

}
;
function eraFame03(t)
{
   var a;


/* Mean longitude of Mercury (IERS Conventions 2003). */
   a = ((4.402608842 + 2608.7903141574 * t) % (ERFA_D2PI));

   return a;

}
;
function eraFave03(t)
{
   var a;


/* Mean longitude of Venus (IERS Conventions 2003). */
   a = ((3.176146697 + 1021.3285546211 * t) % (ERFA_D2PI));

   return a;

}
;
function eraFae03(t)
{
   var a;


/* Mean longitude of Earth (IERS Conventions 2003). */
   a = ((1.753470314 + 628.3075849991 * t) % (ERFA_D2PI));

   return a;

}
;
function eraFama03(t)
{
   var a;


/* Mean longitude of Mars (IERS Conventions 2003). */
   a = ((6.203480913 + 334.0612426700 * t) % (ERFA_D2PI));

   return a;

}
;
function eraFaju03(t)
{
   var a;


/* Mean longitude of Jupiter (IERS Conventions 2003). */
   a = ((0.599546497 + 52.9690962641 * t) % (ERFA_D2PI));

   return a;

}
;
function eraFasa03(t)
{
   var a;


/* Mean longitude of Saturn (IERS Conventions 2003). */
   a = ((0.874016757 + 21.3299104960 * t) % (ERFA_D2PI));

   return a;

}
;
function eraFaur03(t)
{
   var a;


/* Mean longitude of Uranus (IERS Conventions 2003). */
   a = ((5.481293872 + 7.4781598567 * t) % (ERFA_D2PI));

   return a;

}
;
function eraEe00(date1, date2, epsa, dpsi)
{
   var ee;


/* Equation of the equinoxes. */
   ee = dpsi * Math.cos(epsa) + eraEect00(date1, date2);

   return ee;

}
;
function eraEect00(date1, date2)
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

function eraNut00a(date1, date2)
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


// export functions
era.ab      = eraAb;
era.anp     = eraAnp;
era.cal2jd  = eraCal2jd;
era.cp      = eraCp;
era.cr      = eraCr;
era.dat     = eraDat;
era.eform   = eraEform;
era.era00   = eraEra00;
era.fw2m    = eraFw2m;
era.gd2gc   = eraGd2gc;
era.gd2gce  = eraGd2gce;
era.gmst00  = eraGmst00;
era.gst00a  = eraGst00a;
era.ir      = eraIr;
era.jd2cal  = eraJd2cal;
era.num00a  = eraNum00a;
era.obl06   = eraObl06;
era.pdp     = eraPdp;
era.pfw06   = eraPfw06;
era.pnm06a  = eraPnm06a;
era.pom00   = eraPom00;
era.pvtob   = eraPvtob;
era.rx      = eraRx;
era.rxp     = eraRxp;
era.ry      = eraRy;
era.rz      = eraRz;
era.taitt   = eraTaitt;
era.taiut1  = eraTaiut1;
era.tr      = eraTr;
era.trxp    = eraTrxp;
era.utctai  = eraUtctai;
era.utcut1  = eraUtcut1;
era.zp      = eraZp;
			;
})(ERFA)
