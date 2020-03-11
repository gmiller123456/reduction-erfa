"use strict";

const au= 149597870.691;
const T0 = 2451545.00000000;
const ASEC360 = 1296000.0;
const ASEC2RAD = 4.848136811095359935899141e-6;
const TWOPI=Math.PI*2;
const DEG2RAD=Math.PI/180.0;
const ERAD = 6378136.6;
const F = 0.003352819697896;
const ANGVEL = 7.2921150e-5;
const AU_KM = 1.4959787069098932e+8;
const C_AUDAY = 173.1446326846693;
const AU = 1.4959787069098932e+11;
const RAD2DEG = 57.295779513082321;
const C = 299792458.0;

function getTopoCentric(jd_utc,observer,bodyNumber){
	const jd_tai=ERFA.utctai(jd_utc[1],jd_utc[2]);
	const jd_tt=ERFA.taitt(jd_tai[1],jd_tai[2]);
	const jd_ut1=ERFA.utcut1(jd_utc[1],jd_utc[2],0.0);
	//const jd_ut1=IAU.utcut1(jd_utc[1],jd_utc[2],0);

	//TODO: compute TDB
	const TDB_TT=0;
	//const jd_tdb=IAU.tttdb(jd_tt[1],jd_tt[2],TDB_TT);
	//const jd_tdb=[0,2454580.942629462,0];
	const jd_tdb=jd_tt;

	const earth=getBodyPV(2,jd_tdb[1]+jd_tdb[2]);
	const gast=ERFA.gst00a(jd_ut1[1],jd_ut1[2],jd_tt[1],jd_tt[2]);
	const geopv=ERFA.pvtob(observer.longitude*DEG2RAD,observer.latitude*DEG2RAD,0,0,0,0,gast); //meters, m/s

	const rnpb=ERFA.pnm06a(jd_tt[1], jd_tt[2]);
	const rnpb0=ERFA.tr(rnpb);
	const geop=ERFA.rxp(rnpb0,geopv[0]).map(x=>x/AU);
	const geov=ERFA.rxp(rnpb0,geopv[1]).map(x=>x/AU*60*60*24);

	for(let i=0;i<3;i++){
		geop[i]=earth[i]+geop[i];
		geov[i]=earth[i+3]+geov[i];
	}

	const pv1=getBodyPV(bodyNumber,jd_tdb[1]+jd_tdb[2]);

	const pos2=new Array();
	for(let i=0;i<3;i++){
		pos2[i]=pv1[i]-geop[i];
	}

   let lighttime = Math.sqrt (pos2[0] * pos2[0] + pos2[1] * pos2[1] + pos2[2] * pos2[2]) / C_AUDAY ;

	let pos3=getBodyPV(bodyNumber,jd_tdb[1]+jd_tdb[2]-lighttime);
	for(let i=0;i<3;i++){
		pos3[i]=pos3[i]-geop[i];
	}

	//TODO: gavitational deflection
	//const temp=limb_angle(pos3,geopv);
	//const frlimb=temp[1];
	//let loc=1;
	//if (frlimb < 0.8) loc = 0;

	//const pos4=grav_def (jd_tdb,loc,pos3,pvb);
	const pos4=pos3;

	const pos5=aberration (pos4,geop,geov,jd_tdb);

	const j2000=toRaDec(pos3);
	const pos8=ERFA.rxp(rnpb,pos5);

	//const radec=vector2radec(pos8);
	const radec=toRaDec(pos8);
	//console.log(radec);
	//console.log(radec[0]*15);

	const altaz=convertRaDecToAltAz(gast,observer.latitude*DEG2RAD,observer.longitude*DEG2RAD,radec[0]*15.0*DEG2RAD,radec[1]*DEG2RAD);
	//console.log(`Az:${altaz[0]*RAD2DEG} Alt:${altaz[1]*RAD2DEG} Zd:${90-altaz[1]*RAD2DEG}`);

	return [radec[0],radec[1],altaz[0]*RAD2DEG,altaz[1]*RAD2DEG,j2000[0],j2000[1],gast,jd_tt[1]+jd_tt[2]-jd_utc[1]-jd_utc[2]];
}

function aberration(body,observerPosition,observerVelocity,jd_tdb){
	const bodyDist=Math.sqrt(body[0]*body[0]+body[1]*body[1]+body[2]*body[2]);
	const bodyUnit=body.map(x=>x/bodyDist);

	let sun=getBodyPV(10,jd_tdb[1]+jd_tdb[2]);
	for(let i=0;i<3;i++){
		sun[i]=sun[i]-observerPosition[i];
	}
	const sunDistance=Math.sqrt(sun[0]*sun[0]+sun[1]*sun[1]+sun[2]*sun[2]);
	const vel2=observerVelocity.map(x=>x/C_AUDAY);

	let vlen=Math.sqrt(observerVelocity[0]*observerVelocity[0]+observerVelocity[1]*observerVelocity[1]+observerVelocity[2]*observerVelocity[2]);
	vlen=vlen/C_AUDAY;
	const bm1=Math.sqrt(1-vlen*vlen);

	return ERFA.ab(bodyUnit,vel2,sunDistance,bm1);
}

function toRaDec(xyz){
	let t = new Array();
	t[2] = Math.sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);
	t[1] = Math.acos(xyz[2] / t[2]);
	t[0] = Math.atan2(xyz[1], xyz[0]);

	if(t[0]<0){t[0]+=2*Math.PI;}
	if(t[1]<0){t[1]+=2*Math.PI;}

	t[1]=90-t[1]*180/Math.PI; //Return Declination in degrees
	t[0]=t[0]*180/Math.PI/15; //Return Right Ascension in hours

	return t;
}

function convertRaDecToAltAz(gast,lat,lon,ra,dec){
	let h=gast + lon - ra;
	
	const sina=Math.sin(dec)*Math.sin(lat)+Math.cos(dec)*Math.cos(h)*Math.cos(lat);
	const a=Math.asin(sina);

	const cosAz=(Math.sin(dec)*Math.cos(lat)-Math.cos(dec)*Math.cos(h)*Math.sin(lat))/Math.cos(a);
	let Az=Math.acos(cosAz);

	if(Math.sin(h)>0){Az=2.0*Math.PI-Az;}

	let t=new Array();
	t[0]=Az;
	t[1]=a;
	
	return t;
}

//TODO: Actually implement
function tt2tdb(tt) {
	return tt;
}

//Performs the rotation from ecliptic coordinates to J2000 coordinates for the given vector x
function rotvsop2J2000(x){
	/* From VSOP87.doc
	  X        +1.000000000000  +0.000000440360  -0.000000190919   X
	  Y     =  -0.000000479966  +0.917482137087  -0.397776982902   Y
	  Z FK5     0.000000000000  +0.397776982902  +0.917482137087   Z VSOP87A
	*/
	let t = new Array();
	t[0] = x[0] + x[1] * 0.000000440360 + x[2] * -0.000000190919;
	t[1] = x[0] * -0.000000479966 + x[1] * 0.917482137087 + x[2] * -0.397776982902;
	t[2] = x[1] * 0.397776982902 + x[2] * 0.917482137087;

	return t;
}

function getBodyPV(body,jd_tdb){
	const AU_KM = 1.4959787069098932e+8;
	const t=(jd_tdb-2451545.0)/365250.0;

	let b;
	let v;

	switch (body){
		case 0: //Mercury
			b=vsop87a_full.getMercury(t);
			v=vsop87a_full_velocities.getMercury(t);
			break;
		case 1:
			b=vsop87a_full.getVenus(t);
			v=vsop87a_full_velocities.getVenus(t);
			break;
		case 2: //Earth
			b=vsop87a_full.getEarth(t);
			v=vsop87a_full_velocities.getEarth(t);
			break;
		case 3:
			b=vsop87a_full.getMars(t);
			v=vsop87a_full_velocities.getMars(t);
			break;
		case 4:
			b=vsop87a_full.getJupiter(t);
			v=vsop87a_full_velocities.getJupiter(t);
			break;
		case 5:
			b=vsop87a_full.getSaturn(t);
			v=vsop87a_full_velocities.getSaturn(t);
			break;
		case 6:
			b=vsop87a_full.getUranus(t);
			v=vsop87a_full_velocities.getUranus(t);
			break;
		case 7:
			b=vsop87a_full.getNeptune(t);
			v=vsop87a_full_velocities.getNeptune(t);
			break;
		case 8: 
			//Not implemented
			//b=vsop87a_full.getPluto(t);
			//v=vsop87a_full_velocities.getPluto(t);
			break;
		case 9: //Moon
			const e=vsop87a_full.getEarth(t);
			const c=vsop87a_full.getEmb(t);
			b=vsop87a_full.getMoon(e,c);
			const ev=vsop87a_full_velocities.getEarth(t);
			const cv=vsop87a_full_velocities.getEmb(t);
			v=vsop87a_full_velocities.getMoon(ev,cv);
			break;
		case 10: //Sun
			b=[0,0,0];
			v=[0,0,0];
			break;
	}

	b=rotvsop2J2000(b);
	v=rotvsop2J2000(v);

	for(let i=0;i<3;i++){
		b[i+3]=v[i];
	}

	return b;

}
