#include "mdefs.h"
#include "math.h"

int scalecount(double *,int,double,int *,int *,int *);

// CONSTANTS and VARIABLES
//const double pi=3.14159265358979;
//int verb=1;
const double RTerra=6370e5; // in cm!

// GENERAL FUNCTIONS
void boxmuller(double *y1,double *y2) {
	float x1, x2, w;

	do {
		x1 = (2.0 * rand())/RAND_MAX - 1.0;
		x2 = (2.0 * rand())/RAND_MAX - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( (-2.0 * log( w ) ) / w );
	*y1 = x1 * w;
	*y2 = x2 * w;
}

double lon2posx(double lon) {
	double a1=8;
	double a0=42;
	double posx;

	lon=lon;
	if(lon>180.0) lon-=360.0;

	posx=a1*lon+a0;
	return posx;
}

double lat2posy(double lat) {
	double a4=0.000015445556;
	double a3=-0.001187709723;
	double a2=0.070791787735;
	double a1=6.526851079341;
	double a0=-242.193123554951;
	double lat2,posy;

	lat2=lat*lat;

	posy=a4*lat2*lat2+a3*lat2*lat+a2*lat2+a1*lat+a0;
	return posy;
}

double posy2lat(double posy) {
	double a4=0.04035068698576;
	double a3=-0.14932066760361;
	double a2=-1.42895092126193;
	double a1=16.73888922849881;
	double a0=30.24146555579375;
	double lat,posy2;

	posy=posy/155;
	posy2=posy*posy;

	lat=a4*posy2*posy2+a3*posy2*posy+a2*posy2+a1*posy+a0;
	return lat;
}

double mod(double a,double b) {
	int div;
	double r;

	div=int(a/b);

	r=a-(double)div*b;
	if(r<0.0) r+=b;
	return r;
}

float intpl(float *x,double lambdax,double lambday,double lambdat) {
	float ix1,ix2,ix;

	ix1=(1-lambdax)*(1-lambday)*x[0]+lambdax*(1-lambday)*x[1]+lambdax*lambday*x[2]+(1-lambdax)*lambday*x[3];
	ix2=(1-lambdax)*(1-lambday)*x[4]+lambdax*(1-lambday)*x[5]+lambdax*lambday*x[6]+(1-lambdax)*lambday*x[7];

	ix=(1-lambdat)*ix1+lambdat*ix2;
	return ix;
}

class license{
	public:
		license() {
			printf("Lamta 0.4 by Francesco d'Ovidio (francesco.dovidio@iscpif.fr).\n\nThis program is free software: you can redistribute it and/or modify\nit under the terms of the GNU General Public License as published by\nthe Free Software Foundation, either version 3 of the License, or\n(at your option) any later version.\n\nThis program is distributed in the hope that it will be useful,\nbut WITHOUT ANY WARRANTY; without even the implied warranty of\nMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\nGNU General Public License for more details.\n\nYou should have received a copy of the GNU General Public License\nalong with this program.  If not, see <http://www.gnu.org/licenses/>.\n\n");
		}
} License;

// THE FIELD FAMILY 
class field{
	public:
		field() {default_par();}
		int datatype; //1 is U and V in deg./sec, 2 for U and V in cm/sec
		int disttype; //1 is Euclidean, 2 is over the sphere.
		int gridtype; //0 is flat, 1 is sphere, regular spacing in deg., 2 is sphere, Mercator grid
		int noise; //0 no noise, not zero noisy velocity field
		double noiselevel;
		int freeze; //freeze<>0, to freeze the velocity field at time frozentime.
		double frozentime;

		void freezewarning() {
			if(freeze) printf("WARNING: time frozen at t=%lf. Use unfreezetime() to release.\n",frozentime);
			else printf("Time is not frozen.\n");
		}

		void noisewarning() {
			if(noise) printf("WARNING: Noisy velocity field (sigma=%lf). Use field_noise(0,0); to remove noise.\n",noiselevel);
			// 16 Oct 2012 ======================================
			//else printf("Time is not frozen.\n");
			else printf("Velocity field set without noise.\n");
			//===================================================
		}

		int xyp(double t,double x,double y, double *xp) {

			if(disttype==2) {
				while(x<-180.0) x+=360.0;
				while(x>180.0) x-=360.0;
            // Added 2019 07 22
            // To avoid issues at North pole
            while(y>90.0) y=180.0-y;
            while(y<-90.0) y=-180.0-y;
			}

			if(freeze) xypt(frozentime,x,y,xp);
			else xypt(t,x,y,xp);

			return 0;
		}

		virtual int xypt(double t,double x,double y, double *xp) {
			xp[0]=0.0;
			xp[1]=0.0;
			return 0;
		}

		virtual int rawdata_out(long indx, double *U,double *V) {
			*U=0;
			*V=0;
			printf("No data in memory! (Field not initialised or field is a function.)");
			return 0;
		}

		virtual void print_par() {
			printf("Empty field is active.\n");
		}

		virtual void default_par() {
			datatype=1;//1 is U and V in deg./sec, 2 for U and V in cm/sec
			noise=0;
			freeze=0;
			datatype=1;
			disttype=1;
			gridtype=1;

		}

		virtual void set_par(double *) {
		}

} Genericfield, *pField=&Genericfield;

class lut_field: public field {
	public:
		lut_field() {default_par();}
		double xi,xf,yi,yf,ti,tf,latLUT[8192],invdlat;
		float *datax,*datay;
		int numx,numy,numt,memoryisallocated,numlatLUT;

		void default_par() {
			printf("Lut created.\n");
			xi=xf=yi=yf=ti=tf=0.0;
			datax=NULL;
			datay=NULL;
			memoryisallocated=0;
			numx=numy=numt=0;
			datatype=1;
			disttype=2;
			numlatLUT=8192;
			invdlat=1;
		}

		void print_par() {
			printf("LUT field is active:\n xi=%f xf=%f numx=%d yi=%f yf=%f numy=%d ti=%f tf=%f numt=%d memoryisallocated=%d\n",xi,xf,numx,yi,yf,numy,ti,tf,numt,memoryisallocated);
			freezewarning();
			noisewarning();
		}

		void set_par(double *par) {
			int j;
			double xL,c,ddeg;
			xi=par[0];
			xf=par[1];
			numx=(int)par[2];
			yi=par[3];
			yf=par[4];
			numy=(int)par[5];
			ti=par[6];
			tf=par[7];
			numt=(int)par[8];

			if(memoryisallocated) {
				delete datax;
				delete datay;
				printf("U and V freed in LUT field.\n");
				memoryisallocated=0;
			}

			ddeg=(xf-xi)/(numx-1);
			c=1/(ddeg/180.0*Pi);

			for(j=0;j<numlatLUT;j++) {
				xL=((double)j)/numlatLUT*(yf-yi)+yi;
				xL=xL/180.0*Pi;
				latLUT[j]=c*(log(cos(xL/2)+sin(xL/2))-log(cos(xL/2)-sin(xL/2)));
				//printf("j:%d\txL:%lf\t%lf\n",j,xL,latLUT[j]);
			}

			invdlat=1.0/(yf-yi)*(numlatLUT-1);
		}

		double clatLUT(double lat) {
			int latp1;
			latp1=(int)((lat-yi)*invdlat);
			//printf("latp1:%d\n",latp1);
			if(latp1<0) latp1=0;
			if(latp1>numlatLUT) latp1=numlatLUT;
			return latLUT[latp1];
		}


		double mercatorpos(double lat) {
			double Fphi0,Fphii,c,ddeg,ilat;
			Fphi0=clatLUT(yi);
			Fphii=clatLUT(lat);
			ilat=Fphii-Fphi0;
			//printf("Fphi0:%lf Fphii:%lf\n",Fphi0,Fphii);
			return ilat;
		}

		int xypt(double t,double x,double y,double *xyp) {
			double posx,posy,post,lambdax,lambday,lambdat;
			float xpi[8],ypi[8],xp1,yp1;
			int ix,iy,it,ixl,iyl,itl;
			posx=(x-xi)/(xf-xi)*(numx-1);
			ix=(int)posx;
			lambdax=posx-(double)ix;

			if(gridtype==2) {
				posy=mercatorpos(y);
				//printf("mercatorpos:%lf\n",posy);
			}
			else posy=(y-yi)/(yf-yi)*(numy-1);

			iy=(int)posy;
			lambday=posy-(double)iy;

			//if statement to have fixed time option
			//with only one time passed
			if (tf!=ti){
				post=(t-ti)/(tf-ti)*(numt-1);
				it=(int)post;
				lambdat=post-(double)it;
			}
			else it=0;

			if ((datax==NULL)||(datay==NULL)) {
				xyp[0]=0.0;
				xyp[1]=0.0;
				printf("No data!\n");
				return 1;
			}

         /*
         if((ix<0)||(ix>=numx-1)||(iy<0)||(iy>=numy-1)||(it<0)||(it>numt-1)) {
         //if((ix<0)||(ix>=numx-1)||(iy<0)||(iy>=numy-1)||(it<0)||(it>=numt-1)) {
         xyp[0]=0.0;
         xyp[1]=0.0;
         //printf("Out of limits of stored data!\n");
         printf("Something wrong with indexing!\n");
			//printf("x:%f y:%f t:%f\n",x,y,t);
			//printf("ix:%f iy:%f it:%f\n",ix,iy,it);
			//printf("numx:%d numy:%d numt:%d\n",numx,numy,numt);
			printf("\n");
         return 0;
         }
         */
         
         //The code below is needed only for non-preiodic domains
         //(beyond boundary uv are set to boundary value)
         /*
			if(ix<0) ix=0;
			if(iy<0) iy=0;
			if(it<0) it=0;
			if(ix>numx-1) ix=numx-1;
			if(iy>numy-1) iy=numy-1;
			if(it>numt-1) it=numt-1;	
         */
         

			ixl=ix+1;
			iyl=iy+1;
			itl=it+1;

         /*
         if((ixl>numx-1)||(iyl>numy-1)) {
         //if((ixl>numx-1)||(iyl>numy-1)||(itl>numt-1)) {
         xyp[0]=0.0;
         xyp[1]=0.0;
         //printf("Out of limits of stored data!\n");
         printf("Something wrong with indexing (ixl,iyl,itl)!\n");
			printf("ixl:%d iyl:%d itl:%d\n",ixl,iyl,itl);
			printf("\n");
         return 0;
         }
         */
         
         //The code below is needed only for non-preiodic domains
         //(beyond boundary uv are set to boundary value)
         /*
			if(ixl>numx-1) ixl=numx-1;
         */

         // This needed to fix north pole issue
         // (since velocities there are 0)
         // It is triggered only for the particles deployed at North pole
			if(iyl>numy-1) {
            //printf("Particle's at the North pole (ix=%d)\n",ix);
            iyl=numy-1;
         }

         // This is needed to avoid segmentatin fault 
         // on the first step of backward advection
         // (when it=numt-1)
			if(itl>numt-1) itl=numt-1;	
         

			// printf("posx:%lf ix:%d lambdax:%lf posy:%lf iy:%d lambday:%lf post:%lf it:%d	lambdat:%f\n",posx,ix,lambdax,posy,iy,lambday,post,it,lambdat);
			// printf("ixl:%d iyl:%d itl:%d\n",ixl,iyl,itl);

			xpi[0]=datax[ix+iy*numx+it*numx*numy];
			xpi[1]=datax[(ixl)+iy*numx+it*numx*numy];
			xpi[2]=datax[(ixl)+(iyl)*numx+it*numx*numy];
			xpi[3]=datax[ix+(iyl)*numx+it*numx*numy];
			// printf("xpi[0]:%f xpi[1]:%f xpi[2]:%f xpi[3]:%f\n",xpi[0],xpi[1],xpi[2],xpi[3]);

			xpi[4]=datax[ix+iy*numx+(itl)*numx*numy];
			xpi[5]=datax[(ixl)+iy*numx+(itl)*numx*numy];
			xpi[6]=datax[(ixl)+(iyl)*numx+(itl)*numx*numy];
			xpi[7]=datax[ix+(iyl)*numx+(itl)*numx*numy];
			// printf("xpi[4]:%f xpi[5]:%f xpi[6]:%f xpi[7]:%f\n",xpi[4],xpi[5],xpi[6],xpi[7]);
         //printf("Test test test\n");

			ypi[0]=datay[ix+iy*numx+it*numx*numy];
			ypi[1]=datay[(ixl)+iy*numx+it*numx*numy];
			ypi[2]=datay[(ixl)+(iyl)*numx+it*numx*numy];
			ypi[3]=datay[ix+(iyl)*numx+it*numx*numy];
			// printf("ypi[0]:%f ypi[1]:%f ypi[2]:%f ypi[3]:%f\n",ypi[0],ypi[1],ypi[2],ypi[3]);

			ypi[4]=datay[ix+iy*numx+(itl)*numx*numy];
			ypi[5]=datay[(ixl)+iy*numx+(itl)*numx*numy];
			ypi[6]=datay[(ixl)+(iyl)*numx+(itl)*numx*numy];
			ypi[7]=datay[ix+(iyl)*numx+(itl)*numx*numy];
			// printf("ypi[4]:%f ypi[5]:%f ypi[6]:%f ypi[7]:%f\n",ypi[4],ypi[5],ypi[6],ypi[7]);

			xp1=intpl(xpi,lambdax,lambday,lambdat);
			yp1=intpl(ypi,lambdax,lambday,lambdat);

			xyp[0]=(double)xp1;
			xyp[1]=(double)yp1;
         // printf("Velocities: xyp[0]=%lf xyp[1]=%lf\n",xyp[0],xyp[1]);	
			return 0;
		}
} Lut_field;

// FUNCTIONS FOR HANDLING THE OBJECTS, INTEGRATING ETC.
int select_sys(int sys_code) {
	int ret=0;
	switch(sys_code) {
		case 4: pField=&Lut_field; break;
		default: pField=&Genericfield; break;
	}
	return ret;
}

void pol2cart(double lon,double lat,double *x,double *y,double *z) {
	*z=RTerra*sin(lat/360.0*2*Pi);
	*x=RTerra*cos(lat/360.0*2*Pi)*cos(lon/360.0*2*Pi);
	*y=RTerra*cos(lat/360.0*2*Pi)*sin(lon/360.0*2*Pi);
}

void cart2pol(double x,double y,double z,double *lon,double *lat) {
	/*
	   if(z>RTerra) z=RTerra;
	   else if(z<-RTerra) z=-RTerra;
	 *lat=asin(z/RTerra)*360.0/(2*Pi);
	 printf("cart2pol:\n lat:%lf lat2:%lf\n",*lat,atan2(z,sqrt(x*x+y*y))*360/(2*Pi));
	 */
	*lat=atan2(z,sqrt(x*x+y*y))*360.0/(2*Pi);
	*lon=atan2(y,x)*360.0/(2*Pi);
	//printf("cart2pol:\nx:%lf y:%lf z:%lf lon:%lf lat:%lf test:%lf\n",x,y,z,*lon,*lat,atan2(y,x)*360.0/(2*Pi));
}

void Tdist(double x1,double y1,double x2,double y2,double *d) {
	int distswitch;
	double dx,dy,dz,dist1,dist,cx1,cx2,cy1,cy2,cz1,cz2,lat1,lon1,lat2,lon2;
	//printf("Tdist\n");
	if(pField->gridtype==0)
		distswitch=0;
	else
		distswitch=1;

	switch(distswitch) {
		case 1: {
				lon1=x1/360.0*2*Pi;
				lat1=y1/360.0*2*Pi;
				lon2=x2/360.0*2*Pi;
				lat2=y2/360.0*2*Pi;
				dist=acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1))*360.0/(2.0*Pi);
			} break;
			/*
			   case 2: {
	//printf("Tdist case 2 (Euclidean on the sphere)\n");
	pol2cart(x1,y1,&cx1,&cy1,&cz1);
	pol2cart(x2,y2,&cx2,&cy2,&cz2);
	dx=cx2-cx1;
	dy=cy2-cy1;
	dz=cz2-cz1;
	dist1=sqrt(dx*dx+dy*dy+dz*dz);
	dist=2*asin(dist1/(2*RTerra))*360.0/(2*Pi);
	} break;
	*/
		default:
			{
				dx=x2-x1;
				dy=y2-y1;
				dist=sqrt(dx*dx+dy*dy);
			} break;
	}
	*d=dist;
	return;
}

void RK1cartesianstep(double lon,double lat,double h,double u,double v,double *lonn,double *latn)
{
	double x,y,z,xn,yn,zn;
	pol2cart(lon,lat,&x,&y,&z);
	lat=lat/360*2*Pi;
	lon=lon/360*2*Pi;
	//NB: It assumes that u and v are in cm/sec, not deg./sec!
	xn=x+(-u*sin(lon)-v*cos(lon)*sin(lat))*h;
	yn=y+(u*cos(lon)-v*sin(lat)*sin(lon))*h;
	zn=z+(v*cos(lat))*h;
	cart2pol(xn,yn,zn,lonn,latn);
	//printf("RK1cartesianstep:\nxn:%lf yn:%lf zn:%lf lonn:%lf latn:%lf\n",xn,yn,zn,*lonn,*latn);
}

void RK4cartesianstep(double t,double x,double y,double h,double *newstatus) {
	double xyp1[2],xyp2[2],xyp3[2],xyp4[2],thalf,tfull,halfh,xtmp,ytmp;
	//RK step 1
	halfh=h/2.0;

	pField->xyp(t,x,y,xyp1);
	thalf=t+halfh;
	//xtmp=x+xyp1[0]*halfh;
	//ytmp=y+xyp1[1]*halfh;
	RK1cartesianstep(x,y,halfh,xyp1[0],xyp1[1],&xtmp,&ytmp);

	//RK step 2
	pField->xyp(thalf,xtmp,ytmp,xyp2);
	//xtmp=x+xyp2[0]*halfh;
	//ytmp=y+xyp2[1]*halfh;
	RK1cartesianstep(x,y,halfh,xyp2[0],xyp2[1],&xtmp,&ytmp);

	//RK step 3
	pField->xyp(thalf,xtmp,ytmp,xyp3);
	tfull=t+h;
	//xtmp=x+xyp3[0]*h;
	//ytmp=y+xyp3[1]*h;
	RK1cartesianstep(x,y,h,xyp3[0],xyp3[1],&xtmp,&ytmp);

	//RK step 4
	pField->xyp(tfull,xtmp,ytmp,xyp4);

	//newstatus[0]=x+h/6.0*(xyp1[0]+2.0*(xyp2[0]+xyp3[0])+xyp4[0]);
	//newstatus[1]=y+h/6.0*(xyp1[1]+2.0*(xyp2[1]+xyp3[1])+xyp4[1]);
	RK1cartesianstep(x,y,h/6.0,xyp1[0]+2.0*(xyp2[0]+xyp3[0])+xyp4[0],xyp1[1]+2.0*(xyp2[1]+xyp3[1])+xyp4[1],newstatus,newstatus+1);

	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//Euler method, next line only for testing! 
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//RK1cartesianstep(x,y,h,xyp1[0],xyp1[1],newstatus,newstatus+1);

	//done in RK4degstep
	//while(y>90) y=180-y;
	//while(y<-90) y=-180-y;
	//while(x>179) x=x-360;
	//while(x<-180) x=x+360;
}



void RK4degstep(double t,double x,double y,double h,double *newstatus) {
	double xyp1[2],xyp2[2],xyp3[2],xyp4[2],thalf,tfull,halfh,xtmp,ytmp;

	//RK step 1
	halfh=h/2.0;

   /*
    This is not needed, since it is repeated in xyp each time
    Must be moved at the end to correct new status instead
	////TEST
	while(x>180.0) x-=360.0;
	while(x<-180.0) x+=360.0;
	while(y>90.0) y=180.0-y;
	while(y<-90.0) y=-180.0-y;
	//////////////
   */
	//printf("========================================\n");	
	//printf("RK 1 step\n");	
	pField->xyp(t,x,y,xyp1);
	thalf=t+halfh;
	xtmp=x+xyp1[0]*halfh;
	ytmp=y+xyp1[1]*halfh;

	//printf("xyp1[0]=%lf xyp1[1]=%lf\n",xyp1[0],xyp1[1]);	
	//printf("RK 1 step xtmp=%lf ytmp=%lf\n",xtmp,ytmp);	

	//RK step 2
	//printf("RK 2 step\n");	
	pField->xyp(thalf,xtmp,ytmp,xyp2);
	xtmp=x+xyp2[0]*halfh;
	ytmp=y+xyp2[1]*halfh;
	//printf("RK 2 step xtmp=%lf ytmp=%lf\n",xtmp,ytmp);	

	//RK step 3
	//printf("RK 3 step\n");	
	pField->xyp(thalf,xtmp,ytmp,xyp3);
	tfull=t+h;
	xtmp=x+xyp3[0]*h;
	ytmp=y+xyp3[1]*h;
	//printf("RK 3 step xtmp=%lf ytmp=%lf\n",xtmp,ytmp);	

	//RK step 4
	//printf("RK 4 step\n");	
	pField->xyp(tfull,xtmp,ytmp,xyp4);
	//printf("========================================\n");	

	newstatus[0]=x+h/6.0*(xyp1[0]+2.0*(xyp2[0]+xyp3[0])+xyp4[0]);
	newstatus[1]=y+h/6.0*(xyp1[1]+2.0*(xyp2[1]+xyp3[1])+xyp4[1]);
   // Correct new positions if beyond limits
	while(newstatus[0]>180.0) newstatus[0]-=360.0;
	while(newstatus[0]<-180.0) newstatus[0]+=360.0;
	while(newstatus[1]>90.0) newstatus[1]=180.0-newstatus[1];
	while(newstatus[1]<-90.0) newstatus[1]=-180.0-newstatus[1];
}

void RK4flatstep(double t,double x,double y,double h,double *newstatus) {
	double xyp1[2],xyp2[2],xyp3[2],xyp4[2],thalf,tfull,halfh,xtmp,ytmp;

	//RK step 1
	halfh=h/2.0;
	//////////////
	pField->xyp(t,x,y,xyp1);
	thalf=t+halfh;
	xtmp=x+xyp1[0]*halfh;
	ytmp=y+xyp1[1]*halfh;
	//printf("xyp1[0]=%lf xyp1[1]=%lf\n",xyp1[0],xyp1[1]);	

	//RK step 2
	pField->xyp(thalf,xtmp,ytmp,xyp2);
	xtmp=x+xyp2[0]*halfh;
	ytmp=y+xyp2[1]*halfh;

	//RK step 3
	pField->xyp(thalf,xtmp,ytmp,xyp3);
	tfull=t+h;
	xtmp=x+xyp3[0]*h;
	ytmp=y+xyp3[1]*h;

	//RK step 4
	pField->xyp(tfull,xtmp,ytmp,xyp4);

	newstatus[0]=x+h/6.0*(xyp1[0]+2.0*(xyp2[0]+xyp3[0])+xyp4[0]);
	newstatus[1]=y+h/6.0*(xyp1[1]+2.0*(xyp2[1]+xyp3[1])+xyp4[1]);
}

void RK4step(double t,double x,double y,double h,double *newstatus) {
	double noise1,noise2;
	if (pField->gridtype==0) RK4flatstep(t,x,y,h,newstatus);
	else
		if (pField->datatype==2) RK4cartesianstep(t,x,y,h,newstatus);
		else RK4degstep(t,x,y,h,newstatus);
	if(pField->noise) {
		boxmuller(&noise1,&noise2);
		if(pField->datatype==1) noise2=noise2/cos(y/180*Pi);
		newstatus[0]+=noise1*(pField->noiselevel)*sqrt(fabs(h));
		newstatus[1]+=noise2*(pField->noiselevel)*sqrt(fabs(h));
	}
}	

int RK4(double *tspan,double *x0,int numx0,int Nstep,double *trj)
{
	/*
	   Like RK4_tau, but computing the trajectories. The trajectories are
	   stored in trj as (x1(t1)...x1(tmax) y1(1)...y1(tmax) ... xn(t1)...xn(tmax)).
	   */
	int j1,j2,indx,indy;
	double h,halfh,thalf,t,tfull,xyp1[2],xyp2[2],xyp3[2],xyp4[2],x,y,xtmp,ytmp,RKout[2],noise1,noise2;

	h=(tspan[1]-tspan[0])/(double)Nstep;
	halfh=h/2.0;

	t=tspan[0];
	//trj=x1 x2 x3 x4 ... y1 y2 y3 y4 ...
	//printf("numx0=%d\n",numx0);
	for(j2=0;j2<numx0;j2++) {
		trj[2*j2*Nstep]=x0[2*j2];
		trj[(2*j2+1)*Nstep]=x0[2*j2+1];
	}

	for (j1=1;j1<Nstep;j1++) {

		t=tspan[0]+(j1-1)*h;
		//printf("j1=%d Nstep=%d\n",j1,Nstep);
		for(j2=0;j2<numx0;j2++) {
			indx=j1+2*j2*Nstep;
			indy=j1+(2*j2+1)*Nstep;
			x=trj[indx-1];
			y=trj[indy-1];

			RK4step(t,x,y,h,RKout);

			//printf("RKout[0]=%lf RKout[1]=%lf\n",RKout[0],RKout[1]);
			trj[indx]=RKout[0];
			trj[indy]=RKout[1];		
		}

	}
	return 0;
}

int RK4_tau(double *tspan,double *x0,int numx0,int Nstep,double delta,int *elmax,double *tau)
{
	/*
	   This function integrates a number numx0 of 2d points which initial conditions are
	   contained in the vector x0 as: (x1 y1 x2 y2 x3 y3... xn yn). The idea
	   is to put in x1 y1 a point, and in the other variables other points around
	   x1 y1. The function computes the distance of all trajectories from x1 y1 and
	   stops when one has a distance greater than delta. It then provides the time
	   (in tau) and the number of the trajectory (in elmax). (I.e., "3" if the trajectory started)
	   at point x3 y3.) The integrator is a Runge-Kutta of 4th order, with fixed
	   time step (Nstep number of steps in the time window tspan[1], tspan[2]).
	   */
	int j1,j2,indx,indy;
	double h,halfh,thalf,t,tfull,xyp1[2],xyp2[2],xyp3[2],xyp4[2],x,y,xtmp,ytmp,dist,dx,dy,RKout[2];
	double *newstatus,*oldstatus,*swapstatus;

	*tau=0.0;
	*elmax=-1;

	newstatus=new double[numx0*2];
	oldstatus=new double[numx0*2];

	h=(tspan[1]-tspan[0])/(double)Nstep;
	halfh=h/2.0;

	t=tspan[0];
	//trj=x1 x2 x3 x4 ... y1 y2 y3 y4 ...
	//printf("numx0=%d\n",numx0);
	for(j2=0;j2<numx0*2;j2++) oldstatus[j2]=x0[j2];

	for (j1=1;j1<Nstep;j1++) {

		t=tspan[0]+(j1-1)*h;
		tfull=t+h;


		for(j2=0;j2<numx0;j2++) {
			x=oldstatus[j2*2];
			y=oldstatus[j2*2+1];

			RK4step(t,x,y,h,RKout);

			newstatus[2*j2]=RKout[0];
			newstatus[2*j2+1]=RKout[1];

			if(j2>0)
			{
				//dx=newstatus[2*j2]-newstatus[0];
				//dy=newstatus[2*j2+1]-newstatus[1];
				//dist=sqrt(dx*dx+dy*dy);
				Tdist(newstatus[0],newstatus[1],newstatus[2*j2],newstatus[2*j2+1],&dist);
				if(dist>delta) {

					*elmax=j2;
					*tau=tfull-tspan[0]; //ADDED: -tspan[0]
					j2=numx0*2;
					j1=Nstep;
				}


			}


			//printf("x=%lf\t y=%lf\n",trj[j1],trj[j1+Nstep]);
		}

		swapstatus=oldstatus;
		oldstatus=newstatus;
		newstatus=swapstatus;
	}
	delete oldstatus;
	delete newstatus;
	return 0;
}

int RK4_delta(double *tspan,double *x0,int numx0,int Nstep,double *direxp,double *delta,double *dirx0)
{
	/*
	   Like RK4_tau, but returning the maximum distance between the trajectories
	   starting in x2 y2... xn yn and the one starting in x1 y1.
	   Used for the Finite Time Lyapunov Exponent.
	   */
	int j1,j2,indx,indy;
	double h,halfh,thalf,t,tfull,xyp1[2],xyp2[2],xyp3[2],xyp4[2],x,y,xtmp,ytmp,dist,dx,dy,lastf[2];
	double *newstatus,*oldstatus,*swapstatus;



	newstatus=new double[numx0*2];
	oldstatus=new double[numx0*2];

	h=(tspan[1]-tspan[0])/(double)Nstep;
	halfh=h/2.0;

	t=tspan[0];
	//trj=x1 x2 x3 x4 ... y1 y2 y3 y4 ...
	//printf("numx0=%d\n",numx0);
	for(j2=0;j2<numx0*2;j2++) oldstatus[j2]=x0[j2];

	for (j1=1;j1<Nstep;j1++) {
		t=tspan[0]+(j1-1)*h;
		//printf("t=%f\n",t);
		for(j2=0;j2<numx0;j2++) {
			x=oldstatus[j2*2];
			y=oldstatus[j2*2+1];
			RK4step(t,x,y,h,newstatus+2*j2);
		}
		swapstatus=oldstatus;
		oldstatus=newstatus;
		newstatus=swapstatus;
	}
	*delta=0.0;
	*direxp=-1.0;
	pField->xyp(t,newstatus[0],newstatus[1],lastf);
	*dirx0=2*Pi+atan(lastf[1]/lastf[0]);
	for(j2=1;j2<numx0;j2++) {
		dx=newstatus[2*j2]-newstatus[0];
		dy=newstatus[2*j2+1]-newstatus[1];
		//dist=sqrt(dx*dx+dy*dy);
		Tdist(newstatus[0],newstatus[1],newstatus[2*j2],newstatus[2*j2+1],&dist);
		if(dist>(*delta)) {
			*delta=dist;
			*direxp=2*Pi+atan(dy/dx);
		}
	}
	delete oldstatus;
	delete newstatus;
	if (dist==0) {
		printf("RK4_delta:\nError! dist=0 dx:%lf dy:%lf u:%lf v:%lf \n",dx,dy,lastf[0],lastf[1]);
		printf("x0_x%lf x0_y%lf\n",x0[0],x0[1]);
	}
	return 0;
}

int RK4_mix(double *tspan,double *x0,int numx0,int Nstep,double sc,int *flagct,int *ndisk,int *quality)
{
	/*
	   Like RK4_delta, but returning the number of discs of radius sc needed to cover the trajectories
	   starting in x1 y1... xn yn.
	   */
	int j1,j2,indx,indy;
	double h,halfh,thalf,t,tfull,xyp1[2],xyp2[2],xyp3[2],xyp4[2],x,y,xtmp,ytmp,dist,dx,dy,lastf[2];
	double *newstatus,*oldstatus,*swapstatus;

	newstatus=new double[numx0*2];
	oldstatus=new double[numx0*2];

	h=(tspan[1]-tspan[0])/(double)Nstep;
	halfh=h/2.0;

	t=tspan[0];
	//trj=x1 x2 x3 x4 ... y1 y2 y3 y4 ...
	//printf("numx0=%d\n",numx0);
	for(j2=0;j2<numx0*2;j2++) oldstatus[j2]=x0[j2];

	for (j1=1;j1<Nstep;j1++) {
		t=tspan[0]+(j1-1)*h;
		//printf("t=%f\n",t);
		for(j2=0;j2<numx0;j2++) {
			x=oldstatus[j2*2];
			y=oldstatus[j2*2+1];
			RK4step(t,x,y,h,newstatus+2*j2);
		}
		swapstatus=oldstatus;
		oldstatus=newstatus;
		newstatus=swapstatus;
	}
	*ndisk=0;
	scalecount(oldstatus,numx0,sc,flagct,ndisk,quality);	

	delete oldstatus;
	delete newstatus;
	if ((*ndisk<=0)||(*ndisk>numx0)) {
		printf("RK4_mix:\nError!: wrong disk number (ndisk=%d)!\n",*ndisk);
		//printf("dx:%lf dy:%lf u:%lf v:%lf \n",dx,dy,lastf[0],lastf[1]);
		//printf("x0_x%lf x0_y%lf\n",x0[0],x0[1]);
	}
	return 0;
}

int scalecount(double *v,int N,double sc,int *flagct,int *ndisk,int *quality) {

	int ct,ct1,ct2,already_covered;
	double dist;
	*ndisk=0;
	flagct[0]=1;
	for(ct=1;ct<N;ct++) flagct[ct]=0; 

	for(ct1=1;ct1<N;ct1++) {
		//printf("scalecount:\nct1:%d\n",ct1);
		already_covered=0;
		for(ct2=0;ct2<ct1;ct2++) if(flagct[ct2]){
			Tdist(v[2*ct1],v[2*ct1+1],v[2*ct2],v[2*ct2+1],&dist);
			//printf("scalecount:\nct1:%d ct2:%d dist:%lf\n",ct1,ct2,dist);
			if(dist<sc) {
				already_covered=1;
				ct2=ct1;
			}
		}
		if(already_covered==0) {flagct[ct1]=1;*ndisk++;}
	}
	//quality: not yet implemented!
	return 0;
}
