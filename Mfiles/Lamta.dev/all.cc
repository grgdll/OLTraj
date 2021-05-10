#include "oct.h"
#include "field.h"
//#include "inertial.h"
#include "lyapext.h"

//27-july-2005 added RK4_tau_v... 
//all.cc 1.0.8 12/04/2006
//adding AVISO Eastern Mediterranean.....
//20/02/2007 adding the Kerguelen plateau.....

DEFUN_DLD (select, args, ,
		"Octave file select(sys_code).")
{

	int sys_code;
	sys_code=(int)args(0).scalar_value();

	Matrix out(1,1);

	// select_sys defined in field.h
	out(0,0)=select_sys(sys_code);
	// print_par defined in field.h
	pField->print_par();

	return octave_value(out);
}

DEFUN_DLD (lon2posx, args, ,
		"Octave file lon2posx(lon).")
{

	double lon;
	lon=args(0).scalar_value();

	Matrix out(1,1);

	out(0,0)=lon2posx(lon);

	return octave_value(out);
}

DEFUN_DLD (boxmuller, , ,
		"Octave file boxmuller().")
{
	double y1,y2;
	Matrix out(2,1);
	boxmuller(&y1,&y2);
	out(0,0)=y1;
	out(1,0)=y2;
	return octave_value(out);
}

DEFUN_DLD (Tpol2cart, args, ,
		"Octave file pol2cart(lon,lat).")
{

	double lat,lon,x,y,z;
	lon=args(0).scalar_value();
	lat=args(1).scalar_value();
	Matrix out(3,1);

	pol2cart(lon,lat,&x,&y,&z);

	out(0,0)=x;
	out(1,0)=y;
	out(2,0)=z;

	return octave_value(out);
}

DEFUN_DLD (Tcart2pol, args, ,
		"Octave file cart2pol(x,y,z).")
{

	double lat,lon,x,y,z;
	x=args(0).scalar_value();
	y=args(1).scalar_value();
	z=args(2).scalar_value();

	Matrix out(2,1);
	cart2pol(x,y,z,&lon,&lat);

	out(0,0)=lon;
	out(1,0)=lat;

	return octave_value(out);
}

DEFUN_DLD (Tdist, args, ,
		"Octave file d=Tdist(x1,y1,x2,y2). x1,y1 are lon and lat of point 1, same for point 2. Tdist gives the either the Euclidean distance (default) or the distance over the sphere (the length of the arc, in deg.) depending on pField->datatype.")
{

	double x1,y1,x2,y2,d;
	x1=args(0).scalar_value();
	y1=args(1).scalar_value();
	x2=args(2).scalar_value();
	y2=args(3).scalar_value();

	Matrix out(1,1);
	Tdist(x1,y1,x2,y2,&d);

	out(0,0)=d;

	return octave_value(out);
}

DEFUN_DLD (scalecount_, args, ,
		"Octave file n=scalecount_(v,sc). Given a COLUMN vector: [x1 y1 x2 y2..xN yN] returns a N-long flag vector (one flag for each point) to put a disc in order to cover the whole set.")
{

	double sc,*v;
	int N,quality,ndisk,ct,*flagct;
	Matrix Mv=args(0).matrix_value();
	sc=args(1).scalar_value();
	N=(int)Mv.rows()/2;
	v=new double[N*2];
	flagct=new int[N];
	Matrix out(N,1);
	//printf("scalecount_:\nN:%d\n",N);

	for(ct=0;ct<2*N;ct++) v[ct]=Mv(ct,0);
	//printf("scalecount_:first cycle ok.\n");

	scalecount(v,N,sc,flagct,&ndisk,&quality);
	//printf("scalecount_:scalecount ok.\n");
	for(ct=0;ct<N;ct++) out(ct,0)=flagct[ct];

	return octave_value(out);
}





DEFUN_DLD (lat2posy, args, ,
		"Octave file lat2posy(lat).")
{

	double lat;
	lat=args(0).scalar_value();

	Matrix out(1,1);

	out(0,0)=lat2posy(lat);

	return octave_value(out);
}










DEFUN_DLD (posy2lat, args, ,
		"Octave file posy2lat(posy).")
{

	double posy;
	posy=args(0).scalar_value();

	Matrix out(1,1);

	out(0,0)=posy2lat(posy);

	return octave_value(out);
}

DEFUN_DLD (set_par, args, ,
		"Octave file set_par([par1 par2 ...]').")
{
	double *par_list;
	int num_par,j1;
	Matrix out(1,1);
	Matrix Mpar_list=args(0).matrix_value();
	num_par=(int)Mpar_list.rows();

	par_list=new double[num_par];

	for (j1=0;j1<num_par;j1++) par_list[j1]=Mpar_list(j1,0);
	pField->set_par(par_list);
	pField->freezewarning();
	pField->print_par();


	out(0,0)=0;
	delete par_list;
	return octave_value(out);
}

DEFUN_DLD (freezetime, args, ,
		"Octave file freezetime(fr_time). Freezes the velocity field at the state of time fr_time. Use unfreezetime() to cancel.")
{
	double fr_time;


	Matrix out(1,1);
	if (args.length()==1) {
		fr_time=args(0).scalar_value();
		pField->freeze=1;
		pField->frozentime=fr_time;
		out(0,0)=fr_time;
	}
	pField->freezewarning();
	return octave_value(out);
}

DEFUN_DLD (unfreezetime, args, ,
		"Octave file unfreezetime(fr_time). Reset the normal time flow for the velocity field.")
{
	printf("Attempting to unfreeze the time..\n");
	Matrix out(1,1);

	pField->freeze=0;
	out(0,0)=pField->frozentime;

	pField->freezewarning();
	return octave_value(out);
}




DEFUN_DLD (RK4, args, ,
		"Octave file out=RK4(tspan,x0,Nstep). tspan is a 2x1 column with initial and final time, x0 is a (1+N)x1 column (see field.h).")
{
	double tspan[2],*x0,*trj;
	int numx0,Nstep,j;


	Matrix Mtspan=args(0).matrix_value();
	Matrix Mx0=args(1).matrix_value();
	numx0=(int)Mx0.rows()/2;
	Nstep=(int)args(2).scalar_value();
	Matrix out(Nstep*numx0*2,1);

	trj=new double[Nstep*numx0*2];

	x0=new double[numx0*2];
	if (x0==NULL||trj==NULL) {
		printf("Error in the number of initial points or allocating the trajectories.\n");
		return octave_value(out);
	}

	tspan[0]=Mtspan(0,0);
	tspan[1]=Mtspan(1,0);
	for(j=0;j<numx0*2;j++) x0[j]=Mx0(j,0);

	//printf("tspan:%lf %lf numx0:%d Nstep:%d\n",tspan[0],tspan[1],numx0,Nstep);
	//for(j=0;j<numx0;j++) printf("x0_x:%lf x0_y:%lf\n",x0[2*j],x0[2*j+1]);

	RK4(tspan,x0,numx0,Nstep,trj);
	for(j=0;j<numx0*2*Nstep;j++) out(j,0)=trj[j];

	delete x0;
	delete trj;
	return octave_value(out);
}

DEFUN_DLD (RK4_tau, args, ,
		"Octave file out=RK4_tau(tspan,x0,Nstep,delta). tspan is a 2x1 column with initial and final time, x0 is a (1+N)x1 column (see field.h).")
{
	double tspan[2],*x0,delta,tau;
	int numx0,Nstep,elmax,j;

	Matrix out(2,1);
	Matrix Mtspan=args(0).matrix_value();
	Matrix Mx0=args(1).matrix_value();
	numx0=(int)Mx0.rows()/2;
	Nstep=(int)args(2).scalar_value();
	delta=args(3).scalar_value();

	x0=new double[numx0*2];
	if (x0==NULL) {
		printf("Error in the number of initial points.\n");
		return octave_value(out);
	}

	tspan[0]=Mtspan(0,0);
	tspan[1]=Mtspan(1,0);
	for(j=0;j<numx0*2;j++) x0[j]=Mx0(j,0);

	//printf("tspan:%lf %lf numx0:%d Nstep:%d delta:%lf\n",tspan[0],tspan[1],numx0,Nstep,delta);
	//for(j=0;j<numx0;j++) printf("x0_x:%lf x0_y:%lf\n",x0[2*j],x0[2*j+1]);

	RK4_tau(tspan,x0,numx0,Nstep,delta,&elmax,&tau);

	out(0,0)=(double)tau;
	out(1,0)=(double)elmax;
	delete x0;
	return octave_value(out);
}

DEFUN_DLD (FSLEext, args, ,
		"Octave file out=FSLEext(tspan,x00,delta0,Nstep,delta). tspan is a 2x1 column with initial and final time, x00 is a column vector with x and y component, out is a matrix giving: [lambda1,lambda2;theta1,theta2] (see lyapext.h).")
{
	double tspan[2],x00[2],delta0,delta,l1,l2,theta1,theta2;
	int Nstep,elmax,j;

	Matrix out(2,2);
	Matrix Mtspan=args(0).matrix_value();
	Matrix Mx00=args(1).matrix_value();
	delta0=args(2).scalar_value();
	Nstep=(int)args(3).scalar_value();
	delta=args(4).scalar_value();
	x00[0]=Mx00(0,0);
	x00[1]=Mx00(1,0);


	tspan[0]=Mtspan(0,0);
	tspan[1]=Mtspan(1,0);

	//printf("tspan:%lf %lf numx0:%d Nstep:%d delta:%lf\n",tspan[0],tspan[1],numx0,Nstep,delta);
	//for(j=0;j<numx0;j++) printf("x0_x:%lf x0_y:%lf\n",x0[2*j],x0[2*j+1]);
	FSLEext_(tspan,x00,delta0,Nstep,delta,&l1,&l2,&theta1,&theta2);

	out(0,0)=(double)l1;
	out(1,0)=(double)l2;
	out(0,1)=(double)theta1;
	out(1,1)=(double)theta2;


	return octave_value(out);
}

DEFUN_DLD (FTLEext, args, ,
		"Octave file out=FTLEext(tspan,x00,delta0,Nstep). tspan is a 2x1 column with initial and final time, x00 is a column vector with x and y component, out is a matrix giving: [lambda1,lambda2;theta1,theta2] (see lyapext.h).")
{
	double tspan[2],x00[2],delta0,l1,l2,theta1,theta2;
	int Nstep,elmax,j;

	Matrix out(2,2);
	Matrix Mtspan=args(0).matrix_value();
	Matrix Mx00=args(1).matrix_value();
	delta0=args(2).scalar_value();
	Nstep=(int)args(3).scalar_value();

	x00[0]=Mx00(0,0);
	x00[1]=Mx00(1,0);


	tspan[0]=Mtspan(0,0);
	tspan[1]=Mtspan(1,0);

	//printf("tspan:%lf %lf numx0:%d Nstep:%d delta:%lf\n",tspan[0],tspan[1],numx0,Nstep,delta);
	//for(j=0;j<numx0;j++) printf("x0_x:%lf x0_y:%lf\n",x0[2*j],x0[2*j+1]);
	FTLEext_(tspan,x00,delta0,Nstep,&l1,&l2,&theta1,&theta2);

	out(0,0)=(double)l1;
	out(1,0)=(double)l2;
	out(0,1)=(double)theta1;
	out(1,1)=(double)theta2;


	return octave_value(out);
}

DEFUN_DLD (FSLEn, args, ,
		"Octave file out=FSLEn(tspan,x00,numpt,delta0,Nstep,delta). tspan is a 2x1 column with initial and final time, x00 is a column vector with x and y component, out is a vector giving: [lambda1;theta1] (see lyapext.h).")
{
	double tspan[2],x00[2],delta0,delta,l1,theta1;
	int Nstep,numpt;

	Matrix out(2,1);
	Matrix Mtspan=args(0).matrix_value();
	Matrix Mx00=args(1).matrix_value();
	numpt=(int)args(2).scalar_value();
	delta0=args(3).scalar_value();
	Nstep=(int)args(4).scalar_value();
	delta=args(5).scalar_value();
	x00[0]=Mx00(0,0);
	x00[1]=Mx00(1,0);


	tspan[0]=Mtspan(0,0);
	tspan[1]=Mtspan(1,0);

	//printf("tspan:%lf %lf numx0:%d Nstep:%d delta:%lf\n",tspan[0],tspan[1],numx0,Nstep,delta);
	//for(j=0;j<numx0;j++) printf("x0_x:%lf x0_y:%lf\n",x0[2*j],x0[2*j+1]);
	FSLEn_(tspan,x00,numpt,delta0,Nstep,delta,&l1,&theta1);

	out(0,0)=(double)l1;
	out(1,0)=(double)theta1;


	return octave_value(out);
}


DEFUN_DLD (RK4_delta, args, ,
		"Octave file out=RK4_delta(tspan,x0,Nstep). tspan is a 2x1 column with initial and final time, x0 is a (1+N)x1 column (see field.h).")
{
	double tspan[2],*x0,delta,tau,direxp,dirx0;
	int numx0,Nstep,elmax,j;

	Matrix out(3,1);
	Matrix Mtspan=args(0).matrix_value();
	Matrix Mx0=args(1).matrix_value();
	numx0=(int)Mx0.rows()/2;
	Nstep=(int)args(2).scalar_value();

	x0=new double[numx0*2];
	if (x0==NULL) {
		printf("Error in the number of initial points.\n");
		return octave_value(out);
	}

	tspan[0]=Mtspan(0,0);
	tspan[1]=Mtspan(1,0);
	for(j=0;j<numx0*2;j++) x0[j]=Mx0(j,0);

	//printf("tspan:%lf %lf numx0:%d Nstep:%d delta:%lf\n",tspan[0],tspan[1],numx0,Nstep,delta);
	//for(j=0;j<numx0;j++) printf("x0_x:%lf x0_y:%lf\n",x0[2*j],x0[2*j+1]);

	RK4_delta(tspan,x0,numx0,Nstep,&direxp,&delta,&dirx0);

	out(0,0)=(double)delta;
	out(1,0)=(double)direxp;
	out(2,0)=(double)dirx0;
	delete x0;
	return octave_value(out);
}


DEFUN_DLD (rawdata_out, args, ,
		"Octave file vel=rawdata_out(indx).")
{

	long indx;
	double U,V;


	indx=long(args(0).scalar_value());

	Matrix out(2,1);
	pField->rawdata_out(indx,&U,&V);

	out(0,0)=U;
	out(1,0)=V;

	return octave_value(out);
}

DEFUN_DLD (xyp, args, ,
		"Octave file xyp(t,x,y).")
{

	double t,x,y,xyp[2];


	t=args(0).scalar_value();
	x=args(1).scalar_value();
	y=args(2).scalar_value();

	Matrix out(2,1);
	pField->xyp(t,x,y,xyp);

	out(0,0)=xyp[0];
	out(1,0)=xyp[1];

	return octave_value(out);
}
DEFUN_DLD (field_geometry, args, ,
		"Octave file field_geometry(disttype,datatype,gridtype).\ndisttype=1 (Euclidean), 2 (sphere).\ndatatype=1 (deg./sec.), 2 (cm/sec.)\ngridtype=0 (flat), 1 (sphere, regular), 2 (sphere Mercator).")
{

	int disttype,datatype,gridtype;
	Matrix out(3,1);

	disttype=(int)args(0).scalar_value();
	datatype=(int)args(1).scalar_value();
	gridtype=(int)args(2).scalar_value();

	pField->datatype=datatype;
	pField->disttype=disttype;
	pField->gridtype=gridtype;

	out(0,0)=disttype;
	out(1,0)=datatype;
	out(2,0)=gridtype;

	return octave_value(out);
}

DEFUN_DLD (field_noise, args, ,
		"Octave file field_noise(noise_on_off,noiselevel).\n noise_on_off=0 for no noise.")
{

	int noise_on_off;
	double noiselevel;
	Matrix out(2,1);

	noise_on_off=(int)args(0).scalar_value();
	noiselevel=(double)args(1).scalar_value();

	pField->noise=noise_on_off;
	pField->noiselevel=noiselevel;

	out(0,0)=noise_on_off;
	out(1,0)=noiselevel;

	// 16 Oct 2012 ======================
	pField->noisewarning();
	//===================================

	return octave_value(out);
}




DEFUN_DLD (field_geometry_read, args, ,
		"Octave file out=field_geometry_read().\nout(1)=disttype, out(2)=datatype out(3)=gridtype (see field_geometry for help).")
{

	Matrix out(3,1);


	out(0,0)=(double)pField->disttype;
	out(1,0)=(double)pField->datatype;
	out(2,0)=(double)pField->gridtype;

	return octave_value(out);
}

DEFUN_DLD (print_par, args, ,
		"Octave file print_par(). Print the parameters of the active field.")
{
	Matrix out(1,1);

	pField->print_par();

	out(0,0)=0;
	return octave_value(out);
}

DEFUN_DLD (LUT_fill, args, ,
		"Octave file LUT_fill(). Fill the LUT field with values from the active field.")
{
	unsigned long sz,jx,jy,jt,indx;
	double x,y,t,stx,sty,stt,xyp[2];
	Matrix out(1,1);
	out(0,0)=0;

	sz=(unsigned long)Lut_field.numx*(unsigned long)Lut_field.numy*(unsigned long)Lut_field.numt;

	if(Lut_field.memoryisallocated) {
		delete Lut_field.datax;
		delete Lut_field.datay;
	}


	if (sz<=0) {
		printf("Lut_field not initialised. Use set_par.\n");
		return octave_value(out);
	}
	printf("Trying to reserve two %ld float arrays...\n",sz);

	Lut_field.datax=new float[sz];
	Lut_field.datay=new float[sz];
	Lut_field.memoryisallocated=1;

	if ((Lut_field.datax==NULL)||(Lut_field.datay==NULL)) {
		printf("Not enough memory! Try with a smaller size.\n");
		return octave_value(out);
	}

	if (pField==&Lut_field) {
		printf("The LUT cannot be filled from itself. Choose another active field.\n");
		return octave_value(out);
	}

	printf("Done.\n\n");
	pField->print_par();
	printf("Filling the LUT with data from the active field...\n",sz);

	stx=(Lut_field.xf-Lut_field.xi)/(double)(Lut_field.numx-1);
	sty=(Lut_field.yf-Lut_field.yi)/(double)(Lut_field.numy-1);
	stt=(Lut_field.tf-Lut_field.ti)/(double)(Lut_field.numt-1);


	for(jt=0;jt<Lut_field.numt;jt++) for(jy=0;jy<Lut_field.numy;jy++) for(jx=0;jx<Lut_field.numx;jx++) {

		x=Lut_field.xi+stx*(double)jx;
		y=Lut_field.yi+sty*(double)jy;
		t=Lut_field.ti+stt*(double)jt;

		pField->xyp(t,x,y,xyp);
		indx=jx+jy*Lut_field.numx+jt*Lut_field.numx*Lut_field.numy;
		//printf("x:%lf y:%lf t:%lf u:%lf v:%lf indx:%ld\n",x,y,t,xyp[0],xyp[1],indx);
		Lut_field.datax[indx]=(float)xyp[0];
		Lut_field.datay[indx]=(float)xyp[1];
	}


	printf("Done.\n\n");



	out(0,0)=0;
	return octave_value(out);
}

DEFUN_DLD (LUT_frame_fill, args, ,
		"Octave file LUT_frame_fill(U,V,framenumber). Fill the LUT frame framenumber with matrix M.")
{
	unsigned long sz,jx,jy,jt,indx,framenumber;
	Matrix out(1,1);
	Matrix U=args(0).matrix_value();
	Matrix V=args(1).matrix_value();

	framenumber=(unsigned long)args(2).scalar_value();
	out(0,0)=0;

	sz=(unsigned long)Lut_field.numx*(unsigned long)Lut_field.numy*(unsigned long)Lut_field.numt;

	if (sz<=0) {
		printf("Lut_field not initialised. Use set_par.\n");
		return octave_value(out);
	}

	if((framenumber<0)||(framenumber>Lut_field.numt)) {
		printf("Invalid frame number.\n");
		return octave_value(out);
	}

	if((U.rows()!=Lut_field.numx)||(U.columns()!=Lut_field.numy)) {
		printf("Frame size mismatch.\n");
		return octave_value(out);
	}

	if((V.rows()!=Lut_field.numx)||(V.columns()!=Lut_field.numy)) {
		printf("Frame size mismatch.\n");
		return octave_value(out);
	}

	jt=framenumber;

	if((jt>(Lut_field.numt-1))||(jt<0)) {
		printf("Wrong framenumber (use 0<=framenumber<=%d, or reset the parameters with set_par).\n",Lut_field.numt-1);
		return octave_value(out);
	}

	if(Lut_field.memoryisallocated==0) {
		printf("Trying to reserve two %ld float arrays...\n",sz);

		Lut_field.datax=new float[sz];
		Lut_field.datay=new float[sz];
		Lut_field.memoryisallocated=1;

		if ((Lut_field.datax==NULL)||(Lut_field.datay==NULL)) {
			printf("Not enough memory! Try with a smaller size.\n");
			Lut_field.memoryisallocated=0;

			return octave_value(out);
		}
		for(indx=0;indx<sz;indx++) {Lut_field.datax[indx]=0;Lut_field.datay[indx]=0;}
	}

	for(jy=0;jy<Lut_field.numy;jy++) for(jx=0;jx<Lut_field.numx;jx++) {

		indx=jx+jy*Lut_field.numx+jt*Lut_field.numx*Lut_field.numy;
		//printf("jx:%d jy:%d u:%lf v:%lf indx:%ld\n",jx,jy,U(jx,jy),V(jx,jy),indx);
		Lut_field.datax[indx]=(float)U(jx,jy);
		Lut_field.datay[indx]=(float)V(jx,jy);
	}

	printf("Loaded frame %ld.\n",jt);

	out(0,0)=0;
	return octave_value(out);
}

DEFUN_DLD (RK4_tau_v, args, ,
		"Octave file out=RK4_tau_v(tspan,x_v,x0,Nstep,delta). tspan is a 2x1 column with initial and final time, x0 is a (1+N)x1 column (see field.h).")
{
	/*
	   x_v is a vector with all the points where to compute the FSLEs. x0 is a vector with N points around the origin.
	   */
	double tspan[2],*x0,delta,tau,xc[2];
	int numx0,Nstep,elmax,j,jv,numx_v;

	Matrix Mtspan=args(0).matrix_value();
	Matrix Mx_v=args(1).matrix_value();
	Matrix Mx0=args(2).matrix_value();
	numx0=(int)Mx0.rows()/2;
	numx_v=(int)Mx_v.rows()/2;
	Matrix out(2,numx_v);
	Nstep=(int)args(3).scalar_value();
	delta=args(4).scalar_value();

	x0=new double[numx0*2+2];
	if (x0==NULL) {
		printf("Error in the number of initial points.\n");
		return octave_value(out);
	}

	tspan[0]=Mtspan(0,0);
	tspan[1]=Mtspan(1,0);
	//printf("tspan:%lf %lf numx0:%d numx_v:%d Nstep:%d delta:%lf\n",tspan[0],tspan[1],numx0,numx_v,Nstep,delta);
	//printf("tspan:%lf %lf numx0:%d Nstep:%d delta:%lf\n",tspan[0],tspan[1],numx0,Nstep,delta);
	//for(j=0;j<numx0;j++) printf("x0_x:%lf x0_y:%lf\n",x0[2*j],x0[2*j+1]);
	for(jv=0;jv<numx_v;jv++) {
		x0[0]=Mx_v(2*jv,0);
		x0[1]=Mx_v(2*jv+1,0);
		//printf("jv:%d x0:%lf %lf\n",jv,Mx_v(2*jv,0),Mx_v(2*jv+1,0));
		for(j=0;j<numx0;j++) {
			x0[2+2*j]=Mx0(2*j,0)+Mx_v(2*jv,0);
			x0[2+2*j+1]=Mx0(2*j+1,0)+Mx_v(2*jv+1,0);
		}
		//for(j=0;j<(2*numx0+2);j++) printf("x0[%d]:%lf\n",j,x0[j]);
		RK4_tau(tspan,x0,numx0+1,Nstep,delta,&elmax,&tau);
		//printf("elmax:%lf tau:%lf\n",(double)elmax,(double)tau);
		out(0,jv)=(double)tau;
		out(1,jv)=(double)elmax;
	}


	delete x0;
	return octave_value(out);
}

DEFUN_DLD (FSLEext_v, args, ,
		"Octave file out=FSLEext_v(tspan,x_v,delta0,Nstep,delta). tspan is a 2x1 column with initial and final time, x0 is a (1+N)x1 column (see field.h).")
{
	/*
	   x_v is a vector with all the points where to compute the FSLEs. x0 is a vector with N points around the origin.
	   */
	double tspan[2],x0[2],delta,delta0,theta1,theta2,l1,l2;
	int numx0,Nstep,j,jv,numx_v,dpct,mpct,Mpct;

	Matrix Mtspan=args(0).matrix_value();
	Matrix Mx_v=args(1).matrix_value();
	numx_v=(int)Mx_v.rows()/2;
	Matrix out(4,numx_v);
	delta0=args(2).scalar_value();
	Nstep=(int)args(3).scalar_value();
	delta=args(4).scalar_value();

	if (x0==NULL) {
		printf("Error in the number of initial points.\n");
		return octave_value(out);
	}

	tspan[0]=Mtspan(0,0);
	tspan[1]=Mtspan(1,0);
	//printf("tspan:%lf %lf numx0:%d numx_v:%d Nstep:%d delta:%lf\n",tspan[0],tspan[1],numx0,numx_v,Nstep,delta);
	//printf("tspan:%lf %lf numx0:%d Nstep:%d delta:%lf\n",tspan[0],tspan[1],numx0,Nstep,delta);
	//for(j=0;j<numx0;j++) printf("x0_x:%lf x0_y:%lf\n",x0[2*j],x0[2*j+1]);
	dpct=5;
	mpct=5;
	Mpct=mpct+dpct;
	printf("\n=========== START Computing FSLE ===========\n");
	for(jv=0;jv<numx_v;jv++) {
		x0[0]=Mx_v(2*jv,0);
		x0[1]=Mx_v(2*jv+1,0);
		// FSLEext_ defined in lyapexp.h
		FSLEext_(tspan,x0,delta0,Nstep,delta,&l1,&l2,&theta1,&theta2);

		//printf("jv:%d\n",jv);
		//printf("jv/numx_v*100 = %f%%.\n",(double)jv/numx_v*100);
		if ((double)jv/numx_v*100>mpct && (double)jv/numx_v*100<Mpct){
			printf("Computing FSLE: %d%% done.\n",mpct);
			mpct=Mpct;
			Mpct=mpct+dpct;
		}
		//	if (jv%100==0){
		//		printf("Computed FSLE for point %d of %d \n",jv,numx_v);
		//	}
		out(0,jv)=(double)l1;
		out(1,jv)=(double)l2;
		out(2,jv)=(double)theta1;
		out(3,jv)=(double)theta2;

	}

	return octave_value(out);
}

DEFUN_DLD (FTLEext_v, args, ,
		"Octave file out=FTLEext_v(tspan,x_v,delta0,Nstep). tspan is a 2x1 column with initial and final time, x0 is a (1+N)x1 column (see field.h).")
{
	/*
	   x_v is a vector with all the points where to compute the FSLEs. x0 is a vector with N points around the origin.
	   */
	double tspan[2],x0[2],delta0,theta1,theta2,l1,l2;
	int numx0,Nstep,j,jv,numx_v;

	Matrix Mtspan=args(0).matrix_value();
	Matrix Mx_v=args(1).matrix_value();
	numx_v=(int)Mx_v.rows()/2;
	Matrix out(4,numx_v);
	delta0=args(2).scalar_value();
	Nstep=(int)args(3).scalar_value();

	if (x0==NULL) {
		printf("Error in the number of initial points.\n");
		return octave_value(out);
	}

	tspan[0]=Mtspan(0,0);
	tspan[1]=Mtspan(1,0);
	//printf("tspan:%lf %lf numx0:%d numx_v:%d Nstep:%d delta:%lf\n",tspan[0],tspan[1],numx0,numx_v,Nstep,delta);
	//printf("tspan:%lf %lf numx0:%d Nstep:%d delta:%lf\n",tspan[0],tspan[1],numx0,Nstep,delta);
	//for(j=0;j<numx0;j++) printf("x0_x:%lf x0_y:%lf\n",x0[2*j],x0[2*j+1]);
	for(jv=0;jv<numx_v;jv++) {
		x0[0]=Mx_v(2*jv,0);
		x0[1]=Mx_v(2*jv+1,0);
		FTLEext_(tspan,x0,delta0,Nstep,&l1,&l2,&theta1,&theta2);

		//printf("jv:%d\n",jv);
		out(0,jv)=(double)l1;
		out(1,jv)=(double)l2;
		out(2,jv)=(double)theta1;
		out(3,jv)=(double)theta2;

	}

	return octave_value(out);
}
DEFUN_DLD (FSLEn_v, args, ,
		"Octave file out=FSLEn_v(tspan,x_v,numpt,delta0,Nstep,delta). tspan is a 2x1 column with initial and final time, x0 is a (1+N)x1 column (see field.h).")
{
	/*
	   x_v is a vector with all the points where to compute the FSLEs. x0 is a vector with N points around the origin.
	   */
	double tspan[2],x0[2],delta,delta0,theta1,theta2,l1,l2;
	int numx0,Nstep,j,jv,numx_v,numpt;

	Matrix Mtspan=args(0).matrix_value();
	Matrix Mx_v=args(1).matrix_value();
	numx_v=(int)Mx_v.rows()/2;
	Matrix out(2,numx_v);
	numpt=(int)args(2).scalar_value();
	delta0=args(3).scalar_value();
	Nstep=(int)args(4).scalar_value();
	delta=args(5).scalar_value();

	if (x0==NULL) {
		printf("Error in the number of initial points.\n");
		return octave_value(out);
	}

	tspan[0]=Mtspan(0,0);
	tspan[1]=Mtspan(1,0);
	//printf("tspan:%lf %lf numx0:%d numx_v:%d Nstep:%d delta:%lf\n",tspan[0],tspan[1],numx0,numx_v,Nstep,delta);
	//printf("tspan:%lf %lf numx0:%d Nstep:%d delta:%lf\n",tspan[0],tspan[1],numx0,Nstep,delta);
	//for(j=0;j<numx0;j++) printf("x0_x:%lf x0_y:%lf\n",x0[2*j],x0[2*j+1]);
	for(jv=0;jv<numx_v;jv++) {
		x0[0]=Mx_v(2*jv,0);
		x0[1]=Mx_v(2*jv+1,0);
		// FSLEn_ defined in lypext.h
		FSLEn_(tspan,x0,numpt,delta0,Nstep,delta,&l1,&theta1);

		//printf("jv:%d\n",jv);
		out(0,jv)=(double)l1;
		out(1,jv)=(double)theta1;

	}

	return octave_value(out);
}


DEFUN_DLD (UVext_v, args, ,
		"Octave file out=UVext_v(t,x_v). tspan is time, x0 is a (1+N)x1 column (see field.h).")
{

	double t,xyp[2],x,y;
	int jv,numx_v;

	t=args(0).scalar_value();
	Matrix Mx_v=args(1).matrix_value();
	numx_v=(int)Mx_v.rows()/2;
	Matrix out(2,numx_v);

	for(jv=0;jv<numx_v;jv++) {
		x=Mx_v(2*jv,0);
		y=Mx_v(2*jv+1,0);
		pField->xyp(t,x,y,xyp);

		out(0,jv)=xyp[0];
		out(1,jv)=xyp[1];
	}

	return octave_value(out);
}



DEFUN_DLD (LE2d, args, ,
		"Octave file out=LE2d(x0,deltat,delta0).x0 contains the evoluted of three points, forming at time 0 a L-shaped reference frame, deltat is the time interval, delta0 the initial separation.")
{
	double deltat,l1,l2,theta1,theta2,x0[6],delta0;
	int j;

	Matrix out(4,1);
	Matrix Mx0=args(0).matrix_value();
	deltat=args(1).scalar_value();
	delta0=args(2).scalar_value();

	for(j=0;j<6;j++) x0[j]=Mx0(j,0);

	LE2d_(x0,deltat,delta0,&l1,&l2,&theta1,&theta2);

	out(0,0)=(double)l1;
	out(1,0)=(double)l2;
	out(2,0)=(double)theta1;
	out(3,0)=(double)theta2;

	return octave_value(out);
}
