/*
   Lyapunov extension 0.1, 16 august 2005

   This code extends the field.h objects allowing to compute
   the Lyapunov exponent as in Ott, "Chaos in dynamical systems", Sec. 4.4.
   */
//LE2d modified september 4 2008 for the sphere.

int LE2d_(double *x0,double deltat,double delta0,double *l1,double *l2,double *theta1,double *theta2)
{
	/*
	   Lyapunov exponents from the matrix of the evolutions (see Ott's book).
	   x0 contains the evoluted of three points, forming at time 0 a "L"-shaped reference frame.
	   */

	double a11,a12,a21,a22,sq,a11q,a12q,a21q,a22q,atan1,atan2,l1p,l2p,delta0q,dist1,dist2,norm1,norm2;
	// if((pField->gridtype==1)||(pField->gridtype==2)) { //if check that points are not on the 180 deg line
	// }

	a11=x0[2]-x0[0];
	a12=x0[4]-x0[0];
	a21=x0[3]-x0[1];
	a22=x0[5]-x0[1];

	if((pField->gridtype==1)||(pField->gridtype==2)) { //if check that points are not on the 180 deg line
		if(a11>180) a11-=360;
		if(a11<-180) a11+=360;
		if(a12>180) a12-=360;
		if(a12<-180) a12+=360;
	}

	a11q=a11*a11;
	a12q=a12*a12;
	a21q=a21*a21;
	a22q=a22*a22;

	//The next condition should be activated in the next version of lamta
	/*
	   if((pField->gridtype==1)||(pField->gridtype==2)) { //if on a sphere, normalize the distances
	   Tdist(x0[0],x0[1],x0[2],x0[3],&dist1);
	   Tdist(x0[0],x0[1],x0[4],x0[5],&dist2);

	   norm1=sqrt(dist1*dist1/(a11q+a21q));
	   norm2=sqrt(dist2*dist2/(a12q+a22q));

	   a11=a11*norm1;
	   a21=a21*norm1;

	   a12=a12*norm2;
	   a22=a22*norm2;



	   }
	   */
	a11q=a11*a11;
	a12q=a12*a12;
	a21q=a21*a21;
	a22q=a22*a22;

	delta0q=delta0*delta0;
	//printf("delta0:%e delat0q:%e\n",delta0,delta0q);

	sq=sqrt(((a12+a21)*(a12+a21)+(a11-a22)*(a11-a22))*((a12-a21)*(a12-a21)+(a11+a22)*(a11+a22)));

	l1p=0.5*(a11q+a12q+a21q+a22q+sq);
	l2p=0.5*(a11q+a12q+a21q+a22q-sq);

	//printf("a11:%e a12:%e a21:%e a22:%e l1p:%e l2p:%e sq:%e\n",a11,a12,a21,a22,l1p,l2p,sq);

	*l1=1/(2*deltat)*log(l1p/(delta0q));
	*l2=1/(2*deltat)*log(l2p/(delta0q));

	if((a12==0)&&(a21==0)) {
		if(a11>a22) {*theta1=0.0;*theta2=Pi/2;}
		if(a11<=a22) {*theta2=0.0;*theta1=Pi/2;}
	}
	if((a12!=0)||(a21!=0)) {
		atan1=2*(a11*a12+a21*a22);
		atan2=a11q-a12q+a21q-a22q;
		*theta1=atan(atan1/(atan2+sq));
		*theta2=-atan(atan1/(-atan2+sq));

	}

	//printf("l1:%e l2:%e theta1:%e theta2:%e sq:%e atan1:%e atan2:%e\n",*l1,*l2,*theta1,*theta2,sq,atan1,atan2);



	return 0;

}



int FSLEext_(double *tspan,double *x00,double delta0,int Nstep,double delta,double *l1,double *l2,double *theta1,double *theta2)
{
	//Calculation of finite size Lyapunov expoenents as in Ott's book
	int j1,j2,indx,indy,numx0=3;
	double h,t,x,y,dist1=0,dist2=0,dist=0,dx,dy,RKout[2],olddist=0;

	double *newstatus,*oldstatus,*swapstatus,x0[6];


	x0[0]=x00[0];
	x0[1]=x00[1];

	x0[2]=x00[0]+delta0;
	x0[3]=x00[1];

	x0[4]=x00[0];
	x0[5]=x00[1]+delta0;

	*theta1=-99;
	*theta2=-99;
	*l1=0;
	*l2=0;

	newstatus=new double[numx0*2];
	oldstatus=new double[numx0*2];

	h=(tspan[1]-tspan[0])/(double)Nstep;

	t=tspan[0];
	//trj=x1 x2 x3 x4 ... y1 y2 y3 y4 ...
	//printf("numx0=%d\n",numx0);
	for(j2=0;j2<numx0*2;j2++) oldstatus[j2]=x0[j2];

	for (j1=1;j1<Nstep;j1++) {

		t=tspan[0]+(j1-1)*h;



		for(j2=0;j2<numx0;j2++) {
			x=oldstatus[j2*2];
			y=oldstatus[j2*2+1];

			RK4step(t,x,y,h,RKout);

			newstatus[2*j2]=RKout[0];
			newstatus[2*j2+1]=RKout[1];
		}
		Tdist(newstatus[0],newstatus[1],newstatus[2],newstatus[3],&dist1);
		Tdist(newstatus[0],newstatus[1],newstatus[4],newstatus[5],&dist2);
		if(dist1>dist2) dist=dist1;
		else dist=dist2;
		//printf("dist1:%e dist2:%e dist:%e\n",dist1,dist2,dist);
		if(dist>delta) {

			LE2d_(newstatus,t-tspan[0],delta0,l1,l2,theta1,theta2);
			j2=numx0*2;
			j1=Nstep;
		}		

		swapstatus=oldstatus;
		oldstatus=newstatus;
		newstatus=swapstatus;
	}
	delete oldstatus;
	delete newstatus;
	return 0;
}


int FTLEext_(double *tspan,double *x00,double delta0,int Nstep,double *l1,double *l2,double *theta1,double *theta2)
{
	//Calculation of finite size Lyapunov expoenents as in Ott's book
	int j1,j2,indx,indy,numx0=3;
	double h,t,x,y,dist1=0,dist2=0,dist=0,dx,dy,RKout[2],olddist=0;

	double *newstatus,*oldstatus,*swapstatus,x0[6];


	x0[0]=x00[0];
	x0[1]=x00[1];

	x0[2]=x00[0]+delta0;
	x0[3]=x00[1];

	x0[4]=x00[0];
	x0[5]=x00[1]+delta0;

	*theta1=-99;
	*theta2=-99;
	*l1=0;
	*l2=0;

	newstatus=new double[numx0*2];
	oldstatus=new double[numx0*2];

	h=(tspan[1]-tspan[0])/(double)Nstep;

	t=tspan[0];
	//trj=x1 x2 x3 x4 ... y1 y2 y3 y4 ...
	//printf("numx0=%d\n",numx0);
	for(j2=0;j2<numx0*2;j2++) oldstatus[j2]=x0[j2];

	for (j1=1;j1<Nstep;j1++) {

		t=tspan[0]+(j1-1)*h;



		for(j2=0;j2<numx0;j2++) {
			x=oldstatus[j2*2];
			y=oldstatus[j2*2+1];

			RK4step(t,x,y,h,RKout);

			newstatus[2*j2]=RKout[0];
			newstatus[2*j2+1]=RKout[1];
		}

		swapstatus=oldstatus;
		oldstatus=newstatus;
		newstatus=swapstatus;
	}

	LE2d_(newstatus,t-tspan[0],delta0,l1,l2,theta1,theta2);	
	delete oldstatus;
	delete newstatus;
	return 0;
}


int FSLEn_(double *tspan,double *x00,int numpt,double delta0,int Nstep,double delta,double *lambda,double *thetam)
{
	//Calculation of finite size Lyapunov exponenents from separations

	int elmax,ct;
	double *x,theta,U,V,tau,PI=3.14159265358979;

	x=new double[(numpt+1)*2];
	x[0]=x00[0];
	x[1]=x00[1];

	for(ct=0;ct<2*numpt;ct+=2) {
		theta=2*PI/numpt*ct;
		U=cos(theta);
		V=sin(theta);
		x[2+ct]=x[0]+delta0*U;
		x[3+ct]=x[1]+delta0*V;
		//printf("x:%lf y:%lf\n",x[2+ct],x[3+ct]);
	}


	RK4_tau(tspan,x,numpt+1,Nstep,delta,&elmax,&tau);
	*lambda=1.0/tau*log(delta/delta0);
	*thetam=2*PI/numpt*elmax;
	delete x;
	return 0;
}





