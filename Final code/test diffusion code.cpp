#include "output.cpp"
FILE *fp1,*fp2,*fp3;
/*
float gauss(vec P,float t,float &S)//square diffusion
{
	const float Lx=1.0,Ly=1.0;
	float X,Y,T,x=P.i,y=P.j,tg=pi*pi*kf*(1/(Lx*Lx)+1/(Ly*Ly));

	X=sin(pi*x/Lx);				//X2=-pi*pi*X/(L*L);
	Y=sin(pi*y/Ly);				//Y2=-pi*pi*Y/(L*L);
	T=100*exp(-t*tg);	//T1=-2*kf*pi*pi*T/(L*L);			
	S=0;
	return X*Y*T;
}	

float gauss(vec P,float t,float &S)//square test function
{
	S=0;
	float To=10,Tx=5,Ty=7,Txx=10,Tyy=15,Txy=20,Txxx=0,Tyyy=0,Txxy=0,Txyy=0,x=P.i,y=P.j;
	return To+Tx*x+Ty*y+Txx*x*x/2+Tyy*y*y/2+Txy*x*y+Txxx*x*x*x/6+Tyyy*y*y*y/6+Txxy*x*x*y/2+Txyy*x*y*y/2;
}

float gauss(vec P,float t,float &S){S=1;return 0;}


float gauss(vec P,float t,float &S)//circular flutter
{
	const float Ro=1,Ri=.1,Bg=(Ro-Ri),Tg=0.1,c=Bg/Tg,Co=50;
	float s=.2,rp,r1p,E,E1,E2,R,R1,R2,St,Sr,Srr,x=P.i,y=P.j,r=sqrt(x*x+y*y);

Bg=Bg-k*cos(2*pi*(r-Ri-c*t)*cos(2*pi*(r-Ri-c*t)/Bg
	E=sin(2*pi*(r-Ri-c*t)/Bg);	E1=cos(2*pi*(r-Ri-c*t)/Bg)*2*pi/Bg;	E2=-4*pi*pi/(Bg*Bg)*E;
	R=(r/Ri-1)*(1-r/Ro);		R1=(r/Ri-1)/Ri-(1-r/Ro)/Ro;	R2=-2/(Ri*Ro);
	Srr=E2*R+E*R2+2*E1*R1;		Sr=E1*R+E*R1;
	
	St=-c*E1*R*Co;				Srr=Srr+Sr/r;			S=St-kf*Srr*Co;
	return E*R*Co;
}
float gauss(vec P,float t,float &S)//circular gaussian
{
	const float Ro=1,Ri=.1,Ag=0.5*(Ro+Ri),Bg=0.25*(Ro-Ri),Tg=0.1,Co=50;
	float s=.1,rp,r1p,E,E1,E2,R,R1,R2,St,Sr,Srr,x=P.i,y=P.j,r=sqrt(x*x+y*y);

	rp=Ag+Bg*sin(2*pi*t/Tg);	r1p=Bg*cos(2*pi*t/Tg)*2*pi/Tg;

	E=exp((rp-r)*(r-rp)/(s*s));	E1=-2*(r-rp)*E/(s*s);	E2=-2*(E+(r-rp)*E1)/(s*s);
	R=(r/Ri-1)*(1-r/Ro);		R1=1/Ri+1/Ro-2*r/(Ri*Ro);	R2=-2/(Ri*Ro);
	Srr=E2*R+E*R2+2*E1*R1;		Sr=E1*R+E*R1;
	
	St=-R*E1*r1p*Co;				Srr=Srr+Sr/r;			S=St-kf*Srr*Co;
	return E*R*Co;
}
float gauss(vec P,float t,float &S)//square gaussian
{
	const float L=1.0,Ag=.5*L,Bg=.25*L,Tg=0.1,Co=400.0,C1=0;//1/sqrt(kf);
	float s=.1,xp,x1p,E,E1,E2,X,X1,X2,Y,Y2,T,St,Sx,Sy,x=P.i,y=P.j;

	xp=Ag+Bg*sin(2*pi*t/Tg);	x1p=Bg*cos(2*pi*t/Tg)*2*pi/Tg;

	E=exp((xp-x)*(x-xp)/(s*s));	E1=-2*(x-xp)*E/(s*s);	E2=-2*(E+(x-xp)*E1)/(s*s);
	X=x*(L-x);					X1=L-2*x;				X2=-2;
	Y=sin(pi*y/L);				Y2=-pi*pi*Y/(L*L);
	T=Co+C1*t;					
	
	St=X*Y*(E*C1-T*E1*x1p);	Sx=(E2*X+E*X2+2*E1*X1)*Y*T;	Sy=E*X*Y2*T;	S=St-kf*(Sx+Sy);
	return E*X*Y*T;
}

*/
float gauss(vec P,float t,float &S)//square gaussian
{
	const float L=1.0,Ag=.5*L,Bg=.25*L,Tg=0.1,Co=400.0,C1=0;//1/sqrt(kf);
	float s,smax=.2,smin=.05,s1,xp,x1p,E,E1,E2,Es,X,X1,X2,Y,Y2,T,St,Sx,Sy,x=P.i,y=P.j;

	xp=Ag+Bg*sin(2*pi*t/Tg);	x1p=Bg*cos(2*pi*t/Tg)*2*pi/Tg;
	s=(smax+smin)/2+(smax-smin)/2*cos(4*pi*t/Tg);	s1=-(smax-smin)/2*sin(4*pi*t/Tg)*4*pi/Tg;

	E=exp((xp-x)*(x-xp)/(s*s));	
	E1=-2*(x-xp)*E/(s*s);		Es=2*(x-xp)*(x-xp)*E/(s*s*s);
	E2=-2*(E+(x-xp)*E1)/(s*s);	
	X=x*(L-x);					X1=L-2*x;				X2=-2;
	Y=sin(pi*y/L);				Y2=-pi*pi*Y/(L*L);
	
	St=X*Y*(Es*s1-E1*x1p);	Sx=(E2*X+E*X2+2*E1*X1)*Y;	Sy=E*X*Y2;	S=St*Co-kf*(Sx+Sy)*Co;
	return E*X*Y*Co;
}


//float gauss(vec P,float t,float &S){float x=P.i,y=P.j;S=-2*(x*x+y*y-y-x)/6.25e-4;return x*(1-x)*y*(1-y)/6.25e-4;}
//float gauss(vec P,float t,float &S){float x=P.i,y=P.j;S=0;return x*x/2;}
//float gauss(vec P,float t,float &S){float x=P.i,y=P.j;S=-80-600*x-600*y;return 1+12*x+20*x*x+100*x*x*x+1+12*y+20*y*y+100*y*y*y;}
//float gauss(vec P,float t,float &S){float x=P.i,y=P.j;S=-80-600*x-600*y-1200*x*x-1200*y*y;
//return 1+12*x+20*x*x+100*x*x*x+100*x*x*x*x+1+12*y+20*y*y+100*y*y*y+100*y*y*y*y;}
//float gauss(vec P,float t,float &S){float x=P.i,y=P.j;S=4;S=200*S;return 200*x*(1-x)+200*y*(1-y);}

//x2(1-2x+x2)+y2(1-2y+y2)2-12x-24x+
void inputpoints(FILE *fp)
{	
	for(int i=1;i<=pmax;i++)
	{
//		fscanf (fp,"%f %f ",&P[i].i,&P[i].j);if(i<=bmax)P[i].BN=true;	else P[i].BN=false;
		fscanf (fp,"%f %f %i",&P[i].i,&P[i].j,&P[i].BN);

//		if(i<=32)P[i].Tc=100;if(i>32&&i<=bmax)P[i].Tc=400;if(i>bmax)P[i].Tc=250;
//		if(i<=bqmax/4+1)P[i].Tc=100;else P[i].Tc=0;
		P[i].Tc=gauss(P[i],0,P[i].S);
		P[i].Tf=P[i].Tc;
	}
}	

void triangle::inputcell(FILE *fp,int jj)
{
	fscanf (fp,"%i %i %i",&p[0],&p[1],&p[2]);
	initialize(p[0],p[1],p[2],jj);
	fscanf (fp,"%i %i %i",&N[0],&N[1],&N[2]);

	if(volmin>cell[j].Vol)volmin=cell[j].Vol;//C.Tc=250;
	C.Tc=gauss(C,0,S);
	Tf=C.Tc;Tc=C.Tc;//To=C.Tc;
}
//	X=P[p[0]].level>=P[p[1]].level;	Y=P[p[1]].level>=P[p[2]].level;Z=P[p[2]].level>=P[p[0]].level;
//	if(X&&!Z) l=P[p[0]].level;	if(Y&&!X) l=P[p[1]].level;	if(Z&&!Y) l=P[p[2]].level;

void triangle::timestep()
{
	if(fabs(E1)*dt2>500) dt2=500/fabs(E1);
	if(1/dtl<fabs(Flux+S)) dtl=1/fabs(Flux+S);
}

void triangle::gridtaging()
{
 	float w;
	if(P[p[0]].level>=P[p[1]].level&&P[p[0]].level>=P[p[2]].level) level=P[p[0]].level;	
	if(P[p[1]].level>=P[p[2]].level&&P[p[1]].level>=P[p[0]].level) level=P[p[1]].level;	
	if(P[p[2]].level>=P[p[0]].level&&P[p[2]].level>=P[p[1]].level) level=P[p[2]].level;
	if(lmax<level)lmax=level;

	w=(E21*Vol+cell[N[0]].E21*cell[N[0]].Vol+cell[N[1]].E21*cell[N[1]].Vol+cell[N[2]].E21*cell[N[2]].Vol);
	w=w/(Vol+cell[N[0]].Vol+cell[N[1]].Vol+cell[N[2]].Vol);
	w=E2;
	if(fabs(w)>5&&level!=5) split=-1;
	if(fabs(w)<1&&level!=0) coarse=1;
}

int main()
{
fp1=fopen("sqper.txt","r"),*fp2;float w,Dcoef,rem=-.05;
//fscanf (fp1,"%i ",&bmax);
fscanf (fp1,"%i ",&pmax);
inputpoints(fp1);volmin=100.0;
fscanf (fp1,"%i",&nmax);for(int j=1;j<=nmax;j++) cell[j].inputcell(fp1,j);fclose(fp1);
Dcoef=.2;dt=Dcoef*volmin*rho/kf;
//dt=1e-4;
//dt=.1*dt;Dcoef=kf*dt/(rho*volmin);
//.28 exlpoded 20x20 grid.
//.22 exlpoded CC grid.

string p="y";
printf("max boundry node=%i \t nodes=%i \t triangles=%i\n",bmax,pmax,nmax);
printf("Diffusion coef=%1.3f\ttime step = %e\n\n",Dcoef,dt);getch();

reinitialize();Nmax=nmax;
for(int j=1;j<=nmax;j++) cell[j].parameters();
errors(1);derivatives(1);
checkmesh();
pass=0;
//for(int j=1;j<=nmax;j++) cell[j].split=-1;celloperations();
//for(int j=1;j<=nmax;j++) cell[j].split=-1;celloperations();

for(int j=1;j<=nmax;j++) cell[j].Tc=gauss(cell[j].C,0,cell[j].S);
for(int i=1;i<=pmax;i++) if(!P[i].BN) P[i].Tc=gauss(P[i],0,P[i].S);
errors(1);derivatives(1);

output();printf("initial\n");	p=getch();if(p=="q")exit(0);

fp1=fopen("emax.txt","w");
fp2 = fopen("temptec.dat","w");fprintf(fp2,"VARIABLES=x,y,Temp\n");

//----------------**********************----------------------------------------	
Count=-1;
for(t=0;true;t=t+dt)//time loop begins.
{
 	emax=0;Count++;	  
	printf("\nNo=%i\n",Count);
 	pass=0;R=0;
	do
	{	
//		errors(1);for(int j=1;j<=nmax;j++) cell[j].trucationerrors(1);
//		for(int j=1;j<=nmax;j++) gauss(cell[j].C,t,cell[j].S);
//		for(int j=1;j<=nmax;j++) cell[j].heatflux(1);
//		for(int i=1;i<=pmax;i++) {if(!P[i].BN) P[i].nodetemp();else P[i].Flux=0;}

		errors(1);
		pass++;lmax=0;
		
		for(int j=1;j<=nmax;j++) cell[j].gridtaging();
	 	printf("level=%i\tpass=%i\tnmax=%i\tpmax=%i\t",lmax,pass,nmax,pmax);

// 		output();p=getch();if(p=="q")exit(0);

	    R=0;celloperations();
		printf("change=%i\n",R);	p="y";
		//----------------**********************----------------------------------------	
		if(t==0)
		{
			for(int j=1;j<=nmax;j++) cell[j].Tc=gauss(cell[j].C,0,cell[j].S);
			for(int i=1;i<=pmax;i++) if(!P[i].BN) P[i].Tc=gauss(P[i],0,P[i].S);
		}
		//----------------**********************----------------------------------------	
	}while(R==1&&pass<=8);	//rintf("\n");	

	volmin=cell[1].Vol;
	for(int j=1;j<=nmax;j++) if(volmin>cell[j].Vol)volmin=cell[j].Vol;dt=Dcoef*volmin/kf;
	printf("dts=%e\t\t",dt);

	errors(1);//derivatives(1);for(int j=1;j<=nmax;j++) cell[j].trucationerrors(1);

	for(int j=1;j<=nmax;j++) gauss(cell[j].C,t,cell[j].S);
	for(int j=1;j<=nmax;j++) cell[j].heatflux(1);
	
	dtl=100;	for(int j=1;j<=nmax;j++) cell[j].timestep();
	printf("\t\tdt 1limit=%e\tdt\n\n",dtl);

	if(dtl<dt){dt=dtl;for(int j=1;j<=nmax;j++) cell[j].heatflux(1);}
	for(int i=1;i<=pmax;i++) if(!P[i].BN) P[i].nodetemp();
	
//	if(Count%1==0){output();if(p=="q")exit(0);}

	for(int j=1;j<=nmax;j++){Rmax=fabs(gauss(cell[j].C,t+dt,w)-cell[j].Tf);if(emax<Rmax)emax=Rmax;}
	fprintf(fp1,"%f %f\n",t/0.1,emax);fflush(fp1);

	if(t/.1-rem>=.05)
	{
		rem=rem+.05;output();
		fprintf(fp2,"ZONE T =\"t/tg=%f\"\n",t/0.1);
		for(int j=1;j<=nmax;j++)cell[j].plot3D(fp2);
//		for(int i=1;i<=pmax;i++)fprintf(fp2,"%f %f %f\n",P[i].i,P[i].j,P[i].Tf);
		fflush(fp2);fprintf(fp2,"\n");//getch();
	}

	for(int j=1;j<=nmax;j++) {cell[j].Tc=cell[j].Tf;if(cell[j].Tc>200){printf("200>odd\n");getch();}}
	for(int i=1;i<=pmax;i++) {if(!P[i].BN)P[i].Tc=P[i].Tf;if(P[i].Tc>200){printf("200>odd\n");getch();}}
	if(rem>=1)exit(0);
}
return 0;
}	
/*
	if(dt2>10*dt&&0)
	{
		for(int j=1;j<=nmax;j++) cell[j].heatflux(1);
		for(int i=1;i<=pmax;i++) if(!P[i].BN) P[i].nodetemp();

	 	for(int j=1;j<=nmax;j++) gauss(cell[j].C,t+dt,cell[j].S);
		for(int i=1;i<=pmax;i++) gauss(P[i],t+dt,P[i].S);
		for(int cot=0;Rmax>10;cot++)
		{	
			Rmax=0;
			for(int j=1;j<=nmax;j++) cell[j].heatflux(2);
	//		for(int j=1;j<=nmax;j++) cell[j].Tf=cell[j].To;
			for(int i=1;i<=pmax;i++) if(!P[i].BN) P[i].nodetemp();
			printf("Count=%i Rmax=%f\t\n",cot,Rmax);
		}while(Rmax>10);
	 	printf("\n");
	}

void triangle::errormark()
{
float ;	int l,k;static bool found=0; bool X,Y,Z;
X=P[p[0]].level>=P[p[1]].level;	Y=P[p[1]].level>=P[p[2]].level;Z=P[p[2]].level>=P[p[0]].level;
if(X&&!Z) l=P[p[0]].level;	if(Y&&!X) l=P[p[1]].level;	if(Z&&!Y) l=P[p[2]].level;

dtt=sqrt(0.1/fabs(E1));
dtd=.1/fabs(E2);
dts=0.15*Vol/kf;
k=log10(dtd/dts)/log10(16);
//-.5<log10(dtd/dts)/log10(16)+k<.5;

if(k>-0.5&&k<0.5)k=0;
if(k>0.5&&k<1.5)k=1;		if(k<1.5&&k>2.5;)k=2;
if(k<-0.5&&k>-1.5)k=-1;		if(k<-1.5&&k>-2.5)k=-2;
X=k<0&&l==0;Y=k>0&&l==6;Z=dtt<dts&&dtt<dtd;if(X||Y||Z)k=0;

dtopt=min(dts*pow(4,k),dtd*pow(4,-k))
dtmin=min(dtopt,dtt);


split=k;
}*/
