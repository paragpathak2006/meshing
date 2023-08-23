	#include "implicit method.cpp"
//-----------**************--------------output functions--------------*************-----------------
void triangle::plot3D(FILE *fp)
{
 	float w;
	fprintf (fp,"%G %G %G\n",P[p[0]].i,P[p[0]].j,P[p[0]].Tc);
	fprintf (fp,"%G %G %G\n",P[p[1]].i,P[p[1]].j,P[p[1]].Tc);
	fprintf (fp,"%G %G %G\n",P[p[2]].i,P[p[2]].j,P[p[2]].Tc);
	fprintf (fp,"%G %G %G\n",P[p[0]].i,P[p[0]].j,P[p[0]].Tc);
	fprintf (fp,"\n\n");

//	fprintf (fp,"%G %G %G %i\n",P[p[0]].i,P[p[0]].j,P[p[0]].Flux,p[0]);
//	fprintf (fp,"%G %G %G %i\n",P[p[1]].i,P[p[1]].j,P[p[1]].Flux,p[1]);
//	fprintf (fp,"%G %G %G %i\n",P[p[2]].i,P[p[2]].j,P[p[2]].Flux,p[2]);
//	fprintf (fp,"%G %G %G\n",P[p[0]].i,P[p[0]].j,P[p[0]].Flux);
//	fprintf (fp,"\n\n");

//	fprintf (fp,"%G %G %G %i\n",P[p[0]].i,P[p[0]].j,(gauss(P[p[0]],t+dt,S)-P[p[0]].Tf),p[0]);
//	fprintf (fp,"%G %G %G %i\n",P[p[1]].i,P[p[1]].j,(gauss(P[p[1]],t+dt,S)-P[p[1]].Tf),p[1]);
//	fprintf (fp,"%G %G %G %i\n",P[p[2]].i,P[p[2]].j,(gauss(P[p[2]],t+dt,S)-P[p[2]].Tf),p[2]);
//	fprintf (fp,"%G %G %G %i\n",P[p[0]].i,P[p[0]].j,(gauss(P[p[0]],t+dt,S)-P[p[0]].Tf),p[0]);
//	fprintf (fp,"\n\n");

//w=(E21*Vol+cell[N[0]].E21*cell[N[0]].Vol+cell[N[1]].E21*cell[N[1]].Vol+cell[N[2]].E21*cell[N[2]].Vol);
//w=w/(Vol+cell[N[0]].Vol+cell[N[1]].Vol+cell[N[2]].Vol);

//	fprintf(fp,"%f %f %f\n\n",C.i,C.j,Tf);
//	fprintf(fp,"%f %f %e %i\n\n",C.i,C.j,(gauss(C,t+dt,S)-Tf),j);
//	fprintf(fp,"%f %f %e %i\n\n",C.i,C.j,E2,j);
//	fprintf(fp,"%f %f %e %i\n\n",C.i,C.j,(gauss(C,t+dt,w)-gauss(C,t,w))/dt-Flux-S,j);
}
void triangle::plot3DN(FILE *fp)
{
/*	plot3D(fp);int nk;
	for(int k=0;k<3;k++)for(int kk=0;kk<3;kk++)
	{
		if(j==cell[N[k]].N[kk]) nk=N[k];
		else nk=cell[N[k]].N[kk];
		if(nk!=0)cell[nk].plot3D(fp);
	}
*/
	fprintf (fp,"%G %G %G %i\n",P[p[0]].i,P[p[0]].j,gauss(P[p[0]],t+dt,S),p[0]);
	fprintf (fp,"%G %G %G %i\n",P[p[1]].i,P[p[1]].j,gauss(P[p[1]],t+dt,S),p[1]);
	fprintf (fp,"%G %G %G %i\n",P[p[2]].i,P[p[2]].j,gauss(P[p[2]],t+dt,S),p[2]);
	fprintf (fp,"%G %G %G\n",P[p[0]].i,P[p[0]].j,gauss(P[p[0]],t+dt,S));
	fprintf (fp,"%f %f %f\n\n",C.i,C.j,gauss(C,t+dt,S));
	fprintf (fp,"\n\n");
}
void triangle::draw2D(FILE *fp)
{
	fprintf (fp,"%G %G 0\n",P[p[0]].i,P[p[0]].j);
	fprintf (fp,"%G %G 0\n",P[p[1]].i,P[p[1]].j);
	fprintf (fp,"%G %G 0\n",P[p[2]].i,P[p[2]].j);
	fprintf (fp,"%G %G 0\n",P[p[0]].i,P[p[0]].j);
	fprintf (fp,"\n\n");
}
//-----------**************--------------output functions--------------*************-----------------
void triangle::err(FILE *fp,int Ar)
{
 	float w1,w2,w3,w0,wt,w,cx=.00,cy=0;
	if(Ar==1){w=0;w1=-E3c;cx=0;}
//	if(Ar==1){w=0;if(E2!=0)w1=log10(dtd/dts)/log10(4);cx=0;}
//	if(Ar==1){w=0;w1=Ec;cx=0;}
	if(Ar==1){w=0;w1=0;cx=-2*cx;}

	if(Ar==2){w=0;w1=Tf-Tc;cx=0;}
//	if(Ar==2){w=0;w1=dts;cx=0;}
	if(Ar==2){w=0;w1=0;cx=-cx;}

	if(Ar==3){w=0;w1=E1*dt+E3c;cx=0;}
//	if(Ar==3){w=0;w1=(gauss(C,t+dt,w)-Tf)/dt;cx=0;}
	if(Ar==3){w=0;w1=0;cx=0;}

	if(Ar==4){w=0;w1=E1/dt;cx=0;}
	if(Ar==4){w=0;w1=0;cx=cx;}

	if(Ar==0){w=0;w1=(gauss(C,t+dt,w)-Tf)/1;}
//	if(Ar==0){w=0;w1=gauss(C,t+dt,w);}
//	if(Ar==0){w=0;w1=(gauss(C,t+dt,w)-Tf)/dt-(gauss(C,t,w)-Tc)/dt;}
//	if(Ar==0){w=0;w1=Flux/Vol+0*C.i;cx=0;}
//	if(Ar==0){w=0;w1=0;cx=0;}

//	fprintf(fp,"%f %f %f\n",C.i+cx,C.j+cy,w);
	if(w1!=0)fprintf(fp,"%f %f %f\n",C.i+cx,C.j+cy,w1);
}
//-----------**************--------------output functions--------------*************-----------------
void output()
{
 	FILE *fp;int nk=0,jj=691;
	fp=fopen ("temp.txt","w");for(int j=1;j<=nmax;j++)cell[j].plot3D(fp);fclose(fp);
//	fp=fopen ("temp1.txt","w");for(int jj=1;jj<=nmax;jj++)cell[jj].plot3DN(fp);fclose(fp);

	fp=fopen ("model.txt","w");for(int j=1;j<=nmax;j++)cell[j].err(fp,0);fclose(fp);
	fp=fopen ("err1.txt","w"); for(int j=1;j<=nmax;j++)cell[j].err(fp,1);fclose(fp);
	fp=fopen ("err2.txt","w"); for(int j=1;j<=nmax;j++)cell[j].err(fp,2);fclose(fp);
	fp=fopen ("err3.txt","w"); for(int j=1;j<=nmax;j++)cell[j].err(fp,3);fclose(fp);

	fp=fopen ("err.txt","w");for(int j=1;j<=nmax;j++)cell[j].err(fp,4);fclose(fp);
	fp=fopen ("draw.txt","w");for(int j=1;j<=nmax;j++)cell[j].draw2D(fp);
	fclose(fp);
//	fp=fopen ("label.txt","w");for(int j=1;j<=nmax;j++)cell[j].plot3D(fp,3);fclose(fp);

}
void triangle::map3D(FILE *fp)
{
 	vec a,b,c;float dx=.004,dy=.004,bv;
	a=P[p[1]]-P[p[0]];b=P[p[2]]-P[p[0]];
	bv=a%b;
	for(int k=-20;k<=20;k++)
	{for(int kk=-20;kk<=20;kk++)
	{//	c=P[p[0]]+kk*a/20+k*b/20;
		c.i=C.i+k*dx;c.j=C.j+kk*dy;
	
//		if(((c-P[p[1]])%(c-P[p[2]]))*bv>0)
		fprintf(fp,"%f %f %f\n",c.i,c.j,map(c-C,Tc,j));
	}fprintf(fp,"\n");
	}
}
//-----------**************--------------output functions--------------*************-----------------
//-----------**************--------------output functions--------------*************-----------------
void meshop(FILE *op)
{
	fprintf (op,"%i\n",pmax);
	for(int i=1;i<=pmax;i++)
	{
		fprintf (op,"%f %f ",P[i].i,P[i].j);
		fprintf (op,"%i %i ",P[i].p1,P[i].p2);
		for(int j=0;j<=P[i].t[0];j++)fprintf (op,"%i ",P[i].t[j]);
		fprintf (op,"\n");
	}
	fprintf (op,"%i\n",nmax);
	for(int j=1;j<=nmax;j++)
	{
		fprintf (op,"%i %i %i ",cell[j].type,cell[j].split,cell[j].coarse);
		fprintf (op,"%i %i %i ",cell[j].p[0],cell[j].p[1],cell[j].p[2]);
		fprintf (op,"%i %i %i\n",cell[j].N[0],cell[j].N[1],cell[j].N[2]);
	}
	fclose(op);
}
//-----------**************--------------output functions--------------*************-----------------

void celloperations()
{
	string p;

	do{change=0;for(int j=1;j<=nmax;j++) cell[j].splitmarking();}while(change==1);

	for(int i=1;i<=pmax;i++) P[i].pointmarking();
	for(int j=1;j<=nmax;j++) cell[j].colapsebility();

	do{change=0;for(int j=1;j<=nmax;j++) cell[j].coarsemarking();}while(change==1);

	for(int j=1;j<=nmax;j++) cell[j].midpts();
	for(int j=1;j<=Nmax;j++) cell[j].cellmodification();renumber();

//	for(int j=1;j<=nmax;j++) cell[j].setneighbour();reinitialize();
//	for(int j=1;j<=nmax;j++) cell[j].parameters();
	for(int j=1;j<=nmax;j++) if(cell[j].mod)cell[j].setneighbour();
	for(int j=1;j<=nmax;j++) if(cell[j].mod)cell[j].parameters();reinitialize();
}
void checkmesh()
{
	int pt,nk;bool found1,found2,exit1=0;vec mp;
	for(int i=1;i<=pmax;i++)
	{
		while(P[i].deleted){i++;}
		for(int j=1;j<=P[i].t[0];j++)
		{
		 	found1=0;	
			for(int k=0;k<3;k++)if(cell[P[i].t[j]].p[k]==i){found1=1;break;}
			if(!found1){printf("point %i -->triangle %i found. Reverse not true.error\n",i,P[i].t[j]);exit1=1;}
		}
		mp=(P[P[i].p1]+P[P[i].p2])/2-P[i];
		if(d(mp)>1e-5&&P[i].p1!=0)
		{printf("point %i parents %i %i d=%f pmax=%i error\n",i,P[i].p1,P[i].p2,d(mp),pmax);exit1=1;getch();}
	}
	for(int j=1;j<=nmax;j++)
	{if(cell[j].Vol<0){printf("negetive volume cell %i\n",j);exit1=1;}
	for(int k=0;k<3;k++)
	{
		pt=cell[j].p[k];found2=0;
		if(P[pt].deleted){printf("triangle %i wearing deleted point %i.error\n",j,pt);exit1=1;}
		for(int jj=1;jj<=P[pt].t[0];jj++) if(P[pt].t[jj]==j){found2=true;break;}
		if(!found2){printf("triangle %i -->point %i found. Reverse not true.error\n",j,pt);exit1=1;}
	}
	}
	for(int j=1;j<=nmax;j++)
	{
		for(int k=0;k<3;k++)
		{
			found1=0;
			for(int kk=0;kk<3;kk++)
			{
	 			nk=cell[j].N[k];
				if(cell[nk].N[kk]==j)found1=1;
				if(found1)break;
			}
			if(!found1&&nk!=0){printf("triangle %i Bn=%i,%i,%i-->tringle %i found. Reverse not true.error\n",j,cell[j].p[0],cell[j].p[1],cell[j].p[2],nk);exit1=1;getch();}
		}
	}
if(exit1)getch();
}

//-----------**************--------------output functions--------------*************-----------------
//-----------**************--------------output functions--------------*************-----------------

/*	
//fprintf(fp2,"ZONE T =\"t=%i\"\n",t);
//for(int i=1;i<=pmax;i++) fprintf (fp2,"%f %f %f\n",P[i].i,P[i].j,P[i].Tc);fflush(fp2);
//fp2 = fopen ("temptec.dat","w");fprintf (fp2,"VARIABLES=X,Y,T\n");
//for(int i=1;i<=pmax;i++) fprintf (fp2,"%f %f %f\n",P[i].i,P[i].j,P[i].Tc);fflush(fp2);

vec fk;
	for(k=0;k<3;k++)
	{	
		if(k==0) {r=1;s=2;}    if(k==1) {r=2;s=0;}    if(k==2) {r=0;s=1;}
		Ao=P[p[k]].i>=P[p[r]].i&&P[p[k]].i<=P[p[s]].i;
		Bo=P[p[k]].i<=P[p[r]].i&&P[p[k]].i>=P[p[s]].i;
		if(Ao||Bo) {if(Bo){nk=r;r=s;s=nk;}break;}
	}if(k==3) k--;	

	x0=P[p[r]].i-C.i;		x1=P[p[k]].i-C.i;		x2=P[p[s]].i-C.i;
	y0=P[p[r]].j-C.j;		y1=P[p[k]].j-C.j;		y2=P[p[s]].j-C.j;
	yi=y0*(x2-x1)/(x2-x0)+y2*(x1-x0)/(x2-x0);if(fabs(x2-x0)<1e-8) yi=y0/2+y1/2; 
	
	X=x1-x0;Y=y1-yi;	I0=Y*X*(X*X/4+2*x0*X/3+x0*x0/2);//x1>x0
	X=x1-x2;Y=y1-yi;	I2=Y*X*(X*X/4+2*x2*X/3+x2*x2/2);//x1>x2
	Ixx=fabs(I0-I2);	

	X=x1-x0;Y=y1-yi;	I0=y0*Y*X*(X/3+x0)+Y*(y1+yi)*X*(X/4+x0)/2;
	X=x1-x2;Y=y1-yi;	I2=y2*Y*X*(X/3+x2)+Y*(y1+yi)*X*(X/4+x2)/2;
	Ixy=I0-I2;if(Y<0) Ixy=-Ixy;

	for(k=0;k<3;k++)
	{	
		if(k==0) {r=1;s=2;}    if(k==1) {r=2;s=0;}    if(k==2) {r=0;s=1;}
		Ao=P[p[k]].j>=P[p[r]].j&&P[p[k]].j<=P[p[s]].j;
		Bo=P[p[k]].j<=P[p[r]].j&&P[p[k]].j>=P[p[s]].j;
		if(Ao||Bo) {if(Bo){nk=r;r=s;s=nk;}break;}
	}if(k==3) k--;	
 
	x0=P[p[r]].i-C.i;		x1=P[p[k]].i-C.i;		x2=P[p[s]].i-C.i;
	y0=P[p[r]].j-C.j;		y1=P[p[k]].j-C.j;		y2=P[p[s]].j-C.j;
	xi=x0*(y2-y1)/(y2-y0)+x2*(y1-y0)/(y2-y0);if(fabs(y2-y0)<1e-8) xi=x0/2+x1/2; 

	X=x1-xi;Y=y1-y0;	I0=Y*X*(Y*Y/4+2*y0*Y/3+y0*y0/2);//x1>x0
	X=x1-xi;Y=y1-y2;	I2=Y*X*(Y*Y/4+2*y2*Y/3+y2*y2/2);//x1>x2
	Iyy=fabs(I0-I2);

//	I0=(y1-yi)*((x1*x1+x0*x0)*(x1+x0)/4-x0*(x1*x1+x0*x0+x1*x0)/3);
//	printf("Ixx=%G Iyy=%G Ixy=%G\n",Ixx,Iyy,Ixy);getch();

void triangle::derivatives2J()
{
float fx[3],fy[3];vec a,b;
fx[0]=P[p[1]].Tx/2+P[p[2]].Tx/2;
fy[0]=P[p[1]].Ty/2+P[p[2]].Ty/2;

fx[1]=P[p[2]].Tx/2+P[p[0]].Tx/2;
fy[1]=P[p[2]].Ty/2+P[p[0]].Ty/2;

fx[2]=P[p[0]].Tx/2+P[p[1]].Tx/2;
fy[2]=P[p[0]].Ty/2+P[p[1]].Ty/2;

a.i=(fx[0]*A[0].i+fx[1]*A[1].i+fx[2]*A[2].i)/Vol;
a.j=(fx[0]*A[0].j+fx[1]*A[1].j+fx[2]*A[2].j)/Vol;

b.i=(fy[0]*A[0].i+fy[1]*A[1].i+fy[2]*A[2].i)/Vol;
b.j=(fy[0]*A[0].j+fy[1]*A[1].j+fy[2]*A[2].j)/Vol;

Txx=a.i;Tyy=b.j;Txy=a.j/2+b.i/2;Txy=a.j/2+b.i/2;
}
void triangle::derivatives1J()
{
	float fx[3];vec a;
	fx[0]=P[p[1]].Tc/2+P[p[2]].Tc/2;
	fx[1]=P[p[2]].Tc/2+P[p[0]].Tc/2;
	fx[2]=P[p[0]].Tc/2+P[p[1]].Tc/2;
	
	Tx=(fx[0]*A[0].i+fx[1]*A[1].i+fx[2]*A[2].i)/Vol;
	Ty=(fx[0]*A[0].j+fx[1]*A[1].j+fx[2]*A[2].j)/Vol;
}

*/
/*void triangle::solution()
{
	int nkk,nkk1=-1,nkk2=-1,n=5;float A[5][6],x[5],p;
if(N[0]!=0&&N[1]!=0&&N[2]!=0)
{
	for(int k=0;k<3;k++)
	{
		for(int kk=0;kk<3;kk++)
		{
			nkk=cell[N[k]].N[kk];
			if(nkk!=j&&nkk!=0&&N[k]!=0&&nkk1<0)nkk1=nkk;
			if(nkk!=j&&nkk!=0&&N[k]!=0&&nkk2<0)nkk2=nkk;
			if(nkk1>=0&&nkk2>=0)break;
		}
		if(nkk1>=0&&nkk2>=0)break;
	}
	if(nkk1<0||nkk2<0)
	{printf("Not enough neibours, in mapping function for cell %i.\n",j);getch();exit(0);}
	
	Tc;
	cell[N[0]].Tc,cell[N[1]].Tc,cell[N[2]].Tc;
	cell[nkk1].Tc,cell[nkk2].Tc;
	for(int i=0;i<5;i++)
	{
		if(i<3)nkk=N[i];if(i==3)nkk=nkk1;if(i==4)nkk=nkk2;
		A[i][0]=cell[nkk].C.i-C.i;			A[i][1]=cell[nkk].C.j-C.j;		A[i][5]=cell[nkk].Tc-Tc;
		A[i][2]=A[i][0]*A[i][0]/2;			A[i][3]=A[i][1]*A[i][1]/2;		A[i][4]=A[i][0]*A[i][1];
	}
	
	for(int i=0;i<n;i++)//row change down
	{
		p=A[i][i]; if(fabs(p)<1E-15){printf("pivot %e too small in row %i\n",p,i);getch();exit(0);}
	
		for(int j=i;j<=n;j++) A[i][j]=A[i][j]/p;//aii=1
		for(int ii=i+1;ii<n;ii++)//row change down
		{
			p=A[ii][i];
			for(int j=i;j<=n;j++) A[ii][j]=A[ii][j]-p*A[i][j]; 
		}
	}
	
	x[n-1]=A[n-1][n]/A[n-1][n-1];//backward sweep on [U]
	for(int i=n-2;i>=0;i--)
	{
		p=0;
		for(int j=i+1;j<n;j++) p=p+x[j]*A[i][j];
		x[i] = (A[i][n] - p)/A[i][i];
	}Tx=x[0];	Ty=x[1];	Txx=x[2];	Tyy=x[3];	Txy=x[4];
} else {Tx=0;	Ty=0;	Txx=0;	Tyy=0;	Txy=0;}
}
*/


//-----------**************---------------operator descriptions----------*************-----------------
	vec vec::operator+(vec b)
	{
		vec r;
		//r.To=To+b.To;	
		r.Tc=Tc+b.Tc;	r.i=i+b.i;	r.Tx=Tx+b.Tx;	
		r.Tf=Tf+b.Tf;	r.j=j+b.j;	r.Ty=Ty+b.Ty;
		return r;			  
	}
//-----------**************---------------operator descriptions----------*************-----------------
	vec vec::operator-(vec b)
	{
		vec r;
		//r.To=To-b.To;	
		r.Tc=Tc-b.Tc;	r.i=i-b.i;	r.Tx=Tx-b.Tx;	
		r.Tf=Tf-b.Tf;	r.j=j-b.j;	r.Ty=Ty-b.Ty;
				
		return r;			  
	}
//-----------**************---------------operator descriptions----------*************-----------------
	float vec::operator*(vec b)
	{
		float r;
		r=i*b.i+j*b.j;
		return r;			  
	}
//-----------**************---------------operator descriptions----------*************-----------------
	float vec::operator%(vec b)
	{
		float r;
		r=i*b.j-j*b.i;
		return r;			  
	}
//-----------**************---------------operator descriptions----------*************-----------------
	vec vec::operator*(float b)
	{
		vec r;
		//r.To=To*b;
		r.Tc=Tc*b;	r.i=i*b;	r.Tx=Tx*b;	
		r.Tf=Tf*b;	r.j=j*b;	r.Ty=Ty*b;	
		return r;			  
	}
//-----------**************---------------operator descriptions----------*************-----------------
	vec vec::operator/(float b)
	{
		vec r;
		//r.To=To/b;		
		r.Tc=Tc/b;	r.Tx=Tx/b;	r.i=i/b;	
		r.Tf=Tf/b;	r.Ty=Ty/b;r.j=j/b;
		return r;			  
	}
//-----------**************---------------operator descriptions----------*************-----------------
	vec operator*(float b,vec a){return a*b;}
//-----------**************---------------operator descriptions----------*************-----------------
//-----------**************---------------operator descriptions----------*************-----------------
//-----------**************---------------operator descriptions----------*************-----------------
//-----------**************---------------operator descriptions----------*************-----------------
