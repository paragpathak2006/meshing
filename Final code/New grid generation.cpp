#include "output.cpp"
FILE *fp1,*fp2,*fp3;
float gauss(vec P,float t,float &S)
{
	const float L=1.0,Ag=.5*L,Bg=.25*L,Tg=0.1,Co=4.0,C1=0;//1/sqrt(kf);
	float s=.2,xp,x1p,E,E1,E2,X,X1,X2,Y,Y2,T,St,Sx,Sy,x=P.i,y=P.j;

	xp=Ag+Bg*sin(2*pi*t/Tg);	x1p=Bg*cos(2*pi*t/Tg)*2*pi/Tg;

	E=exp((xp-x)*(x-xp)/(s*s));	E1=-2*(x-xp)*E/(s*s);	E2=-2*(E+(x-xp)*E1)/(s*s);
	X=x*(L-x);					X1=L-2*x;				X2=-2;
	Y=sin(pi*y/L);				Y2=-pi*pi*Y/(L*L);
	T=Co+C1*t;					
	
	St=X*Y*(E*C1-T*E1*x1p);	Sx=(E2*X+E*X2+2*E1*X1)*Y*T;	Sy=E*X*Y2*T;	S=St-kf*(Sx+Sy);
	return E*X*Y*T;
}
void inputpoints(FILE *fp)
{	
	for(int i=1;i<=pmax;i++)
	{
		if(i<=bmax)P[i].BN=true;	else P[i].BN=false;
	}
}	

void triangle::inputcell(FILE *fp,int jj)
{
	initialize(1,2,3,jj);
	N[0]=0;N[1]=0;N[2]=0;

	if(volmin>cell[j].Vol)volmin=cell[j].Vol;C.Tc=250;
//	C.Tc=gauss(C,0,S);
	Tf=C.Tc;Tc=C.Tc;//To=C.Tc;
}

 
int main()
{pass=0;
bmax=3;pmax=3;
P[1].i=0;		P[2].i=1;	P[3].i=0.5;		
P[1].j=0;		P[2].j=0;	P[3].j=sqrt(3)/2;

inputpoints(fp1);

nmax=1;cell[1].inputcell(fp1,1);

string p="y";
printf("max boundry node=%i \t nodes=%i \t triangles=%i\n",bmax,pmax,nmax);getch();
reinitialize();Nmax=nmax;
//output();getch();

printf("nmax=%i pmax=%i\n\n",nmax,pmax);

cell[1].split=-1;celloperations();
printf("nmax=%i pmax=%i\n\n",nmax,pmax); 
//output();getch();

for(int j=1;j<=nmax;j++)cell[j].split=-1;celloperations();
printf("nmax=%i pmax=%i\n\n",nmax,pmax); 
//output();getch();

for(int j=1;j<=nmax;j++)cell[j].split=-1;celloperations();
printf("nmax=%i pmax=%i\n\n",nmax,pmax); 
//output();getch();

for(int j=1;j<=nmax;j++)cell[j].split=-1;celloperations();
printf("nmax=%i pmax=%i\n\n",nmax,pmax); 
//output();getch();

for(int j=1;j<=nmax;j++)cell[j].split=-1;celloperations();
printf("nmax=%i pmax=%i\n\n",nmax,pmax); 
//output();getch();

float L=P[153].i-P[151].i,dx=L/10;
for(int i=1;i<=pmax;i++)if(P[i].i<P[151].i&&P[i].i>P[151].i-dx*2/3)P[i].i=P[151].i;
for(int i=1;i<=pmax;i++)if(P[i].i>P[153].i&&P[i].i<P[153].i+dx*2/3)P[i].i=P[153].i;

for(int i=1;i<=pmax;i++)
{if(P[i].i==P[151].i||P[i].i==P[153].i||P[i].j==P[153].j||P[i].j==P[252].j)
P[i].BN=1;else P[i].BN=0;}

for(int i=1;i<=pmax;i++){
if(P[i].j<=P[151].j&&P[i].j>=P[151].j-L&&P[i].i>=P[151].i&&P[i].i<=P[153].i)
P[i].coarse=0;else P[i].coarse=1;}

for(int j=1;j<=nmax;j++)
{
if(P[cell[j].p[0]].coarse==1||P[cell[j].p[1]].coarse==1||P[cell[j].p[2]].coarse==1)
cell[j].del();else cell[j].deleted=0;
}

printf("nmax=%i pmax=%i\n\n",nmax,pmax); 
output();getch();


int nk;
while(P[pmax].coarse) pmax--;
for(int i=1;i<pmax;i++)
if(P[i].coarse)
{
 	P[i]=P[pmax];
	for(int j=1;j<=P[i].t[0];j++)
	{
		nk=P[i].t[j];
		for(int k=0;k<3;k++)if(cell[nk].p[k]==pmax)	cell[nk].p[k]=i;
	}
	P[pmax].coarse=1;
	while(P[pmax].coarse) {pmax--;}
}

//for(int j=1;j<nmax;j++)
//if(cell[j].deleted){while(cell[nmax].deleted) nmax--;cell[j]=cell[nmax];cell[j].j=j;nmax--;}

renumber();
int r,s;pass=1;
for(int j=1;j<=nmax;j++) cell[j].setneighbour();
for(int j=1;j<=nmax;j++) 
{
for(int k=0;k<3;k++) 
{
if(k==0) {r=1;s=2;}if(k==1) {r=2;s=0;}if(k==2) {r=0;s=1;}
if(!P[cell[j].p[k]].BN&&P[cell[j].p[r]].BN&&P[cell[j].p[s]].BN)cell[j].N[k]=0;
}
}

vec a;float fd=1;
while(fd>1e-7)
{
	fd=0;
	for(int i=1;i<=pmax;i++)
	{
//		printf("\ni=%i BN=%i\n",i,P[i].BN);
		if(!P[i].BN)
		{
//		printf("\n\ni=%i\n",i);
			 a.i=0;a.j=0;
			for(int j=1;j<=P[i].t[0];j++)
			{
				nk=P[i].t[j];
//				printf("nk=%i ",nk);
				for(int k=0;k<3;k++)
				{
				if(cell[nk].p[k]!=i) {a.i=a.i+P[cell[nk].p[k]].i-P[i].i;
			a.j=a.j+P[cell[nk].p[k]].j-P[i].j;	}
				}
			}
			a.i=a.i/(2*P[i].t[0]);	a.j=a.j/(2*P[i].t[0]);
			P[i].i=P[i].i+a.i;	P[i].j=P[i].j+a.j;

			fd=fd+a.i*a.i+a.j*a.j;
		}
	}
	printf("fd=%f\n\n",fd); 

}

printf("nmax=%i pmax=%i\n\n",nmax,pmax); 
output();getch();

float fdx=P[42].i,fdy=P[95].j;
for(int i=1;i<=pmax;i++){P[i].i=P[i].i-fdx;P[i].j=P[i].j-fdy;}

fdx=P[52].i,fdy=P[138].j;
for(int i=1;i<=pmax;i++)
{P[i].i=P[i].i/fdx;P[i].j=P[i].j/fdy;}
cell[51].setneighbour();cell[58].setneighbour();

FILE *fp=fopen("sqper.txt","w");
fprintf(fp,"%i\n",pmax);
for(int i=1;i<=pmax;i++)fprintf(fp,"%f %f %i\n",P[i].i,P[i].j,P[i].BN);

fprintf(fp,"%i\n",nmax);
for(int j=1;j<=nmax;j++) 
{
for(int k=0;k<3;k++) if(cell[j].N[k]>nmax)cell[j].N[k]=0;
fprintf(fp,"%i %i %i ",cell[j].p[0],cell[j].p[1],cell[j].p[2]);
fprintf(fp,"%i %i %i\n",cell[j].N[0],cell[j].N[1],cell[j].N[2]);
}
fclose(fp);

printf("%i\n",pmax);
for(int i=1;i<=pmax;i++)printf("i=%i %f %f %i\n",i,P[i].i,P[i].j,P[i].BN);

printf("%i\n",nmax);
for(int j=1;j<=nmax;j++) 
{
printf("j=%i\tp0=%i  p1=%i  p2=%i\t\t",j,cell[j].p[0],cell[j].p[1],cell[j].p[2]);
printf("N0=%i  N1=%i  N2=%i\n",cell[j].N[0],cell[j].N[1],cell[j].N[2]);
if(j%100==0)getch();
}


getch();
return 0;
}	
