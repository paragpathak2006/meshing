#include <iostream>
#include<math.h>
#include<stdlib.h>
#include<conio.h>
#include<string.h>
using namespace std;
float max(float a,float b){if(a>b)return a;else return b;}

int main()
{
FILE *fp;
const int X=10,Y=10;string p;
float Re=100,R1=1,R2=1,L=1,H=1,dx=L/X,dy=H/Y,gamma=L/Re,ae,aw,an,as,aP;
float u[X+2][Y+2],v[X+2][Y+2],apu[X+2][Y+2],apv[X+2][Y+2],bp[X+2][Y+2],P[X+2][Y+2],Pcorr[X+2][Y+2];

// BC  for velocity
for (int i=0;i<=X+1;i++)for (int j=0;j<=Y+1;j++){u[i][j]=0;v[i][j]=0;P[i][j]=0;Pcorr[i][j]=0;bp[i][j]=0;}	
for (int i=0;i<=X+1;i++)u[i][Y+1]=1;

fp=fopen("temp.txt","w");
for (int i=0;i<=X+1;i++){for (int j=0;j<=Y+1;j++)fprintf(fp,"%i %i %f\n",i,j,u[i][j]);fprintf(fp,"\n");}fprintf(fp,"\n");
//for (int i=0;i<=X+1;i++){for (int j=0;j<=Y+1;j++)fprintf(fp,"%i %i %f\n",i,j,v[i][j]);fprintf(fp,"\n");}
fclose(fp);p=getch();if(p=="q")exit(0);
//main iteration loop
for (int r=0;R1>1e-6;r++)
{
	for (int i=1;i<=X;i++)for (int j=1;j<=Y;j++)
	{
		//momentum equation
		if(i!=1)
		{
			aw=gamma/dx + 0.5*(u[i-1][j]+u[i][j]);	
			ae=gamma/dx - 0.5*(u[i+1][j]+u[i][j]);	

			as=gamma/dy + 0.5*(v[i-1][j]+v[i][j]);
			an=gamma/dy - 0.5*(v[i-1][j+1]+v[i][j+1]);

/*			aw=gamma/dx + max(u[i-1][j],0);	
			ae=gamma/dx - max(u[i+1][j],0);	

			as=gamma/dy + max((v[i-1][j]+v[i][j])/2,0);
			an=gamma/dy - max((v[i-1][j+1]+v[i][j+1])/2,0);
*/
			if(j==1)as=gamma*2/dx;if(j==Y) an=gamma*2/dx; 
//			ae=gamma/dx;aw=gamma/dx;	an=gamma/dx;	as=gamma/dx;		
			apu[i][j]=ae+aw+an+as;

			u[i][j]=(ae*u[i+1][j]+aw*u[i-1][j]+an*u[i][j+1]+as*u[i][j-1]+(P[i-1][j]-P[i][j])*dy)/apu[i][j];

		}
		
		if(j!=1)
		{
			aw=gamma/dx + 0.5*(u[i][j-1]+u[i][j]);		
			ae=gamma/dx - 0.5*(u[i+1][j-1]+u[i+1][j]);	

			as=gamma/dy + 0.5*(v[i][j-1]+v[i][j]);
			an=gamma/dy - 0.5*(v[i][j+1]+v[i][j]);

/*			aw=gamma/dx + max((u[i][j-1]/2+u[i][j]/2),0);		
			ae=gamma/dx - max((u[i+1][j-1]+u[i+1][j])/2,0);	

			as=gamma/dy + max(v[i][j-1],0);
			an=gamma/dy - max(v[i][j+1],0);
	
*/			if(i==1)aw=gamma*2/dx;if(i==X) ae=gamma*2/dx;	
			apv[i][j]=ae+aw+an+as;
			v[i][j]=(ae*v[i+1][j]+aw*v[i-1][j]+an*v[i][j+1]+as*v[i][j-1]+(P[i][j-1]-P[i][j])*dx)/apv[i][j]; 
		}
	}


	for (int i=1;i<=X;i++)for (int j=1;j<=Y;j++)
	bp[i][j]=(u[i][j]-u[i+1][j])*dy+(v[i][j]-v[i][j+1])*dx;
	
	for (int k=0;k<5;k++){
	for (int i=1;i<=X;i++)for (int j=1;j<=Y;j++)
	{
		ae=dy/apu[i+1][j]; 	aw=dy/apu[i][j]; 		
		an=dx/apv[i][j+1]; 	as=dx/apv[i][j];
		if(i==1)aw=0;if(i==X)ae=0;if(j==1)as=0;if(j==Y)an=0;
		aP=ae+aw+as+an; 		
		Pcorr[i][j]=(aw*Pcorr[i-1][j]+ae*Pcorr[i+1][j]+as*Pcorr[i][j-1]+an*Pcorr[i][j+1]+bp[i][j])/aP;
	}}
	
	//velocity correction:
	for (int i=1;i<=X;i++)for (int j=1;j<=Y;j++)
	{
		if(i!=1) u[i][j]=u[i][j]+dy*(Pcorr[i-1][j]-Pcorr[i][j])/apu[i][j];
		if(j!=1) v[i][j]=v[i][j]+dx*(Pcorr[i][j-1]-Pcorr[i][j])/apv[i][j];	
		P[i][j]=P[i][j]+.5*Pcorr[i][j];
	}
	R1=0;	for (int i=1;i<=X;i++)for (int j=1;j<=Y;j++)if(R1<fabs(bp[i][j])) R1=fabs(bp[i][j]);
	printf("r=%i R1=%e\n",r,R1);

}
	fp=fopen("u.txt","w");
	for (int i=1;i<=X;i++){for (int j=1;j<=Y;j++)fprintf(fp,"%f %f %f %f\n",i*dx-dx/2,j*dy-dy/2,u[i][j]/2+u[i+1][j]/2,v[i][j]/2+v[i][j+1]/2);fprintf(fp,"\n");}fprintf(fp,"\n");

	for (int i=1;i<=X;i++)fprintf(fp,"%f %f %f %f\n",i*dx-dx/2,0,0.0,0.0);fprintf(fp,"\n");
	for (int i=1;i<=X;i++)fprintf(fp,"%f %f %f %f\n",i*dx-dx/2,H,1.0,0.0);fprintf(fp,"\n");

	for (int j=1;j<=Y;j++)fprintf(fp,"%f %f %f %f\n",0.0,j*dy-dy/2,0.0,0.0);fprintf(fp,"\n");
	for (int j=1;j<=Y;j++)fprintf(fp,"%f %f %f %f\n",L,j*dy-dy/2,0.0,0.0);fprintf(fp,"\n");

	fclose(fp);p=getch();if(p=="q")exit(0);

return 0;
}
