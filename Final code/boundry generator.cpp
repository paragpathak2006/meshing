#include <iostream>
#include<math.h>
#include<stdlib.h>
#include<conio.h>
#include<string.h>
using namespace std;
FILE *fp;
const float pi=3.14159;
void Rec(float xo,float yo,float L,float H,int X,int Y)
{
	xo=xo-L/2;
	yo=yo-H/2;
	fprintf(fp,"%i\n",(X+Y)*2);
	for(int j=0;j<=Y-1;j++) fprintf(fp,"%f %f\n",xo,yo+H*j/Y);
	for(int i=0;i<=X-1;i++) fprintf(fp,"%f %f\n",xo+L*i/X,yo+H);
	for(int j=Y;j>=1;j--) fprintf(fp,"%f %f\n",xo+L,yo+H*j/Y);
	for(int i=X;i>=1;i--) fprintf(fp,"%f %f\n",xo+L*i/X,yo);
}
void Sq(float xo,float yo,float L,int X){Rec(xo,yo,L,L,X,X);}
void Cr(float xo,float yo,float r,int O,float z)
{
	fprintf(fp,"%i\n",O);
	for(int j=0;j<=O-1;j++)
	{
		fprintf(fp,"%f %f\n",xo+r*sin(2*pi*j/O+z*2*pi/64),yo+r*cos(2*pi*j/O+z*2*pi/64));
	}
}
void tri(float xo,float yo,float L,int X)
{
	fprintf(fp,"%i\n",3*X);
	for(int j=0;j<X;j++)fprintf(fp,"%f %f\n",xo+L*j/X,yo);
	for(int j=0;j<X;j++)fprintf(fp,"%f %f\n",xo+L-L*j*.5/X,yo+L*j*sqrt(3)*.5/X);
	for(int j=0;j<X;j++)fprintf(fp,"%f %f\n",xo+L*.5-L*j*.5/X,yo+L*sqrt(3)*.5-L*j*sqrt(3)*.5/X);
}

int main()
{
float r,T,x;
fp = fopen ("boundry.txt","w");
float xo=0,yo=0,L=1,H=1;int X=10,Y=10;
//Rec(xo,yo,L,H,X,Y);
Cr(xo,yo,1,32,0);
Cr(xo,yo,.1,32,0.5);
//Sq(xo,yo,L,X);
//tri(xo,yo,L,X);
//Cr(xo,yo,1,32,0);

fprintf(fp,"0\n");
fclose(fp);
exit(0);
}
