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
int main()
{
float r,T,x;
fp = fopen ("boundry.txt","w");
float xo=5,yo=5,L=10,H=5;int X=10,Y=5;
Rec(xo,yo,L,H,X,Y);

//Cr(xo,yo,r,O,z);
//Cr(xo,yo,.5,32,0.5);
//Sq(xo,yo,10,8);
//Cr(xo,yo,1,32,0);

fprintf(fp,"0\n");
fclose(fp);
exit(0);
}
