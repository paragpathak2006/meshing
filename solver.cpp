#include <iostream>
#include<math.h>
#include<stdlib.h>
#include<conio.h>
#include<string.h>
using namespace std;
FILE *fp1,*fp2,*fp3;

class point
{public :float i,j;};

class traingle
{public:int N[3],p[3];};

int nmax=2,pmax,bmin,bmax,cmin,cmax;
const int cellmax=1000,pnmax=2200;
float pi=3.14159;
traingle cell[cellmax+10];point P[pnmax+10];

float d(point a,point b)//distace between pt a and b.
{
 	  return sqrt((a.i-b.i)*(a.i-b.i)+(a.j-b.j)*(a.j-b.j));
}
void seg(point a,point b,FILE *fp)//draw seg a to b
{
	fprintf(fp,"%f %f\n",a.i,a.j);
	fprintf(fp,"%f %f\n\n",b.i,b.j);
}
void tri(traingle t,FILE *fp)//draw triangle t.
{
	seg(P[t.p[0]],P[t.p[1]],fp);
	seg(P[t.p[1]],P[t.p[2]],fp);
	seg(P[t.p[2]],P[t.p[0]],fp);
}

int main()
{
	fp1 = fopen ("rectangle.txt","r");

	fscanf (fp1,"%i %i",&bmax,&pmax);
	for(int i=1;i<=pmax;i++)fscanf (fp1,"%f %f",&P[i].i,&P[i].j);

	fscanf (fp1,"%i",&nmax);

	if(pmax>pnmax){cout<<"points exeed limit";getch();exit(0);}	
	if(nmax>cellmax){cout<<"triangles exeed limit";getch();exit(0);}	

	for(int j=1;j<=nmax;j++)
	{
	 	fscanf (fp1,"%i %i %i",&cell[j].p[0],&cell[j].p[1],&cell[j].p[2]);
		fscanf (fp1,"%i %i %i",&cell[j].N[0],&cell[j].N[1],&cell[j].N[2]);
	}fclose(fp1);


printf("Points exterior-->\tbmin=%i\tbmax=%i\n",1,bmax);
printf("Points interior-->\tpmin=%i\tpmax=%i\n",1+bmax,pmax);
printf("triangles-->\t\tnmin=%i\tnmax=%i\n",1,nmax);

fp2 = fopen ("temp.txt","w");
for(int j=1;j<=nmax;j++) tri(cell[j],fp2);fclose(fp2);
getch();
return 0;
}
