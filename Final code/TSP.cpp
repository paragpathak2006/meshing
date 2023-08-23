#include <iostream>
#include<math.h>
#include<stdlib.h>
#include<conio.h>
#include<string.h>
using namespace std;
class point
{public :float i,j;int N;};
FILE *sq,*bd,*op;
const int pnmax=20;
point P[pnmax+10],a,b,c,e;

float d(point a,point b)//distace between pt a and b.
{
 	  return sqrt((a.i-b.i)*(a.i-b.i)+(a.j-b.j)*(a.j-b.j));
}
point M(point a,point b)//point a-b.
{
 	  point d;
 	  d.i=a.i-b.i;
 	  d.j=a.j-b.j;
 	  return d;
}
point CC(point a,point b,point c)//calulates circumcentre of point a,b,c
{
	point d;
	float a1,a2,b1,b2,c1,c2,x1,x2,x3,y1,y2,y3;
	x1=a.i;x2=b.i;x3=c.i;
	y1=a.j;y2=b.j;y3=c.j;
	a1=2*(x1-x2);
	a2=2*(x3-x2);
	b1=2*(y1-y2);
	b2=2*(y3-y2);
	c1=-(x1*x1 - x2*x2 + y1*y1 - y2*y2);
	c2=-(x3*x3 - x2*x2 + y3*y3 - y2*y2);
	d.i=-(c1*b2-c2*b1)/(a1*b2-a2*b1);
	d.j=-(c1*a2-c2*a1)/(b1*a2-b2*a1);
	return d;
}

int main()
{
int ii;float min,l;
P[0].i=0;P[0].j=0;
P[1].i=1;P[1].j=2;
P[2].i=6;P[2].j=4;
P[3].i=4;P[3].j=2;
P[4].i=9;P[4].j=7;
P[5].i=11;P[5].j=0;
P[6].i=0;P[6].j=4;
P[7].i=3;P[7].j=6;

P[0].N=1;P[1].N=2;P[2].N=0;

for(int j=3;j<6;j++)
{
	min=1e5;
	for(int i=0;i<j;i++)
	{
		l=d(P[j],P[i])+d(P[j],P[P[i].N])-d(P[i],P[P[i].N]);
		if(l<min){min=l;ii=i;}
		printf("j=%i i=%i in=%i Pji=%f Pjin=%f Piin=%f l=%f\n",j,i,P[i].N,d(P[j],P[i]),d(P[j],P[P[i].N]),d(P[i],P[P[i].N]),l);
	}printf("\n\n");
	P[j].N=P[ii].N;P[ii].N=j;

	FILE *fp=fopen("TSP.txt","w");
	ii=0;
	do
	{
		fprintf(fp,"%f %f %i\n",P[ii].i,P[ii].j,ii);
		ii=P[ii].N;
	}while(ii!=0);
	fprintf(fp,"%f %f %i\n",P[ii].i,P[ii].j,ii);
	fclose(fp);getch();
}
FILE *fp=fopen("TSP.txt","w");

ii=0;
do
{
	fprintf(fp,"%f %f %i\n",P[ii].i,P[ii].j,ii);
	ii=P[ii].N;
}while(ii!=0);
fprintf(fp,"%f %f %i\n",P[ii].i,P[ii].j,ii);

l=0;ii=0;
do
{
	l=l+d(P[ii],P[P[ii].N]);
	ii=P[ii].N;
}while(ii!=0);
l=l+d(P[ii],P[P[ii].N]);

printf("dmin_algo=%f\n",l);

P[0].N=1;
P[1].N=6;
P[6].N=2;
P[2].N=4;
P[4].N=5;
P[5].N=3;
P[3].N=0;
l=0;
l=d(P[0],P[1]) + d(P[1],P[6]) + d(P[6],P[2]) + d(P[2],P[4]) + d(P[4],P[5]) + d(P[5],P[3]) + d(P[3],P[0]);
l=d(P[0],P[1]) + d(P[1],P[6]) + d(P[6],P[5]) + d(P[5],P[4]) + d(P[4],P[2]) + d(P[2],P[3]) + d(P[3],P[0]);
printf("dmin_test=%f\n",l);



//printf("\n",10.2+7.23-);//1-4/
//printf("\n",6.35+7.23);//2-4
//printf("\n",6.35+7.27);//2-3

getch();
return 0;
}
