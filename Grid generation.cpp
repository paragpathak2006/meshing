#include <iostream>
#include<math.h>
#include<stdlib.h>
#include<conio.h>
#include<string.h>
using namespace std;
class point
{public :float x,y,f;};
FILE *sq,*bd,*op;
class traingle
{public :
	point V,C;//V-circumcentre,C-
	int N[3],p[3],e,m,a;//e--existence,m-cell number modification,a-active triangles
	float r;
};
int deleted,nn=2,nmax=2,pmax,s,nk,r,convergence=0,N,bmax,cmin,cmax;
const int cellmax=10000,pnmax=22000,bmin=4;
float ab,sx,sy,pi=3.14159,tolerence=-.5e-4;
traingle cell[cellmax+10];
point P[pnmax+10],a,b,c,e;

float d(point a,point b)//distace between pt a and b.
{
 	  return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
float C(point a,point b)//cross product of point a and b.
{
 	  return (a.x*b.y-a.y*b.x);
}
point M(point a,point b)//point a-b.
{
 	  point d;
 	  d.x=a.x-b.x;
 	  d.y=a.y-b.y;
 	  return d;
}
void seg(point a,point b)//draw seg a to b
{
	fprintf(op,"%f %f\n",a.x,a.y);
	fprintf(op,"%f %f\n\n",b.x,b.y);
}
void tri(traingle t)//draw triangle t.
{
seg(P[t.p[0]],P[t.p[1]]);
seg(P[t.p[1]],P[t.p[2]]);
seg(P[t.p[2]],P[t.p[0]]);
}
point CC(point a,point b,point c)//calulates circumcentre of point a,b,c
{
	point d;
	float a1,a2,b1,b2,c1,c2,x1,x2,x3,y1,y2,y3;
	x1=a.x;x2=b.x;x3=c.x;
	y1=a.y;y2=b.y;y3=c.y;
	a1=2*(x1-x2);
	a2=2*(x3-x2);
	b1=2*(y1-y2);
	b2=2*(y3-y2);
	c1=-(x1*x1 - x2*x2 + y1*y1 - y2*y2);
	c2=-(x3*x3 - x2*x2 + y3*y3 - y2*y2);
	d.x=-(c1*b2-c2*b1)/(a1*b2-a2*b1);
	d.y=-(c1*a2-c2*a1)/(b1*a2-b2*a1);
	return d;
}
void bayeralgo(int min,int max)//implemented the bayer algorithm .Given an existing voronoi diagram 
//from points 0 to min-1 .Put points from min to max.
{
	for(int i=min;i<=max;i++)//additon of boundry and cavity points.
	{
		cell[0].e=1;
		cell[0].m=0;
	
		for(int j=1;j<=nmax;j++)//reinitializing cell parameters
		{
			cell[j].V = CC(P[cell[j].p[0]],P[cell[j].p[2]],P[cell[j].p[1]]);
			cell[j].r = d(cell[j].V,P[cell[j].p[0]]);
			cell[j].e = 1;
			cell[j].m = 0;
		}
		//-------------------------------------------------------------------------------------------------------------------------
		for(int j=1;j<=nmax;j++)//deleting unrequired cells
		{
			if (d(P[i],cell[j].V)+tolerence<=cell[j].r)
			{
				cell[j].e=0;
				cell[j].m=1;
	
				cell[cell[j].N[0]].m=1;
				cell[cell[j].N[1]].m=1;
				cell[cell[j].N[2]].m=1;
			}
		}
		//-------------------------------------------------------------------------------------------------------------------------
	/*	for(int j=1;j<=nmax;j++)//deletion refined for resoving circular set of points 
		{
			if (cell[j].e==0)
			{
			for(int k=0;k<=2;k++)
			{
			nk=cell[j].N[k];
	
			if(cell[nk].e==1)
			{
				if(k==0) {r=1;s=2;}
				if(k==1) {r=0;s=2;}
				if(k==2) {r=0;s=1;}
				
				a=M(P[cell[j].p[r]],P[i]);
				b=M(P[cell[j].p[s]],P[i]);
				c=M(P[cell[j].p[r]],P[cell[j].p[k]]);
				e=M(P[cell[j].p[s]],P[cell[j].p[k]]);
				
				if(C(a,b)*C(c,e)<=0)
				{
	   			cell[nk].e=0;
				cell[nk].m=1;
				cell[cell[nk].N[0]].m=1;
				cell[cell[nk].N[1]].m=1;
				cell[cell[nk].N[2]].m=1;
				}
			}
			}
			}
		}*/	
		nn=nmax;
		//-------------------------------------------------------------------------------------------------------------------------
		for(int j=1;j<=nn;j++)//new cells creation
		{
			if (cell[j].e==0)//looking at deleted cells
			{
				for(int k=0;k<=2;k++)
				{nk=cell[j].N[k];
				if (cell[nk].e==1)//seeing undeleted neighbours of deleted cells
				{
					nmax++;
					if(nmax==cellmax){cout<<"traigles exeed limit";getch();exit(0);}	
				    
					cell[nmax]=cell[j];
					cell[nmax].p[k]=i;
					cell[nmax].e=1;
					cell[nmax].m=1;
				}
				}
			}
		}
		//-------------------------------------------------------------------------------------------------------------------------
		for(int j=1;j<=nmax;j++)//cells number reassginment
		{
			if (cell[j].e==0)
			{ 
			cell[j]=cell[nmax];
			nmax--;
			cell[j].m=1;
			}
			if(cell[j].m==1)
			{
				cell[j].N[0]=0;
				cell[j].N[1]=0;
				cell[j].N[2]=0;
			}
		}nn=nmax;	
		//-------------------------------------------------------------------------------------------------------------------------
		for(int j=1;j<=nn;j++)//neighbour realocation
		{
		if(cell[j].m==1)
		{
			for(int ii=0;ii<=nmax;ii++) cell[ii].e=1;
			for(int k=0;k<=2;k++)
			{
				for(int ii=1;ii<=nmax;ii++)//consider another triangle score asgingment 
				{
					if(j==ii) ii++;
					for(int kk=0;kk<=2;kk++) {if(cell[j].p[k]==cell[ii].p[kk]){cell[ii].e=cell[ii].e*(k+1)*10;}}
				}
			}
			for(int ii=1;ii<=nmax;ii++) 
			{
				if(cell[ii].e==600) cell[j].N[0]=ii;
				if(cell[ii].e==300) cell[j].N[1]=ii;
				if(cell[ii].e==200) cell[j].N[2]=ii;
			}		
		}
		}
	
	}
}
void polygon(int &jp,int &j,int i)//tracing a polygonal line around point i.
{
	int fbj,nr,ns,k,r,s;
	for(k=0;k<=2;k++)
	{
	if(cell[j].p[k]==i) break;
	}
	
	if(k==0) {r=1;s=2;}
	if(k==1) {r=0;s=2;}
	if(k==2) {r=0;s=1;}
	
	nr=cell[j].N[r];
	ns=cell[j].N[s];
	fbj=j;
	if(nr!=jp) j=nr;
	if(ns!=jp) j=ns;
	jp=fbj;
}

int main()
{
bd = fopen ("boundry.txt","r+");
int fmax,fmin;

fscanf (bd,"%i",&bmax);bmax=bmax+3;
for(int i=bmin;i<=bmax;i++) fscanf (bd,"%f %f",&P[i].x,&P[i].y);

fscanf (bd,"%i",&cmax);
if(cmax==0)  {pmax=bmax+1;cmax=bmax;cmin=bmin;}
else 
{
cmax=cmax+bmax;cmin=bmax+1;
for(int i=cmin;i<=cmax;i++) fscanf (bd,"%f %f",&P[i].x,&P[i].y);
pmax=cmax+1;
}

float xmax,xmin,ymax,ymin,dx,dy;
xmax=P[bmin].x;
xmin=P[bmin].x;
ymax=P[bmin].y;
ymin=P[bmin].y;
for(int i=bmin;i<=cmax;i++)
{
	if(xmax<P[i].x) xmax=P[i].x;
	if(xmin>P[i].x) xmin=P[i].x;
	if(ymax<P[i].y) ymax=P[i].y;
	if(ymin>P[i].y) ymin=P[i].y;
} 
a.x=(xmax+xmin)/2;	dx=(xmax-xmin);
a.y=(ymax+ymin)/2;	dy=(ymax-ymin);
P[0].x=a.x-dx;	P[1].x=a.x+dx;	P[2].x=a.x-dx;	P[3].x=a.x+1.2*dx;
P[0].y=a.y-dy;	P[1].y=a.y-dy;	P[2].y=a.y+dy;	P[3].y=a.y+1.2*dy;

if(pmax>pnmax){cout<<"points exeed limit";getch();exit(0);}	
cout<<"bmin="<<bmin<<"  bmax="<<bmax<<"\ncmin="<<cmin<<"  cmax="<<cmax<<"\npmax="<<pmax<<endl;getch();

cell[0].e=1;

cell[1].N[0]=2;	cell[1].p[0]=0;
cell[1].N[1]=0;	cell[1].p[1]=1;
cell[1].N[2]=0;	cell[1].p[2]=2;

cell[2].N[0]=1;	cell[2].p[0]=3;
cell[2].N[1]=0;	cell[2].p[1]=2;
cell[2].N[2]=0;	cell[2].p[2]=1;

//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
bayeralgo(bmin,cmax);//boundry trigulated
op = fopen ("temp.txt","w");
for(int j=1;j<=nmax;j++)//consider a triangle 
{
cell[j].e=1;
cell[j].m=1;
cell[j].a=1;
tri(cell[j]);
}
fclose(op);
nn=nmax;
printf("boundry triangulation compleat\n"); 
getch();

for(int j=1;j<=nn;j++)//exterior cells deletion
{
cell[j].e=1;
cell[j].m=0;
	for(int i=0;i<=3;i++) 
	{
	for(int k=0;k<=2;k++)
	{
		if(cell[j].p[k]==i&&cell[j].e==1)
		{ 
			cell[j].e=0;
		}
	}
	}
}
//-------------------------------------------------------------------------------------------------------------------------
if(cmax!=bmax){for(int j=1;j<=nn;j++)//cavity cells deletion
{
	for(int k=0;k<=2;k++)
	{
	for(int i=cmin;i<=cmax;i++) 
	{
		if(cell[j].p[k]==i&&cell[j].e!=0)
		{ 
			cell[j].e=cell[j].e*10;
			if(cell[j].e==1000) 
			{
				cell[j].e=0;
			}
		}
	}
	}
}}
//-------------------------------------------------------------------------------------------------------------------------
for(int j=1;j<=nmax;j++)//cells number reassginment
{
	if (cell[j].e==0)
	{ 
		if(cell[nmax].e!=0) cell[j]=cell[nmax];
		if(cell[nmax].e==0&&j!=nmax) j--;
		nmax--;
	}
cell[j].N[0]=0;
cell[j].N[1]=0;
cell[j].N[2]=0;

}
//-------------------------------------------------------------------------------------------------------------------------
for(int j=1;j<=nmax;j++)//neighbour re-allocation
{
	for(int i=0;i<=nmax;i++) cell[i].e=1;

	for(int k=0;k<=2;k++)
	{
		for(int i=1;i<=nmax;i++)//consider another triangle score asgingment 
		{
			if(j==i) i++;
			for(int kk=0;kk<=2;kk++) {if(cell[j].p[k]==cell[i].p[kk]){cell[i].e=cell[i].e*(k+1)*10;}}
		}
	}
	for(int i=1;i<=nmax;i++)//neibour asg 
	{
		if(cell[i].e==600) cell[j].N[0]=i; 
		if(cell[i].e==300) cell[j].N[1]=i;
		if(cell[i].e==200) cell[j].N[2]=i;  
	}
}
op = fopen ("temp.txt","w");
for(int j=1;j<=nmax;j++)//output 
{
cell[j].e=1;
cell[j].m=1;
cell[j].a=1;
tri(cell[j]);
}
fclose(op);
printf("boundry triangulation cleaned up\n"); 
getch();
//-------------------------------------------------------------------------------------------------------------------------
//           Automatic point insertion
//-------------------------------------------------------------------------------------------------------------------------
for(int i=bmin;i<=cmax;i++)//difining PDF for boundry nodes
{
	if(i!=bmin&&i!=bmax&&i!=cmin&&i!=cmax) P[i].f=(d(P[i],P[i-1])+d(P[i],P[i+1]))/2;

	if(i==bmin) P[bmin].f=(d(P[bmin],P[bmin+1])+d(P[bmin],P[bmax]))/2;
	if(i==bmax) P[bmax].f=(d(P[bmax],P[bmax-1])+d(P[bmax],P[bmin]))/2;

	if(i==cmin&&cmax!=cmin) P[cmin].f=(d(P[cmin],P[cmin+1])+d(P[cmin],P[cmax]))/2;
	if(i==cmax&&cmax!=cmin) P[cmax].f=(d(P[cmax],P[cmax-1])+d(P[cmax],P[cmin]))/2;

	if(cmax==cmin) P[cmax].f=P[bmax].f*1;
}
fmax=cmax;
//-------------------------------------------------------------------------------------------------------------------------
int i,count=0;point b;
float dd,O;
while(convergence==0)
{
	count++;
	fmin=fmax+1;
	i=fmin-1;
	for(int j=1;j<=nmax;j++)//points generator
	{
	 	if(cell[j].a==1)//if triangle is active.
		{
			i++;
			if(i==pnmax){cout<<"points exeed limit";getch();exit(0);}	
			P[i].x=(P[cell[j].p[0]].x+P[cell[j].p[1]].x+P[cell[j].p[2]].x)/3;	
			P[i].y=(P[cell[j].p[0]].y+P[cell[j].p[1]].y+P[cell[j].p[2]].y)/3;	
			P[i].f=(P[cell[j].p[0]].f+P[cell[j].p[1]].f+P[cell[j].p[2]].f)/3;	

			s=1;
			if(d(P[i],P[cell[j].p[0]])<.7*P[i].f)s=0;
			if(d(P[i],P[cell[j].p[1]])<.7*P[i].f)s=0;
			if(d(P[i],P[cell[j].p[2]])<.7*P[i].f)s=0;
			for(int ii=4;ii<=fmax;ii++)
			{
				if(d(P[i],P[ii])<.7*P[i].f&&P[i].f>P[ii].f) s=0;
			}

			for(int ii=fmax+1;ii<=i-1;ii++)
			{
				if(d(P[i],P[ii])<.7*P[i].f&&P[i].f>P[ii].f) s=0;
				if(d(P[i],P[ii])<.7*P[ii].f&&P[i].f<P[ii].f) {s=0;P[ii]=P[i];}
			}

			if(s==0) i--;	
		}

	}
	fmax=i;
	if(fmax==fmin-1){printf("No more prospective points left.Solution converged\n");break;} 
	bayeralgo(fmin,fmax);//using bayer algorthm for inserting generated points into mesh.

	convergence=1;
	for(int j=1;j<=nmax;j++)//cell activation critaria. 
	{
		ab=(P[cell[j].p[0]].f+P[cell[j].p[1]].f+P[cell[j].p[2]].f)/3;
		cell[j].a=0;
		if(d(P[cell[j].p[0]],P[cell[j].p[1]])>1.5*ab) {cell[j].a=1;convergence=0;} 
		if(d(P[cell[j].p[1]],P[cell[j].p[2]])>1.5*ab) {cell[j].a=1;convergence=0;}
		if(d(P[cell[j].p[2]],P[cell[j].p[0]])>1.5*ab) {cell[j].a=1;convergence=0;}
	}

	op = fopen ("temp.txt","w");
	for(int j=1;j<=nmax;j++)//output
	{
	tri(cell[j]);
	}
	fclose(op);
	printf("iteration=%i compleat\n",count);getch();
	if(convergence==1) {printf("All traingles inactive.Solution converged\n%i triangles ,%i nodes generated\n",nmax,fmax-3);break;}

}//printf("nmax=%i fmax=%i\n",nmax,fmax);

//getch();
//-------------------------------------------------------------------------------------------------------------------------
//           Smoothning
//-------------------------------------------------------------------------------------------------------------------------

int jp,J,k,nr,ns,ss,j;
float sx,sy,px=0,py=0,A=.9,S; 
ss=1;
while(ss==1)
{
	ss=0;S=0;
	for(int i=pmax;i<=fmax;i++)//smoothning
	{

		N=0;sx=0;sy=0;
		for(j=1;j<=nmax;j++) 
		{
		for(k=0;k<=2;k++)
		{
			if(i==cell[j].p[k]) break;
		}
			if(k!=3) break;
		}
		N=0;J=j;
		
		if(k==0) {r=1;s=2;}
		if(k==1) {r=2;s=0;}
		if(k==2) {r=0;s=1;}

		nr=cell[j].N[r];
		jp=nr;
		do{
			polygon(jp,j,i);
			for(k=0;k<=2;k++) {if(i==cell[j].p[k]) {break;}}

			if(k==0) {r=1;s=2;}
			if(k==1) {r=2;s=0;}
			if(k==2) {r=0;s=1;}
			
			N++;
			sx=sx+(P[cell[j].p[r]].x+P[cell[j].p[s]].x)/2;
			sy=sy+(P[cell[j].p[r]].y+P[cell[j].p[s]].y)/2;
		}while(j!=J);

		px=P[i].x;
		py=P[i].y;
		P[i].x=P[i].x*(1-A)+A*sx/N;
		P[i].y=P[i].y*(1-A)+A*sy/N;

		if(fabs(P[i].x-px)>1e-6*P[i].f) ss=1;
		if(fabs(P[i].y-py)>1e-6*P[i].f) ss=1;

		S=S+fabs(P[i].x-px);
		S=S+fabs(P[i].y-py);//printf("Sum=%f\n",S);
//getch();

	}
//printf("Sum=%e\n",S);
}
printf("mesh with %i triangles  & %i nodes created\n",nmax,fmax-3);
getch();

op = fopen ("temp.txt","w");
for(int j=1;j<=nmax;j++)//output
{
tri(cell[j]);
}
fclose(op);

op = fopen ("rectangle.txt","w");

fprintf (op,"%i %i\n",bmax-3,fmax-3);
for(int i=4;i<=fmax;i++)	fprintf (op,"%f %f\n",P[i].x,P[i].y);

fprintf (op,"%i\n",nmax);
for(int j=1;j<=nmax;j++)
{
	fprintf (op,"%i %i %i ",cell[j].p[0]-3,cell[j].p[1]-3,cell[j].p[2]-3);
	fprintf (op,"%i %i %i\n",cell[j].N[0],cell[j].N[1],cell[j].N[2]);
}
fclose(op);
exit(0);
}
