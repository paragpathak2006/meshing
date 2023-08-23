#include <iostream>
#include<math.h>
#include<stdlib.h>
#include<conio.h>
#include<string.h>
using namespace std;
int const Gnmax=4000,Gpmax=4000,ptlimit=15;const float pi=3.14159;
float kf=1,rho=1.0,dt,dto,volmin=1.0e5,t=0,Norm=1;
int pmax=3,nmax,Nmax,nn,bmax,change=1,pass=0,R=1,test=0,lmax=0,Count=0;
class vec
{
	public:
	float i,j,Tc,Tf,Tx,Ty,Txx,Tyy,Txy,S;
	bool coarse,deleted,BN;
	int p1,p2,t[ptlimit+1],level;
	
	vec(){deleted=false;BN=false;level=0;coarse=false;p1=0;p2=0;for(int i=0;i<=9;i++){t[i]=0;}}
	
	void pointmarking();
	void nodetemp();
	void nodederivatives();
	void nodederivatives2();
	void add(int i);
	void del(int i);
	void del();

	vec operator+(vec b);
	vec operator-(vec b);
	float operator*(vec b);
	float operator%(vec b);
	vec operator*(float b);
	vec operator/(float b);
	friend vec operator*(float b,vec a);
};

class triangle
{
	public:
	int N[3],p[3],midpt[3],type,split,j,level;
	bool midptallocation,mod,coarse,deleted,served,dM;

	vec V,C,A[3],e[3],n[3],f[3];
	float S,Vol,Tf,Tc,To,E1,E2,Ec,Ef,E3c,E3f,Sdc,Sdf;
	float Ixx,Iyy,Ixy,Sxx,Syy,Sxy,Tx,Ty,Txx,Tyy,Txy,Txxx,Txxy,Txyy,Tyyy,Ttc,Ttf;
	float de[3],dn[3],dF[3],dc[3],dr[3],ds[3],Flux,Ai[9][9],ai[5],dtt,dts,dtd;
	
	triangle(){type=0;deleted=false;split=0;coarse=false;mod=false;
	Tx=0;Txx=0;Txxx=0; Ty=0;Tyy=0;Tyyy=0; Txxy=0;Txyy=0;Txy=0;};
	void initialize(int p0,int p1,int p2,int j);

//---------triangle functions-------------
	void splitmarking();//iterative
	void colapsebility();
	void coarsemarking();//iterative
	void midpts();
	void cellmodification();
	void setneighbour();
	void parameters();
	void del();

//---------triangle modifications--------
	void halfsplit();
	void fullsplit();
	void boundrysplit();

	void halfcolapse(int i);
	void fullcolapse();
	void boundrycolapse();
	void modified();
//---------diffusion functions------------
	void errors();
	void trucationerrors(int Ar);
	void errormark();
	
	void heatflux(int A);
	void Interiorderivatives(int A);
	void Boundryderivatives(int A);
	void smothner();
//	void solution(float A[][4]);
	void inputcell(FILE *fp,int jj);
	void timestep();
	void gridtaging();

//-----------output functions--------------
	void plot3D(FILE *fp);
	void plot3DN(FILE *fp);
	void map3D(FILE *fp);
	
	void draw2D(FILE *fp);
	void err(FILE *fp,int Ar);
};
triangle cell[Gnmax+1];
vec P[Gpmax+1];

int newpt();
float gauss(vec P,float t,float &S);
float gauss(vec P,float t);

float map(vec P,float T,int j)
{
 	float x=P.i,y=P.j,Trec=0;
	Trec=cell[j].Txxx*x*x*x/6+cell[j].Tyyy*y*y*y/6+Trec;
	Trec=cell[j].Txxy*x*x*y/2+cell[j].Txyy*x*y*y/2+Trec;

	Trec=cell[j].Txx*x*x/2+cell[j].Tyy*y*y/2+cell[j].Txy*x*y+Trec;
	Trec=T+cell[j].Tx*x+cell[j].Ty*y+Trec;
 	return Trec;
}

void triangle::modified(){mod=1;cell[N[0]].mod=1;cell[N[1]].mod=1;cell[N[2]].mod=1;}
float d(vec a){float f;f=sqrt(a.i*a.i+a.j*a.j);return f;}
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::initialize(int p0,int p1,int p2,int jj)
{
 	int swap;
	p[0]=p0;	p[1]=p1;	p[2]=p2;
	N[0]=0;		N[1]=0;		N[2]=0;

	C=(P[p[0]]+P[p[1]]+P[p[2]])/3;
	mod=true;	coarse=false;	served=false;	deleted=false;	j=jj;	
	Vol=(P[p1]-P[p0])%(P[p2]-P[p0]);Vol=Vol/2;
	Tc=C.Tc;
	if(Vol<0) {cout<<"-ve volume cell"<<j<<endl;getch();swap=p[1];p[1]=p[2]; p[2]=swap;Vol=-Vol;}  
	P[p0].add(j);	P[p1].add(j);	P[p2].add(j);
}
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::splitmarking()//iterative
{
	int nk,kk,r,s,pa,c1,c2;cell[0].split=0;
	bool Ar,As,Br,Bs;
	if(split>0&&type==0&&0)//full cell addtitional gaurd
	{ 
		kk=split-1;
		if(kk==0) {r=1;s=2;}	if(kk==1) {r=2;s=0;}	if(kk==2) {r=0;s=1;}	
	
		Ar=cell[N[r]].split>0&&cell[N[r]].type==0;
		As=cell[N[s]].split>0&&cell[N[s]].type==0;
		Br=cell[N[r]].type>0&&cell[N[r]].N[0]!=j;
		Bs=cell[N[s]].type>0&&cell[N[s]].N[0]!=j;

		if(Ar||As||Br||Bs) split=-1;
	}
	if(!served&&split<0&&type==0)//full cell
	{ 
		change=1;served=true;
		for(int k=0;k<3;k++)
		{
			if(cell[N[k]].split>0) cell[N[k]].split=-1;
			if(cell[N[k]].split==0){for(kk=0;kk<3;kk++) if(cell[N[k]].N[kk]==j) cell[N[k]].split=kk+1;}
			if(cell[N[k]].type>0&&cell[N[k]].N[0]==j) cell[N[k]].dM=true;
		}	
	}
	if(!served&&split!=0&&type>0)//half cell
	{ 
		change=1;served=true;
		for(int k=0;k<3;k++)
		{
			if(k!=0)
			{

			if(cell[N[k]].split>0) cell[N[k]].split=-1;
			if(cell[N[k]].split==0){for(kk=0;kk<3;kk++) if(cell[N[k]].N[kk]==j) cell[N[k]].split=kk+1;}
			if(cell[N[k]].type>0&&cell[N[k]].N[0]==j&&k!=0) cell[N[k]].dM=true;
			}
		}	
	}cell[0].split=0;
}
//-----------**************--------------triangle functions--------------*************-----------------
void vec::pointmarking()
{
	bool A,B;
	coarse=true;
	if(p1==0||p2==0) coarse=false;
	if(coarse)
	{
	 	for(int j=1;j<=t[0];j++)
		{
			A=cell[t[j]].coarse==false&&cell[t[j]].type==0;
			B=cell[t[j]].split!=0;
			if(A||B){coarse=false;break;}
		}
	}
}
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::colapsebility()
{
	bool A,B,C,D;
	int p0A,p0B,p1A,p1B,p2A,p2B,E;

		p0A=P[p[0]].p1;		p0B=P[p[0]].p2;
		p1A=P[p[1]].p1;		p1B=P[p[1]].p2;
		p2A=P[p[2]].p1;		p2B=P[p[2]].p2;

		D = P[p[0]].coarse||P[p[1]].coarse||P[p[2]].coarse;
		if(!D) coarse=false;

	if(type==0&&D)
	{
		A = p0A==p1A || p0A==p1B || p0B==p1A || p0B==p1B;
		B = p1A==p2A || p1A==p2B || p1B==p2A || p1B==p2B;
		C = p2A==p0A || p2A==p0B || p2B==p0A || p2B==p0B;
		E = p0A*p1A*p2A;
		
		if(A&&B&&C&&E!=0) coarse=true;else coarse=false;

		if(!coarse)
		{
			if(!A&&!B) P[p[1]].coarse=false;
			if(!B&&!C) P[p[2]].coarse=false;
			if(!C&&!A) P[p[0]].coarse=false;
		}
	}
	if(type==1&&D)	
	{
		if(P[p[2]].coarse) coarse=true;
		else coarse=false;
		P[p[0]].coarse=false;
		P[p[1]].coarse=false;
	}
	if(type==2&&D) 
	{
		coarse=false;
		P[p[0]].coarse=false;
		P[p[2]].coarse=false;
	}
}
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::coarsemarking()//iterative
{
	bool A,B,C;
	int jj,E;
	if(coarse&&type==0)
	{
		E=P[p[0]].coarse+P[p[1]].coarse+P[p[2]].coarse;
		if(E<2){P[p[0]].coarse=false;P[p[1]].coarse=false;P[p[2]].coarse=false;coarse=false;change=1;}

		if(E==2&&0)//additional gaurd
		{
			for(int k=0;k<3;k++)
			{
		 	if(P[p[k]].coarse) 
			{
			for(int i=1;i<=P[p[k]].t[0];i++)
			{
			 	jj=P[p[k]].t[i];
				E=P[cell[jj].p[0]].coarse+P[cell[jj].p[1]].coarse+P[cell[jj].p[2]].coarse;
				if(cell[jj].coarse&&cell[jj].type==0&&jj!=j&&E<3)
				{ 
				P[p[0]].coarse=false;
				P[p[1]].coarse=false;
				P[p[2]].coarse=false;
				coarse=false;
				change=1;break;
				}
			}
			}
			}
		}
	}
	if(coarse&&type==1&&!P[p[2]].coarse){coarse=false;change=1;}
}
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::midpts()
{
	int r,s,nk,kk,pa,c1,c2,n0;
	bool A,B,CC;
	float deno;
	if(split!=0)
	{
	midptallocation=true;
	for(int k=0;k<3;k++)
	{
		if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}	

		A=(split==k+1||split==-1)&&type==0;
		B=k!=type&&k!=0&&type>0;
		CC=dM&&k!=type;

		if(cell[N[k]].midptallocation&&(A||B||CC))//grab
		{for(kk=0;kk<3;kk++) {if(cell[N[k]].N[kk]==j) midpt[k]=cell[N[k]].midpt[kk];}}

		if((A||B||CC)&&!cell[N[k]].midptallocation)//produce
		{
			pa=newpt();if(pmax>Gpmax) {printf("points exeed limit of %i",Gpmax);getch();exit(0);}
			P[pa]=(P[p[r]]+P[p[s]])/2;

			deno=1/dF[k]+1/dc[k];
//			if(N[k]!=0)P[pa].Tc=(4*P[pa].Tc/dn[k]+cell[N[k]].Tc/dF[k]+Tc/dc[k])/deno;
//			if(N[k]!=0) P[pa].Tc==(map(P[pa]-C,Tc,j)/dF[k]+map(P[pa]-cell[N[k]].C,cell[N[k]].Tc,N[k])/dF[k])/deno;
			if(N[k]!=0) P[pa].Tc==map(P[pa]-C,Tc,j)/2+map(P[pa]-cell[N[k]].C,cell[N[k]].Tc,N[k])/2;		
			P[pa].p1=p[r];P[pa].p2=p[s];
			
			midpt[k]=pa;
			if(P[p[r]].BN&&P[p[s]].BN&&N[k]==0) P[pa].BN=true;	else P[pa].BN=false;

			if(P[p[r]].level>P[p[s]].level)P[pa].level=P[p[r]].level+1;
			else P[pa].level=P[p[s]].level+1;

		}
	}	
	}
}	
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::cellmodification()
{
	bool A,B,C;

	if(split==-1&&type==0) fullsplit();
	if(split>0&&type==0) halfsplit();
	if(split!=0&&type==1) boundrysplit(); 
	if(coarse&&type==0)
	{
		A=P[p[0]].coarse&&P[p[1]].coarse;
		B=P[p[1]].coarse&&P[p[2]].coarse;
		C=P[p[2]].coarse&&P[p[0]].coarse;
		
	 	if(A&&!B) halfcolapse(2);
	 	if(B&&!C) halfcolapse(0);
	 	if(C&&!A) halfcolapse(1);

	 	if(A&&B&&C)fullcolapse();
	}
	if(coarse&&type==1) boundrycolapse();
}

//-----------**************--------------triangle functions--------------*************-----------------
//-----------@@@@@@@@@@@@@@-*************@@@@@@@@@@@@@@@@@@@@***********@@@@@@@@@@@@@@-----------------
//-----------@@@@@@@@@@@@@@-************triangle modifications***********@@@@@@@@@@@@@-----------------
void triangle::halfsplit()
{
 	modified();vec Co=C;
	R=1;float Trec=Tc;
	int k,r,s,pk,pr,ps,mk,Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	if(!dM) k=split-1;
    if(dM) k=0;
	if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}	
	
	pk=p[k];	pr=p[r];	ps=p[s];	mk=midpt[k];
	P[ps].del(j);

	nmax++;
	if(nmax>Gnmax) {printf("triangles exeed limit");getch();exit(0);}
	cell[nmax].initialize(pk,mk,ps,nmax);
	cell[nmax].type=2;
	cell[nmax].split=0;

//	if(Bn==0) cell[nmax].Tc=map(cell[nmax].C-C,Tc,j);	
//	if(Bn!=0)	cell[nmax].Tc=Trec/2+cell[nmax].Tc/2;

	initialize(pk,pr,mk,j);
//	if(Bn==0) Tc=map(C-Co,Trec,j);	
	//Tc=Trec/2+Tc/2;
	type=1;split=0;
}
//-------------@@@@@@@@@@@@-------------triangle modifications-----------@@@@@@@@@@--------------------
void triangle::fullsplit()
{
 	R=1;
 	modified();
	int Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	int k,r,s,m0=midpt[0],m1=midpt[1],m2=midpt[2];float Trec=Tc;
	P[p[0]].del(j);		P[p[1]].del(j);		P[p[2]].del(j);
	for(k=0;k<3;k++)
	{
		if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}	

		nmax++;   
		if(nmax>Gnmax) {printf("triangles exeed limit");getch();exit(0);}

		cell[nmax].initialize(p[k],midpt[s],midpt[r],nmax);	
		cell[nmax].split=0;	
//		if(Bn==0)cell[nmax].Tc=map(cell[nmax].C-C,Tc,j);
	}
	initialize(m0,m1,m2,j);Tc=Trec;
	type=0;split=0;
}
//-------------@@@@@@@@@@@@-------------triangle modifications-----------@@@@@@@@@@--------------------
void triangle::boundrysplit()
{	
	R=1;
 	modified();cell[N[1]].modified();
	int Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	int c2=N[1],p0=p[0],ms=midpt[2],mr=cell[c2].midpt[1],mk=p[2];	
	float Tr=Tc,Ts=cell[c2].Tc;
	vec Cr=C,Cs=cell[c2].C;

	P[p[0]].del(j);			P[p[0]].del(c2);
	initialize(ms,p[1],p[2],j);type=0;
//	Tc=map(C-Cr,Tr,j);	
//	Tc=Tr/2+Tc/2;
	
	cell[c2].initialize(mr,cell[c2].p[1],cell[c2].p[2],c2);
//	if(Bn==0) cell[c2].Tc=map(cell[c2].C-Cs,Ts,c2);	
//	cell[c2].Tc=Ts/2+cell[c2].Tc/2;
	cell[c2].type=0;
	//----------------**********************----------------------------------------	
	nmax++;if(nmax>Gnmax) {printf("triangles exeed limit");getch();exit(0);}
	cell[nmax].initialize(mk,mr,ms,nmax);
	cell[nmax].type=0;
	cell[nmax].split=0;
//	if(Bn==0) cell[nmax].Tc=map(cell[nmax].C-Cr,Tr,j)/2+map(cell[nmax].C-Cs,Ts,c2)/2;
//	cell[nmax].Tc=Tr/2+Ts/2;

	nmax++;if(nmax>Gnmax) {printf("triangles exeed limit");getch();exit(0);}
	cell[nmax].initialize(p0,ms,mr,nmax);
	cell[nmax].type=0;
	cell[nmax].split=0;
//	if(Bn==0) cell[nmax].Tc=map(cell[nmax].C-Cr,Tr,j)/2+map(cell[nmax].C-Cs,Ts,c2)/2;
//	cell[nmax].Tc=cell[nmax-1].Tc/2+cell[nmax].Tc/2;

	//----------------**********************----------------------------------------	
	if(dM) 
	{
		halfsplit();
//		if(Bn==0) Tc=map(C-Cr,Tr,j);
//		if(Bn==0) cell[nmax].Tc=map(cell[nmax].C-Cr,Tr,j);
	}
	if(cell[c2].dM) 
	{
		cell[c2].halfsplit();
//		if(Bn==0) cell[c2].Tc=map(cell[c2].C-Cs,Ts,c2);
//		if(Bn==0) cell[nmax].Tc=map(cell[nmax].C-Cs,Ts,c2);
	} 
	split=0;cell[c2].split=0;
}
//-------------@@@@@@@@@@@@-------------triangle modifications-----------@@@@@@@@@@--------------------
void triangle::halfcolapse(int k)
{
	R=1;modified();cell[N[0]].modified();cell[N[1]].modified();cell[N[2]].modified();
	int r,s,pt[3],m,n;
	int Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	if(k==0){r=1;s=2;}	if(k==1){r=2;s=0;}	if(k==2){r=0;s=1;}	


	float Tr=cell[N[r]].Tc,Ts=cell[N[s]].Tc;
	for(m=0;m<3;m++)
	{for(n=0;n<3;n++){if(cell[N[m]].N[n]==j) pt[m]=cell[N[m]].p[n];}}

	cell[N[r]].initialize(pt[k],pt[r],p[k],N[r]);cell[N[r]].type=1;	
	cell[N[s]].initialize(pt[k],p[k],pt[s],N[s]);cell[N[s]].type=2;	

//	if(Bn==0)cell[N[r]].Tc=map(cell[N[r]].C-C,Tc,j);
//	if(Bn==0)cell[N[s]].Tc=map(cell[N[s]].C-C,Tc,j);

//	cell[N[r]].Tc=(2*Tr+Tc+cell[N[k]].Tc)/4;
//	cell[N[s]].Tc=(2*Ts+Tc+cell[N[k]].Tc)/4;

	P[pt[k]].del(N[k]);
	P[p[k]].del(j);
	cell[N[k]].del();
	del();
}
//-------------@@@@@@@@@@@@-------------triangle modifications-----------@@@@@@@@@@--------------------
void triangle::fullcolapse()
{
 	R=1;modified();cell[N[0]].modified();cell[N[1]].modified();cell[N[2]].modified();
	int Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	float Trec=Tc;
	int pt[3],k,kk;
	for(k=0;k<3;k++) 
	{
		for(kk=0;kk<3;kk++) if(cell[N[k]].N[kk]==j) pt[k]=cell[N[k]].p[kk];
		P[pt[k]].del(N[k]);
		cell[N[k]].del();
	}
	initialize(pt[0],pt[1],pt[2],j);Tc=Trec;
}
//-------------@@@@@@@@@@@@-------------triangle modifications-----------@@@@@@@@@@--------------------
void triangle::boundrycolapse()
{
	R=1;modified();cell[N[1]].modified();
	int n1=N[1],p2=cell[N[1]].p[2];
	int Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	float Tr=Tc,Ts=cell[n1].Tc;
	vec Cr=C,Cs=cell[n1].C;
	
	P[p2].del(n1);
	P[p[0]].del(n1);
	initialize(p[0],p[1],p2,j);type=0;
//	if(Bn==0)Tc=map(C-Cr,Tr,j)/2+map(C-Cs,Ts,n1)/2;
//	else 
	Tc=Tr/2+Ts/2;
	cell[n1].del();
}
//-------------@@@@@@@@@@@@-------------triangle modifications-----------@@@@@@@@@@@@@-----------------
//-------------@@@@@@@@@@@@@*************@@@@@@@@@@@@@@@@@@@@***********@@@@@@@@@@@@@@-----------------
//-----------**************--------------triangle functions--------------*************-----------------


//-----------**************--------------triangle functions--------------*************-----------------
void triangle::setneighbour()
{
 	int k,kk,r,s,i,l;
 	bool found;
	if(!deleted)
	{
		for(k=0;k<3;k++)
		{
			if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}	
	 
			for(i=1;i<=P[p[r]].t[0];i++)
			{for(l=1;l<=P[p[s]].t[0];l++)
			{found=P[p[r]].t[i]==P[p[s]].t[l]&&P[p[s]].t[l]!=j;
			if(found) break;}
			if(found) break;}
			if(found) N[k]=P[p[r]].t[i];
			if(!found&&P[p[r]].BN&&P[p[s]].BN) N[k]=0;
			
			if(!found&&(!P[p[r]].BN||!P[p[s]].BN)) 
			{printf("triangle %i neibour error at N[%i]=%i %i po=%i",j,k,N[k],deleted,P[p[0]].deleted);
			printf(" p1=%i p2=%i.\n",P[p[1]].deleted,P[p[2]].deleted);getch();};
		}
	}
}	
//-----------**************--------------triangle functions--------------*************-----------------
void triangle::parameters()
{
	/*int r,s,k,nk,l;
	bool Ao,Bo;
	float X,Y,x0,y0,x1,y1,x2,y2,xi,yi,I0,I2;
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
*/
	float x,y;int nk,i,ii=0,l,k,Bn,r,s,ni=9;
	double A1[9][18],pp=0,x1[9],S1[18];

//	for(int k=0;k<3;k++) if(P[p[k]].t[0]==4) ii=p[k];
	for(int i=0;i<9;i++)for(int ii=0;ii<9;ii++){if(i==ii) A1[i][9+ii]=1;else A1[i][9+ii]=0;}

	Bn = P[p[0]].BN + P[p[1]].BN + P[p[2]].BN;
//if(Bn==2){printf("j=%i\n",j);}

	if(Bn==0)
	{
	
		for(int k=0;k<3;k++)
		{
		for(int kk=0;kk<3;kk++)
		{
			if(j==cell[N[k]].N[kk]) nk=N[k];	else nk=cell[N[k]].N[kk];
	
			if(nk!=N[k]&&ii>0) for(int h=0;h<3;h++)
			if(cell[nk].p[h]==ii){nk=cell[nk].N[h];ii=0;}
		
			x=cell[nk].C.i-C.i;	y=cell[nk].C.j-C.j;
	
			i=3*k+kk;
	
			A1[i][0]=x;			A1[i][1]=y;
			A1[i][2]=x*x/2;		A1[i][3]=y*y/2;		A1[i][4]=x*y;
			A1[i][5]=x*x*x/6;	A1[i][6]=y*y*y/6;	A1[i][7]=x*x*y/2;	A1[i][8]=x*y*y/2;		
		}
		}
	
		for(int jj=0;jj<2*ni;jj++) {S1[jj]=0;for(int i=0;i<ni;i++) S1[jj]=S1[jj]+A1[i][jj];}
		for(int i=0;i<ni;i++)for(int jj=0;jj<2*ni;jj++) A1[i][jj]=S1[jj]-A1[i][jj];
	
		for(int i=0;i<ni;i++)//row change down
		{
			pp=A1[i][i]; 
	//		if(fabs(pp)<1E-15)	break;printf("pivot %e too small in row %i\n",pp,i);getch();exit(0);}
			for(int jj=i;jj<2*ni;jj++) A1[i][jj]=A1[i][jj]/pp;//aii=1
			for(int ii=0;ii<ni;ii++)//row change down
			{	
				if(ii==i)ii++; if(ii==ni)break;	
				pp=A1[ii][i];
				for(int jj=i;jj<2*ni;jj++) A1[ii][jj]=A1[ii][jj]-pp*A1[i][jj]; 
			}
		}
	for(int k=0;k<9;k++)for(int kk=9;kk<18;kk++) Ai[k][kk-9]=A1[k][kk];

if(j==180||1)
{
//	for(int k=0;k<9;k++){for(int kk=0;kk<9;kk++) printf("%1.0f ",A1[k][kk]);printf("\n");}printf("\n");
//	for(int k=0;k<9;k++){for(int kk=9;kk<18;kk++) printf("%1.2e ",A1[k][kk]);printf("j=%i ",j);}
//printf("j=%i\t",j);
//	for(int k=9;k<18;k++){printf("D%i=%1.2f\t",k-8,A1[2][k]+A1[3][k]);}printf("\n");getch();
};

	}

	vec a,c,b,B,Xo,Yo;float tt;
	for(int k=0;k<=2;k++)
	{
	    if(k==0){r=1;s=2;}	if(k==1) {r=2;s=0;}	if(k==2){r=0;s=1;}
		nk=N[k];if(nk==0) cell[nk].C = (P[p[s]]+P[p[r]])/2;
		
		A[k].i = P[p[s]].j - P[p[r]].j;
		A[k].j = P[p[r]].i - P[p[s]].i;
	
		de[k]= d(cell[nk].C-C);
		dn[k]= d(P[p[s]]-P[p[r]]);

		e[k] = (cell[nk].C - C)/de[k];
		n[k] = (P[p[s]]-P[p[r]])/dn[k];
	
		tt=((P[p[s]]-C)%n[k])/(e[k]%n[k]);
		f[k]=C+tt*e[k];

		dc[k]=d(C-f[k]);		dF[k]=de[k]-dc[k];
		dr[k]=d(P[p[r]]-f[k]);	ds[k]=dn[k]-dr[k];

	}
}

//-----------**************--------------triangle functions--------------*************-----------------
//-----------**************--------------triangle functions--------------*************-----------------
//-----------**************--------------triangle functions--------------*************-----------------
//-----------**************--------------diffusion functions-------------*************-----------------
//-----------**************--------------diffusion functions-------------*************-----------------
void triangle::heatflux(int Ar)
{		
	int r,s,nk,kk,Bn;
	float Te,Tn,ap,w;Flux=0;
	Bn = P[p[0]].BN + P[p[1]].BN + P[p[2]].BN;
	for(kk=0;kk<3;kk++) {if(P[p[kk]].BN&&Bn==1) break;if(!P[p[kk]].BN&&Bn==2) break;if(!N[kk]!=0&&Bn==3) break;}

	for(int k=0;k<=2;k++)//Jf.Af calculations.
	{
		if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}
		w=e[k]*n[k];		nk=N[k];

		if(Ar==1)
		{
			if(nk==0) cell[nk].Tc=(P[p[s]].Tc+P[p[r]].Tc)/2;
			Te=(cell[nk].Tc-Tc)/de[k];
			Tn=(P[p[s]].Tc-P[p[r]].Tc)/dn[k];Tn=0;
		}

		if(Ar==2)
		{
			if(nk==0) cell[nk].Tf=(P[p[s]].Tf+P[p[r]].Tf)/2;
			Te=(cell[nk].Tf-Tf)/de[k];
			Tn=(P[p[s]].Tf-P[p[r]].Tf)/dn[k];Tn=0;
		}
		Flux=Flux-(Te-w*Tn)*kf*e[k]*A[k]/(1-w*w);
	}
//	ap=Vol/dt;
//Tf=gauss(C,dt,S);

	if(Ar==1){Tf=Tc+(Txx+Tyy)*dt+S*dt;if(Bn==3)Tf=Tc-Flux*dt/Vol+S*dt;Ttc=Txx+Tyy+S;}
	if(Ar==2){Ttf=Txx+Tyy+S;E1=(Ttf-Ttc)*.5/dt;
	}
}

//-----------**************--------------derivatives functions-------------*************-----------------
void triangle::Interiorderivatives(int Ar)
{
	float x,y,ff=0;int nk,i,ii=0,k,Bn,r,s,ni=9;
	double A1[9][10],pp=0,x1[9],S1[10];bool A=0;

	for(int k=0;k<3;k++) if(P[p[k]].t[0]==4) ii=p[k];

	Bn = P[p[0]].BN + P[p[1]].BN + P[p[2]].BN;
//if(Bn==2){printf("j=%i\n",j);}

	if(Bn==0)
	{for(int r=0;r<5;r++)ai[r]=0;
		ff=0;
		for(int k=0;k<3;k++)
		{
		for(int kk=0;kk<3;kk++)
		{
			if(j==cell[N[k]].N[kk]) nk=N[k];	else nk=cell[N[k]].N[kk];
	
			if(nk!=N[k]&&ii>0) for(int h=0;h<3;h++)
			if(cell[nk].p[h]==ii){nk=cell[nk].N[h];ii=0;}
		
			x=cell[nk].C.i-C.i;	y=cell[nk].C.j-C.j;
	
			i=3*k+kk;
			if(Ar==1) A1[i][9]=cell[nk].Tc-Tc;
			if(Ar==2) A1[i][9]=cell[nk].Tf-Tf;
	
			A1[i][0]=x;			A1[i][1]=y;
			A1[i][2]=x*x/2;		A1[i][3]=y*y/2;		A1[i][4]=x*y;
			A1[i][5]=x*x*x/6;	A1[i][6]=y*y*y/6;	A1[i][7]=x*x*y/2;	A1[i][8]=x*y*y/2;		
		for(int r=0;r<5;r++)ai[r]=ai[r]+(Ai[i][2]+Ai[i][3])*pow(x,r)*pow(y,4-r);
//		ff=ff+Ai[i][2]*A1[i][9];
		}
		}
//printf("a0=%e\n",ai[0]);getch();
		if(j==26&&pass==3&&nmax==687&&A){printf("\nj=%i initial matrix\n",j);for(int k=0;k<9;k++){for(int kk=0;kk<10;kk++)printf("%e ",A1[k][kk]); printf("\n");}printf("\n");getch();}
	
//		for(int jj=0;jj<=ni;jj++) {S1[jj]=0;for(int i=0;i<ni;i++) S1[jj]=S1[jj]+A1[i][jj];}
//		for(int i=0;i<ni;i++)for(int jj=0;jj<=ni;jj++) A1[i][jj]=S1[jj]-A1[i][jj];

	
		for(int i=0;i<ni;i++)//row change down
		{
//	if(j==26&&pass==3&&nmax==687&&A){printf("\nj=%i before row change\n",j);for(int k=0;k<9;k++){for(int kk=0;kk<10;kk++)printf("%e ",A1[k][kk]); printf("\n");}printf("\n");getch();};
			pp=A1[i][i]; 

		if(i<ni-1&&fabs(pp)<1e-13)
		{
		 	pp=0;nk=0;
			for(ii=i+1;ii<ni-1;ii++)if(fabs(A1[ii][i])>fabs(pp)){pp=A1[ii][i];nk=ii;}
 			for(int jj=i;jj<=ni;jj++) {pp=A1[i][jj];A1[i][jj]=A1[nk][jj];A1[nk][jj]=pp;}
		}	
		if(j==26&&pass==3&&nmax==687&&A){printf("\nj=%i after row change %i-%i\n",j,i,nk);for(int k=0;k<9;k++){for(int kk=0;kk<10;kk++)printf("%e ",A1[k][kk]); printf("\n");}printf("\n");getch();};

			pp=A1[i][i]; 
	//		if(fabs(pp)<1E-15)	break;printf("pivot %e too small in row %i\n",pp,i);getch();exit(0);}
			for(int jj=i;jj<=ni;jj++) A1[i][jj]=A1[i][jj]/pp;//aii=1
			for(int ii=0;ii<ni;ii++)//row change down
			{	
				if(ii==i)ii++; if(ii==ni)break;	
				pp=A1[ii][i];
				for(int jj=i;jj<=ni;jj++) A1[ii][jj]=A1[ii][jj]-pp*A1[i][jj]; 
			}
			if(j==401&&nmax==687&&A){printf("\nj=%i after Normalisation\n",j);for(int k=0;k<9;k++){for(int kk=0;kk<10;kk++)printf("%e ",A1[k][kk]); printf("\n");}printf("\n");getch();};
		}
		if(j==401&&nmax==687&&A){printf("\nj=%i final matrix\n",j);for(int k=0;k<9;k++){for(int kk=0;kk<10;kk++)printf("%e ",A1[k][kk]); printf("\n");}printf("\n");getch();};

		/*
	
		x1[ni-1]=A1[ni-1][ni]/A1[ni-1][ni-1];//backward sweep on [U]
		for(int i=ni-2;i>=0;i--)
		{
			pp=0;
			for(int jj=i+1;jj<ni;jj++) pp=pp+x1[jj]*A1[i][jj];
			x1[i] = (A1[i][ni] - pp)/A1[i][i];
		}
		
		Tx=x1[7];	Ty=x1[8];	Txx=x1[4];	Tyy=x1[5];	Txy=x1[6];
		if(nk==0){	Tx=0;	Ty=0;	Txx=0;	Tyy=0;	Txy=0;}*/

		Tx=A1[0][9];	Ty=A1[1][9];	
		Txx=A1[2][9];	Tyy=A1[3][9];	Txy=A1[4][9];
		Txxx=A1[5][9];	Tyyy=A1[6][9];	Txxy=A1[7][9];	Txyy=A1[8][9];
	}

	ii=P[p[0]].t[0]==4+P[p[1]].t[0]==4+P[p[2]].t[0]==4;
	if(ii>1){printf("mesh has four square problem");getch();exit(0);}
	
	if(Bn==3||ii>1)
	{
		Tx=0;	Ty=0;	
		Txx=0;	Tyy=0;	Txy=0;
		Txxx=0;	Tyyy=0;	Txxy=0;	Txyy=0;
	}
}

//-----------**************--------------derivatives functions-------------*************-----------------
void triangle::Boundryderivatives(int Ar)
{
	float x,y;int nk,ii=0,i,k,Bn,r,s,ni=5;
	double A1[5][6],pp=0,x1[5],S1[5];

	Bn = P[p[0]].BN + P[p[1]].BN + P[p[2]].BN;
//	if(Bn==1) printf("j=%i\n",j);
	for(k=0;k<3;k++) {if(P[p[k]].BN&&Bn==1) break;if(!P[p[k]].BN&&Bn==2) break;if(!N[k]!=0&&Bn==3) break;}
	if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}

	if(Bn==1)
	{
		for(int i=0;i<5;i++)
		{

			if(i==3)nk=N[r];if(i==4) nk=N[s];

		 	if(i<3){if(j!=cell[N[k]].N[i]) nk=cell[N[k]].N[i];else nk=p[k];}

		 	if(i<3&&j==cell[N[k]].N[i])
			{
			  	x=P[nk].i-C.i;	y=P[nk].j-C.j;	
			  	if(Ar==1)A1[i][5]=P[nk].Tc-Tc;
				if(Ar==2)A1[i][5]=P[nk].Tf-Tf;
			}

			else
			{
				x=cell[nk].C.i-C.i;	y=cell[nk].C.j-C.j;	
				if(Ar==1)A1[i][5]=cell[nk].Tc-Tc;
				if(Ar==2)A1[i][5]=cell[nk].Tf-Tf;
			}

			A1[i][0]=x;			A1[i][1]=y;
			A1[i][2]=x*x/2;		A1[i][3]=y*y/2;		A1[i][4]=x*y;
//			if(j==163){for(ii=0;ii<=5;ii++) printf("%e ",A1[i][ii]);printf("\n");}
		}
	}

	if(Bn==2)
	{
		for(int i=0;i<5;i++)
		{
			if(i==2)nk=N[r];if(i==3) nk=N[s];
		 	if(i==0) nk=p[r];if(i==1) nk=p[s];
			if(i==4) 
			{
			 	for(int kk=0;kk<3;kk++)
				{
			 	if(cell[N[r]].p[kk]==p[s]&&cell[N[r]].N[kk]!=0){nk=cell[N[r]].N[kk];break;}
			 	if(cell[N[s]].p[kk]==p[r]&&cell[N[s]].N[kk]!=0){nk=cell[N[s]].N[kk];break;}
				}
			}
		 	if(i>1)
			{
				x=cell[nk].C.i-C.i;	y=cell[nk].C.j-C.j;	
				if(Ar==1)A1[i][5]=cell[nk].Tc-Tc;
				if(Ar==2)A1[i][5]=cell[nk].Tf-Tf;
			}
		 	else
			{
			  	x=P[nk].i-C.i;	y=P[nk].j-C.j;	
			  	if(Ar==1)A1[i][5]=P[nk].Tc-Tc;
				if(Ar==2)A1[i][5]=P[nk].Tf-Tf;
			}

			A1[i][0]=x;			A1[i][1]=y;
			A1[i][2]=x*x/2;		A1[i][3]=y*y/2;		A1[i][4]=x*y;
//			if(j==163){for(ii=0;ii<=5;ii++) printf("%e ",A1[i][ii]);printf("\n");}
		}
	}
/*	if(Bn==3)
	{
		for(int i=0;i<5;i++)
		{
		if(i<3){for(int kk=0;kk<3;kk++){if(j==cell[N[k]].N[i]) nk=N[k];else nk=cell[N[k]].N[kk];}}		if(i==0) nk=p[r];if(i==1) nk=p[s];
		 	if(i==3) nk=p[r];if(i==4) nk=p[s];
		 	if(i<3)
			{
				x=cell[nk].C.i-P[p[k]].i;	y=cell[nk].C.j-P[p[k]].j;	
				if(Ar==1)A1[i][5]=cell[nk].Tc-P[p[k]].Tc;
				if(Ar==2)A1[i][5]=cell[nk].Tf-P[p[k]].Tc;
			}
		 	else
			{
			  	x=P[nk].i-P[p[k]].i;	y=P[nk].j-P[p[k]].j;	
			  	if(Ar==1)A1[i][5]=P[nk].Tc-P[p[k]].Tc;
				if(Ar==2)A1[i][5]=P[nk].Tf-P[p[k]].Tc;
			}

			A1[i][0]=x;			A1[i][1]=y;
			A1[i][2]=x*x/2;		A1[i][3]=y*y/2;		A1[i][4]=x*y;
//			if(j==163){for(ii=0;ii<=5;ii++) printf("%e ",A1[i][ii]);printf("\n");}
		}
	}*/

	if(Bn!=0&&Bn!=3)
	{
		for(int jj=0;jj<=ni;jj++) {S1[jj]=0;for(int i=0;i<ni;i++) S1[jj]=S1[jj]+A1[i][jj];}
		for(int i=0;i<ni;i++)for(int jj=0;jj<=ni;jj++) A1[i][jj]=S1[jj]-A1[i][jj];
	
//		if(j==163){printf("\n");for(int i=0;i<ni;i++){for(ii=0;ii<=ni;ii++) printf("%e ",A1[i][ii]);printf("\n");}}

		for(int i=0;i<ni;i++)//row change down
		{
			pp=A1[i][i]; 
			for(int jj=i;jj<=ni;jj++) A1[i][jj]=A1[i][jj]/pp;//aii=1
			for(int ii=0;ii<ni;ii++)//row change down
			{	
				if(ii==i)ii++; if(ii==ni)break;	
				pp=A1[ii][i];
				for(int jj=i;jj<=ni;jj++) A1[ii][jj]=A1[ii][jj]-pp*A1[i][jj]; 
			}
//		if(j==163){printf("\ni=%i\n",i);for(int k=0;k<ni;k++){for(int kk=0;kk<=ni;kk++) printf("%e ",A1[k][kk]);printf("\n");}getch();}
		}
//		if(j==163){printf("\n");for(int k=0;k<ni;k++){for(int kk=0;kk<=ni;kk++) printf("%e ",A1[k][kk]);printf("\n");}}
		
		/*x1[ni-1]=A1[ni-1][ni]/A1[ni-1][ni-1];//backward sweep on [U]
		for(int i=ni-2;i>=0;i--)
		{
			pp=0;
			for(int jj=i+1;jj<ni;jj++) pp=pp+x1[jj]*A1[i][jj];
			x1[i] = (A1[i][ni] - pp)/A1[i][i];
		}
		
		Tx=x1[7];	Ty=x1[8];	Txx=x1[4];	Tyy=x1[5];	Txy=x1[6];
		if(nk==0){	Tx=0;	Ty=0;	Txx=0;	Tyy=0;	Txy=0;}
		*/
		Tx=A1[0][5];	Ty=A1[1][5];	
		Txx=A1[2][5];	Tyy=A1[3][5];	Txy=A1[4][5];
		Txxx=0;	Tyyy=0;	Txxy=0;	Txyy=0;
	}
	if(Bn==3)
	{
		x=cell[N[k]].C.i-C.i;	y=cell[N[k]].C.j-C.j;	
		Tx=cell[N[k]].Tx-x*cell[N[k]].Txx-y*cell[N[k]].Txy;		Ty=cell[N[k]].Ty-x*cell[N[k]].Txy-y*cell[N[k]].Tyy;	
		Txx=cell[N[k]].Txx-x*cell[N[k]].Txxx-y*cell[N[k]].Txxy;	
		Tyy=cell[N[k]].Tyy-x*cell[N[k]].Txyy-y*cell[N[k]].Tyyy;
		Txy=cell[N[k]].Txy-x*cell[N[k]].Txxy-y*cell[N[k]].Txyy;
		Txxx=cell[N[k]].Txxx;	Tyyy=cell[N[k]].Tyyy;	
		Txxy=cell[N[k]].Txxy;	Txyy=cell[N[k]].Txyy;
	}
/*
		if(Bn==1)
		{
		 	x=C.i-cell[N[k]].C.i;y=C.j-cell[N[k]].C.j;
			Txx=cell[N[k]].Txx+cell[N[k]].Txxx*x+cell[N[k]].Txxy*y;
			Tyy=cell[N[k]].Tyy+cell[N[k]].Tyyy*y+cell[N[k]].Txyy*x;
			Txy=cell[N[k]].Txy+cell[N[k]].Txxy*x+cell[N[k]].Txyy*y;


//			A1[2][5];	Tyy=A1[3][5];	Txy=A1[4][5];
			Txxx=cell[N[k]].Txxx;	Tyyy=cell[N[k]].Tyyy;	Txxy=cell[N[k]].Txxy;	Txyy=cell[N[k]].Txyy;
		}
		if(Bn==2)
		{
		 	x=C.i-cell[N[r]].C.i;y=C.j-cell[N[r]].C.j;
			Txx=(cell[N[r]].Txx+cell[N[r]].Txxx*x+cell[N[r]].Txxy*y)/de[r];
			Tyy=(cell[N[r]].Tyy+cell[N[r]].Tyyy*y+cell[N[r]].Txyy*x)/de[r];
			Txy=(cell[N[r]].Txy+cell[N[r]].Txxy*x+cell[N[r]].Txyy*y)/de[r];

		 	x=C.i-cell[N[s]].C.i;y=C.j-cell[N[s]].C.j;
			Txx=Txx+(cell[N[s]].Txx+cell[N[s]].Txxx*x+cell[N[s]].Txxy*y)/de[s];
			Tyy=Tyy+(cell[N[s]].Tyy+cell[N[s]].Tyyy*y+cell[N[s]].Txyy*x)/de[s];
			Txy=Txy+(cell[N[s]].Txy+cell[N[s]].Txxy*x+cell[N[s]].Txyy*y)/de[s];

			Txx=Txx/(1/de[s]+1/de[r]);
			Tyy=Tyy/(1/de[s]+1/de[r]);
			Txy=Txy/(1/de[s]+1/de[r]);

			Txx=(cell[N[r]].Txx/de[s]+cell[N[s]].Txx/de[s])/(1/de[s]+1/de[r]);
			Tyy=(cell[N[r]].Tyy/de[s]+cell[N[s]].Tyy/de[s])/(1/de[s]+1/de[r]);
			Txy=(cell[N[r]].Txy/de[s]+cell[N[s]].Txy/de[s])/(1/de[s]+1/de[r]);


//			A1[2][5];	Tyy=A1[3][5];	Txy=A1[4][5];
			Txxx=cell[N[k]].Txxx;	Tyyy=cell[N[k]].Tyyy;	Txxy=cell[N[k]].Txxy;	Txyy=cell[N[k]].Txyy;
		}*/

}
//-----------**************--------------error functions-------------*************-----------------
/*void triangle::Boundryderivatives(int Ar)
{
	float x,y;int nk,ii=0,i,k,Bn,r,s,ni=9;
	double A1[9][10],pp=0,x1[9],S1[10];

	for(int k=0;k<3;k++) if(P[p[k]].t[0]==4) ii=p[k];

	Bn = P[p[0]].BN + P[p[1]].BN + P[p[2]].BN;
//	if(Bn==2) printf("j=%i\n",j);
	if(Bn==1)
	{
		for(k=0;k<3;k++) if(P[p[k]].BN) break;
		if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}
		for(int i=0;i<9;i++)
		{
		 	if(i<3){if(j==cell[N[k]].N[i]) nk=N[k];	else nk=cell[N[k]].N[i];}
			if(i==3)nk=N[r];if(i==4) nk=N[s];

			x=cell[nk].C.i-C.i;	y=cell[nk].C.j-C.j;	A1[i][9]=cell[nk].Tc-Tc;
		 	if(i>4)
			{	
			 	for(int kk=0;kk<3;kk++) if(cell[N[k]].N[kk]==j&&i==5) nk=cell[N[k]].p[kk];
		 		if(i>5)nk=p[8-i];
				x=P[nk].i-C.i;	y=P[nk].j-C.j;	A1[i][9]=P[nk].Tc-Tc;
			}

			A1[i][0]=x;			A1[i][1]=y;
			A1[i][2]=x*x/2;		A1[i][3]=y*y/2;		A1[i][4]=x*y;
			A1[i][5]=x*x*x/6;	A1[i][6]=y*y*y/6;	A1[i][7]=x*x*y/2;	A1[i][8]=x*y*y/2;		
		}
	}	

	if(Bn==2)
	{
		for(k=0;k<3;k++) if(!P[p[k]].BN) break;
		if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}
//		if(j==163){printf("k=%i p0=%i Bn=%i p1=%i Bn=%i p2=%i Bn=%i\n",k,p[0],P[p[0]].BN,p[1],P[p[1]].BN,p[2],P[p[2]].BN);}
	
		for(int i=0;i<9;i++)
		{

		 	if(i==5)for(int kk=0;kk<3;kk++)	if(cell[N[r]].p[kk]==p[s]) nk=cell[N[r]].N[kk];	
			if(i==6)for(int kk=0;kk<3;kk++)	if(cell[N[s]].p[kk]==p[r]) nk=cell[N[s]].N[kk];
			if(i==7)nk=N[r];if(i==8) nk=N[s];

			x=cell[nk].C.i-C.i;	y=cell[nk].C.j-C.j;	A1[i][9]=cell[nk].Tc-Tc;
	
		 	if(i<5)
			{	
		 		if(i<3)nk=p[i];
			 	if(i==3)for(int kk=0;kk<3;kk++)	if(cell[N[r]].N[kk]==j) nk=cell[N[r]].p[kk];	
			 	if(i==4)for(int kk=0;kk<3;kk++)	if(cell[N[s]].N[kk]==j) nk=cell[N[s]].p[kk];	

				x=P[nk].i-C.i;	y=P[nk].j-C.j;	A1[i][9]=P[nk].Tc-Tc;
			}

			A1[i][7]=x;			A1[i][8]=y;
			A1[i][4]=x*x/2;		A1[i][5]=y*y/2;		A1[i][6]=x*y;
			A1[i][3]=x*x*x/6;	A1[i][2]=y*y*y/6;	A1[i][1]=x*x*y/2;	A1[i][0]=x*y*y/2;		
			if(j==163){for(ii=0;ii<=9;ii++) printf("%e ",A1[i][ii]);printf("\n");}

			Tx=A1[7][9];	Ty=A1[8][9];	
			Txx=A1[4][9];	Tyy=A1[5][9];	Txy=A1[6][9];
			Txxx=A1[3][9];	Tyyy=A1[2][9];	Txxy=A1[1][9];	Txyy=A1[0][9];
		}
	}
	if(Bn==1||Bn==2)
	{
		for(int jj=0;jj<=ni;jj++) {S1[jj]=0;for(int i=0;i<ni;i++) S1[jj]=S1[jj]+A1[i][jj];}
		for(int i=0;i<ni;i++)for(int jj=0;jj<=ni;jj++) A1[i][jj]=S1[jj]-A1[i][jj];
	
		if(j==163){printf("\n");for(int i=0;i<ni;i++){for(ii=0;ii<=ni;ii++) printf("%e ",A1[i][ii]);printf("\n");}}

		for(int i=0;i<ni;i++)//row change down
		{
			pp=A1[i][i]; 
			for(int jj=i;jj<=ni;jj++) A1[i][jj]=A1[i][jj]/pp;//aii=1
			for(int ii=0;ii<ni;ii++)//row change down
			{	
				if(ii==i)ii++; if(ii==ni)break;	
				pp=A1[ii][i];
				for(int jj=i;jj<=ni;jj++) A1[ii][jj]=A1[ii][jj]-pp*A1[i][jj]; 
			}
//		if(j==163){printf("\ni=%i\n",i);for(int k=0;k<ni;k++){for(int kk=0;kk<=ni;kk++) printf("%e ",A1[k][kk]);printf("\n");}getch();}

		}
		if(j==163){printf("\n");for(int k=0;k<ni;k++){for(int kk=0;kk<=ni;kk++) printf("%e ",A1[k][kk]);printf("\n");}}
		
		x1[ni-1]=A1[ni-1][ni]/A1[ni-1][ni-1];//backward sweep on [U]
		for(int i=ni-2;i>=0;i--)
		{
			pp=0;
			for(int jj=i+1;jj<ni;jj++) pp=pp+x1[jj]*A1[i][jj];
			x1[i] = (A1[i][ni] - pp)/A1[i][i];
		}
		
		Tx=x1[7];	Ty=x1[8];	Txx=x1[4];	Tyy=x1[5];	Txy=x1[6];
		if(nk==0){	Tx=0;	Ty=0;	Txx=0;	Tyy=0;	Txy=0;}

	}

	if(Bn==1)
	{
		Tx=A1[0][9];	Ty=A1[1][9];	
		Txx=A1[2][9];	Tyy=A1[3][9];	Txy=A1[4][9];
		Txxx=A1[5][9];	Tyyy=A1[6][9];	Txxy=A1[7][9];	Txyy=A1[8][9];
	}	

}*/
void solution(float A[][4],float Sol[][3],int k)
{
	float D0=0,D1=0,D2=0,D=0;
	
	D=D + A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
	D=D + A[0][1]*(A[1][2]*A[2][0]-A[1][0]*A[2][2]);
	D=D + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
	
	D0=D0+A[0][3]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
	D0=D0+A[0][1]*(A[1][2]*A[2][3]-A[1][3]*A[2][2]);
	D0=D0+A[0][2]*(A[1][3]*A[2][1]-A[1][1]*A[2][3]);
	
	D1=D1+A[0][0]*(A[1][3]*A[2][2]-A[1][2]*A[2][3]);
	D1=D1+A[0][3]*(A[1][2]*A[2][0]-A[1][0]*A[2][2]);
	D1=D1+A[0][2]*(A[1][0]*A[2][3]-A[1][3]*A[2][0]);
	
	D2=D2+A[0][0]*(A[1][1]*A[2][3]-A[1][3]*A[2][1]);
	D2=D2+A[0][1]*(A[1][3]*A[2][0]-A[1][0]*A[2][3]);
	D2=D2+A[0][3]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
	
	Sol[k][0]=D0/D;Sol[k][1]=D1/D;Sol[k][2]=D2/D;
}

void triangle::errors()
{
	float Af[3][5],Ac[3][5],err,w,xe,ye,xn,yn,Te,Tn,Te2,Tn2,axc[3],ayc[3];
	E1=fabs(Tf+To-2*Tc)/4;

	int nk,r,s;
	for(int k=0;k<3;k++)
	{
		if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}	
		w=e[k]*n[k];		nk=N[k];

		if(nk==0) cell[nk].Tc=(P[p[s]].Tc+P[p[r]].Tc)/2;
		Te=(cell[nk].Tc-Tc)/de[k];
		Tn=(P[p[s]].Tc-P[p[r]].Tc)/dn[k];

		Ac[k][0]=1;		Ac[k][1]=f[k].i;		Ac[k][2]=f[k].j;	
		Ac[k][3]=((Te-w*Tn)*e[k].i+(Tn-w*Te)*n[k].i)/(1-w*w);
		Ac[k][4]=((Te-w*Tn)*e[k].j+(Tn-w*Te)*n[k].j)/(1-w*w);

		if(nk==0) cell[nk].Tf=(P[p[s]].Tf+P[p[r]].Tf)/2;
		Te=(cell[nk].Tf-Tf)/de[k];
		Tn=(P[p[s]].Tf-P[p[r]].Tf)/dn[k];

		Af[k][0]=1;		Af[k][1]=f[k].i;		Af[k][2]=f[k].j;	
		Af[k][3]=((Te-w*Tn)*e[k].i+(Tn-w*Te)*n[k].i)/(1-w*w);
		Af[k][4]=((Te-w*Tn)*e[k].j+(Tn-w*Te)*n[k].j)/(1-w*w);

	}solution(Ac,axc,ayc);solution(Af,axf,ayf);

	Ec = Ixx*axc[1]+ Iyy*ayc[2]+Ixy*(ayc[1]+ayc[2]);		Ec = fabs(Ec)*0.5/Vol;
	Ef = Ixx*axf[1]+ Iyy*ayf[2]+Ixy*(ayf[1]+ayf[2]);		Ef = fabs(Ef)*0.5/Vol;
	
/*	for(int k=0;k<3;k++)
	{
	 	xe=e[k].i;	ye=e[k].j;
		xn=n[k].i;	yn=n[k].j;

	 	Te2=xe*xe*axf[1]+ye*ye*ayf[2]+xe*ye*(ayf[1]+axf[2])/2;
		Tn2=xn*xn*axf[1]+yn*yn*ayf[2]+xn*yn*(ayf[1]+axf[2])/2;
	
		E3=-Ce[k]*(d(C-f[k])-d(cell[nk].C-f[k]))*Te2/2;
		E3= E3-Cn[k]*(d(P[r]-f[k])-d(P[s]-f[k]))*Tn2/2;
	}
	E3 = fabs(E3)*dt/Vol;*/
}

//-----------**************--------------error functions-------------*************-----------------

void triangle::smothner()
{
	float a[9];int nk,i;
	
	i=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN-1;if(i==-1)i=0;
	if(fabs(Tyy)>1e4&&i!=0)
	{
//printf("j=%i\n",j);
		Tx=(cell[N[0]].Tx+cell[N[1]].Tx+cell[N[2]].Tx)/(3-i);
		Ty=(cell[N[0]].Ty+cell[N[1]].Ty+cell[N[2]].Ty)/(3-i);

		Txx=(cell[N[0]].Txx+cell[N[1]].Txx+cell[N[2]].Txx)/(3-i);
		Tyy=(cell[N[0]].Tyy+cell[N[1]].Tyy+cell[N[2]].Tyy)/(3-i);
		Txy=(cell[N[0]].Txy+cell[N[1]].Txy+cell[N[2]].Txy)/(3-i);

		Txxx=(cell[N[0]].Txxx+cell[N[1]].Txxx+cell[N[2]].Txxx)/(3-i);
		Txxy=(cell[N[0]].Txxy+cell[N[1]].Txxy+cell[N[2]].Txxy)/(3-i);
		Txyy=(cell[N[0]].Txyy+cell[N[1]].Txyy+cell[N[2]].Txyy)/(3-i);
		Tyyy=(cell[N[0]].Tyyy+cell[N[1]].Tyyy+cell[N[2]].Tyyy)/(3-i);
	}
}
void Derivatives(int A)
{
	for(int j=1;j<=nmax;j++) cell[j].Interiorderivatives(A);
	for(int j=1;j<=nmax;j++) cell[j].Boundryderivatives(A);
//	for(int j=1;j<=nmax;j++) cell[j].smothner();
}
//-----------**************--------------diffusion functions-------------*************----
/*void triangle::trucationerrors(int Ar)
{
	int nk,r,s,Bn;
	float x,y,xe,ye,xn,yn,deno,Cxx,Cyy,Cxy,Fxx,Fyy,Fxy,Pxx,Pyy,Pxy,Tee,Tnn,Ten,w,Ce,E3=0;//E3c=0;E3f=0;
 	if(P[p[0]].BN+P[p[1]].BN+P[p[2]].BN==0)
 	{
		for(int k=0;k<3;k++)
		{
			if(k==0) {r=1;s=2;}		if(k==1) {r=2;s=0;}		if(k==2) {r=0;s=1;}	
	
			xe=e[k].i;	ye=e[k].j;	xn=n[k].i;	yn=n[k].j;
	
			deno=(1/dc[k]+1/dF[k]);
			x=f[k].i-C.i;	y=f[k].j-C.j;

			Cxx=Txx;//+Txxx*x+Txxy*y;
			Cyy=Tyy;//+Tyyy*y+Txyy*x;
			Cxy=Txy;//+Txxy*x+Txyy*y;

			x=f[k].i-cell[N[k]].C.i;	y=f[k].j-cell[N[k]].C.j;

			Fxx=cell[N[k]].Txx;//+cell[N[k]].Txxx*x+cell[N[k]].Txxy*y;
			Fyy=cell[N[k]].Tyy;//+cell[N[k]].Tyyy*y+cell[N[k]].Txyy*x;
			Fxy=cell[N[k]].Txy;//+cell[N[k]].Txxy*x+cell[N[k]].Txyy*y;

			if(P[cell[N[k]].p[0]].BN+P[cell[N[k]].p[1]].BN+P[cell[N[k]].p[2]].BN>0)
			{Fxx=Cxx;Fyy=Cyy;Fxy=Cxy;}
			
			Pxx=(Cxx/dc[k]+Fxx/dF[k])/deno;
			Pyy=(Cyy/dc[k]+Fyy/dF[k])/deno;
			Pxy=(Cxy/dc[k]+Fxy/dF[k])/deno;
	
			//Te=xeTx+yeTy;
		 	Tee=xe*xe*Pxx+ye*ye*Pyy+2*xe*ye*Pxy;
			Tnn=xn*xn*Pxx+yn*yn*Pyy+2*xn*yn*Pxy;
			Ten=xe*xn*Pxx+ye*yn*Pyy+(xe*yn+ye*xn)*Pxy;
	
			w=e[k]*n[k];Ce=e[k]*A[k]/(1-w*w);
			E3=E3+((dF[k]-dc[k])*Tee/2-(ds[k]-dr[k])*w*Tnn/2)*Ce;
			E3=E3+((ds[k]-dr[k])*(Ten-w*Tnn)/2)*Ce;
		}
	}
	E3 = -kf*dt*E3/Vol;
	if(Ar==1)E3c=E3;if(Ar==2)E3f=E3;
	
	vec a,dv;float EE,tt,Sd,Sxf,Sxb,Syf,Syb,Sff,Sbb,Sfb,Sbf;
	if(Ar==1)tt=t;if(Ar==2)tt=t+dt;

	dv.i=.01;dv.j=0;a=C+dv;gauss(a,tt,Sxf);		a=C-dv;gauss(a,tt,Sxb);
	dv.i=0;dv.j=.01;a=C+dv;gauss(a,tt,Syf);		a=C-dv;gauss(a,tt,Syb);

	dv.i=.01;dv.j=.01;a=C+dv;gauss(a,tt,Sff);	a=C-dv;gauss(a,tt,Sbb);
	dv.i=.01;dv.j=-.01;a=C+dv;gauss(a,tt,Sfb);	a=C-dv;gauss(a,tt,Sbf);

	Sxx=(Sxf+Sxb-2*S)/1e-4;	
	Syy=(Syf+Syb-2*S)/1e-4;	
	Sxy=(Sff+Sbb-Sbf-Sfb)/4e-4;

	Sd = Ixx*Sxx/2 + Iyy*Syy/2 + Ixy*Sxy;		Sd = Sd/Vol;
	EE = Ixx*Txx/2 + Iyy*Tyy/2 + Ixy*Txy;		EE = EE/Vol;

	Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
	if(Ar==1) {Ec=EE;Sdc=Sd;}
	if(Ar==2){Ef=EE;E2=Ef-Ec;Sdf=Sd;	if(Bn>0)E2=0;}
}
//-----------**************--------------diffusion functions-------------*************----
void solution(float A[][4],float Sol[][3],int k)
{
	float D0=0,D1=0,D2=0,D=0;
	
	D=D + A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
	D=D + A[0][1]*(A[1][2]*A[2][0]-A[1][0]*A[2][2]);
	D=D + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
	
	D0=D0+A[0][3]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
	D0=D0+A[0][1]*(A[1][2]*A[2][3]-A[1][3]*A[2][2]);
	D0=D0+A[0][2]*(A[1][3]*A[2][1]-A[1][1]*A[2][3]);
	
	D1=D1+A[0][0]*(A[1][3]*A[2][2]-A[1][2]*A[2][3]);
	D1=D1+A[0][3]*(A[1][2]*A[2][0]-A[1][0]*A[2][2]);
	D1=D1+A[0][2]*(A[1][0]*A[2][3]-A[1][3]*A[2][0]);
	
	D2=D2+A[0][0]*(A[1][1]*A[2][3]-A[1][3]*A[2][1]);
	D2=D2+A[0][1]*(A[1][3]*A[2][0]-A[1][0]*A[2][3]);
	D2=D2+A[0][3]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
	
	Sol[k][0]=D0/D;Sol[k][1]=D1/D;Sol[k][2]=D2/D;
}*/
int fact(int r){if(r==0||r==4)return 24;if(r==1||r==3)return 6;if(r==2) return 4;}
void triangle::trucationerrors(int Ar)
{

	int nk,r,s,Bn,ni=9,i,ii;float x,y,A2[3][4],Sol[3][3],Ti[5];bool A=0;
	double A1[9][10],pp=0,x1[9],S1[10];

	Bn=P[p[0]].BN+P[p[1]].BN+P[p[2]].BN;
 	if(Bn==0)
 	{

//	Txx=Txx+Txxx*x+Txxy*y+ Txxxx*x*x/2 + Txxyy*y*y/2 + Txxxy*x*y;
//	Tyy=Tyy+Tyyy*y+Txyy*x+ Txxyy*x*x/2 + Tyyyy*y*y/2 + Txyyy*x*y;
//	Txy=Txy+Txxy*x+Txyy*y+ Txxxy*x*x/2 + Txyyy*y*y/2 + Txxyy*x*y;

for(int k=0;k<3;k++){x=cell[N[k]].C.i-C.i;y=cell[N[k]].C.j-C.j;A2[k][0]=x*x/2;A2[k][1]=y*y/2;A2[k][2]=x*y;}
for(int k=0;k<3;k++){x=cell[N[k]].C.i-C.i;y=cell[N[k]].C.j-C.j;A2[k][3]=cell[N[k]].Txx-Txx-Txxx*x-Txxy*y;}solution(A2,Sol,0);
for(int k=0;k<3;k++){x=cell[N[k]].C.i-C.i;y=cell[N[k]].C.j-C.j;A2[k][3]=cell[N[k]].Tyy-Tyy-Txyy*x-Tyyy*y;}solution(A2,Sol,1);
for(int k=0;k<3;k++){x=cell[N[k]].C.i-C.i;y=cell[N[k]].C.j-C.j;A2[k][3]=cell[N[k]].Txy-Txy-Txxy*x-Txyy*y;}solution(A2,Sol,2);
Ti[4]=Sol[0][0];Ti[0]=Sol[1][1];Ti[2]=(Sol[1][0]+Sol[0][1]+Sol[2][2])/3;
Ti[3]=(Sol[2][0]+Sol[0][2])/2;Ti[1]=(Sol[1][2]+Sol[2][1])/2;
E2=Ti[0]*ai[0]/24+Ti[1]*ai[1]/6+Ti[2]*ai[2]/4+Ti[3]*ai[3]/6+Ti[4]*ai[4]/24;
E2=Ti[2]*ai[2]/4;

	for(int k=0;k<3;k++)
	{
	for(int kk=0;kk<3;kk++)
	{
		if(j==cell[N[k]].N[kk]) nk=N[k];	else nk=cell[N[k]].N[kk];

//		if(nk!=N[k]&&ii>0) for(int h=0;h<3;h++)
//		if(cell[nk].p[h]==ii){nk=cell[nk].N[h];ii=0;}
	
		x=cell[nk].C.i-C.i;	y=cell[nk].C.j-C.j;

		i=3*k+kk;
		A1[i][9]=0;
	for(int r=0;r<5;r++)A1[i][9]=A1[i][9]-Ti[r]*pow(x,r)*pow(y,4-r)/fact(r);

		A1[i][0]=x;			A1[i][1]=y;
		A1[i][2]=x*x/2;		A1[i][3]=y*y/2;		A1[i][4]=x*y;
		A1[i][5]=x*x*x/6;	A1[i][6]=y*y*y/6;	A1[i][7]=x*x*y/2;	A1[i][8]=x*y*y/2;		
//		ff=ff+Ai[i][2]*A1[i][9];

	}
	}
//printf("a0=%e\n",ai[0]);getch();
		if(j==403&&pass==3&&nmax==687&&A){printf("\nj=%i initial matrix\n",j);for(int k=0;k<9;k++){for(int kk=0;kk<10;kk++)printf("%e ",A1[k][kk]); printf("\n");}printf("\n");getch();}
//		for(int jj=0;jj<=ni;jj++) {S1[jj]=0;for(int i=0;i<ni;i++) S1[jj]=S1[jj]+A1[i][jj];}
//		for(int i=0;i<ni;i++)for(int jj=0;jj<=ni;jj++) A1[i][jj]=S1[jj]-A1[i][jj];

		for(int i=0;i<ni;i++)//row change down
		{
//	if(j==403&&pass==3&&nmax==687&&A){printf("\nj=%i before row change\n",j);for(int k=0;k<9;k++){for(int kk=0;kk<10;kk++)printf("%e ",A1[k][kk]); printf("\n");}printf("\n");getch();};
		pp=A1[i][i]; 
		if(i<ni-1&&fabs(pp)<1e-13)
		{
		 	pp=0;nk=0;
			for(ii=i+1;ii<ni-1;ii++)if(fabs(A1[ii][i])>fabs(pp)){pp=A1[ii][i];nk=ii;}
 			for(int jj=i;jj<=ni;jj++) {pp=A1[i][jj];A1[i][jj]=A1[nk][jj];A1[nk][jj]=pp;}
		}	
		if(j==403&&pass==3&&nmax==687&&A){printf("\nj=%i after row change %i-%i\n",j,i,nk);for(int k=0;k<9;k++){for(int kk=0;kk<10;kk++)printf("%e ",A1[k][kk]); printf("\n");}printf("\n");getch();};

			pp=A1[i][i]; 
	//		if(fabs(pp)<1E-15)	break;printf("pivot %e too small in row %i\n",pp,i);getch();exit(0);}
			for(int jj=i;jj<=ni;jj++) A1[i][jj]=A1[i][jj]/pp;//aii=1
			for(int ii=0;ii<ni;ii++)//row change down
			{	
				if(ii==i)ii++; if(ii==ni)break;	
				pp=A1[ii][i];
				for(int jj=i;jj<=ni;jj++) A1[ii][jj]=A1[ii][jj]-pp*A1[i][jj]; 
			}
			if(j==403&&nmax==687&&A){printf("\nj=%i after Normalisation\n",j);for(int k=0;k<9;k++){for(int kk=0;kk<10;kk++)printf("%e ",A1[k][kk]); printf("\n");}printf("\n");getch();};
		}
		if(j==403&&nmax==687&&A){printf("\nj=%i final matrix\n",j);for(int k=0;k<9;k++){for(int kk=0;kk<10;kk++)printf("%e ",A1[k][kk]); printf("\n");}printf("\n");getch();};
	/*
	if(j==218){for(int k=0;k<9;k++){for(int kk=0;kk<10;kk++)printf("%e ",A1[k][kk]); printf("\n");}};

	x1[ni-1]=A1[ni-1][ni]/A1[ni-1][ni-1];//backward sweep on [U]
	for(int i=ni-2;i>=0;i--)
	{
		pp=0;
		for(int jj=i+1;jj<ni;jj++) pp=pp+x1[jj]*A1[i][jj];
		x1[i] = (A1[i][ni] - pp)/A1[i][i];
	}
	
	Tx=x1[7];	Ty=x1[8];	Txx=x1[4];	Tyy=x1[5];	Txy=x1[6];
	if(nk==0){	Tx=0;	Ty=0;	Txx=0;	Tyy=0;	Txy=0;}*/

	E2=A1[2][9]+A1[3][9];if(j==403&&nmax==687&&A){printf("\nE2=%f\n",E2);getch();}
	}else E2=0;
	E3c=Ti[4];

//	if(Ar==1) {Ec=EE;Sdc=Sd;}
//	if(Ar==2){Ef=EE;E2=Ef-Ec;Sdf=Sd;	if(Bn>0)E2=0;}
}


//-----------**************--------------diffusion functions-------------*************---
//-------------@@@@@@@@@@@@@*************@@@@@@@@@@@@@@@@@@@@***********@@@@@@@@@@@@@@---
//-----------**************---------------general functions--------------*************---
//-----------**************---------------general functions--------------*************-----------------
//-----------**************---------------general functions--------------*************-----------------
void reinitialize()
{
	for(int j=0;j<=Nmax;j++) 
	{
		cell[j].split=0;
		cell[j].served=false;
		cell[j].coarse=false;
		cell[j].mod=false;
		cell[j].midptallocation=false;
		cell[j].dM=false;
	}Nmax=nmax;
}
//-----------**************---------------general functions--------------*************-----------------
void renumber()
{
	int p,nk,r,s;
	while(cell[nmax].deleted) {nmax--;}
	for(int j=1;j<nmax;j++) 
	{
		if(cell[j].deleted)
		{
			cell[j]=cell[nmax];
			cell[j].modified();
			cell[j].j=j;
			cell[j].deleted=false;

			for(int k=0;k<3;k++)
			{
				p=cell[j].p[k];
				P[p].del(nmax);
				P[p].add(j);
			}cell[nmax].del();
			while(cell[nmax].deleted) {nmax--;}
		}
	}
	/*while(P[pmax].coarse) {P[pmax].del();pmax--;}
	for(int i=1;i<=pmax;i++) 
	{
		if(P[i].coarse)
		{
			P[i]=P[pmax];
			for(int j=1;j<=P[i].t[0];j++)
			{
		 	nk=P[i].t[j];
 			for(int k=0;k<3;k++)
			{
				if(cell[nk].p[k]==pmax) 
				{
					if(k==0){r=1;s=2;}	if(k==1){r=2;s=0;}	if(k==2){r=0;s=1;}	
		    		cell[nk].p[k]=i;
					if(P[cell[nk].p[r]].p1==pmax)P[cell[nk].p[r]].p1=i;
					if(P[cell[nk].p[r]].p2==pmax)P[cell[nk].p[r]].p2=i;

					if(P[cell[nk].p[s]].p1==pmax)P[cell[nk].p[s]].p1=i;
					if(P[cell[nk].p[s]].p2==pmax)P[cell[nk].p[s]].p2=i;

				}
			}
			}
			P[pmax].coarse=true;
			while(P[pmax].coarse) {P[pmax].del();pmax--;}
			P[i].del();
		}
	}*/
}


//-----------**************---------------general functions--------------*************-----------------
//-------------@@@@@@@@@@@@@*************@@@@@@@@@@@@@@@@@@@@***********@@@@@@@@@@@@@@-----------------
//-----------**************---------------vector functions---------------*************-----------------
void vec::nodetemp()
{
	float S=0,Sf=0,Bn;int pt;
	for(int j=1;j<=t[0];j++) 
	{
		pt=t[j];Bn=P[cell[pt].p[0]].BN+P[cell[pt].p[1]].BN+P[cell[pt].p[2]].BN;
  		Sf=Sf+cell[pt].Tf/d(*this-cell[pt].C);
		S=S+1/d(*this-cell[pt].C);
	}Tf=Sf/S;
	
/*	for(int j=1;j<=t[0];j++) 
	{
		pt=t[j];Bn=P[cell[pt].p[0]].BN+P[cell[pt].p[1]].BN+P[cell[pt].p[2]].BN;
		if(Bn==0)
		{
	  		Sf=Sf+map(*this-cell[pt].C,cell[pt].Tf,pt)/d(*this-cell[pt].C);
			S=S+1/d(*this-cell[pt].C);
		}
	}Tf=Sf/S;*/

}
void vec::add(int i)
{
	bool found=false;
	if(t[0]==ptlimit)
	{
		printf("error .Triangle %i exceded point limit of 12.\n",i);
		getch();exit(1);
	}
	for(int j=1;j<=t[0];j++)if(i==t[j]) found=true;
	if(i==0) found=true;
	if(!found){t[0]++;t[t[0]]=i;}

}
//-----------**************---------------vector functions---------------*************-----------------
void vec::del(int i)
{
	bool found=false;
	for(int j=1;j<=t[0];j++)
	{
		if(t[j]==i) found=true;
		if(found&&j!=t[0]) t[j]=t[j+1];
	}
	if(!found)
	{
	for(int j=1;j<=t[0];j++){printf("\n%i ",t[j]);}

		printf("\nerror .Triangle %i not found in triangle removal from point\n",i);
		printf("\nmax=%i\tpmax=%i",nmax,pmax);
		
		getch();
		exit(1);
	}
	t[t[0]]=0;t[0]--;
}
void triangle::del()
{
	deleted=true;coarse=false;
	split=0;type=0;
	mod=false;dM=false;
	midptallocation=false;
	served=true;
	modified();
	
//	P[p[0]].del(j);P[p[1]].del(j);P[p[2]].del(j);

	N[0]=0;N[1]=0;N[2]=0;
	p[0]=0;p[1]=0;p[2]=0;
}
	
//-----------**************---------------vector functions---------------*************-----------------
void vec::del()
{
	deleted=true;
	coarse=false;
	BN=true;
	p1=0;p2=0;i=0;j=0;level=0;
	for(int i=0;i<=9;i++){t[i]=0;}
}
//-----------**************---------------general functions--------------*************-----------------
int newtri()
{
	int j;
	for(j=nmax;j>0;j--) if(cell[j].deleted) break;

	if(j==0) {nmax++;return nmax;}
	else{cell[j].deleted=false;return j;}
}
//-----------**************---------------general functions--------------*************-----------------
int newpt()
{
	int i;
	for(i=pmax;i>0;i--) if(P[i].deleted) break;

	if(i==0) {pmax++;return pmax;}
	else{P[i].deleted=false;return i;}
}

//-----------**************---------------general functions--------------*************-----------------
//-----------**************---------------general functions--------------*************-----------------
//-----------**************---------------vector functions---------------*************-----------------
//-------------@@@@@@@@@@@@@*************@@@@@@@@@@@@@@@@@@@@***********@@@@@@@@@@@@@@-----------------
