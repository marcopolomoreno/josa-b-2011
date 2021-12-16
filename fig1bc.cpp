#include <stdio.h>
#include <math.h>

double soma,q33,q12,q13,gama,phi,x1,x2,x3,z1,z2;
double d,w,w1,w2,w3,Omega,h,t,area;
double Tp,Tr,fr,a10,a20,a30,alpha;
int Pulsos,i,j,k,m,n,g,N,PassoRKFento;
double a[9],b[9],c[9],k1[9],k2[9],k3[9],k4[9];
double Pi=3.141592654;

double f1(double a1, double a2, double a3,
          double a12,double b12,double a13,
          double b13,double a23,double b23,int i)  //sistema Lambda
{
    if (i==1) return  2*Omega*(cos(alpha)*b13-sin(alpha)*a13)                +0.5*q33*a3     -(a1-a10)*gama;
    if (i==2) return  2*Omega*(cos(alpha)*b23-sin(alpha)*a23)                +0.5*q33*a3     -(a2-a20)*gama;
    if (i==3) return  2*Omega*(sin(alpha)*(a13+a23)-cos(alpha)*(b13+b23))        -q33*a3       -(a3-0)*gama;
    
    if (i==4) return -(w2-w1)*b12+Omega*(cos(alpha)*(b13+b23)-sin(alpha)*(a13+a23))-a12*q12-a12*gama;
    if (i==5) return  (w2-w1)*a12+Omega*(cos(alpha)*(a23-a13)+sin(alpha)*(b23-b13))-b12*q12-b12*gama;
   
    if (i==6) return  -(d-w1)*b13+Omega*cos(alpha)*b12+Omega*sin(alpha)*(a12+a1-a3)-0.5*a13*q33-a13*gama;
    if (i==7) return   (d-w1)*a13+Omega*sin(alpha)*b12+Omega*cos(alpha)*(a3-a1-a12)-0.5*b13*q33-b13*gama;
    if (i==8) return  -(d-w2)*b23-Omega*cos(alpha)*b12+Omega*sin(alpha)*(a12+a2-a3)-0.5*a23*q33-a23*gama;
    if (i==9) return   (d-w2)*a23-Omega*sin(alpha)*b12+Omega*cos(alpha)*(a3-a2-a12)-0.5*b23*q33-b23*gama;
    
}

main() 
{
    FILE *arquivo;
    arquivo=fopen("dados.dat","w");
    
    //31-05-2010
    //modifica��es em 28/07/2010
    //Nome da pasta: tempo - fases - semianalitico - Omega
    //introdu��o das fases w*Tr e phi
    //calculos numericos na presen�a do campo
    //calculos analiticos no decaimento
    //Plota Popula��es e coer�ncias vs tempo
    //Para sistemas em lambda
    
    q33=(2*Pi)*4e6;
    q13=0.5*q33;
    q12=0*q13;
    fr=50*q13/(2*Pi);
    Tp=100e-15;
    w1=0;
    w2=2*Pi*70*fr;
    w3=(2*Pi)*4e6*fr+0;
    w= (2*Pi)*4e6*fr-0.5*w2; 
    PassoRKFento=10;
    phi=(2*Pi)*0;
    
    Omega=2*q13/(fr*Tp);
    Pulsos=200;
    gama=0; 
    a10=0.5;                     //popula��o inicial do estado 1
    a20=0.5;                     //popula��o inicial do estado 2
    a30=0;                       //popula��o inicial do estado 3

    a[1]=a10;a[2]=a20;a[3]=a30;
    a[4]=0;a[5]=0;a[6]=0;
    a[7]=0;a[8]=0;a[9]=0;
    t=0;N=-1;
    d=w3-w;
    Tr=1/fr;
    t=-Tp/2;
    
       for (i=1;i<=2*Pulsos+1;i++)
       {
           printf("*");                   
           if (i % 2 == 0)         
           {
                 N=N+1;
                 g=PassoRKFento;
                 h=(1/double(g))*Tp;
                 
                 alpha=-N*w*Tr+N*phi;  
                                  
                         for (k=1;k<=g;k++)
                         {
                         t=t+h;
                                                                         
                         for (j=1;j<=9;j++)
                         k1[j]=f1(a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],j);
                                                           
                         for (j=1;j<=9;j++)
                         k2[j]=f1(a[1]+k1[1]*h/2,a[2]+k1[2]*h/2,a[3]+k1[3]*h/2,
                         a[4]+k1[4]*h/2,a[5]+k1[5]*h/2,a[6]+k1[6]*h/2,a[7]+k1[7]*h/2,
                         a[8]+k1[8]*h/2,a[9]+k1[9]*h/2,j);
                  
                         for (j=1;j<=9;j++)
                         k3[j]=f1(a[1]+k2[1]*h/2,a[2]+k2[2]*h/2,a[3]+k2[3]*h/2,
                         a[4]+k2[4]*h/2,a[5]+k2[5]*h/2,a[6]+k2[6]*h/2,a[7]+k2[7]*h/2,
                         a[8]+k2[8]*h/2,a[9]+k2[9]*h/2,j);
                  
                         for (j=1;j<=9;j++)           
                         k4[j]=f1(a[1]+k3[1]*h,a[2]+k3[2]*h,a[3]+k3[3]*h,a[4]+k3[4]*h,
                         a[5]+k3[5]*h,a[6]+k3[6]*h,a[7]+k3[7]*h,a[8]+k3[8]*h,
                         a[9]+k3[9]*h,j);

                         for (j=1;j<=9;j++)
                         b[j]=a[j]+h*(k1[j]/6+k2[j]/3+k3[j]/3+k4[j]/6);   
                  
                         for (m=1;m<=9;m++)
                         a[m]=b[m];     
                         } 
                             
                 }
               
           if (i % 2 == 1)
           {
                t=t+(Tr-Tp);
                x1=(w2-w1)*(Tr-Tp);
                x2=(d -w1)*(Tr-Tp);
                x3=(d -w2)*(Tr-Tp);
                z1=q33*(Tr-Tp);
                z2=q12*(Tr-Tp);
                                         
              b[1]=a[1]+0.5*a[3]*(1-exp(-z1));
              b[2]=a[2]+0.5*a[3]*(1-exp(-z1));
              b[3]=a[3]*exp(-z1);
              b[4]=(a[4]*cos(x1)-a[5]*sin(x1))*exp(-z2);
              b[5]=(a[4]*sin(x1)+a[5]*cos(x1))*exp(-z2);
              b[6]=(a[6]*cos(x2)-a[7]*sin(x2))*exp(-0.5*z1);
              b[7]=(a[6]*sin(x2)+a[7]*cos(x2))*exp(-0.5*z1);
              b[8]=(a[8]*cos(x3)-a[9]*sin(x3))*exp(-0.5*z1);
              b[9]=(a[8]*sin(x3)+a[9]*cos(x3))*exp(-0.5*z1);
              
              for (m=1;m<=9;m++)
              a[m]=b[m];     
                 }      
          
        soma=b[1]+b[2]+b[3];     
        printf("%d %12.10f %12.10f %12.10f %12.10f %12.10f", 
               N,b[1],b[2],b[3],soma,sqrt(b[4]*b[4]+b[5]*b[5]));
        printf("\n");       
        fprintf(arquivo,"%12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f",
               t/Tr-0.9999,b[1],b[2],b[3],sqrt(b[4]*b[4]+b[5]*b[5]),sqrt(b[6]*b[6]+b[7]*b[7]),
               sqrt(b[8]*b[8]+b[9]*b[9]),0.5*b[1]+0.5*b[2]+b[4],
                0.5*b[1]+0.5*b[2]-b[4],b[4],b[5],b[6],b[7],b[8],b[9]);   
        fprintf(arquivo,"\n");  
           
        }
        
     fclose(arquivo);
     printf("\a");
}
