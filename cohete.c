////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//      El problema de los tres cuerpos - Algoritmo Runge-Kutta
//
//      El ejercicio trata de mandar un cohete a la luna. Tenemos que calcular la trayectoria que sigue
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define TMAX 500000
#define h 50     //1 minuto aprox

#define PI 3.14159265359
#define dTL 1 
#define w 2.66E-6
#define delta 7.01E-12
#define mu 1.23E-2

#define G 6.67e-11
#define MT 5.9736e24
#define ML 0.07349e24
#define RT 6.378160e6
#define RL 1.7374e6

double r0,phi0,Pr0,Pphi0;           //Condiciones iniciales
double r,phi,Pr,Pphi,rprima;        //Ecuaciones de movimiento
double k1,k2,k3,k4;                 //Coeficientes para el algoritmo 
double rluna[2], rtierra[2], rcohete[2];
double auxr, auxphi, auxPr, auxPphi;
double n,g;                           //contador que imprima 1 de cada x iteraciones
double H, H1, r_luna;

FILE *f1, *f2;

int EcMovimiento(int t); //Funcion que teniendo r0,phi0,Pr0,Pphi0 calcula las r,phi,Pr,Pphi (teniendo r calcula r'...)
int PrintfEcIniciales();    
int CalcCoeficientesr(int t);
int CalcCoeficientesphi(int t);
int CalcCoeficientesPr(int t);
int CalcCoeficientesPphi(int t);


///////////////////////////////////////////////////FUNCIÓN PRINCIPAL///////////////////////////////////////////////////////////////////

int main()
{
    f1=fopen("Posiciones.dat","w");
    f2=fopen("Hamiltoniano.dat","w");
    int t;

    //Damos valores de las condiciones iniciales
    r0=0.0166;
    phi0=PI/3;
    Pr0=0.000029205;
    Pphi0=0.000000;

    rtierra[0]=rtierra[1]=0.0;
    t=0;
    n=0;
    g=0;

    while (t<TMAX)
    {
        rluna[0]=dTL*cos(w*t);  //Calculamos la posición de la luna 
        rluna[1]=dTL*sin(w*t);

        rcohete[0]=r0*cos(phi0);  //Calculamos la posición del cohete
        rcohete[1]=r0*sin(phi0);

        //printf("Condiciones iniciales para t=%i\n\n",t);
        //PrintfEcIniciales();

        

        if(n==50) //Cuando coincida con 1 de cada 1000 valores, guardamos el dato de la posición en un fichero
        {
            fprintf(f1,"%E\t%E\t",rtierra[0],rtierra[1]);
            fprintf(f1,"%E\t%E\t",rcohete[0],rcohete[1]);
            fprintf(f1,"%E\t%E\t",rluna[0],rluna[1]);   

            fprintf(f1,"\n"); 
            n=0;

        }

        


        EcMovimiento(t);    //Calculamos r' teniendo r etc.

        //printf("Condiciones despues de calcular las ec de mvto para t=%i\n\n",t);
        //PrintfEcIniciales();

        //Calculamos los coeficientes para cada variable: r,phi,Pr,Pphi y calculamos y(t+h)


        CalcCoeficientesr(t);
        CalcCoeficientesphi(t);
        CalcCoeficientesPr(t);
        CalcCoeficientesPphi(t);

        r_luna=sqrt((r*r)+(dTL*dTL)-(2.0*r*dTL*cos(phi-(w*t))));

        H=Pr*Pr/(2.0)+Pphi*Pphi/(2.0*r*r)-G*MT/r-G*ML/(r_luna);

        H1=H-Pphi*w;

        fprintf(f2,"%E\t%E\t",g, H1);
        fprintf(f2,"\n");

        //printf("Condiciones despues de calcular los coeficientes para t=%i\n\n",t);
        //PrintfEcIniciales();

        //Tras esto tenemos el valor de r(t+h), phi(t+h),... en las variables aux. Las pasamos a las variables originales

        r0=auxr;
        phi0=auxphi;
        Pr0=auxPr;
        Pphi0=auxPphi;

        printf("Condiciones finales para t=%i\n\n",t);
        PrintfEcIniciales();

        t=t+h;
        n=n+1;
        g=g+1;
    }

    fclose(f1);
    fclose(f2);

    return 0;
}






///////////////////////////////////////////////////////FUNCIONES////////////////////////////////////////////////////////////////

int EcMovimiento(int t)
{
    r=Pr0;
    phi=Pphi0/(pow(r0,2.0));
    rprima=sqrt(1.0+pow(r0,2.0)-2*r0*cos(phi0-w*t));
    Pr=((pow(Pphi0,2.0))/(pow(r0,3.0)))-delta*(pow(r0,-2.0)+mu*pow(rprima,-3.0)*(r0-cos(phi0-w*t)));
    Pphi=-delta*mu*r0*sin(phi0-w*t)/pow(rprima,3.0);


    return 0;
}


int PrintfEcIniciales()
{
    printf("Valores de las Ec de movimiento:\n\n");
    printf("r=%E\n\n",r);
    printf("phi=%E\n\n",phi);
    printf("Pr=%E\n\n",Pr);
    printf("Pphi=%E\n\n",Pphi);
    printf("\n\nValores de las constantes\n\n");
    printf("r0=%E\n\n",r0);
    printf("phi0=%E\n\n",phi0);
    printf("Pr0=%E\n\n",Pr0);
    printf("Pphi0=%E\n\n",Pphi0);


    return 0;
}


int CalcCoeficientesr(int t)
{
    k1=h*r;

    r0=r0+(k1/2.0);

    EcMovimiento((t+(h/2)));

    k2=h*r;

    r0=r0-(k1/2.0)+(k2/2.0);

    EcMovimiento((t+(h/2)));

    k3=h*r;

    r0=r0-(k2/2.0)+(k3);

    EcMovimiento((t+h));

    k4=h*r;

    r0=r0-k3;   //Volvemos al mismo valor de r0 que teníamos antes de calcular los coeficientes

    //Calculamos r(t+h)

    auxr=r0+((1.0/6.0)*(k1+(2.0*k2)+(2.0*k3)+k4));  //Guardamos en una variable aux el valor de r(t+h)  

    //printf("\n\nValores de las k y aux de r para t=%i\n\n",t);
    //printf("k1=%E\nk2=%E\nk3=%E\nk4=%E\nauxr=%E\n\n",k1,k2,k3,k4,auxr);


    return 0;
}


int CalcCoeficientesphi(int t)
{
    k1=h*phi;

    phi0=phi0+(k1/2.0);

    EcMovimiento((t+(h/2)));

    k2=h*phi;

    phi0=phi0-(k1/2.0)+(k2/2.0);

    EcMovimiento((t+(h/2)));

    k3=h*phi;

    phi0=phi0-(k2/2.0)+(k3);

    EcMovimiento((t+h));

    k4=h*phi;

    phi0=phi0-k3;   //Volvemos al mismo valor de phi que teníamos antes de calcular los coeficientes

    //Calculamos phi(t+h)

    auxphi=phi0+((1.0/6.0)*(k1+(2.0*k2)+(2.0*k3)+k4));  //Guardamos en una variable aux el valor de phi(t+h)  

    //printf("\n\nValores de las k y aux de Pr para t=%i\n\n",t);
    //printf("k1=%E\nk2=%E\nk3=%E\nk4=%E\nauxphi=%E\n\n",k1,k2,k3,k4,auxphi);


    return 0;
}


int CalcCoeficientesPr(int t)
{
    k1=h*Pr;

    Pr0=Pr0+(k1/2.0);

    EcMovimiento((t+(h/2)));

    k2=h*Pr;

    Pr0=Pr0-(k1/2.0)+(k2/2.0);

    EcMovimiento((t+(h/2)));

    k3=h*Pr;

    Pr0=Pr0-(k2/2.0)+(k3);

    EcMovimiento((t+h));

    k4=h*Pr;

    Pr0=Pr0-k3;   //Volvemos al mismo valor de Pr que teníamos antes de calcular los coeficientes

    //Calculamos Pr(t+h)

    auxPr=Pr0+((1.0/6.0)*(k1+(2.0*k2)+(2.0*k3)+k4));  //Guardamos en una variable aux el valor de Pr(t+h)  

    //printf("\n\nValores de las k y aux de Pr para t=%i\n\n",t);
    //printf("k1=%E\nk2=%E\nk3=%E\nk4=%E\nauxPr=%E\n\n",k1,k2,k3,k4,auxPr);


    return 0;
}


int CalcCoeficientesPphi(int t)
{
    k1=h*Pphi;

    Pphi0=Pphi0+(k1/2.0);

    EcMovimiento((t+(h/2)));

    k2=h*Pphi;

    Pphi0=Pphi0-(k1/2.0)+(k2/2.0);

    EcMovimiento((t+(h/2)));

    k3=h*Pphi;

    Pphi0=Pphi0-(k2/2.0)+(k3);

    EcMovimiento((t+h));

    k4=h*Pphi;

    Pphi0=Pphi0-k3;   //Volvemos al mismo valor de Pphi que teníamos antes de calcular los coeficientes

    //Calculamos Pphi(t+h)

    auxPphi=Pphi0+((1.0/6.0)*(k1+(2.0*k2)+(2.0*k3)+k4));  //Guardamos en una variable aux el valor de Pphi(t+h)  

    //printf("\n\nValores de las k y aux de Pr para t=%i\n\n",t);
    //printf("k1=%E\nk2=%E\nk3=%E\nk4=%E\nauxPphi=%E\n\n",k1,k2,k3,k4,auxPphi);


    return 0;
}