#include <iostream>
#include <cmath>

using namespace std;
int I = 10;
double M[] = { 20, 75,100,105,75,0 };
double V[] = {0, 75, 150, 200, 250, 300};
double Tp = 110;
double Hm = 0.01;
double Hv = 0.0001;
double C = 0.1;
int size_of = 5;

//dT/dt = Vн + Vc
//M(V) = k*V + m на каждом из кусочно-линейных участков
//dV/dt = M(V)/I линейное ДУ из условия на ускорение
//V(t) = c*exp(k/I * t)-m/k - решение ДУ, на каждом из участков с, k, m - разные


double k(double a, double b) { //функция, вычисляющая тангенс угла наклона на каждом участке
    if (b != 0) {
        double r = a/b;
        return r;
    }
}

double m(int x, int y, double k) { //коэффициент m
    double r = y - k * x;
    return r;
}

void init_koef(double moment[], double speed[], double *K, double *M) { //вычисляет коэффициенты на каждом участке
    
    for (int i = 0; i < size_of; i++) {
        K[i] = k(moment[i + 1] - moment[i], speed[i + 1] - speed[i]);
        M[i] = m(speed[i], moment[i], K[i]);
    }
}

void init_time(double* m, double* k, double* time, double* coef) {

    for (int i = 1; i < size_of; i++) {
        coef[i - 1] = (V[i - 1] + m[i - 1] / k[i - 1]) / (exp((k[i - 1] / I) * time[i - 1]));
        time[i] = I * ((log((V[i] + m[i - 1] / k[i - 1]) / coef[i - 1])) / k[i - 1]);
    }
    coef[size_of - 1] = (V[size_of - 1] + m[size_of - 1] / k[size_of - 1]) / (exp(k[size_of - 1] / I * time[size_of - 1]));
    }

double function(double t, double Temp, double k, double m, double c, double parametr) { 
    double velocity = c* exp(t * (k / I)) - m/k;
    double moment = k * velocity + m;
    return moment * Hm + velocity*velocity*Hv + C*parametr - C*Temp;
}

double runge_kutta(double* time, double k, double m, double c, double parametr, double* T_nach) { //решатель ДУ для каждого участка 
    double a = time[0];
    double b;
    if (time[1] < 0) {
        b = time[0] + 10;
        time[1] = b;
    }
    else {
        b = time[1];
    }
    double h = (b - a) / 100;
    const int n = 99;
    double t[100];
    double V1[100];
    double V2[100];
    double V3[100];
    double V4[100];
    double T[100];
    t[0] = time[0];
    T[0] = T_nach[0];
       
    for (int i = 1; i <= n; i++) {
        t[i] = a + i * h;
        V1[i] = h * function(t[i - 1], T[i - 1], k, m, c, parametr);
        V2[i] = h * function(t[i - 1] + h / 2.0, T[i - 1] + V1[i] / 2.0, k, m, c, parametr);
        V3[i] = h * function(t[i - 1] + h / 2, T[i - 1] + V2[i] / 2, k, m, c, parametr);
        V4[i] = h * function(t[i - 1] + h, T[i - 1] + V3[i], k, m, c, parametr);
        T[i] = T[i - 1] + (V1[i] + 2 * V2[i] + 2 * V3[i] + V4[i]) / 6;
    }
    T_nach[1] = T[99];
    
    for (int i = 0; i < 100; i++) {
        if (T[i] > Tp || T[i] == Tp) {
            return t[i];
        }
    }
    return -1;
}
void model() {

    double* T_nach = new double[size_of];
    double* timer = new double[size_of];
    double* coef = new double[size_of];
    double* k = new double[size_of];
    double* m = new double[size_of];
    timer[0] = 0;

    init_koef(M, V, k, m);
    init_time(m, k, timer, coef);

    timer[0] = 0;
    double param;
    cin >> param;
    T_nach[0] = param;
    int i;
    for (i = 0; i < size_of - 1; i++) {
        double r = runge_kutta(timer + i, k[i], m[i], coef[i], param, T_nach + i);
        if (r != -1) {
            cout << "Time is " << r << endl;
            break;
        }
    }

    double r = -1;
    int lenght = size_of;
    if (i == size_of - 1) {
        int j = i;
        while (r == -1) {
            if (j == 5000) {
                break;
            }
            lenght++;

            double* temp = new double[lenght];
            memcpy(temp, T_nach, (lenght - 1) * sizeof(double));
            delete[] T_nach;
            T_nach = temp;

            double* q = new double[lenght];
            memcpy(q, timer, (lenght - 1) * sizeof(double));
            delete[] timer;
            timer = q;

            r = runge_kutta(timer + j, k[i], m[i], coef[i], param, T_nach + j);
            if (r != -1) {
                cout << "Time is " << r << endl;
                return;
            }
            j++;
        }
        cout << "Ne peregeetsya" << endl;
        return;
    }
}

int main(){

    model();

    return 0;
}