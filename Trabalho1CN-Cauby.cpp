#include <iostream>
#include <fstream>
#include <cmath>
#include<iomanip>

using namespace std;

double f(double x) {
// função que define a equação a ser resolvida
return x*x*x - x - 1;
}
double df(double x) {
// Defina a derivada da função que deseja utilizar aqui
return 3*x*x - 1;
}

void bisseccao(double a, double b, double eps, int max_it){
    double x, f_a, f_b, f_x;
    int i = 1;
    ofstream file;
    file.open("bisseccao.txt");
    file << "Iteracao\tRaiz\t\tErro\n";
    file << "--------------------------------------------------\n";

    while(i <= max_it){
        x = (a + b) / 2;
        f_a = f(a);
        f_b = f(b);
        f_x = f(x);

        file << i << "\t\t" << x << "\t\t" << abs(b - a) << "\n";

        if(abs(b - a) < eps){
            cout << "Solucao encontrada: x = " << x << "\n";
            cout << i << " iteracoes.\n";
            file.close();
            return;
        }

        if(f_x * f_a < 0){
            b = x;
        } else {
            a = x;
        }

        i++;
    }

    cout << "O metodo da bisseccao falhou apos " << max_it << " iteracoes.\n";
    file.close();
}

void mil(double a, double b, double eps, int max_iter) {
    double x0 = a;
    double x1 = b;
    int i = 0;
    ofstream file("mil.txt"); // arquivo de saída
    file << "Iteracao\tRaiz\t\tErro\n";
     file << "--------------------------------------------------\n";
    while (i < max_iter) {
        double f0 = f(x0);
        double f1 = f(x1);
        double x = x1 - f1*(x1 - x0)/(f1 - f0);
        double fx = f(x);
        file << i+1 << "\t\t" << x << "\t\t" << abs(fx) << endl; // Salva a interação e a solução no arquivo
        if (abs(fx) < eps) {
            cout << "A raiz encontrada e: " << x << endl;
            cout << i << " iteracoes.\n";
            file.close();
            return;
        }
            if (f0*fx < 0) {
            x1 = x;
        } else {
            x0 = x;
        }
            i++;
    }
    cout << "O Metodo nao convergiu." << endl;
    file.close();
}


void newton(double x0, double eps, int max_it){
    double x, f_x, df_x;
    int i = 1;
    ofstream file;
    file.open("newton.txt");
    file << "Iteracao\tRaiz\t\tErro\n";
     file << "--------------------------------------------------\n";

    while(i <= max_it){
        f_x = f(x0);
        df_x = df(x0);

        x = x0 - f_x / df_x;

        file << i << "\t\t" << x << "\t\t" << abs(x-x0) << "\n";

        if(abs(x - x0) < eps){
            cout << "Solucao encontrada: x = " << x << "\n";
            cout << i << " iteracoes.\n";
            file.close();
            return;
        }

        x0 = x;
        i++;
    }

    cout << "O metodo de Newton falhou apos " << max_it << " iteracoes.\n";
    file.close();
}

void secante(double x0, double x1, double eps, int max_it){
    double x, f_x0, f_x1, f_x;
    int i = 1;
    ofstream file;
    file.open("secante.txt");
    file << "Iteracao\tRaiz\t\tErro\n";
     file << "--------------------------------------------------\n";

    f_x0 = f(x0);
    f_x1 = f(x1);

    while(i <= max_it){
        x = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);

        f_x = f(x);
        file << i << "\t\t" << x << "\t\t" << abs(x - x1) << "\n";

        if(abs(x - x1) < eps){
            cout << "Solucao encontrada: x = " << x << "\n";
            cout << i << " iteracoes.\n";
            file.close();
            return;
        }

        x0 = x1;
        f_x0 = f_x1;

        x1 = x;
        f_x1 = f_x;

        i++;
    }

    cout << "O metodo da Secante falhou apos " << max_it << " iteracoes.\n";
    file.close();
}

void regulaFalsi(double a, double b, double eps, int maxIter) {
    double fa = f(a);
    double fb = f(b);
    double cp = b - a;
    int k = 0;
    ofstream file;
    file.open("regula_falsi.txt");
    file << "Iteracao\tRaiz\t\tErro\n";
     file << "--------------------------------------------------\n";
    double x, fx, fx_ant = 0.0;
    while (cp > eps && k < maxIter) {
        x = (a * fb - b * fa) / (fb - fa);
        fx = f(x);
        file << k << "\t\t" << x << "\t\t" <<  cp << endl;
        if (fa * fx < 0) {
            b = x;
            fb = fx;
        } else {
            a = x;
            fa = fx;
        }
        cp = abs(x - fx_ant);
        fx_ant = x;
        k++;
    }
    file.close();
    cout << "Raiz: " << x << endl;
    cout << "N de Iteracao: " << k << endl;
}

int main() {
    int num;
    double a, b,eps;
    int maxIter;

    do{
        cout << "Escolha o metodo:" << endl;
        cout << "1. Bisseccao" << endl;
        cout << "2. MIL" << endl;
        cout << "3. Newton" << endl;
        cout << "4. Secante" << endl;
        cout << "5. Regula Falsi" << endl;
        cout << "6. Sair" << endl;
        cin >> num;
        if(num != 6){
            cout << "Digite o numero maximo de iteracoes: ";
            cin >> maxIter;
            cout << "Digite a precisao: ";
            cin >> eps;
        }
        switch (num) {
            case 1:
                cout << "Digite o intervalo [a, b]: ";
                cin >> a >> b;
                bisseccao(a, b, eps, maxIter);
            break;
            case 2:
                cout << "Digite o intervalo [a, b]: ";
                cin >> a >> b;
                mil(a, b, eps, maxIter);
            break;
            case 3:
                cout << "Digite o valor inicial x0: ";
                cin >> a;
                newton(a, eps, maxIter);
            break;
            case 4:
                cout << "Digite o intervalo [a, b]: ";
                cin >> a >> b;
                secante(a, b, eps, maxIter);
            break;
            case 5:
                cout << "Digite o intervalo [a, b]: ";
                cin >> a >> b;
                regulaFalsi(a, b, eps, maxIter);
            break;
            case 6:
            break;
            default:
                cout << "Numero invalido!" << endl;
            break;
        }
    }while (num != 6);
}