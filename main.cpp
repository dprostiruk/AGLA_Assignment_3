#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <cstdio>

#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"

using namespace std;

double v(double v0, double k0, double a1, double b1, double a2, double b2, double t) {
    return a2 / b2 + (v0 - a2 / b2) * cos(sqrt(a1 * a2) * t) - (k0 - a1 / b1) * ((sqrt(a2) * b1) / (b2 * sqrt(a1))) * sin(sqrt(a1 * a2) * t);
}

double k(double v0, double k0, double a1, double b1, double a2, double b2, double t) {
    return a1 / b1 + (k0 - a1 / b1) * cos(sqrt(a1 * a2) * t) + (v0 - a2 / b2) * ((sqrt(a1) * b2) / (b1 * sqrt(a2))) * sin(sqrt(a1 * a2) * t);
}

int main() {
#ifdef WIN32
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
#else
    FILE* pipe = popen(GNUPLOT_NAME, "w");
#endif

    double v0, k0;
    double a1, b1, a2, b2, T, N;
    cin >> v0 >> k0 >> a1 >> b1 >> a2 >> b2 >> T >> N;

    double t = 0;
    double h = T / N;

    cout << "t:" << endl;
    for (int i = 0; i <= N; i++) {
        cout << fixed << setprecision(2) << t << " ";
        t += h;
    }
    cout << endl;

    t = 0;
    cout << "v:" << endl;
    for (int i = 0; i <= N; i++) {
        cout << fixed << setprecision(2) << v(v0, k0, a1, b1, a2, b2, t) << " ";
        t += h;
    }
    cout << endl;

    t = 0;
    cout << "k:" << endl;
    for (int i = 0; i <= N; i++) {
        cout << fixed << setprecision(2) << k(v0, k0, a1, b1, a2, b2, t) << " ";
        t += h;
    }
    cout << endl;

    //v(t) & k(t)

    fprintf(pipe, "v0=%lf\nk0=%lf\na1=%lf\nb1=%lf\na2=%lf\nb2=%lf\n", v0, k0, a1, b1, a2, b2);
    fprintf(pipe, "v(t) = a2 / b2 + (v0 - a2 / b2) * cos(sqrt(a1 * a2) * t) - (k0 - a1 / b1) * ((sqrt(a2) * b1) / (b2 * sqrt(a1))) * sin(sqrt(a1 * a2) * t)\n");
    fprintf(pipe, "k(t) = a1 / b1 + (k0 - a1 / b1) * cos(sqrt(a1 * a2) * t) + (v0 - a2 / b2) * ((sqrt(a1) * b2) / (b1 * sqrt(a2))) * sin(sqrt(a1 * a2) * t)\n");
    fprintf(pipe, "plot [0 : 200] [0 : 20] v(x) title \"v(t)\", k(x) title \"k(t)\"\n");

    //v(k)

//    fprintf(pipe, "plot [0 : 200] [0 : 20] '-' with lines\n");
//    for (int i = 1; i < N; ++i) {
//        fprintf(pipe, "%f %f\n", v(v0, k0, a1, b1, a2, b2, i * h), k(v0, k0, a1, b1, a2, b2, i * h));
//    }
//    _pclose(pipe);


    return 0;
}
