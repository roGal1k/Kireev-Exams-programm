#include "widget.h"
#include "ui_widget.h"
// #include <iostream>
#include <boost/algorithm/algorithm.hpp>
#include <cmath>

#define PI 3.14159265358979323846
using namespace std;

//объявляем вспомогательные функции
double check_ft(double ae, double fi_c);             // получение значения вспомогательной функции f_t(V/V_C, phi_c)
double check_fw(double ae, double fi_c, double nu);  // получение значения вспомогательной функции f_w(ae, phi_c, nu)
double comp_a2(double R, double fi_c2, double alpha);    // получение ширины просмотра полосы в режиме поиска по "наклонным строкам"
double comp_a1(double R, double fi_c1, double alpha);    // получение ширины просмотра полосы в режиме перехода ЛА с одной "наклонной строки" на другую

///расчет по таблицам сформированных А.П.Гришневым

//Формирование работы консоли и вывод результатов выполнения
void Widget::table(double &W_graph, double &T_graph)
{
    double v_c = 956; // Воздушная цель летит прямолинейно с известной скоростью [км/ч]
    double b = 300; // Зона неопределенности представляет собой подвижный прямоугольник со стороной [км]
    double l = 200; // И второй стороной [км]

    // Летный данные истребителя-перехватчика и его поисковые возможности характеризуются
    //	крейсерской скоростью полета
    double v = 1100; // [км/ч]
    // радиусом обнаружения цели
    double R = 60; // [км]
    // и полураствором сектора обнаружения цели
    double fi_c = 60; // [градусов]
    // перевод из градусов в радианы
    double fi_c_rad = fi_c * PI/ 180;

    double T = ui->lineEdit_2->text().toDouble()/60;
    ui->lineEdit_2->clear();

   // std::cout << "Исходные данные:" << std::endl
   // << "Воздушная цель летит прямолинейно с известной скоростью " << v_c << " км/ч" << std::endl
   // << "Зона неопределенности представляет собой подвижный прямоугольник со сторонами " << b << " на " << l << " км" << std::endl
   // << "Летный данные истребителя-перехватчика и его поисковые возможности характеризуются " << std::endl
   // << "крейсерской скоростью полета " << v << " км/ч" << std::endl
   // << "и полураствором сектора обнаружения цели " << fi_c << " градусов" << std::endl
   // << std::endl << "Необходимо при исходных данных этой задачи определить вероятность обнаружения цели" << std::endl
   // << "отведенное время Т " << T*60 << " мин при оптимальном способе поиска, отвечающем этому времени" << std::endl;

    // Определяем величину каппа
    double kappa = v / v_c;
    // Округляем каппа до двух знаков после запятой
    kappa = round(kappa * 100.0) / 100.0;
    // По величине каппа и скорости цели находим f_T
    double f_T = check_ft(kappa, fi_c);
    // Тогда гарантированное время обнаружения цели будет
    double T_g = ((b*l) / (v*R)) *(f_T + (R / b)*((v / (v + v_c)) - 2 * sin(fi_c_rad)*f_T) - (((R*R) / (b*l*sqrt(1 - pow(v_c / v, 2))))*((v / ((v + v_c)*f_T)) - (2 * sin(fi_c_rad)))));
    // Вероятность обнаружения цели методом "строка", за найденное гарантированное время обнаружения цели, составляет 1.


    // Если время Т, отводимое на поиск дели, задано, ограничено и меньше гарантированного,
    // то ориентация поиска на гарантированное обнаружение дели неправомерна.
    // В этом случае мы должны найти оптимальный способ поиска цели, максимизирующий значение вероятности обнаружения цели
    // и определить получающееся при этом значение Для упрощения решения задачи ограничим класс рассматриваемых способов поиска
    // просмотром зоны поиска по "наклонным строкам".

    // Найдем вероятность обнаружения цели за 40 минут

    //возьмем значения (l^*)/(v*T)
    double l_r[5] = {0.5, 0.2, 0.1, 0.05, 0.02};

    double l1 = 0;
    double l2 = 0;
    double u11 = 0;
    double u22 = 0;
    int i = 0;
    double mu_10 = 0;
    double s11_1 = 0;
    double s11_2 = 0;
    double s15 = 0;
    double s14;

    // Найдем интервал пересечения кривой u(x) с прямой u = x
    while (!((u11 < l1) && (u22 > l2)))
    {
        l1 = l_r[i];
        s11_1 = check_fw(kappa, fi_c, 1 / l1);
        mu_10 = (v_c / v) - l1;
        mu_10 = sqrt(1 - pow(mu_10, 2));
        s14 = (v*T*mu_10) / (b - 2 * R*sin(fi_c_rad));
        s15 = int(s14) / s14;
        u11 = l / (v*T) - (R / (b - 2 * R*sin(fi_c_rad))) * s11_1*s15;

        l2 = l_r[i + 1];
        s11_2 = check_fw(kappa, fi_c, 1 / l2);
        mu_10 = (v_c / v) - l2;
        mu_10 = sqrt(1 - pow(mu_10, 2));
        s14 = (v*T*mu_10) / (b - 2 * R*sin(fi_c_rad));
        s15 = int(s14) / s14;
        u22 = l / (v*T) - (R / (b - 2 * R*sin(fi_c_rad))) * s11_2*s15;

        i++;
        if (i + 1 == 5) break;
    }

    // по интерполяционной формуле находим точку пересечения
    double xi = (l1 - u11) / (u22 + l1 - u11 - l2);
    // откуда nu будет
    double nu = 1 / (l1 - l2 * xi);

    // дальше ищем вероятность обнаружения цели
    double s7 = s11_1 + xi*(s11_2 - s11_1);
    double s10 = sqrt(1-pow((v_c / v) - 1 / nu,2));
    int s12 = int(v * T * s10) / (b - 2 * R*sin(fi_c_rad));
    double s151 = (2 * sin(fi_c_rad)*(1 + v_c / v) - s7 )*(R / (v*T*(1 + v_c / v))*(s12*s7 / s10));
    double s16 = s7 + s151;
    double W_T = (v*R*T) / (b*l) * s16;

   // std::cout << std::endl << "В результате решения было получено: уменьшение времени поиска от " << T_g << " ч "
   //     << "до " << T*60 << " мин " << std::endl << "приводит при оптимальном способе поиска, ориентированном "
   //     << "на отведенное время, к уменьшению" << std::endl << "вероятности обнаружения от W(T*) = 1 до W(T) = " << W_T << std::endl;
    W_graph = W_T;
    T_graph = T;
}

double comp_a1(double R, double fi_c1, double alpha)
{
    double a1;
    if (fi_c1 <= alpha)
        a1 = R * sin(alpha + fi_c1);
    else if ((alpha < fi_c1) and (fi_c1 < ((PI / 2) - alpha)))
        a1 = 2 * R * sin(fi_c1) * cos(alpha);
    else if ((PI / 2 - alpha) < fi_c1 and fi_c1 <= (PI / 2 + alpha))
        a1 = R * (1 + sin(fi_c1 - alpha));
    else if (fi_c1 > (PI / 2 + alpha))
        a1 = 2 * R;
    return a1;
}

double comp_a2(double R, double fi_c2, double alpha)
{
    double a2;
    if (fi_c2 <= (PI / 2 - alpha))
        a2 = R * sin(fi_c2 + alpha);
    else if ((PI / 2 - alpha) < fi_c2 and fi_c2 <= alpha)
        a2 = R;
    else if (alpha < fi_c2 and fi_c2 <= (PI / 2 + alpha))
        a2 = R * (1 + sin(fi_c2 - alpha));
    else if (fi_c2 > (PI / 2 + alpha))
        a2 = 2 * R;
    return a2;
}

double check_fw(double ae, double fi_c, double nu)
{
    if (ae == 1.05)
    {
        if (nu == 2)
        {
            if (fi_c >= 30 && fi_c < 40) return 1.020;
            else if (fi_c >= 40 && fi_c < 50) return 1.022;
            else if (fi_c >= 50 && fi_c < 60) return 1.022;
            else if (fi_c >= 60 && fi_c < 70) return 1.091;
            else if (fi_c >= 70 && fi_c < 80) return 1.267;
            else if (fi_c >= 80 && fi_c < 90) return 1.435;
            else if (fi_c >= 90 && fi_c < 100) return 1.591;
            else if (fi_c >= 100 && fi_c < 120) return 1.730;
            else if (fi_c >= 120 && fi_c < 140) return 1.940;
            else if (fi_c >= 140 && fi_c < 160) return 2.030;
            else if (fi_c >= 160 && fi_c < 180) return 2.045;
            else if (fi_c > 180) return 2.045;
        }
        else if (nu == 5)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.688;
            else if (fi_c >= 40 && fi_c < 50) return 0.688;
            else if (fi_c >= 50 && fi_c < 60) return 0.688;
            else if (fi_c >= 60 && fi_c < 70) return 0.688;
            else if (fi_c >= 70 && fi_c < 80) return 0.740;
            else if (fi_c >= 80 && fi_c < 90) return 0.859;
            else if (fi_c >= 90 && fi_c < 100) return 0.972;
            else if (fi_c >= 100 && fi_c < 120) return 1.076;
            else if (fi_c >= 120 && fi_c < 140) return 1.248;
            else if (fi_c >= 140 && fi_c < 160) return 1.351;
            else if (fi_c >= 160 && fi_c < 180) return 1.377;
            else if (fi_c > 180) return 1.377;
        }
        else if (nu == 10)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.532;
            else if (fi_c >= 40 && fi_c < 50) return 0.532;
            else if (fi_c >= 50 && fi_c < 60) return 0.532;
            else if (fi_c >= 60 && fi_c < 70) return 0.532;
            else if (fi_c >= 70 && fi_c < 80) return 0.539;
            else if (fi_c >= 80 && fi_c < 90) return 0.631;
            else if (fi_c >= 90 && fi_c < 100) return 0.721;
            else if (fi_c >= 100 && fi_c < 120) return 0.804;
            else if (fi_c >= 120 && fi_c < 140) return 0.944;
            else if (fi_c >= 140 && fi_c < 160) return 1.055;
            else if (fi_c >= 160 && fi_c < 180) return 1.065;
            else if (fi_c > 180) return 1.065;
        }
        else if (nu == 20)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.434;
            else if (fi_c >= 40 && fi_c < 50) return 0.434;
            else if (fi_c >= 50 && fi_c < 60) return 0.434;
            else if (fi_c >= 60 && fi_c < 70) return 0.434;
            else if (fi_c >= 70 && fi_c < 80) return 0.434;
            else if (fi_c >= 80 && fi_c < 90) return 0.501;
            else if (fi_c >= 90 && fi_c < 100) return 0.574;
            else if (fi_c >= 100 && fi_c < 120) return 0.643;
            else if (fi_c >= 120 && fi_c < 140) return 0.761;
            else if (fi_c >= 140 && fi_c < 160) return 0.839;
            else if (fi_c >= 160 && fi_c < 180) return 0.888;
            else if (fi_c > 180) return 0.888;
        }
        else if (nu == 50)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.362;
            else if (fi_c >= 40 && fi_c < 50) return 0.362;
            else if (fi_c >= 50 && fi_c < 60) return 0.362;
            else if (fi_c >= 60 && fi_c < 70) return 0.362;
            else if (fi_c >= 70 && fi_c < 80) return 0.362;
            else if (fi_c >= 80 && fi_c < 90) return 0.413;
            else if (fi_c >= 90 && fi_c < 100) return 0.574;
            else if (fi_c >= 100 && fi_c < 120) return 0.643;
            else if (fi_c >= 120 && fi_c < 140) return 0.761;
            else if (fi_c >= 140 && fi_c < 160) return 0.839;
            else if (fi_c >= 160 && fi_c < 180) return 0.868;
            else if (fi_c > 180) return 0.868;
        }
    }
    else if (ae == 1.10)
    {
        if (nu == 2)
        {
            if (fi_c >= 30 && fi_c < 40) return 1.032;
            else if (fi_c >= 40 && fi_c < 50) return 1.040;
            else if (fi_c >= 50 && fi_c < 60) return 1.040;
            else if (fi_c >= 60 && fi_c < 70) return 1.170;
            else if (fi_c >= 70 && fi_c < 80) return 1.347;
            else if (fi_c >= 80 && fi_c < 90) return 1.515;
            else if (fi_c >= 90 && fi_c < 100) return 1.669;
            else if (fi_c >= 100 && fi_c < 120) return 1.803;
            else if (fi_c >= 120 && fi_c < 140) return 1.999;
            else if (fi_c >= 140 && fi_c < 160) return 2.060;
            else if (fi_c >= 160 && fi_c < 180) return 2.081;
            else if (fi_c > 180) return 2.081;
        }
        else if (nu == 5)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.733;
            else if (fi_c >= 40 && fi_c < 50) return 0.733;
            else if (fi_c >= 50 && fi_c < 60) return 0.733;
            else if (fi_c >= 60 && fi_c < 70) return 0.733;
            else if (fi_c >= 70 && fi_c < 80) return 0.848;
            else if (fi_c >= 80 && fi_c < 90) return 0.972;
            else if (fi_c >= 90 && fi_c < 100) return 1.068;
            else if (fi_c >= 100 && fi_c < 120) return 1.192;
            else if (fi_c >= 120 && fi_c < 140) return 1.361;
            else if (fi_c >= 140 && fi_c < 160) return 1.432;
            else if (fi_c >= 160 && fi_c < 180) return 1.486;
            else if (fi_c > 180) return 1.486;
        }
        else if (nu == 10)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.596;
            else if (fi_c >= 40 && fi_c < 50) return 0.596;
            else if (fi_c >= 50 && fi_c < 60) return 0.596;
            else if (fi_c >= 60 && fi_c < 70) return 0.596;
            else if (fi_c >= 70 && fi_c < 80) return 0.662;
            else if (fi_c >= 80 && fi_c < 90) return 0.784;
            else if (fi_c >= 90 && fi_c < 100) return 0.861;
            else if (fi_c >= 100 && fi_c < 120) return 0.949;
            else if (fi_c >= 120 && fi_c < 140) return 1.092;
            else if (fi_c >= 140 && fi_c < 160) return 1.175;
            else if (fi_c >= 160 && fi_c < 180) return 1.192;
            else if (fi_c > 180) return 1.192;
        }
        else if (nu == 20)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.514;
            else if (fi_c >= 40 && fi_c < 50) return 0.514;
            else if (fi_c >= 50 && fi_c < 60) return 0.514;
            else if (fi_c >= 60 && fi_c < 70) return 0.514;
            else if (fi_c >= 70 && fi_c < 80) return 0.561;
            else if (fi_c >= 80 && fi_c < 90) return 0.649;
            else if (fi_c >= 90 && fi_c < 100) return 0.661;
            else if (fi_c >= 100 && fi_c < 120) return 0.949;
            else if (fi_c >= 120 && fi_c < 140) return 1.092;
            else if (fi_c >= 140 && fi_c < 160) return 1.175;
            else if (fi_c >= 160 && fi_c < 180) return 1.192;
            else if (fi_c > 180) return 1.192;
        }
        else if (nu == 50)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.458;
            else if (fi_c >= 40 && fi_c < 50) return 0.458;
            else if (fi_c >= 50 && fi_c < 60) return 0.458;
            else if (fi_c >= 60 && fi_c < 70) return 0.458;
            else if (fi_c >= 70 && fi_c < 80) return 0.496;
            else if (fi_c >= 80 && fi_c < 90) return 0.575;
            else if (fi_c >= 90 && fi_c < 100) return 0.650;
            else if (fi_c >= 100 && fi_c < 120) return 0.719;
            else if (fi_c >= 120 && fi_c < 140) return 0.832;
            else if (fi_c >= 140 && fi_c < 160) return 0.900;
            else if (fi_c >= 160 && fi_c < 180) return 0.916;
            else if (fi_c > 180) return 0.916;
        }
    }
    else if (ae == 1.15)
    {
        if (nu == 2)
        {
            if (fi_c >= 30 && fi_c <40) return 1.039;
            else if (fi_c >= 40 && fi_c < 50) return 1.055;
            else if (fi_c >= 50 && fi_c < 60) return 1.055;
            else if (fi_c >= 60 && fi_c < 70) return 1.239;
            else if (fi_c >= 70 && fi_c < 80) return 1.417;
            else if (fi_c >= 80 && fi_c < 90) return 1.563;
            else if (fi_c >= 90 && fi_c < 100) return 1.734;
            else if (fi_c >= 100 && fi_c < 120) return 1.864;
            else if (fi_c >= 120 && fi_c < 140) return 2.047;
            else if (fi_c >= 140 && fi_c < 160) return 2.110;
            else if (fi_c >= 160 && fi_c < 180) return 2.110;
            else if (fi_c > 180) return 2.110;
        }
        else if (nu == 5)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.768;
            else if (fi_c >= 40 && fi_c < 50) return 0.769;
            else if (fi_c >= 50 && fi_c < 60) return 0.769;
            else if (fi_c >= 60 && fi_c < 70) return 0.808;
            else if (fi_c >= 70 && fi_c < 80) return 0.941;
            else if (fi_c >= 80 && fi_c < 90) return 1.068;
            else if (fi_c >= 90 && fi_c < 100) return 1.187;
            else if (fi_c >= 100 && fi_c < 120) return 1.293;
            else if (fi_c >= 120 && fi_c < 140) return 1.454;
            else if (fi_c >= 140 && fi_c < 160) return 1.523;
            else if (fi_c >= 160 && fi_c < 180) return 1.536;
            else if (fi_c > 180) return 1.536;
        }
        else if (nu == 10)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.646;
            else if (fi_c >= 40 && fi_c < 50) return 0.646;
            else if (fi_c >= 50 && fi_c < 60) return 0.646;
            else if (fi_c >= 60 && fi_c < 70) return 0.655;
            else if (fi_c >= 70 && fi_c < 80) return 0.767;
            else if (fi_c >= 80 && fi_c < 90) return 0.867;
            else if (fi_c >= 90 && fi_c < 100) return 0.877;
            else if (fi_c >= 100 && fi_c < 120) return 1.069;
            else if (fi_c >= 120 && fi_c < 140) return 1.210;
            else if (fi_c >= 140 && fi_c < 160) return 1.284;
            else if (fi_c >= 160 && fi_c < 180) return 1.293;
            else if (fi_c > 180) return 1.293;
        }
        else if (nu == 20)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.575;
            else if (fi_c >= 40 && fi_c < 50) return 0.575;
            else if (fi_c >= 50 && fi_c < 60) return 0.575;
            else if (fi_c >= 60 && fi_c < 70) return 0.575;
            else if (fi_c >= 70 && fi_c < 80) return 0.645;
            else if (fi_c >= 80 && fi_c < 90) return 0.772;
            else if (fi_c >= 90 && fi_c < 100) return 0.882;
            else if (fi_c >= 100 && fi_c < 120) return 0.943;
            else if (fi_c >= 120 && fi_c < 140) return 1.073;
            else if (fi_c >= 140 && fi_c < 160) return 1.142;
            else if (fi_c >= 160 && fi_c < 180) return 1.150;
            else if (fi_c > 180) return 1.150;
        }
        else if (nu == 50)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.528;
            else if (fi_c >= 40 && fi_c < 50) return 0.528;
            else if (fi_c >= 50 && fi_c < 60) return 0.528;
            else if (fi_c >= 60 && fi_c < 70) return 0.528;
            else if (fi_c >= 70 && fi_c < 80) return 0.618;
            else if (fi_c >= 80 && fi_c < 90) return 0.705;
            else if (fi_c >= 90 && fi_c < 100) return 0.739;
            else if (fi_c >= 100 && fi_c < 120) return 0.805;
            else if (fi_c >= 120 && fi_c < 140) return 0.983;
            else if (fi_c >= 140 && fi_c < 160) return 1.047;
            else if (fi_c >= 160 && fi_c < 180) return 1.056;
            else if (fi_c > 180) return 1.056;
        }
    }
    //////////////////////////
    else if (ae == 1.20)
    {
        if (nu == 2)
        {
            if (fi_c >= 30 && fi_c < 40) return 1.042;
            else if (fi_c >= 40 && fi_c < 50) return 1.066;
            else if (fi_c >= 50 && fi_c < 60) return 1.115;
            else if (fi_c >= 60 && fi_c < 70) return 1.300;
            else if (fi_c >= 70 && fi_c < 80) return 1.477;
            else if (fi_c >= 80 && fi_c < 90) return 1.642;
            else if (fi_c >= 90 && fi_c < 100) return 1.789;
            else if (fi_c >= 100 && fi_c < 120) return 1.915;
            else if (fi_c >= 120 && fi_c < 140) return 2.068;
            else if (fi_c >= 140 && fi_c < 160) return 2.134;
            else if (fi_c >= 160 && fi_c < 180) return 2.134;
            else if (fi_c > 180) return 2.134;
        }
        else if (nu == 5)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.795;
            else if (fi_c >= 40 && fi_c < 50) return 0.795;
            else if (fi_c >= 50 && fi_c < 60) return 0.795;
            else if (fi_c >= 60 && fi_c < 70) return 0.868;
            else if (fi_c >= 70 && fi_c < 80) return 1.022;
            else if (fi_c >= 80 && fi_c < 90) return 1.132;
            else if (fi_c >= 90 && fi_c < 100) return 1.272;
            else if (fi_c >= 100 && fi_c < 120) return 1.376;
            else if (fi_c >= 120 && fi_c < 140) return 1.531;
            else if (fi_c >= 140 && fi_c < 160) return 1.637;
            else if (fi_c >= 160 && fi_c < 180) return 1.699;
            else if (fi_c > 180) return 1.699;
        }
        else if (nu == 10)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.337;
            else if (fi_c >= 40 && fi_c < 50) return 0.337;
            else if (fi_c >= 50 && fi_c < 60) return 0.687;
            else if (fi_c >= 60 && fi_c < 70) return 0.741;
            else if (fi_c >= 70 && fi_c < 80) return 0.859;
            else if (fi_c >= 80 && fi_c < 90) return 0.972;
            else if (fi_c >= 90 && fi_c < 100) return 1.076;
            else if (fi_c >= 100 && fi_c < 120) return 1.168;
            else if (fi_c >= 120 && fi_c < 140) return 1.304;
            else if (fi_c >= 140 && fi_c < 160) return 1.371;
            else if (fi_c >= 160 && fi_c < 180) return 1.374;
            else if (fi_c > 180) return 1.374;
        }
        else if (nu == 20)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.622;
            else if (fi_c >= 40 && fi_c < 50) return 0.624;
            else if (fi_c >= 50 && fi_c < 60) return 0.624;
            else if (fi_c >= 60 && fi_c < 70) return 0.663;
            else if (fi_c >= 70 && fi_c < 80) return 0.773;
            else if (fi_c >= 80 && fi_c < 90) return 0.878;
            else if (fi_c >= 90 && fi_c < 100) return 0.971;
            else if (fi_c >= 100 && fi_c < 120) return 1.038;
            else if (fi_c >= 120 && fi_c < 140) return 1.183;
            else if (fi_c >= 140 && fi_c < 160) return 1.344;
            else if (fi_c >= 160 && fi_c < 180) return 1.347;
            else if (fi_c > 180) return 1.347;
        }
        else if (nu == 50)
        {
            if (fi_c >= 30 && fi_c < 40) return 0.551;
            else if (fi_c >= 40 && fi_c < 50) return 0.552;
            else if (fi_c >= 50 && fi_c < 60) return 0.562;
            else if (fi_c >= 60 && fi_c < 70) return 0.619;
            else if (fi_c >= 70 && fi_c < 80) return 0.719;
            else if (fi_c >= 80 && fi_c < 90) return 0.815;
            else if (fi_c >= 90 && fi_c < 100) return 0.994;
            else if (fi_c >= 100 && fi_c < 120) return 0.994;
            else if (fi_c >= 120 && fi_c < 140) return 1.115;
            else if (fi_c >= 140 && fi_c < 160) return 1.294;
            else if (fi_c >= 160 && fi_c < 180) return 1.329;
            else if (fi_c > 180) return 1.329;
        }
    }
}

double check_ft(double ae, double fi_c)
{
    if (ae > 0 && ae <= 1.05)
    {
        if (fi_c >= 30 && fi_c <= 70) return 3.28;
        else if (fi_c == 80) return 2.89;
        else if (fi_c == 90) return 2.51;
        else if (fi_c == 100) return 2.24;
        else if (fi_c == 120) return 1.885;
        else if (fi_c == 140) return 1.703;
        else if (fi_c == 160) return 1.640;
        else if (fi_c == 180) return 1.640;
    }

    else if (ae > 1.05 && ae <= 1.10)
    {
        if (fi_c >= 30 && fi_c <= 60) return 2.40;
        else if (fi_c == 70) return 2.22;
        else if (fi_c == 80) return 1.917;
        else if (fi_c == 90) return 1.694;
        else if (fi_c == 100) return 1.531;
        else if (fi_c == 120) return 1.322;
        else if (fi_c == 140) return 1.222;
        else if (fi_c == 160) return 1.200;
        else if (fi_c == 180) return 1.200;
    }

    else if (ae > 1.10 && ae <= 1.15)
    {
        if (fi_c >= 30 && fi_c <= 60) return 2.02;
        else if (fi_c == 70) return 1.735;
        else if (fi_c == 80) return 1.516;
        else if (fi_c == 90) return 1.356;
        else if (fi_c == 100) return 1.237;
        else if (fi_c == 120) return 1.087;
        else if (fi_c == 140) return 1.021;
        else if (fi_c == 160) return 1.012;
        else if (fi_c == 180) return 1.012;
    }

    else if (ae > 1.15 && ae <= 1.20)
    {
        if (fi_c == 30) return 1.813;
        else if (fi_c == 40) return 1.809;
        else if (fi_c == 50) return 1.809;
        else if (fi_c == 60) return 1.703;
        else if (fi_c == 70) return 1.465;
        else if (fi_c == 80) return 1.292;
        else if (fi_c == 90) return 1.165;
        else if (fi_c == 100) return 1.071;
        else if (fi_c == 120) return 0.954;
        else if (fi_c == 140) return 0.907;
        else if (fi_c == 160) return 0.904;
        else if (fi_c == 180) return 0.904;
    }

    else if (ae > 1.20 && ae <= 1.30)
    {
        if (fi_c == 30) return 1.588;
        else if (fi_c == 40) return 1.565;
        else if (fi_c == 50) return 1.565;
        else if (fi_c == 60) return 1.339;
        else if (fi_c == 70) return 1.170;
        else if (fi_c == 80) return 1.046;
        else if (fi_c == 90) return 0.955;
        else if (fi_c == 100) return 0.888;
        else if (fi_c == 120) return 0.807;
        else if (fi_c == 140) return 0.782;
        else if (fi_c == 160) return 0.782;
        else if (fi_c == 180) return 0.782;
    }

    else if (ae > 1.30 && ae <= 1.50)
    {
        if (fi_c == 30) return 1.412;
        else if (fi_c == 40) return 1.355;
        else if (fi_c == 50) return 1.174;
        else if (fi_c == 60) return 1.022;
        else if (fi_c == 70) return 0.911;
        else if (fi_c == 80) return 0.829;
        else if (fi_c == 90) return 0.769;
        else if (fi_c == 100) return 0.725;
        else if (fi_c == 120) return 0.678;
        else if (fi_c == 140) return 0.671;
        else if (fi_c == 160) return 0.671;
        else if (fi_c == 180) return 0.671;
    }

    else if (ae > 1.50 && ae <= 2.00)
    {
        if (fi_c == 30) return 1.333;
        else if (fi_c == 40) return 1.037;
        else if (fi_c == 50) return 0.870;
        else if (fi_c == 60) return 0.770;
        else if (fi_c == 70) return 0.703;
        else if (fi_c == 80) return 0.654;
        else if (fi_c == 90) return 0.619;
        else if (fi_c == 100) return 0.595;
        else if (fi_c == 120) return 0.577;
        else if (fi_c == 140) return 0.577;
        else if (fi_c == 160) return 0.577;
        else if (fi_c == 180) return 0.577;
    }

    else if (ae > 2.00 && ae <= 5.00)
    {
        if (fi_c == 30) return 1.042;
        else if (fi_c == 40) return 0.810;
        else if (fi_c == 50) return 0.680;
        else if (fi_c == 60) return 0.801;
        else if (fi_c == 70) return 0.544;
        else if (fi_c == 80) return 0.529;
        else if (fi_c == 90) return 0.516;
        else if (fi_c == 100) return 0.510;
        else if (fi_c == 120) return 0.510;
        else if (fi_c == 140) return 0.510;
        else if (fi_c == 160) return 0.510;
        else if (fi_c == 180) return 0.510;
    }

    else if (ae > 5.00 && ae <= 10.00)
    {
        if (fi_c == 30) return 1.010;
        else if (fi_c == 40) return 0.786;
        else if (fi_c == 50) return 0.659;
        else if (fi_c == 60) return 0.583;
        else if (fi_c == 70) return 0.537;
        else if (fi_c == 80) return 0.513;
        else if (fi_c == 90) return 0.504;
        else if (fi_c == 100) return 0.502;
        else if (fi_c == 120) return 0.502;
        else if (fi_c == 140) return 0.502;
        else if (fi_c == 160) return 0.502;
        else if (fi_c == 180) return 0.502;
    }

    else if (ae > 10.00 && ae <= 20.00)
    {
        if (fi_c == 30) return 1.002;
        else if (fi_c == 40) return 0.780;
        else if (fi_c == 50) return 0.634;
        else if (fi_c == 60) return 0.579;
        else if (fi_c == 70) return 0.533;
        else if (fi_c == 80) return 0.509;
        else if (fi_c == 90) return 0.501;
        else if (fi_c == 100) return 0.501;
        else if (fi_c == 120) return 0.501;
        else if (fi_c == 140) return 0.501;
        else if (fi_c == 160) return 0.501;
        else if (fi_c == 180) return 0.501;
    }

    else if (ae > 20.00 && ae <= 50.00)
    {
        if (fi_c == 30) return 1.000;
        else if (fi_c == 40) return 0.778;
        else if (fi_c == 50) return 0.653;
        else if (fi_c == 60) return 0.578;
        else if (fi_c == 70) return 0.532;
        else if (fi_c == 80) return 0.508;
        else if (fi_c == 90) return 0.500;
        else if (fi_c == 100) return 0.500;
        else if (fi_c == 120) return 0.500;
        else if (fi_c == 140) return 0.500;
        else if (fi_c == 160) return 0.500;
        else if (fi_c == 180) return 0.500;
    }

    else if (ae > 50.00 && ae <= 100.00)
    {
        if (fi_c == 30) return 1.000;
        else if (fi_c == 40) return 0.778;
        else if (fi_c == 50) return 0.653;
        else if (fi_c == 60) return 0.577;
        else if (fi_c == 70) return 0.532;
        else if (fi_c == 80) return 0.508;
        else if (fi_c == 90) return 0.500;
        else if (fi_c == 100) return 0.500;
        else if (fi_c == 120) return 0.500;
        else if (fi_c == 140) return 0.500;
        else if (fi_c == 160) return 0.500;
        else if (fi_c == 180) return 0.500;
    }
    else return 0.0;
}


///Расчет согласно формулам

double f_t(double R, double a1, double V_c, double V)
{
    double f_t = (R/(a1 * sqrt(1-pow((V_c/V),2))));
    return f_t;
}

double f_w(double a_1O,double u_1O,double V, double R)
{
    double f_w = ((a_1O *u_1O)/(V*R));
    return f_w;
}

double optimum(double &a_1O,double &u_1O, double V_c, double V)
{
    double T,T_2,l;
    double l_new = (l-(V+V_c)*T_2);
    double nu_1O = acos((( V_c/V) - (l_new/V*T)));
    u_1O = (V*sqrt(1+pow((V_c/V),2) -2*(V_c/V)*cos(nu_1O)));
    return nu_1O;
}

void calculate()
{
    double v_c = 956; // Воздушная цель летит прямолинейно с известной скоростью [км/ч]
    double b = 300; // Зона неопределенности представляет собой подвижный прямоугольник со стороной [км]
    double l = 200; // И второй стороной [км]
    double v = 1100; // [км/ч]
    double R = 60; // [км]
    double fi_c = 60; // [градусов]
    double fi_c_rad = fi_c * PI/ 180;

    double T = 2.0 / 3.0;

  //  std::cout << "Исходные данные:" << std::endl
  //  << "Воздушная цель летит прямолинейно с известной скоростью " << v_c << " км/ч" << std::endl
  //  << "Зона неопределенности представляет собой подвижный прямоугольник со сторонами " << b << " на " << l << " км" << std::endl
  //  << "Летный данные истребителя-перехватчика и его поисковые возможности характеризуются " << std::endl
  //  << "крейсерской скоростью полета " << v << " км/ч" << std::endl
  //  << "и полураствором сектора обнаружения цели " << fi_c << " градусов" << std::endl
  //  << std::endl << "Необходимо при исходных данных этой задачи определить вероятность обнаружения цели" << std::endl
  //  << "отведенное время Т " << T*60 << " мин при оптимальном способе поиска, отвечающем этому времени" << std::endl;

    double ae= v / v_c;
    ae = round(ae * 100.0) / 100.0;
//    double f_T = check_ft(ae, fi_c);
    double alpha =ae;
    double a1 = comp_a1(R,fi_c_rad,alpha);
    double f_T = f_t(R,a1,v_c,v);
    double T_g = ((b*l) / (v*R)) *(f_T + (R / b)*((v / (v + v_c)) - 2 * sin(fi_c_rad)*f_T) - (((R*R) / (b*l*sqrt(1 - pow(v_c / v, 2))))*((v / ((v + v_c)*f_T)) - (2 * sin(fi_c_rad)))));


    double l_r[5] = {0.5, 0.2, 0.1, 0.05, 0.02};

    double l1 = 0;
    double l2 = 0;
    double u11 = 0;
    double u22 = 0;
    int i = 0;
    double mu_10 = 0;
    double s11_1 = 0;
    double s11_2 = 0;
    double s15 = 0;
    double s14;

    // Найдем интервал пересечения кривой u(x) с прямой u = x
    while (!((u11 < l1) && (u22 > l2)))
    {
        l1 = l_r[i];
        double a_1O,u_1O;
        mu_10=optimum(a_1O,u_1O, v_c, v);
        s11_1 = f_w(a_1O,u_1O,v, R);
        mu_10 = (v_c / v) - l1;
        mu_10 = sqrt(1 - pow(mu_10, 2));
        s14 = (v*T*mu_10) / (b - 2 * R*sin(fi_c_rad));
        s15 = int(s14) / s14;
        u11 = l / (v*T) - (R / (b - 2 * R*sin(fi_c_rad))) * s11_1*s15;

        l2 = l_r[i + 1];
        s11_2 = check_fw(ae, fi_c, 1 / l2);
        mu_10 = (v_c / v) - l2;
        mu_10 = sqrt(1 - pow(mu_10, 2));
        s14 = (v*T*mu_10) / (b - 2 * R*sin(fi_c_rad));
        s15 = int(s14) / s14;
        u22 = l / (v*T) - (R / (b - 2 * R*sin(fi_c_rad))) * s11_2*s15;

        i++;
        if (i + 1 == 5) break;
    }

    double xi = (l1 - u11) / (u22 + l1 - u11 - l2);
    double nu = 1 / (l1 - l2 * xi);

    // дальше ищем вероятность обнаружения цели
    double s7 = s11_1 + xi*(s11_2 - s11_1);
    double s10 = sqrt(1-pow((v_c / v) - 1 / nu,2));
    int s12 = int(v * T * s10) / (b - 2 * R*sin(fi_c_rad));
    double s151 = (2 * sin(fi_c_rad)*(1 + v_c / v) - s7 )*(R / (v*T*(1 + v_c / v))*(s12*s7 / s10));
    double s16 = s7 + s151;
    double W_T = (v*R*T) / (b*l) * s16;

   // std::cout << std::endl << "В результате решения было получено: уменьшение времени поиска от " << T_g << " ч "
   //     << "до " << T*60 << " мин " << std::endl << "приводит при оптимальном способе поиска, ориентированном "
   //     << "на отведенное время, к уменьшению" << std::endl << "вероятности обнаружения от W(T*) = 1 до W(T) = " << W_T << std::endl;

}

