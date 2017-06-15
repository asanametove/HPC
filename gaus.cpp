#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstring>
//#include <conio.h>
#include <time.h>
#include <omp.h>
//#include <windows.h>
#include <vector>
#include <math.h>
#include <iomanip>
#include <cstdlib>
 
using namespace std;

void Task4MatrixInverseParallel()
{
    int i=0,j=0,k=0;
    int n;
    int temp;
    double det = 1;

    const double EPS = 1E-9;
    int timein, timeout, timeres = 0;
    ifstream fin;
    fin.open("input.txt");
    fin >> n;
    cout << "N = " << n << endl;

    //n=2;
    srand(1);
    vector <vector<double>> Matrix (n, vector<double> (n));
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            fin >> Matrix[i][j];
	    //cout << Matrix[i][j] << "\t";
        }
	//cout << "\n";
    }

    //Matrix[0][0] = 1;// = -1
    //Matrix[0][1] = 1;
    //Matrix[1][0] = 1;
    //Matrix[1][1] = 0;

    vector <vector<double>> E (n, vector<double> (n));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
        {
            E[i][j] = 0.0;

            if (i == j)
                E[i][j] = 1.0;
        }

    timein = time(NULL);
    //cout << "timein = " << timein << endl;
        // Проходим по строкам матрицы (назовём их исходными)
    // сверху вниз. На данном этапе происходит прямой ход
    // и исходная матрица превращается в верхнюю треугольную
    for (int k = 0; k < n; ++k)
    {
        // Если элемент на главной диагонали в исходной
        // строке - нуль, то ищем строку, где элемент
        // того же столбца не нулевой, и меняем строки
        // местами
        if (fabs(Matrix[k][k]) < 1e-8)
        {
            // Ключ, говорязий о том, что был произведён обмен строк
            bool changed = false;

            // Идём по строкам, расположенным ниже исходной
            for (int i = k + 1; i < n; ++i)
            {
                // Если нашли строку, где в том же столбце
                // имеется ненулевой элемент
                if (fabs(Matrix[i][k]) > 1e-8)
                {
                    // Меняем найденную и исходную строки местами
                    // как в исходной матрице, так и в единичной
                    swap(Matrix[k],   Matrix[i]);
                    swap(E[k], E[i]);

                    // Взводим ключ - сообщаем о произведённом обмене строк
                    changed = true;

                    break;
                }
            }

            // Если обмен строк произведён не был - матрица не может быть
            // обращена
            if (!changed)
            {
                cout<<"Matrix cannot inverse!!!"<<endl;
                // Сообщаем о неудаче обращения
            }
        }

        // Запоминаем делитель - диагональный элемент
        double div = Matrix[k][k];

        // Все элементы исходной строки делим на диагональный
        // элемент как в исходной матрице, так и в единичной
        for (int j = 0; j < n; ++j)
        {
            Matrix[k][j]   /= div;
            E[k][j] /= div;
        }

        // Идём по строкам, которые расположены ниже исходной
#pragma omp parallel for
        for (int i = k + 1; i < n; ++i)
        {
            // Запоминаем множитель - элемент очередной строки,
            // расположенный под диагональным элементом исходной
            // строки
            double multi = Matrix[i][k];

            // Отнимаем от очередной строки исходную, умноженную
            // на сохранённый ранее множитель как в исходной,
            // так и в единичной матрице
            for (int j = 0; j < n; ++j)
            {
                Matrix[i][j]   -= multi * Matrix[k][j];
                E[i][j] -= multi * E[k][j];
            }
        }
    }

    // Проходим по вернхней треугольной матрице, полученной
    // на прямом ходе, снизу вверх
    // На данном этапе происходит обратный ход, и из исходной
    // матрицы окончательно формируется единичная, а из единичной -
    // обратная
#pragma omp parallel for
    for (int k = n - 1; k > 0; --k)
    {
        // Идём по строкам, которые расположены выше исходной
        for (int i = k - 1; i + 1 > 0; --i)
        {
            // Запоминаем множитель - элемент очередной строки,
            // расположенный над диагональным элементом исходной
            // строки
            double multi = Matrix[i][k];

            // Отнимаем от очередной строки исходную, умноженную
            // на сохранённый ранее множитель как в исходной,
            // так и в единичной матрице
            for (int j = 0; j < n; ++j)
            {
                Matrix[i][j]   -= multi * Matrix[k][j];
                E[i][j] -= multi * E[k][j];
            }
        }
    }

    /*
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
            cout<<E[i][j] << " ";
        cout<<endl;
    }
    */

    timeout = time(NULL);
    timeres = timeout - timein;
//    cout<<E[n-1][n-1]<<endl;
    cout<<"time: "<<timeres<<endl;
}




int main()
{   
    for(int i = 1; i <= 8; i++) {
    	omp_set_num_threads(i);
	cout << "Threads: " << i << " ";
    	Task4MatrixInverseParallel();
    }

    system("pause");

    return 0;
}
