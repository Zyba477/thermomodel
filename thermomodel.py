import matplotlib.pyplot as plt
import numpy as np

N1 = 20 # кол-во промежутков на пластину
N2 = N1
t_end = 7200 # окончание по времени

lambda1 = 12 # теплопроводности матереала А
lambda2 = 24 # теплопроводности матереала B
lambda3 = 12

ro1 = 1000 # плотность матереала А
ro2 = 7800 # плотность матереала B
ro3 = 1000

c1 = 46 # теплоемкость матереала А
c2 = 460 # теплоемкость матереала В
c3 = 46

tl = 40.06 # температура на границе х = 0
tr = 22 # температура на границе х = L 
t0 = 22 # температура начальная

l = 1.5 # толщина пластина

# общее число узов:
N = N1 + N2 + N2 + 1

# расчетный шаг сетки по пространственной координате:
h = l / (N - 1)

# определяем коэфициент температуропроводности:
a1 = lambda1/(ro1 * c1)
a2 = lambda2/(ro2 * c2)
a3 = lambda2/(ro3 * c3)
# определяем расчетный шаг сетки по времени:
tau = t_end / 100.0

# определяем поле температуры в начальный момент времени:
t = [None] * N  
for i in range(N ):
    t[i] = t0

# проводим интегрирование нестационарного уравнения теплопроводности
time = 0
alpha = [None] * N
beta = [None] * N
alpha1 = alpha.copy()
beta1 = beta.copy()
# alpha1 = 
while time < t_end:
    time = time + tau

    # определяем начальные прогоночные коэфициэнты на основе левого граничного условия
    alpha[0] = 0.0
    beta[0] =  tl
    
    # цикл с параметром для определения прогоночных коэфициентов по формуле 8 в первой части пластины
    for i in range(1, N1):
        # ai, bi, ci, fi коэфициент канонического представления СЛАУ с трехдиагональной матрицей
        ai = lambda1 / pow(h, 2)
        bi = 2.0 * lambda1 / pow(h, 2) + (ro1 * c1 / tau)
        ci = lambda1 / pow(h, 2)
        fi = (-ro1) * c1 * t[i] / tau

        # alfa[i], beta[i] – прогоночные коэффициенты
        alpha[i] = ai / (bi - ci * alpha[i-1])
        beta[i] = (ci * beta[i-1]-fi) / (bi -ci * alpha[i-1])

    #  определяем прогоночные коэффициенты на границе раздела двух
    # частей, используем соотношения (28)
    # print (n1)

    al1 =  (2.0 * a1 * a2 * tau * lambda2)
    al2 =  2.0 * a1 * a2 * tau 
    al3 = lambda2 + (lambda1  * (1-alpha[N1-1]))
    al4 = (a1 * lambda2) + (a2 *lambda1)
    alpha[N1] = al1 / ((al2 * al3) + (pow(h,2) * al4 ))
    
    bt1 = 2.0 * a1 * a2 * tau * lambda1 * beta[N1-1]
    bt2 = (a1 * lambda2) + (a2 *lambda1)
    bt3 = (pow(h,2) * bt2 * t[N1])
    bt4 = 2.0 * a1 * a2 * tau
    bt5 = lambda2 + (lambda1 * (1-alpha[N1-1]))
    beta[N1] = (bt1 +  bt3) / ((bt4 * bt5) + (pow(h,2) * bt2))

    # alpha[N1]=2.0*a1*a2*tau*lambda2/(2.0*a1*a2*tau*(lambda2+lambda1 
    # *(1-alpha[N1-1]))+pow(h,2)*(a1*lambda2+a2*lambda1))

 
    # beta[N1]=(2.0*a1*a2*tau*lambda1*beta[N1-1]+pow(h,2)*(a1*lambda2+a2 
    # *lambda1)*t[N1])/(2.0*a1*a2*tau*(lambda2+lambda1 
    # *(1-alpha[N1-1]))+pow(h,2)*(a1*lambda2+a2*lambda1))
    
    # цикл с параметром для определения прогоночных коэффициентов по
    # формуле (8) во второй части пластины
    for i in range(N1+1,N-1):
        # {ai, bi, ci, fi – коэффициенты канонического представления СЛАУ с трехдиагональной матрицей
        ai = lambda2/pow(h,2);
        bi = 2.0*lambda2/pow(h,2)+ro2*c2/tau;
        ci = lambda2/pow(h,2);
        fi = -ro2*c2*t[i]/tau;
        # alfa[i], beta[i] – прогоночные коэффициенты
        alpha[i] = ai/(bi-ci*alpha[i-1]);
        beta[i] = (ci*beta[i-1]-fi)/(bi-ci*alpha[i-1]);



    alpha1[0] = 0.0
    beta1[0] =  tl
    for i in range(1, N1):
        # ai, bi, ci, fi коэфициент канонического представления СЛАУ с трехдиагональной матрицей
        ai = lambda2 / pow(h, 2)
        bi = 2.0 * lambda2 / pow(h, 2) + (ro2 * c2 / tau)
        ci = lambda2 / pow(h, 2)
        fi = (-ro2) * c2 * t[i] / tau

        # alfa[i], beta[i] – прогоночные коэффициенты
        alpha1[i] = ai / (bi - ci * alpha1[i-1])
        beta1[i] = (ci * beta1[i-1]-fi) / (bi -ci * alpha1[i-1])


    al1 =  (2.0 * a2 * a3 * tau * lambda3)
    al2 =  2.0 * a2 * a3 * tau 
    al3 = lambda3 + (lambda2  * (1-alpha1[N1-1]))
    al4 = (a2 * lambda3) + (a3 *lambda2)
    alpha1[N1] = al1 / ((al2 * al3) + (pow(h,2) * al4 ))
    
    bt1 = 2.0 * a2 * a3 * tau * lambda2 * beta1[N1-1]
    bt2 = (a2 * lambda3) + (a3 *lambda2)
    bt3 = (pow(h,2) * bt2 * t[N1])
    bt4 = 2.0 * a2 * a3 * tau
    bt5 = lambda3 + (lambda2 * (1-alpha1[N1-1]))
    beta1[N1] = (bt1 +  bt3) / ((bt4 * bt5) + (pow(h,2) * bt2))

    for i in range(N1+1,N-1):
        # {ai, bi, ci, fi – коэффициенты канонического представления СЛАУ с трехдиагональной матрицей
        ai = lambda3 / pow(h,2);
        bi = 2.0 * lambda3 / pow(h,2) + ro3 * c3 / tau;
        ci = lambda3 / pow(h,2);
        fi = -ro3 * c3 * t[i] / tau;
        # alfa[i], beta[i] – прогоночные коэффициенты
        alpha[i] = ai/(bi-ci*alpha[i-1]);
        beta[i] = (ci*beta[i-1]-fi)/(bi-ci*alpha[i-1]);









    



    # определяем значение температуры на правой границе
    

   
    alpha = alpha[ : -1] + alpha1[int(len(alpha1)/2) : ]
    beta = beta[ : -1] + beta1[int(len(beta1)/2) : ]
    # используя соотношение (7) определяем неизвестное поле
    # температуры
    # t1 = t.copy()
    t[N-2] = tr
    for i in range(N-2,-1,-1):
        t[i] = alpha[i] * t[i + 1] + beta[i]

    


    # for i in range(N-2,-1,-1):
    #     t1[i] = alpha1[i] * t1[i + 1] + beta1[i]


# print(alpha)
# print(alpha1)
# if alpha == alpha1:
#     print("AAAA")
# if beta == beta1:
#     print("BBBBBBBBBBB")
for i in range(len(t)):
     print(" {:.8f} {:.5f}".format(h * (i), t[i]))


plt.plot(t)
plt.show()     
# print(len(t))