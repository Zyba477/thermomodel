
Const mf = 500;

Type 
  vector = array[1..mf] of real;

Var
  i, N, N1, N2 : integer;
  T, alfa, beta : vector;
  ai, bi, ci, fi : real;
  a1, lamda1, ro1, c1 : real;
  a2, lamda2, ro2, c2 : real;
  h, tau, t_end, time : real;
  T0, Tl, Tr, L : real;
Begin
  N1 := 20;
  N2 := N1;
  t_end := 160;

  lamda1 := 48;
  lamda2 := 384;

  ro1 := 7800;
  ro2 := 8800;

  c1 := 460;
  c2 := 381;

  Tl := 100;
  Tr := 50;
  T0 := 10;
  L := 0.3;

  // определяем общее число узлов в пластине
  N := N1+N2+1;

  // определяем расчетный шаг сетки по пространственной координате
  h := L/(N-1);

  // определяем коэффициенты температуропроводности
  a1 := lamda1/(ro1*c1);
  a2 := lamda2/(ro2*c2);

  // определяем расчетный шаг сетки по времени
  tau := t_end/100.0;

  // определяем поле температуры в начальный момент времени
  For i:= 1 To N Do
    T[i] := T0;

  // проводим интегрирование нестационарного уравнения теплопроводности
  time := 0;
  While time < t_end Do
    Begin
      // увеличиваем переменную времени на шаг τ
      time := time + tau;

      // определяем начальные прогоночные коэффициенты на основе левого граничного условия
      alfa[1] := 0.0;
      beta[1] := Tl;

      // цикл с параметром для определения прогоночных коэффициентов по формуле (8) в первой части пластины
      For i:= 2 To N1 Do
        Begin
          // ai, bi, ci, fi – коэффициенты канонического представления СЛАУ с трехдиагональной матрицей
          ai := lamda1/sqr(h);
          bi := 2.0*lamda1/sqr(h)+ro1*c1/tau;
          ci := lamda1/sqr(h);
          fi := -ro1*c1*T[i]/tau;
          // alfa[i], beta[i] – прогоночные коэффициенты
          alfa[i] := ai/(bi-ci*alfa[i-1]);
          beta[i] := (ci*beta[i-1]-fi)/(bi-ci*alfa[i-1]);

          
        End;

     

      // определяем прогоночные коэффициенты на границе раздела двух частей, используем соотношения (28)
      alfa[N1+1] := 2.0*a1*a2*tau*lamda2/(2.0*a1*a2*tau*(lamda2+lamda1
                    *(1-alfa[N1]))+sqr(h)*(a1*lamda2+a2*lamda1));
      beta[N1+1] := (2.0*a1*a2*tau*lamda1*beta[N1]+sqr(h)*(a1*lamda2+a2
                    *lamda1)*T[N1+1])/(2.0*a1*a2*tau*(lamda2+lamda1
                    *(1-alfa[N1]))+sqr(h)*(a1*lamda2+a2*lamda1));
      

      // цикл с параметром для определения прогоночных коэффициентов по формуле (8) во второй части пластины
      For i:= N1+2 To N-1 Do
        Begin
          // ai, bi, ci, fi – коэффициенты канонического представления СЛАУ с трехдиагональной матрицей
          ai := lamda2/sqr(h);
          bi := 2.0*lamda2/sqr(h)+ro2*c2/tau;
          ci := lamda2/sqr(h);
          fi := -ro2*c2*T[i]/tau;
          // alfa[i], beta[i] – прогоночные коэффициенты
          alfa[i] := ai/(bi-ci*alfa[i-1]);
          beta[i] := (ci*beta[i-1]-fi)/(bi-ci*alfa[i-1]);
        End;

      // определяем значение температуры на правой границе
      T[N] := Tr;

      // For i := 1 to N Do
      //   Begin
      //     writeln(alfa[i]);
      //   End;
      // writeln();
      // For i := 1 to N Do
      //   Begin
      //     writeln(beta[i]);
      //   End;
    
      // Break;
      
      // используя соотношение (7) определяем неизвестное поле температуры
      For i:= N-1 Downto 1 Do
        T[i] := alfa[i]*T[i+1]+beta[i];
    End;

  for i:=1 to N do 
   writeln(' ',h*(i-1):10:8,' ',T[i]:8:5);
End.