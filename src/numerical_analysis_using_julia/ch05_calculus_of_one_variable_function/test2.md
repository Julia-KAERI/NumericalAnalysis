# Test

## Test2

</br>
 
### Simpson 1/3 적분

$[x_{i-1},\,x_{i+1}]$ 구간에서의 적분을 생각하자. 테일러 정리에 의해

$$
f(x) = f(x_i) _+ f'(x_i)(x-x_i) + \dfrac{f''(x_i)}{2}(x-x_i)^2 + \dfrac{f^{(3)}(x_i)}{6}(x-x_i)^3 + \dfrac{f^{(4)}(\xi_i)}{24}(x-x_i)^4
$$

를 만족하는 $\xi_i \in [x_{i-1},\, x_{i+1}]$ 이 존재한다. 양변을 적분하면

$$
\begin{aligned}
\int_{x_{i-1}}^{x_{i+1}}f(x)\,dx &= f(x_i)(2h) + \dfrac{f''(x_i)}{3}h^3 + \dfrac{f^{(4)}(\xi_i)}{60}h^5 
\end{aligned}
$$

이다. 또한 @eq-difference_3points 로부터

$$
f''(x) =  \dfrac{f(x+h) - 2f(x) + f(x-h)}{h^2} - \dfrac{f^{(4)}(\xi_i')}{36}h^5 
$$

를 만족하는 $\xi_i'\in [x_{i-1},\,x_{i+1}]$ 가 존재한다는 것을 알고 있다. 따라서

$$
\begin{aligned}
\int_{x_{i-1}}^{x_{i+1}}f(x)\,dx &= f(x_i)(2h) + \dfrac{h}{3}\left( f(x_{i+1}) -2f(x_i) + f(x_{x-i}) \right)  + \left(\dfrac{f^{(4)}(\xi_i)}{60} - \dfrac{f^{(4)}(\xi_i')}{36} \right)h^5 \\
&=\dfrac{1}{3} \left[f(x_{i-1}) + 4f(x_{i}) + f(x_{i+1})\right] -  \left( \dfrac{f^{(4)}(\xi_i')}{36} - \dfrac{f^{(4)}(\xi_i)}{60} \right)h^5
\end{aligned}
$$

따라서 우리는 Simplson 1/3 적분에서의 오차가 $O(h^5)$ 임을 알 수 있다. 여기서는 $\xi_i$ 와 $\xi_i'$ 에서의 4차 도함수값이 필요했지만 실제로는 

$$
\begin{aligned}
\int_{x_{i-1}}^{x_{i+1}}f(x)\,dx =\dfrac{h}{3} \left[f(x_{i-1}) + 4f(x_{i}) + f(x_{i+1})\right] -  \dfrac{f^{(4)}(\overline{\xi}_i)}{90}h^5
\end{aligned}
$$

를 만족하는 $\overline{\xi}_i \in [x_{i-1},\,x_{i+1}]$ 가 존재한다(Atkinson, Kendall E. (1989). An Introduction to Numerical Analysis (2nd ed.). John Wiley & Sons. 을 참고하라)

이제 $[a,\,b]$ 구간을 $a=x_1<x_2<\cdots<x_{n-1}<x_n=b$ 의 $n-1$ 개의 구간으로 나누어 적분하여 합치는 것을 생각하자. Simpson 1/3 적분은 2개의 구간을 한꺼번에 적분하므로 $n$ 이 짝수일 때와 홀수일 때 식이 다르다. $n$ 이 짝수이면 $2h$ 의 너비를 갖는 $n/2$ 개의 구간으로 나누어 적분을 수행 할 수 있다. 에러도 $n/2$ points 에 대한 중간값 정리를 사용하면 다음을 만족하는 $\xi \in [x_1,\,x_n]$ 이 존재한다.

$$
\begin{aligned}
\int_{a}^b f(x)\, dx &= \dfrac{h}{3} \left[f(x_1) + 4f(x_2) + 2f(x_3) + 4f(x_4) + 2f(x_5)+\cdots \right.\\
&\qquad \cdots \left.+ 4f(x_{n-3})+2f(x_{n-2}) + 4f(x_{n-1})+f(x_n)\right] - \dfrac{n}{2}\dfrac{f^{(4)}(\xi)}{90}h^5
\end{aligned}
$$

</br>

### Simpson 3/8 적분

Simpson 1/3 적분이 전체 $n$ 개의 구간을 2개씩 묶어서 적분하여 합쳤다면 3/8 적분은 3개씩 묶어서 합친다. 라그랑쥬 다항식을 이용하면,

$$
L_4(x) = f(x_{i-1})l_{i-1}(x) + f(x_i)l_i(x) + f(x_{i+1})l_{i+1}(x) + f(x_{i+2})l_{i+2}(x)
$$

에 대해

$$
f(x) = L_4 (x) + \dfrac{f^{(4)}(\xi_i)}{4!}\prod_{j=1}^4(x-x_{j-2})
$$

를 만족하는 $\xi_i\in [x_{i-1},\, x_{i+2}]$ 가 존재한다.

$$
\begin{aligned}
\int_{x_{i-1}}^{x_{i+2}} f(x)\, dx = \dfrac{3h}{8} \left( f(x_{i-1}) + 3 f(x_{i}) + 3f(x_{i+1}) + f(x_{i+2})\right) - \dfrac{3}{80}f^{(4)}(\xi) h^5
\end{aligned}
$$

를 얻는다. 이것을 $[a,\,b]$ 구간을 $a=x_1<\cdots <x_{3n+1}=b$ 이며 $x_{i+1}-x_i = h=\text{const.}$ 라면 

$$
\begin{aligned}
\int_a^b f(x)\, dx &= \dfrac{3h}{8} \left[f(x_1) + f(x_{3n+1}) + 3\sum_{i=1}^n \left(f(x_{3i+2}) + f(x_{3i})\right)  + 2\sum_{i=1}^{n-1} f(x_{3i+1})\right] - \dfrac{3n}{80}f^{(4)} (\xi) h^5
\end{aligned}
$$

를 만족하는 $\xi \in [a,\,b]$ 가 존재한다. 여기서 

$$
\mathcal{I}_{3/8}[f,\,a,\,b,\, 3n]=\dfrac{3h}{8} \left[f(x_1) + f(x_{3n+1}) + 3\sum_{i=1}^n \left(f(x_{3i+2}) + f(x_{3i})\right)  + 2\sum_{i=1}^{n-1} f(x_{3i+1})\right] 
$$

는 Simpson 3/8 적분값이며 그 오차는 $- \dfrac{3n}{80}f^{(4)} (\xi) h^5$ 이다.