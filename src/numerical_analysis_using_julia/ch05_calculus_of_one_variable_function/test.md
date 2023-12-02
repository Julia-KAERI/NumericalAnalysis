#

## Richardson's Extrapolation

Richardson extrapolation 은 어떤 값에 대한 추정값이 다항식의 truncation error 를 가질 때 근접한 추정값을 이용하여 truncation error 다항식의 차수를 높이는 방법이다. 예를 들어 보자.

### 1차 truncation error

구간을 $h$ 간격으로 나누었을 때 어떤 값 $M$ 에 대한 근사가 $N_1(h)$ 이며 그 오차가 $h$ 에 대한 다항식으로 다음과 같다고 하자.

$$
M = N_1(h) + k_1 h  + k_2 h^2 + k_3 h^3 + \cdots  \approx N_1 (h)+k_1 h. \tag{1}
$$

이 때의 truncation error 는 $O(h)$ 이다. 즉 어떤 $M$ 값을 계산하는데 구간 간격 $h$ 에 의존하는 계산값과 $h$ 에 대한 다항식으로 표현되는 truncation error 가 예상된다고 하자. 그렇다면,

$$
M = N_1\left(\dfrac{h}{2}\right) + k_1 \dfrac{h}{2}  + k_2 \dfrac{h^2}{4}+ k_3 \dfrac{h^3}{8} \cdots  \tag{2}
$$

이며, 위 두 식에 대해 $2\times (2)-(1)$ 을 수행하면 1차항이 사라지고

$$
M = 2 N_1\left(\dfrac{h}{2}\right) -N_1(h) - \dfrac{k_2}{2} h^2 - \dfrac{3k_3}{4} h^3 - \cdots
$$

를 얻는다. 즉 $N_1(h)$ 와 $N_1 \left(\dfrac{h}{2}\right)$ 값만으로 $O(h^2)$ 의 truncation 에러를 갖는 근사값을 구한 것이다.

</br>

## 2차 truncation error

Truncation error 가 $O(h^2)$ 라고 하자. 즉,

$$
M = N_1(h) + k_2 h^2 + k_3 h^3 + k_4 h^4 \cdots \tag{3}
$$

일 때

$$
M = N_1\left(\dfrac{h}{2}\right) +  k_2 \dfrac{h^2}{4}+ k_3 \dfrac{h^3}{8} + k_4 \dfrac{h^4}{16} + \cdots  \tag{4}
$$

라면 $(4\times (4) - (3)) \div 3$ 을 통해

$$
M = \dfrac{4}{3}N_1\left(\dfrac{h}{2}\right)-\dfrac{1}{3}N_1(h) - \dfrac{k_3}{6} h^3 - \dfrac{k_4}{4} h^4- \cdots
$$

을 얻는다. 즉 $N_2(h)= \dfrac{4}{3}N_1\left(\dfrac{h}{2}\right)-\dfrac{1}{3}N_1(h)$ 라 하면, $M$ 에 대한 근사는 $O(h^3)$ 의 truncation error 를 갖는다. 많은 경우 $k_3=0$ 인데 이 경우라면 $O(h^4)$ 의 truncation error 를 갖는다.

</br>

## $n$ 차 truncation error

$$
M = N_1(h) + k_n h^n + k_{n+1} h^{n+1} + k_{n+2} h^{n+2} \cdots \tag{5}
$$

일 때

$$
M = N_1\left(\dfrac{h}{2}\right) +  k_n \dfrac{h^n}{2^n}+ k_{n+1} \dfrac{h^{n+1}}{2^{n+1}} + k_{n+2} \dfrac{h^{n+2}}{2^{n+2}} + \cdots  \tag{6}
$$

이다. $(2^n \times (6)-(5))\div (2^{n}-1)$ 을 통해

$$
M = \dfrac{1}{2^n-1}\left(2^n N_1 \left(\dfrac{h}{2}\right) - N_1(h)\right) - \dfrac{k_{n+1}}{2(2^n-1)}h^{n+1} - \dfrac{3k_{n+2}}{4(2^n-1)} h^{n+2} - \cdots
$$

를 얻는다.
