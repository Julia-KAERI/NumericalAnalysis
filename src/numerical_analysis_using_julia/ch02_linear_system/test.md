### 하우스홀더 행렬과 QR 분해


$\boldsymbol{v}\in \mathbb{R}^n$ 에 대해 아래와 같이 정의된 $\boldsymbol{P}_v$ 를 **하우스홀더 행렬 (Householder matrix)**(혹은 하우스홀더 행렬, 하우스홀더 변환) 이라 한다.

$$
\begin{aligned}
\boldsymbol{P}_{\boldsymbol{v}} := I_n- \dfrac{2\boldsymbol{v}\boldsymbol{v}^{\ast}}{\|\boldsymbol{v}\|^2}, \qquad \text{i. e. }\quad 
(\boldsymbol{P}_{\boldsymbol{v}})_{ij} := \delta_{ij} - \dfrac{2 v_i \overline{v_j}}{\|\boldsymbol{v}\|^2}.
\end{aligned}
$$

여기서 $\boldsymbol{v}\boldsymbol{v}^{\ast}$ 는 벡터의 내적이 아니라 $n\times 1$ 행렬 $\boldsymbol{v}$ 와 $1 \times n$ 행렬 $\boldsymbol{v}^{\ast}$ 가 곱해진 $n \times n$ 행렬을 의미한다. 벡터 $\boldsymbol{x}\in \mathbb{C}^n$ 에 대해 $\boldsymbol{P}_{\boldsymbol{v}} \boldsymbol{x}\in \mathbb{C}^n$ 을 **하우스홀더 변환** 이라 한다. $\boldsymbol{v}$ 가 단위벡터일 경우, 즉 $\|\boldsymbol{v}\|=1$ 라면 좀 더 간단하게 쓸 수 있다.

$$
\boldsymbol{P}_{\boldsymbol{v}} = I_n- 2\boldsymbol{v}\boldsymbol{v}^{\ast},\qquad \text{where } \|\boldsymbol{v}\|=1.
$$

</br>


::: {#exr-perperties_of_householder_matrix}

위와 같이 정의된 하우스홀더 행렬은 다음의 특징을 가진다.

1. Hermitian 이다. 즉  $\boldsymbol{P}_{\boldsymbol{v}}=  \boldsymbol{P}_{\boldsymbol{v}}^\ast$.

2. 직교행렬이다. 즉 $\boldsymbol{P}_{\boldsymbol{v}} (\boldsymbol{P}_{\boldsymbol{v}})^\ast = I$.

3. $(\boldsymbol{P}_{\boldsymbol{v}})^2=I$ 이다.

:::

::: {.solution}
$\boldsymbol{v}$ 가 단위행렬일 경우에만 보여도 된다. $\boldsymbol{P} = \boldsymbol{P}_{\boldsymbol{v}}$ 라 하면, 

$$
(\boldsymbol{P}^\ast)_{ij} = \overline{P_{ji}}= \overline{\delta_{ij}-2 v_j \overline{v_i}} = \delta_{ij}-2 v_i \overline{v_j} = P_{ji}
$$

이므로 $\boldsymbol{P}$ 는 hermitian matrix 이다. 또한,

$$
\begin{aligned}
\left(\boldsymbol{P} (\boldsymbol{P}^\ast)\right)_{ij} &= \left((\boldsymbol{P})^2\right)_{ij} = \sum_{k}(\delta_{ik}-2 v_i \overline{v_k})(\delta_{kj} -2 v_k \overline{v_j}) \\
&= \sum_k \delta_{ik}\delta_{kj} - 2 \sum_k \delta_{ik}v_k \overline{v_j} - 2 \sum_k \delta_{kj} v_i \overline{v_k} + 4 \sum_{k} v_i\overline{v_j}v_k \overline{v_k} \\
&= \delta_{ij} - 2 v_i \overline{v_{j}} - 2 v_i \overline{v_j} + 4 v_i \overline{v_j} = \delta_{ij}
\end{aligned}
$$

이다. 즉 $\boldsymbol{P}(\boldsymbol{P})^\ast = I$ 이므로 orthogonal 하며, $(\boldsymbol{P}_{\boldsymbol{v}})^2=I$ 이다.

:::

</br>

#### 하우스홀더 변환과 반사

하우스 홀더 변환이 리플렉션(reflection, 반사) 라고 불리는 이유는 다음과 같다. 벡터 $\boldsymbol{v}$ 만으로 $\boldsymbol{v}$ 와 수직이며 원점을 지나는 평면이 유일하게 정의된다. 임의의 벡터 $\boldsymbol{x}$ 에 대해,

$$
\boldsymbol{P}_{\boldsymbol{v}}\boldsymbol{x} = \boldsymbol{x} - \dfrac{2 \langle \boldsymbol{x},\,  \boldsymbol{v} \rangle \boldsymbol{v}}{\|\boldsymbol{v}\|^2}
$$

이다. 이로부터,

$$
\begin{aligned}
\dfrac{1}{2 }(\boldsymbol{x}+\boldsymbol{P_vx}) &= \boldsymbol{x} - \dfrac{\langle \boldsymbol{v,\,  x}\rangle \boldsymbol{v}}{\|\boldsymbol{v}\|^2}= \boldsymbol{x} - \text{Proj}_{\boldsymbol{v}} \boldsymbol{x} \\
\boldsymbol{x}-\boldsymbol{Px} &= \dfrac{2 \langle \boldsymbol{v ,\, x} \rangle  \boldsymbol{v}}{\|\boldsymbol{v}\|^2} = 2 \, \text{Proj}_{\boldsymbol{v}}\boldsymbol{x}
\end{aligned}
$$

이다. 즉 $\boldsymbol{x}$ 와 $\boldsymbol{Px}$ 는 $\boldsymbol{v}$ 에 의해 정의되는 평면에 대해 대칭이다.

$\boldsymbol{x}$ 를 $\boldsymbol{v}$ 와 평행한 부분과 수직한 부분으로 분리하자. 즉 $\boldsymbol{x}_{\|} = \text{Proj}_{\boldsymbol{v}}\boldsymbol{x}$, $\boldsymbol{x}_{\perp} = \boldsymbol{x}-\text{Proj}_{\boldsymbol{v}}\boldsymbol{x}$ 라 하면, $\boldsymbol{x} = \boldsymbol{x}_{\|} + \boldsymbol{x}_{\perp}$ 이며 $\boldsymbol{x}\cdot \boldsymbol{v} = \boldsymbol{x}_{\|}\cdot \boldsymbol{v}$ 이다.

$$
\boldsymbol{P}_{\boldsymbol{v}}\boldsymbol{x} = \boldsymbol{x} - 2 \boldsymbol{x}_{\|} = \boldsymbol{x}_{\perp} - \boldsymbol{x}_{\|}
$$

이다.


</br>

#### 하우스홀더 변환과 QR 분해

$\boldsymbol{x} \in \mathbb{F}^m$ 를 행렬 $\boldsymbol{A}\in \mathbb{F}^{m \times n}$ 의 첫번째 열이라 하자. $\{\boldsymbol{e}_1,\ldots,\,\boldsymbol{e}_m\}$ 은 표준 기저벡터이다. $\alpha$ 를 다음과 같이 정의한다.

$$
\alpha = \left\{ \begin{array}{ll} -\|\boldsymbol{x}\|  & \text{where } \mathbb{F}=\mathbb{R}, \\-\exp (i \,\text{arg} (x_1))\|\boldsymbol{x}\| \qquad & \text{where } \mathbb{F} =\mathbb{C}. \end{array} \right.
$$
$\mathbb{F} = \mathbb{C}$ 이면 $\boldsymbol{x}$ 의 첫번째 성분 $x_1 = a+ib,\, a,\,b\in \mathbb{R}$ 에 대해 $\text{arg} (x_{1})= \arctan(b/a)$ 이다. 다음과 같이 정하자.

$$
\begin{aligned}
\boldsymbol{u} &= \boldsymbol{x} + \alpha \boldsymbol{e}_1, \\
\boldsymbol{P} &= I - 2\dfrac{\boldsymbol{u}\boldsymbol{u}^\ast }{\|\boldsymbol{u}\|^2}
\end{aligned}
$$

이므로,

$$
\boldsymbol{Px} = \boldsymbol{x} - 2 \dfrac{\langle \boldsymbol{u},\,\boldsymbol{x} \rangle}{\|\boldsymbol{u}\|^2}\boldsymbol{u}
$$

이다. 그런데 $|\alpha|^2 = \|\boldsymbol{x}\|^2$ 이므로,

$$
\begin{aligned}
\langle \boldsymbol{u},\,\boldsymbol{x}\rangle & = \langle \boldsymbol{x}-\alpha \boldsymbol{e}_1,\, \boldsymbol{x} \rangle = \|\boldsymbol{x}\|^2-\overline{\alpha} x_1 = \|\boldsymbol{x}\|^2 - \text{Re} (x_1)\|\boldsymbol{x}\|\\
\langle \boldsymbol{u},\,\boldsymbol{u}\rangle &= \langle \boldsymbol{x} -\alpha \boldsymbol{e}_1,\, \boldsymbol{x} -\alpha\boldsymbol{e}_1 \rangle  = \|\boldsymbol{x}\|^2 - \overline{\alpha} x_1 -\alpha \overline{x_1} + |\alpha|^2 \\
&= 2\|\boldsymbol{x}\|^2 - 2\text{Re} ( \alpha x_1) = 2\langle \boldsymbol{u},\, \boldsymbol{x} \rangle
\end{aligned}
$$

이다. 즉,

$$
\boldsymbol{Px} = \boldsymbol{x} - \boldsymbol{u} = -\alpha\boldsymbol{e}_1
$$

이다. $\mathbb{F} = \mathbb{R}$ 이면 $\boldsymbol{Px} = \|\boldsymbol{x}\|\boldsymbol{e}_1$ 이며 $\mathbb{F}=\mathbb{C}$ 이면 $\boldsymbol{Px} = \exp(i \arg(x_1))\|\boldsymbol{x}\| \boldsymbol{e}_1$ 이 된다.


이제 $\boldsymbol{P}_1= \boldsymbol{P},\, \boldsymbol{A}=\boldsymbol{A}_1$ 이라 놓으면, 

$$
\boldsymbol{P}_1\boldsymbol{A}_1 = \begin{bmatrix}\alpha _{1} & \ast &\cdots &\ast \\0 & & &\\ \vdots & & \boldsymbol{A}_2 & \\ 0 & & & \end{bmatrix}
$$

꼴이 된다. 이제 행렬 $m \times n$ 행렬 $\boldsymbol{A}_k$ 가 첫번째 행부터 $k$ 번째 행까지는 상삼각 행렬의 모양을 따른다고 하자. 즉 $\boldsymbol{A}_k$ 가 $k \times k$ 상삼각 행렬 행렬 $\boldsymbol{B}_k$, $(m-k) \times (n-k)$ 행렬 $\boldsymbol{A}_k'$ 로 다음과 같이 표현된다고 하자.

$$
\boldsymbol{A}_k = \begin{bmatrix} \boldsymbol{B}_k & \boldsymbol{C}_k \\ \boldsymbol{0}  & \boldsymbol{A}_k'\end{bmatrix}.
$$

이제 $\boldsymbol{x}$ 에 대해 $\boldsymbol{x}^{(k)}$ 를

$$
\boldsymbol{x}^{(j)} := \begin{bmatrix} \underbrace{\begin{array}{lll} 0 & 0 & 0 \end{array}}_{k-1} & x_k & x_{k+1}  \end{bmatrix}
$$


$$
\boldsymbol{P}_k = \begin{bmatrix} \boldsymbol{I}_{k-1} & 0 \\ 0 & \boldsymbol{P}_k'\end{bmatrix}
$$

이라 하면, $(k-1)\times (k-1)$ 단위행렬 $I_{k-1}$ 과 $\boldsymbol{A}_k$ 를 $k$ 행 과 $k$ 열 부터 잘라 $\boldsymbol{A}_{k} = \begin{bmatrix} B_{k} & C_{k} \\ 0 &\boldsymbol{A}'_{k}\end{bmatrix}$  로 만들자. $B_k$ 는 $(k-1) \times (k-1)$ 행렬이며 $\boldsymbol{A}'_{k}$ 는 $(m-k+1)\times (n-k+1)$ 행렬이다. $\boldsymbol{A}_k$ 가 $k$ 열까지 상삼각 행렬 모양이므로 $B_k$ 아래는 $0$ 행렬이다. 두 행렬의 곱은
$$
\boldsymbol{Q}_k \boldsymbol{A}_k = \begin{bmatrix} \boldsymbol{I}_{k-1} & 0 \\ 0 & \boldsymbol{Q}_k' \end{bmatrix} \begin{bmatrix} B_{k} & C_{k} \\ 0 &\boldsymbol{A}'_{k}\end{bmatrix} = \begin{bmatrix} B_k & C_k \\ 0 & \boldsymbol{Q}'_k \boldsymbol{A'}_k \end{bmatrix}
$$

이 되고 $\boldsymbol{Q}'_k \boldsymbol{A}'_k$ 의 첫번째 열은 첫번째 행을 제외하면 모두 $0$ 이므로 $\boldsymbol{Q}_k \boldsymbol{A}_k$ 는 $k$ 열까지 상삼각 행렬 꼴이 된다.

$L = \min\{m,\,n\}$ 이라 하면 $\boldsymbol{Q}_L \boldsymbol{Q}_{L-1} \cdots \boldsymbol{Q}_1 \boldsymbol{A}$ 는 상삼각행렬꼴이 된다. 이를 $\boldsymbol{R}$ 이라 하자. $\boldsymbol{Q}'_k$ 가 하우스홀더 행렬이므로

$$
\boldsymbol{Q}_k \boldsymbol{Q}_k^\ast = \begin{bmatrix} \boldsymbol{I}_{k-1} & 0 \\ 0 & \boldsymbol{Q}_k'\end{bmatrix} \begin{bmatrix} \boldsymbol{I}_{k-1} & 0 \\ 0 & (\boldsymbol{Q}_k')^\ast\end{bmatrix} = \begin{bmatrix} I_{k-1} & 0 \\0 & \boldsymbol{Q}_k' 
(\boldsymbol{Q}_k')^\ast\end{bmatrix} = I
$$

이다. 즉 $\boldsymbol{Q}_k$ 도 직교행렬이다. $\boldsymbol{Q}_k$ 가 에르미트 행렬임은 쉽게 보일 수 있다. 이제,

$$
\boldsymbol{Q}_L \cdots \boldsymbol{Q}_1 \boldsymbol{A} = \boldsymbol{R} \implies \boldsymbol{A} = \boldsymbol{Q}_1^\ast \cdots \boldsymbol{Q}_L^\ast \boldsymbol{R}
$$

임은 쉽게 보일 수 있다. 직교행렬의 곱은 직교행렬이므로 $\boldsymbol{Q}_1^\ast \cdots \boldsymbol{Q}_L^\ast$ 도 직교행렬이다. 따라서 QR 분해를 할 수 있다. 

