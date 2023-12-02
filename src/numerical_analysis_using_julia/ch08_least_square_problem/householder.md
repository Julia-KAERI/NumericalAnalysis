# Householder

하우스홀더 변환을 통해 변환된 벡터가 원래 백터 $\boldsymbol{x}$ 의 특정 성분이 $0$ 이 되도록 $\boldsymbol{v}$ 를 정할 수 있다. $k$ 번째 성분을 $0$ 으로 하려고 할 때

$$
\boldsymbol{v} = \boldsymbol{x} + \alpha \boldsymbol{e}_k
$$

라 하자. 이 때 $\alpha$ 는 실수이다. 그렇다면, 

$$
\begin{aligned}
(\boldsymbol{e}_k)_{ij} & = \delta_{ik}\delta_{j1}, \\
(\boldsymbol{e}^\ast_k)_{ij} & = \delta_{jk}\delta_{i1}, \\
(\boldsymbol{x})_{ij} & = x_i  \delta_{j1}, \\
(\boldsymbol{x})^\ast_{ij} & = \overline{x_j}  \delta_{i1}, \\
\end{aligned}
$$

을 얻는다. 이로부터 
$$
\begin{aligned}
(\boldsymbol{e}_k\boldsymbol{x}^\ast) _{ij} &=\sum_m (\boldsymbol{e}_k)_{im} (\boldsymbol{x}^\ast)_{mj} = \sum_m \delta_{ik}\delta_{m1}\delta_{m1}\overline{x_j} = \delta_{ik}\overline{x_j}, \\
(\boldsymbol{x}\boldsymbol{e}^\ast_k )_{ij} &= \sum_m(\boldsymbol{x})_{im} (\boldsymbol{e}^\ast_k)_{mj} = \sum_m \delta_{m1}x_{i}\delta_{jk}\delta_{m1} =x_i\delta_{jk} \\
(\boldsymbol{e}_k\boldsymbol{e}^\ast_k)_{ij} &= \sum_m (\boldsymbol{e}_k)_{im} (\boldsymbol{e}_k^\ast)_{mj} = \sum_m \delta_{ik}\delta_{m1}\delta_{m1} \delta_{jk} = \delta_{ik}\delta_{jk}
\end{aligned}
$$

를 얻는다.

$$
\begin{aligned}
\boldsymbol{v}\boldsymbol{v}^\ast &= \boldsymbol{xx}^\ast + \alpha \boldsymbol{e}_k \boldsymbol{x}^\ast  + \alpha \boldsymbol{x}\boldsymbol{e}_k^\ast + \alpha^2 \boldsymbol{e}_k\boldsymbol{e}_k^\ast, \\
\left[((\boldsymbol{xx}^\ast)\boldsymbol{x}\right] _{i} &=  \sum_m x_i \overline{x_m}x_m =  \|\boldsymbol{x}\|_2^2 x_i              \\
\left[(\boldsymbol{e}_k \boldsymbol{x}^\ast)\boldsymbol{x}\right]_{i} &= \sum_m (\boldsymbol{e}_k\boldsymbol{x}^\ast)_{im} (\boldsymbol{x})_{m} = \sum_m\delta_{ik}\overline{x_m} x_m  =   ( \|\boldsymbol{x}\|_2\boldsymbol{e}_k)_{i}, \\
\left[(\boldsymbol{xe}_k^\ast)\boldsymbol{x}\right]_{i} &=  \sum_m (\boldsymbol{x}\boldsymbol{e}^\ast_k)_{im} (\boldsymbol{x})_{m} = \sum_m x_i \delta_{mk} x_m  = x_kx_i\\ %x_k\boldsymbol{x}, \\
\left[(\boldsymbol{e}_k\boldsymbol{e}_k^\ast) \boldsymbol{x}\right]_{i} &= \sum_m \delta_{ik}\delta_{mk}x_m = x_k\delta_{ik} = (x_k \boldsymbol{e}_k)_i
\end{aligned}
$$

이다. 이를 이용해 $\boldsymbol{x}$ 의 $\boldsymbol{v}$ 에 대한 하우스홀더 변환을 계산하면,

$$
\begin{aligned}
\boldsymbol{H_v} \boldsymbol{x}&= \left(1 - \dfrac{2 \boldsymbol{v}\boldsymbol{v}^\ast}{\boldsymbol{x}^\ast \boldsymbol{x}+\alpha (\overline{x_k} + x_k)+{\alpha}^2 }\right)\boldsymbol{x} \\
&=\dfrac{\|\boldsymbol{x}\|_2^2\boldsymbol{x}+\alpha (\overline{x_k} + x_k)\boldsymbol{x} + \alpha^2\boldsymbol{x} - 2\|\boldsymbol{x}\|_2^2\boldsymbol{x} -2\alpha \|\boldsymbol{x}\|_2^2 \boldsymbol{e}_k -2\alpha x_k \boldsymbol{x} -2\alpha^2 x_k \boldsymbol{e}_k}{\boldsymbol{x}^\ast \boldsymbol{x}+\alpha (\overline{x_k} + x_k)+{\alpha}^2  }\\
&= \dfrac{\alpha^2 - \|\boldsymbol{x}\|_2^2}{\boldsymbol{x}^\ast \boldsymbol{x}+\alpha (\overline{x_k} + x_k)+{\alpha}^2 } \boldsymbol{x} - 2\alpha\dfrac{}{\boldsymbol{x}^\ast \boldsymbol{x}+\alpha (\overline{x_k} + x_k)+{\alpha}^2 }
\end{aligned}
$$
