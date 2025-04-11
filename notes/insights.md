# From thick to thin
笔者的一些insights, 便于理解和记忆数理经济学的体系
___
## 期中
___
*前三章*
- 正定矩阵及其判定
- 点集拓扑, 紧性
- 拉格朗日乘数法

### 1 极值与零点
- **无约束优化问题** :  $f : \mathbb{R}^n \to \mathbb{R}$, 目标是求 $x^*$ 使得 $f(x^*)$ 极小
- **求函数组零点** : $F : \mathbb{R}^n \to \mathbb{R}^n$, 目标是求 $x^*$ 使得 $F(x^*) = 0$

以上两个问题的交集在于 FOC: 
$$ F = \nabla f = 0 $$
换句话说, 优化问题求解的第一步是零点问题, 零点问题的一些特例可以还原成优化问题 (显然不是所有的方程组都是某个函数的 FOC), 假如 $F = \nabla f$ :
- $f$ 的 Hessian 矩阵是 $F$ 的Jacobian 矩阵, 他们都是对称矩阵
$$ \begin{pmatrix}
        \frac{\partial^2 f}{\partial x_1^2} & \cdots & \frac{\partial^2 f}{\partial x_1 \partial x_n} \\
        \vdots & \ddots & \vdots \\
        \frac{\partial^2 f}{\partial x_n \partial x_1} & \cdots & \frac{\partial^2 f}{\partial x_n^2}
    \end{pmatrix} = \begin{pmatrix}
        \frac{\partial F_1}{\partial x_1} & \cdots & \frac{\partial F_1}{\partial x_n} \\
        \vdots & \ddots & \vdots \\
        \frac{\partial F_n}{\partial x_1} & \cdots & \frac{\partial F_n}{\partial x_n}
    \end{pmatrix}$$

对于无约束优化问题, 有两步走 : 第一步求一阶条件的解, 第二部确定这个解是极大值, 极小值还是鞍点
#### 1.1 找零点的方法
我们在*第四章*介绍了五种方法
- **Bisection** : 二分法, 适用于单变量函数, 需要函数在区间上连续, 且两端异号
- **Secant** : 割线法, 适用于单变量函数, 需要函数在区间上连续, 且两端异号, 但不需要函数可微
- **False Position** : 假位法, 适用于单变量函数, 需要函数在区间上连续, 且两端异号, 但不需要函数可微, 但比割线法更快收敛
- **Newton Raphson** : 牛顿法, 适用于单变量函数, 需要函数在区间上连续, 且两端异号, 需要函数可微, 但不需要函数二阶可微
- **Gradient Descent** : 梯度下降法, 适用于多变量函数, 需要函数在区间上连续, 且两端异号, 需要函数可微, 但不需要函数二阶可微
#### 1.2 极值点or鞍点
规范方法 : 判断 Hessian 矩阵的正定性, 本课程用顺序主子式法来判断会比较快
|$H_f$|$H_i$ | $x^*$|
|-|-|-|
|正定|$H_i > 0$|极小值|
|负定|$H_i < 0$|极大值|
|不定|$o.w.$|鞍点|

> 正小负大

Trick1 : 正定矩阵的特征值都是正数, 负定矩阵的特征值都是负数
Trick2 : 将 $x^*$ 代入 $f(x)$ 先做排除
Trick3 : 结合 $f$ 的凹凸性判断, 如果严格凸或者严格凹, 则 $x^*$ 一定是极小值或者极大值
### 2 无约束优化与有约束优化

有约束优化是*第五章*的内容, 即在可微函数**组** $g(x) = 0$ 上求 $f$ 的极值, 沟通二者的桥梁是拉格拉日乘数法与 **NDCQ**
$$ L(x, \lambda) = f(x) + \lambda^\top g(x) $$
然后 $(x^*, \lambda^*)$ 是一阶条件 $\nabla L(x^*, \lambda^*) = 0$ 的解, 我们似乎回到了无约束优化

#### 2.1 非退化约束条件
NDCQ (非退化约束条件) 是让拉格朗日乘数法成立的条件

假设 $x$ 是 $k$ 维向量, $\lambda$ 是 $m$ 维向量, 则 $L$ 是 $k + m$ 维向量, 考虑
$$Dg(x^*) = \begin{pmatrix}
        \frac{\partial g_1}{\partial x_1} & \cdots & \frac{\partial g_1}{\partial x_k} \\
        \vdots & \ddots & \vdots \\
        \frac{\partial g_m}{\partial x_1} & \cdots & \frac{\partial g_m}{\partial x_k}
    \end{pmatrix}$$
- **NDCQ** : $Dg(x^*)$ 的秩为 $m$. 

- 假如 **NDCQ** 得到满足, 则 $(x^*, \lambda^*)$ 可能是极值点
- 假如不满足, 即约束 $g$ 是退化的, 此时可能仍然存在极值点, 但无法由拉格朗日乘数法得到

#### 2.2 加边海色矩阵
加边海色矩阵是判断有约束优化问题极值点的工具:
$$ \bar{H} = \begin{pmatrix}
        0 & Dg(x^*) \\
        Dg(x^*)^\top & H_f(x^*)
    \end{pmatrix} $$
加边后显然不正定了, 判断方法也有变化 :

|$\bar{H}$|$H_{i}, i\in [2m+1,n+m]$ | $x^*$|
|-|-|-|
||$\text{sgn} H_{i} = (-1)^{m}$|极小值|
||$\text{sgn} H_{i} = (-1)^{i-m}$|极大值|
||$\text{sgn} H_{i} \neq 0$, 但不属于以上两种|鞍点|

对加边海色矩阵的一个直观理解 : 假如 $f$ 在 $g = 0$ 这个$n-m$ 维的嵌入 $\mathbb{R}^n$ 的流形上表现出了类似正定负定的性质, 那么可以判断是极值还是鞍点, 比如我们在 3 维空间中的球面上找极值.
> 实际上, 课堂上只讲了 $m = 1$ 的情形, 在此情形下极小值对应主子式恒负, 极大值对应符号交替, 这样就不对称, 和无约束优化问题的海色矩阵对比起来显得不自然
> 经济学中最常见的 $m=1, n=2$ 的情形, 自然只需看 $|\bar{H}|$ 的符号, 正大负小
> 
>[这个直观理解来自于 CMU 的笔记](https://www.math.cmu.edu/~hanifc/NotesOnHessians.pdf)
>[判断方法和完整证明来自 Northwestern University 的笔记](https://sites.math.northwestern.edu/~clark/285/2006-07/handouts/lagrange-2deriv.pdf)

___
## 期末
___