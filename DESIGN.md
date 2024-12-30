# Introduction
Linear programming is a method to achieve the best outcome in a mathematical model whose requirements are represented by linear relationships. It is a special case of mathematical programming (mathematical optimization).

To solve a linear programming problem, we may apply several methods, such as the graphical method, the simplex method, and the interior-point method. I learned the simplex method recently and feels that it is really a intriguing and complex method with a lots of details to consider, just like what CPU and GPU design is. So I decide to implement a simple version of the simplex method to help me understand it better, and to ensure that I have a solid understanding of the method.

## Related Work

### 

## Formulation
### Algorithm Description
Psuedo code:
```Python
Simplex(A, b, c): # The inputs are parameters for simplex method slack form, consisting of minimizing a linear objective function subject to linear inequality (less or equal than) constraints and with all variables non-negative
    (B_I, A, B, c) = InitializeSimplex(A, b, c) # Initialize the simplex method. If the LP is feasible, return the initial basis B_I storing the indices of vectors in the correspoinding basis B; otherwise, report INFEASIBLE
    while True do
        # stoping criterion
        if there is no index e (1 <= e <= n) having c_e < 0 then
            x = CalculateX(B_I, A, B, b) # Calculate the solution x
            return x, c^T x # Return the solution x and the optimal value
        end if
        # optimize
        # find a optimization direction
        Choose an index e having c_e < 0 according to a certain rule
        # choose a step size
        for each index i (1 <= i <= m) do
            if A_{ie} > 0 then
                theta_i = b_i / A_{ie}
            else
                theta_i = +inf
            end if
        end for
        Choose the index l (1 <= l <= m) such that minimize theta_i
        # 这里 l 的含义其实已经是A中的行号了
        if theta_l = +inf then
            return UNBOUNDED
        end if
        # update the basis
        B_I, A, b, c, z = Pivot(B_I, A, B, b, c, z, e, l) # Pivot the basis
    end while
```

And the CalculateX and Pivot functions are as follows:
```Python
# Calculate x is used when the simplex method is terminated. It calculats the solution x according to the basis B and the corresponding b_i
# The terminated situation has a B that indicates the optimal solution vertex.
# B_I是一个长度为m的数组，B_I[i]是自变量序号
CalculateX(B_I, A, B, b):
    # assign non-basic variables to 0, and assign basic variables to the corresponding b_i. The inputs: B_I is the indices of vectors in the basis B, A is the matrix of the LP, B is the basis matrix, and b is the right-hand side of the LP
    for j = 1 to n do
        if j is not in B_I then
            x_j = 0
        else
            # 在 B_I 包含的index中，找到被变成identity的index，然后将对应的b_i赋值给x_j
            # 这里可以看到A应当是被B^{-1}过后的A
            for i = 1 to m do
                if a_{ij} = 1 then # 对于浮点数，这里应当是一个范围，或者用某个函数来判断。（但是对于很相近的情况，可能出现问题）
                    x_j = b_i
                end if
            end for
        end if
    end for
    return x
```

```Python
# Pivot is used to update the basis. 
# The inputs are the current basis B_I, the matrix A, 
# the basis matrix B, the right-hand side b, 
# the objective function c, the current optimal value z, 
# the entering index e, and the leaving index l
Pivot(B_I, A, b, c, z, e, l):
    # Aims to set the e-th column of A into one-hot vector, and the l-th variable row of this vector is set to 1
    # we also wants to preserve the property that
    # scale the l-th row
    # 实际上就是进行初等行变换，把序号为e的那一列变成单位向量，替换掉序号为l的那一列的基向量
    # 操作的方法是：获取l对应的基向量的行号l_row。感觉这里B_I应当是一个map，key是基向量的序号（1-m），value是这个基向量的1在A中的行号
    # 注意这些行号是以1开始的，而不是以0开始的。我们需要预留一个0位置，用来给初始化解的L_aux存储x_0
    # 然后把l_row这一行的所有元素（包括b）都除以A_{l_row, e}，然后把这一行scale并加到所有其他行上，使得A_{some_other_row, e} = 0（包括c\hat和z所在的行）
    # 最后把B_I中的l替换成e

    # scale the l-th row
    # lrow = B_I[l]
    # 对于这个算法来说，lrow刚好是l
    b[lrow] = b[lrow] / A[lrow, e]
    for j = 1 to e-1 do
        A[lrow, j] = A[lrow, j] / A[lrow, e]
    end for
    for j = e + 1 to n do
        A[lrow, j] = A[lrow, j] / A[lrow, e]
    end for
    A[lrow, e] = 1.0

    # scale and add to other rows
    for i = 1 to m but i != e do
        scale_factor = A[i, e]
        for j = 1 to n do
            A[i, j] = A[i, j] - scale_factor * A[lrow, j]
        end for
        b[i] = b[i] - scale_factor * b[lrow]
    end for

    # scale and add to c_hat and z
    scale_factor = c[e]
    for j = 1 to n do
        c[j] = c[j] - scale_factor * A[lrow, j]
    end for
    z = z - scale_factor * b[lrow]

    # update B_I
    newpair = {e: B_I[l]}
    B_I.remove(l)
    B_I.add(newpair)

```

```Python
# InitializeSimplex is used to initialize the simplex method.
# The inputs are the matrix A, the right-hand side b, and the objective function c
InitializeSimplex(A, b, c):
    let l be the index of the minimum b_i
    set B_I to include the indices of slack variables
    if b_l >= 0 then
        return (B_I, A, b, c, 0)
    end if
    construct L_{aux} by adding -x0 to each constraint, and using x_0 as the objective function
    let (new_A, new_b, new_c) = ConstructAuxiliary(A, b, c)
    # perform one step of pivot to make all b_i positive, simply by choosing 0 as the entering index and l as the leaving index
    (B_I, new_A, new_b, new_c, z) = Pivot(B_I, new_A, new_b, new_c, z, l, 0)
    iterate the WHILE loop in the Simplex algorithm until an optimal solution L_{aux} is found
    if z == 0 then # 注意这里可能要用一个范围来判断浮点数
        return the final slack form with x_0 removed and the objective function of L restored
        (方法是恢复原来的目标函数c，设置z为0，清空constraints中x0的系数为0，然后对解基的第一行做高斯行变换到0：类似可以复用pivot的逻辑，但是只需要做一行)
    else
        return INFEASIBLE
    end if
```

感觉要进行一些细节设计了。



上述算法都是在单纯型表格（simplex tabular）中利用初等行变换计算的方法。
实际上，我可以用矩阵的形式来形式化这个算法。

不，你认真读一下PPT，这个算法的实现细节。

最后发现在单纯型表格里做初等行变换是一个非常有效的方法。
不用矩阵的一般求逆算法，而是用初等行变换。
并行上，可以在行变换的过程中并行，而不是在一般的矩阵运算上并行。
一般的矩阵运算计算开销很大，且有很多数据传输，而直接在单纯型表格上做初等行变换，只需要传输少量的数据，且计算开销小。

看懂之后只能说这个PPT上的方法太巧妙了，就应该用单纯型表格来计算——始终维护一组base且让它们始终为单位矩阵。
可能历史上最开始有intuition的人就是简单地拓展成松弛型，然后开始做行变换。

### Implementation
用GPT生成了框架，然后自己写完了。用了几个测试用例测试了各种情况：
1. 有可行解
1. 有可行解且需要使用auxiliary LP初始化解
1. unbound
1. infeasible

## 商业软件
MATLAB可以求解线性规划问题。使用函数`linprog`即可。

# 整数线性规划
整数线性规划是线性规划的一个扩展，要求解的变量是整数。这实际上是一个复杂的问题。我们转而阅读《introduction to linear optimization》的第11章，前人有如下研究成果：
1. 整数线性规划问题，在于给出足够强的限制条件。在某些足够强的限制条件下，线性规划求解的结果就是整数规划求解的结果。
1. 理想情况下，可行域变成了一个convex hull——凸网格，此时线性规划的求解结果就是整数规划的求解结果。
1. 可以使用cutting-plane方法，即在每次求解连续线性规划到极端情况的时候，加入一些约束条件，使得解空间更小，最终收敛到整数解。也可以使用branch-and-bound方法，即在每次求解连续线性规划到最优时，将解空间分为两部分，分别求解。在分为两部分的过程中，把非整数变量所在的那个条带清除掉，然后继续求解。对于infeasible的子问题就丢掉，feasible的子问题就继续求解。

写出这两个方法的伪代码。
```Python
SolveILP_CuttingPlane():
    while True do
        SolveLP()
        CalculateX()
        if x is integer then
            return
        end if
        # add a constraint to cut off the current solution
        AddConstraintCuttingPlane()
    end while
```

在求解的过程中，当$B$矩阵变成单位矩阵之后，有如下一些结构：

$$
x_i + \sum_{j \in N} \lfloor\overline{a}_{ij}\rfloor x_j \leq \lfloor b_i \rfloor
$$

下面讨论操作细节：

由于$b_i$一定是非负的，因此$\lfloor b_i \rfloor $也是非负的，所以依然满足约束条件。在加入这个新的约束条件后，要向加入的所有式子加入松弛变量，形成等式。

加入的约束条件之前要考虑是否和原来的约束条件是同样的约束条件。如果是同样的约束条件，那么就不需要加入。如果不是，那么就需要加入。

在加入了这些约束条件之后，


```Python
AddConstraintCuttingPlane():
    # add a constraint to cut off the current solution
    # the constraint is a linear combination of the current solution

```

自己手动尝试了一下，感觉cut不是很聪明的样子——cut的大小依赖于约束条件和整数之间的差距。简单地把约束条件变成整数，无法保证下一次得到的解依然是整数，据说求解会很慢。

而branch and bound方法看起来是对解的范围进行了限制，利用分类讨论包含了所有范围内的问题，看起来更容易成功一些啊。

于是按照《introduction to linear optimization》一书的描述，实现了branch and bound方法。这里不再赘述。


