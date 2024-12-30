#include <Simplex.hpp>
#include <Config.hpp>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <list>
#include <memory>

const double EPSILON = 1e-9; // Define a small tolerance value

bool isClose(double a, double b, double epsilon = EPSILON) {
    return std::fabs(a - b) < epsilon;
}

bool Simplex::isInfeasible() {
    return infesible;
}

#ifdef SIMPLEX_DEBUG
void Simplex::printStatus() {

    // 输出保留两位小数
    std::cout.precision(2);

    std::cout << "[DEBUG] tabular status: " << std::endl;
    // output c
    std::cout << "  | ";
    for (size_t i = 0; i <= n; i++) {
        std::cout << std::setw(5) << c[i] << " ";
    }
    std::cout << "| z: " << std::setw(5) << z << std::endl;
    // output A
    for (size_t i = 0; i <= m; i++) {
        std::cout << B_I[i] << " | ";
        for (size_t j = 0; j <= n; j++) {
            std::cout << std::setw(5) << A[i][j] << " ";
        }
        std::cout << "| " << b[i] << std::endl;
    }
    // x
    // std::cout << "[DEBUG] x: ";
    // for (size_t i = 0; i <= n; i++) {
    //     std::cout << x[i] << " ";
    // }
    std::cout << std::endl;

    /*
    
    std::cout << "[DEBUG] A matrix (include 0): " << std::endl;
    for(size_t i = 0; i <= m; i++) {
        for(size_t j = 0; j <= n; j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    // output the c
    std::cout << std::endl << "[DEBUG] c (include 0): ";
    for(size_t i = 0; i < c.size(); i++) {
        std::cout << c[i] << " ";
    }
    // output the b
    std::cout << std::endl << "[DEBUG] b (include 0): ";
    for(size_t i = 0; i < b.size(); i++) {
        std::cout << b[i] << " ";
    }
    std::cout << std::endl;
    // output the B_I
    std::cout << "[DEBUG] B_I (include 0): ";
    for(size_t i = 0; i <= m; i++) {
        std::cout << B_I[i] << " ";
    }
    std::cout << std::endl;
    // output the x
    std::cout << "[DEBUG] x (include 0): ";
    for(size_t i = 0; i <= n; i++) {
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;
    // output the z
    std::cout << "[DEBUG] z: " << z << std::endl;
    */

}
#endif

// 复制构造函数
Simplex::Simplex(const Simplex& other)
    : A(other.A), b(other.b), c(other.c), c_old(other.c_old), x(other.x), B_I(other.B_I), infesible(other.infesible), unbound(other.unbound), m(other.m), n(other.n), z(other.z), z_old(other.z_old), original_var_num(other.original_var_num) {}

// 赋值运算符
Simplex& Simplex::operator=(const Simplex& other) {
    if (this != &other) {
        A = other.A;
        b = other.b;
        c = other.c;
        c_old = other.c_old;
        x = other.x;
        B_I = other.B_I;
        infesible = other.infesible;
        unbound = other.unbound;
        m = other.m;
        n = other.n;
        z = other.z;
        z_old = other.z_old;
        original_var_num = other.original_var_num;
    }
    return *this;
}

double Simplex::getOptimalValue() {
    return -z;
}

/**
 * Constructor for the Simplex class.
 * @param A Coefficients of the constraints, in standard form
 * @param b Right-hand side values of the constraints
 * @param c Coefficients of the objective function
 * 
 * Copy and store the input data, initialize the tableau, and set the basic variables.
 * We store the data in the class in slack form
 */
Simplex::Simplex(const std::vector<std::vector<double>>& inputA, const std::vector<double>& inputb, const std::vector<double>& inputc) {
    // 标号从1开始，是为了给x0预留空间
    // b的长度比inputb多1，是为了统一
    b.resize(inputb.size()+1);
    std::copy(inputb.begin(), inputb.end(), b.begin()+1);
    m = inputb.size();
    int csize = inputA[0].size();
    original_var_num = csize;
    n = csize + m;
    // c的长度比n多1，是为了给x0预留空间
    c.resize(n+1);
    std::copy(inputc.begin(), inputc.end(), c.begin()+1);
    // copy data of A
    A.resize(m+1, std::vector<double>(n+1, 0));
    for(size_t i = 0; i < m; i++) {
        std::copy(inputA[i].begin(), inputA[i].end(), A[i+1].begin()+1);
    }
    // initialize the base variables and base matrix
    for(size_t i = 1; i <= m; i++) {
        A[i][csize+i] = 1;
    }
    #ifdef SIMPLEX_DEBUG
    std::cout << "[DEBUG] Simplex initialized. m: " << m << " n: " << n  << " csize: " << csize << std::endl;
    // output the A matrix
    std::cout << "[DEBUG] A matrix (include 0): " << std::endl;
    for(size_t i = 0; i <= m; i++) {
        for(size_t j = 0; j <= n; j++) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    // output the c
    std::cout << std::endl << "[DEBUG] c (include 0): ";
    for(size_t i = 0; i < c.size(); i++) {
        std::cout << c[i] << " ";
    }
    // output the b
    std::cout << std::endl << "[DEBUG] b (include 0): ";
    for(size_t i = 0; i < b.size(); i++) {
        std::cout << b[i] << " ";
    }
    std::cout << std::endl;
    #endif
}

// 新的构造函数，用于添加单个变量约束
Simplex::Simplex(const Simplex& other, int x_index, bool isGreaterEqual, double rhs)
    : A(other.A), b(other.b), c(other.c), c_old(other.c_old), x(other.x), B_I(other.B_I), infesible(other.infesible), unbound(other.unbound), m(other.m), n(other.n), z(other.z), z_old(other.z_old), original_var_num(other.original_var_num) {
    /**
     * 1.
     * 拓宽A列
     * 提升x的个数
     * 拓宽A行
     * 拓宽b
     * 加入constraint和松弛变量
     * 拓宽c
     * 拓宽B_I，存储最后一个变量的值
     * 
     * 2.
     * 对原问题进行初等行变换，保持原来的列基
     * 
     * 3.
     * 如果有负数，则重新计算初始解，方法是拓展x0，构造辅助线性规划问题，利用initialize时的处理方法
     */
    
    // 拓宽A列
    for(size_t i = 0; i <= m; i++) {
        A[i].push_back(0);
    }
    n++;
    // 提升x的个数
    x.push_back(0);
    // 拓宽A行
    A.push_back(std::vector<double>(n+1, 0));
    m++;
    // 拓宽b. 
    b.push_back(0); 
    // 拓宽c
    c.push_back(0);
    // 拓宽B_I
    B_I.push_back(0);
    B_I[m] = n;

    // 加入constraint和松弛变量.如果是x_i >= rhs，则化为-x_i + slack = -rhs，如果是x_i <= rhs，则化为x_i + slack = rhs
    if(isGreaterEqual) {
        A[m][x_index] = -1;
        b[m] = -rhs;
    } else {
        A[m][x_index] = 1;
        b[m] = rhs;
    }
    A[m][n] = 1;

    #ifdef SIMPLEX_DEBUG
    std::cout << "[DEBUG] ####### New sub simplex initialized with new constraint. " << "x" << x_index ;
    if(isGreaterEqual) {
        std::cout << " >= ";
    } else {
        std::cout << " <= ";
    }
    std::cout << rhs << " Now m: " << m << " n: " << n  << std::endl;
    printStatus();
    #endif

    // 对原问题进行初等行变换，保持原来的列基
    for(size_t i = 1; i /*<=*/ < m; i++) { // note that we do not consider the last row
        if(B_I[i] != 0) {
            restoreNewConstraintRow(B_I[i], i);
            #ifdef SIMPLEX_DEBUG
            std::cout << "[DEBUG] After restoring the new constraint row with B_I[" << i << "] = " << B_I[i] << std::endl;
            printStatus();
            #endif
        } else {
            // something that i have not considered
            assert(0);
        }
    }


    #ifdef SIMPLEX_DEBUG
    std::cout << "[DEBUG] ####### After restoring the basis. m: " << m << " n: " << n  << std::endl;
    printStatus();
    // output the A matrix
    // std::cout << "[DEBUG] A matrix (include 0): " << std::endl;
    // for(size_t i = 0; i <= m; i++) {
    //     for(size_t j = 0; j <= n; j++) {
    //         std::cout << A[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // // output the b
    // std::cout << std::endl << "[DEBUG] b (include 0): ";
    // for(size_t i = 0; i < b.size(); i++) {
    //     std::cout << b[i] << " ";
    // }
    // std::cout << std::endl;
    // // output the c
    // std::cout << std::endl << "[DEBUG] c (include 0): ";
    // for(size_t i = 0; i < c.size(); i++) {
    //     std::cout << c[i] << " ";
    // }
    // std::cout << std::endl;
    // // output the B_I
    // std::cout << "[DEBUG] B_I after adding new constraint: ";
    // for(size_t i = 0; i <= m; i++) {
    //     std::cout << B_I[i] << " ";
    // }
    // std::cout << std::endl;
    #endif

    // 如果有负数，则重新计算初始解，方法是拓展x0，构造辅助线性规划问题，利用initialize时的处理方法
    if(b[m] < 0) {
        #ifdef SIMPLEX_DEBUG
        std::cout << "[DEBUG] ####### The new constraint has negtive thing. Starts to handle by solving L_aux" << std::endl;
        #endif
        transformNegtiveB(m);
        
        #ifdef SIMPLEX_DEBUG
        printStatus();
        #endif

    }
    // solve it 直接求解
    solverLoop();
    calculateX();
    #ifdef SIMPLEX_DEBUG
    std::cout << "[DEBUG] ####### After solving the new sub simplex" << std::endl;
    printStatus();
    printSolution();
    #endif


}

// solverLoop
void Simplex::solverLoop() {
    while (!isOptimal()) {
        #ifdef SIMPLEX_DEBUG
        std::cout << "[DEBUG] ########## Start a new iteration #######" << std::endl;
        #endif
        int enter = selectEnteringVariable();
        int leave = selectLeavingVariable(enter);
        if(unbound) {
            return;
        }
        pivot(enter, leave);
    }
}

void Simplex::solve() {
    solveContinous();
}

void Simplex::solveContinous() {
    initializeSimplex();
    if(infesible) {
        return;
    }
    solverLoop();
    if(infesible) {
        return;
    }
    calculateX();
}

bool Simplex::xsAreIntegers() {
    for(size_t i = 1; i <= original_var_num; i++) {
        if(!isClose(x[i], std::round(x[i]))) {
            #ifdef SIMPLEX_DEBUG
            std::cout << "[DEBUG] x[" << i << "] = " << x[i] << " is not an integer." << std::endl;
            #endif
            return false;
        }
    }
    #ifdef SIMPLEX_DEBUG
    std::cout << "[DEBUG] All x are integers." << std::endl;
    #endif
    return true;
}

void Simplex::findFirstFractionalX(int& x_index, double& x_value) {
    // if(!isClose(x[2], std::round(x[2]))) {
    //     x_index = 2;
    //     x_value = x[2];
    //     #ifdef SIMPLEX_DEBUG
    //     std::cout << "[DEBUG] Find the first fractional x: x[" << x_index << "] = " << x_value << std::endl;
    //     #endif
    //     return;
    // }
    for(size_t i = 1; i <= original_var_num; i++) {
        if(!isClose(x[i], std::round(x[i]))) {
            x_index = i;
            x_value = x[i];
            #ifdef SIMPLEX_DEBUG
            std::cout << "[DEBUG] Find the first fractional x: x[" << x_index << "] = " << x_value << std::endl;
            #endif
            return;
        }
    }
}

void Simplex::solveILP() {
    solveILP_BB();
}

void Simplex::solveILP_BB() {
    solveContinous();
    if(infesible) return;
    if(unbound) return;
    calculateX();
    if(xsAreIntegers()) {
        return;
    } else {
        // branch and bound
        // a list contains the Simplex objects
        std::list<std::unique_ptr<Simplex>> simplexList;
        simplexList.emplace_back(std::make_unique<Simplex>(*this));
        // 初始化best_so_far的值为无穷大和best_so_far的解为空
        // 注意，我们的代码中 -z 才表示最优解（最小值）
        double best_so_far = std::numeric_limits<double>::max();
        std::vector<double> best_so_far_x;
        while (!simplexList.empty()) {
            #ifdef SIMPLEX_DEBUG
            std::cout << "[DEBUG] Try to process the problem list with a size of: " << simplexList.size() << std::endl;
            #endif
            // pop the first element
            std::unique_ptr<Simplex> currentSimplex = std::move(simplexList.front());
            simplexList.pop_front();
            #ifdef SIMPLEX_DEBUG
            std::cout << "[DEBUG] After poping, the problem list size is: " << simplexList.size() << std::endl;
            #endif
            // already solved
            // currentSimplex->solverLoop();
            if (currentSimplex->isInfeasible()) {
                #ifdef SIMPLEX_DEBUG
                std::cout << "[DEBUG][BB info] The current Simplex instance is infeasible. Skip it." << std::endl;
                #endif
                continue;
            }
            if (currentSimplex->unbound) {
                #ifdef SIMPLEX_DEBUG
                std::cout << "[DEBUG][BB info] The current Simplex instance is unbounded. Skip it." << std::endl;
                #endif
                continue;
            }
            #ifdef SIMPLEX_DEBUG
            std::cout << "[DEBUG] The current Simplex instance's integer lower bound is " << -currentSimplex->z << std::endl;
            currentSimplex->printSolution();
            #endif
            // 和best_so_far比较，如果比best_so_far差，则不用继续走了
            double newOptimalValue = -currentSimplex->z;
            if(newOptimalValue > best_so_far) {
                #ifdef SIMPLEX_DEBUG
                std::cout << "[DEBUG] The current Simplex instance's optimal value is larger than the best so far. Skip it." << std::endl;
                #endif
                continue;
            }
            if (currentSimplex->xsAreIntegers()) {
                double currentOptimalValue = -currentSimplex->z;
                if (currentOptimalValue < best_so_far) {
                    best_so_far = currentOptimalValue;
                    best_so_far_x = currentSimplex->x;
                }
                continue;
            } else {
                // 分支操作，生成新的Simplex实例并加入列表
                // 寻找第一个不是整数的x
                #ifdef SIMPLEX_DEBUG
                std::cout << "[DEBUG] Starts to branch" << std::endl;
                #endif
                int x_index = -1;
                double x_value;
                currentSimplex->findFirstFractionalX(x_index, x_value);
                #ifdef SIMPLEX_DEBUG
                std::cout << "[DEBUG] Branching on x[" << x_index << "] = " << x_value << std::endl;
                std::cout << "[DEBUG] We add two new Simplex instances with x[" << x_index << "] <= " << std::floor(x_value) << " and x[" << x_index << "] >= " << std::ceil(x_value) << std::endl;
                #endif
                simplexList.emplace_back(std::make_unique<Simplex>(*currentSimplex, x_index, true, std::ceil(x_value)));
                simplexList.emplace_back(std::make_unique<Simplex>(*currentSimplex, x_index, false, std::floor(x_value)));
            }
        }
        if(best_so_far == std::numeric_limits<double>::max()) {
            infesible = true;
            #ifdef SIMPLEX_DEBUG
            std::cout << "[DEBUG] The integer linear programming problem is infeasible because no integer solution is found." << std::endl;
            #endif
            return;
        } else {
            #ifdef SIMPLEX_DEBUG
            std::cout << "[DEBUG] The optimal value of the integer linear programming problem is " << best_so_far << std::endl;
            std::cout << "[DEBUG] The optimal solution of the integer linear programming problem is:" << std::endl;
            for (size_t i = 1; i <= original_var_num; i++) {
                std::cout << "x" << i << " = " << best_so_far_x[i] << std::endl;
            }
            #endif
            z = -best_so_far;
            x = best_so_far_x;
        }

    }
}

void Simplex::printSolution() {
    if (infesible) {
        std::cout << "The linear programming problem is infeasible." << std::endl;
        return;
    }
    if (unbound) {
        std::cout << "The linear programming problem is unbounded." << std::endl;
        return;
    }
    std::cout << "Optimal value: " << -z << std::endl;
    std::cout << "Optimal solution:" << std::endl;
    #ifdef SIMPLEX_DEBUG
    std::cout << "inner x: ";
    for (double value : x) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
    #endif
    for (size_t i = 1; i <= original_var_num; i++) {
        std::cout << "x" << i << " = " << x[i] << std::endl;
    }
}

void Simplex::printOptimalValue() {
    if (infesible) {
        std::cout << "The linear programming problem is infeasible." << std::endl;
        return;
    }
    if (unbound) {
        std::cout << "The linear programming problem is unbounded." << std::endl;
        return;
    }
    std::cout << -z << std::endl;
}

void Simplex::fixBIRemoveWithReallocBase(int toFixIndex) {
    // 此时出现的情况是，L_aux是成立的，x0为0，但是B_I[m] = 0，表示这一列也是基本解。这种情况下，为了获得一个解，我继续加入松弛变量，用于让tabular继续满足原来的结构。
    // 我对于这一行，增加一个松弛变量，让它的系数为1，基础值为对应的b
    // 这种情况下（L_aux最优值为0，原问题是feasible的，x0正好取0），对应的rhs的b的这个分量一定为0
    // 因此，可以在这一行任意选取一个非0的系数（它一定不和其他基的位置重合），让这一行乘以它的倒数，然后pivot它。
    assert(isClose(b[toFixIndex], 0));
    // find the first non-zero element
    size_t nonZeroIndex = 0;
    for(size_t i = 1; i <= n; i++) {
        if(!isClose(A[toFixIndex][i], 0)) {
            nonZeroIndex = i;
            break;
        }
    }
    if(nonZeroIndex == 0) {
        // something wrong
        assert(0);
    }
    // pivot
    pivot(nonZeroIndex, toFixIndex);
    #ifdef SIMPLEX_DEBUG
    std::cout << "[DEBUG] After fixing B_I[" << toFixIndex << "] = 0, we have" << std::endl;
    printStatus();
    #endif
}

void Simplex::fixBIRemoveWithNewSlackVariable(int toFixIndex) {
    
    #ifdef SIMPLEX_DEBUG
    std::cout << "[DEBUG] Fix B_I[" << toFixIndex << "] = 0" << std::endl;
    #endif
    // 向第i行的末尾增加一个系数为1的松弛变量，同时调整A、c、B_I
    for (size_t i = 0; i <= m; i++) {
        A[i].push_back(0);
    }
    n++;
    A[toFixIndex][n] = 1;
    // 提升x的个数
    x.push_back(0);
    // 调整c
    c.push_back(0);
    c_old.push_back(0);
    // 调整B_I
    B_I[toFixIndex] = n;
    #ifdef SIMPLEX_DEBUG
    std::cout << "[DEBUG] After fixing B_I[" << toFixIndex << "] = 0, we have" << std::endl;
    printStatus();
    #endif
}

void Simplex::transformNegtiveB(int theMinimumRowIndex) {
    // we will introduce x0, change the objective function to minimize x0
    // add x0 to the objective function
    c_old = c;
    z_old = z;
    z = 0;
    c = std::vector<double>(n+1, 0);
    c[0] = 1.0;
    // add x0 to the constraints
    for(size_t i = 1; i <= m; i++) {
        A[i][0] = -1;
    }
    #ifdef SIMPLEX_DEBUG
    // output the A matrix
    std::cout << "[DEBUG] ### After adding x0:" << std::endl;
    printStatus();
    #endif
    // pivot to make the basic feasible solution feasible
    #ifdef SIMPLEX_DEBUG
    std::cout << "Pivoting to make the basic feasible solution feasible" << std::endl;
    // print enter and leave
    std::cout << "Enter: " << 0 << " Leave: " << theMinimumRowIndex << std::endl;
    #endif
    pivot(0, theMinimumRowIndex);
    solverLoop();
    #ifdef SIMPLEX_DEBUG
    std::cout << "########### After solving L_aux, we have: " << std::endl;
    printStatus();
    #endif
    if(isClose(z, 0)) {
        // remove x0 from the solution and restore the original objective function
        // calculateX();
        #ifdef SIMPLEX_DEBUG
        std::cout << "The initial basic feasible solution is feasible. Starts to solve original LP" << std::endl;
        #endif
        c = c_old;
        z = z_old;
        // remove x0 from the constraints
        for(size_t i = 1; i <= m; i++) {
            A[i][0] = 0;
        }
        #ifdef SIMPLEX_DEBUG
        // output the c
        std::cout << std::endl << "[DEBUG] c restored to original: ";
        for(size_t i = 0; i < c.size(); i++) {
            std::cout << c[i] << " ";
        }
        std::cout << std::endl;
        printStatus();
        #endif
        // for all i in B_I, if B_I[i] != 0, pivot with restoreCZRow
        for(size_t i = 1; i <= m; i++) {
            if(B_I[i] != 0) {
                restoreCZRow(B_I[i], i);
            } else {
                // something that i have not considered
                // assert(0);
                // this row is dead now
                // 此时出现的情况是，L_aux是成立的，x0为0，但是B_I[m] = 0，表示这一列也是基本解。这种情况下，为了获得一个解，我继续加入松弛变量，用于让tabular继续满足原来的结构。
                // 我对于这一行，增加一个松弛变量，让它的系数为1，基础值为对应的b
                // 这种情况下（L_aux最优值为0，原问题是feasible的，x0正好取0），对应的rhs的b的这个分量一定为0
                // 因此，可以在这一行任意选取一个非0的系数（它一定不和其他基的位置重合），让这一行乘以它的倒数，然后pivot它。
                // fixBIRemoveWithNewSlackVariable(i);
                fixBIRemoveWithReallocBase(i);
                #ifdef SIMPLEX_DEBUG
                std::cout << "[DEBUG] Something wrong with B_I[" << i << "] = 0" << std::endl;
                #endif
            }
        }
        #ifdef SIMPLEX_DEBUG
        // output the A matrix
        std::cout << "[DEBUG] After removing x0 and repivot the first row of c and z " << std::endl;
        printStatus();
        // std::cout << "[DEBUG] A matrix:" << std::endl;
        // for(size_t i = 0; i <= m; i++) {
        //     for(size_t j = 0; j <= n; j++) {
        //         std::cout << A[i][j] << " ";
        //     }
        //     std::cout << std::endl;
        // }
        // // output the c
        // std::cout << std::endl << "[DEBUG] c after removing x0: ";
        // for(size_t i = 0; i < c.size(); i++) {
        //     std::cout << c[i] << " ";
        // }
        // std::cout << std::endl;
        // // output the z
        // std::cout << std::endl << "[DEBUG] z after removing x0: " << z << std::endl;
        // // output the b
        // std::cout << std::endl << "[DEBUG] b after removing x0: ";
        // for(size_t i = 0; i < b.size(); i++) {
        //     std::cout << b[i] << " ";
        // }
        // std::cout << std::endl;
        #endif
    } else {
        infesible = true;
        #ifdef SIMPLEX_DEBUG
        std::cout << "The linear programming problem is infeasible." << std::endl;
        #endif
    }
}

void Simplex::initializeSimplex() {
    B_I.resize(m+1);
    for(size_t i = 1; i <= m; i++) {
        B_I[i] = n - m + i;
    }
    #ifdef SIMPLEX_DEBUG
    // output the B_I
    std::cout << "[DEBUG] B_I (include 0): ";
    for(size_t i = 0; i <= m; i++) {
        std::cout << B_I[i] << " ";
    }
    std::cout << std::endl;
    #endif
    // let l be the index of the minimum b_i
    auto minElementIter = std::min_element(b.begin()+1, b.end());
    int l = std::distance(b.begin(), minElementIter);
    #ifdef SIMPLEX_DEBUG
    // print the minimum element and its index
    std::cout << "[DEBUG] min element: " << *minElementIter << " index: " << l << std::endl;
    #endif
    if(*minElementIter >= 0) {
        // no need to pivot
        return;
    } else {
        // std::cout << "The initial basic feasible solution is not feasible. The algorithm will terminate." << std::endl;
        // assert(0);
        transformNegtiveB(l);
    }
}

void Simplex::calculateX() {
    x.resize(n+1);
    for(size_t j = 1; j <= n; j++) {
        // if j is not in B_I, then x[j] = 0
        auto iter = std::find(B_I.begin() + 1, B_I.end(), j);
        if(iter == B_I.end()) {
            x[j] = 0;
            #ifdef SIMPLEX_DEBUG
            std::cout << "[DEBUG] Set x[" << j << "] = 0" << std::endl;
            #endif
        } else {
            int k = std::distance(B_I.begin(), iter);
            // x[j] = b[i];
            for (size_t i = 1; i <= m; i++) {
                if(isClose(A[i][j], 1)) {
                    x[j] = b[i];
                    #ifdef SIMPLEX_DEBUG
                    std::cout << "[DEBUG] Set x[" << j << "] = " << "b[" << i << "]" << " = " << x[j] << std::endl; 
                    #endif
                    break;
                }
            }
        }
    }

    #ifdef SIMPLEX_DEBUG
    // output the x
    std::cout << "[DEBUG] x (include 0): ";
    for(size_t i = 0; i <= n; i++) {
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;
    #endif
}

void Simplex::restoreCZRow(int enter, int leave) {
    // scale and add to c_hat and z
    double factor = c[enter];
    for(size_t j = 0; j <= n; j++) {
        c[j] -= factor * A[leave][j];
    }
    z -= factor * b[leave];
}

void Simplex::restoreNewConstraintRow(int enter, int leave) {
    // scale and add to the last row
    double factor = A[m][enter];
    for(size_t j = 0; j <= n; j++) {
        A[m][j] -= factor * A[leave][j];
    }
    b[m] -= factor * b[leave];
}

/**
 * enter: the index of the entering variable，正好是列号，取值范围是[1, n]（如果有x0的时候就还包括0）
 * leave: the index of leaving row，正好是行号，取值范围是[1, m]
 * 
 * 进行初等行变换，在A中更改base向量，然后更新B_I矩阵
 */
void Simplex::pivot(int enter, int leave) {
    // 进行初等行变换，使得A[leave][enter] = 1，而其他的A[*][enter] = 0
    // scale leave行
    double scalefactor = A[leave][enter];
    for(size_t j = 0; j <= n; j++) {
        A[leave][j] /= scalefactor;
    }
    b[leave] /= scalefactor;

    // scale其他行
    for(size_t i = 1; i <= m; i++) {
        if(i == leave) {
            continue;
        }
        double factor = A[i][enter];
        for(size_t j = 0; j <= n; j++) {
            A[i][j] -= factor * A[leave][j];
        }
        b[i] -= factor * b[leave];
    }

    // scale and add to c_hat and z
    double factor = c[enter];
    for(size_t j = 0; j <= n; j++) {
        c[j] -= factor * A[leave][j];
    }
    z -= factor * b[leave];

    // update B_I
    B_I[leave] = enter;
    #ifdef SIMPLEX_DEBUG
    // output the A matrix
    printStatus();
    #endif   
}

int Simplex::selectEnteringVariable() {
    // Implementation of selecting entering variable
    // if (c[2] < 0) {
    //     #ifdef SIMPLEX_DEBUG
    //     std::cout << "[DEBUG] Selecting Entering variable: " << 2 << std::endl;
    //     #endif
    //     return 2;
    // }
    for(size_t i = 0; i <= n; i++) {
        if(c[i] < 0) {
            #ifdef SIMPLEX_DEBUG
            std::cout << "[DEBUG] Selecting Entering variable: " << i << std::endl;
            #endif
            return i;
        }
    }
    assert(0);
}

int Simplex::selectLeavingVariable(int enter) {
    // 选择出基变量：在大于0的lambda行中选择最小的theta
    double minTheta = std::numeric_limits<double>::max();
    int leave = -1;
    for(size_t i = 1; i <= m; i++) {
        if(A[i][enter] > 0) {
            double theta = b[i] / A[i][enter];
            if(theta < minTheta) {
                minTheta = theta;
                leave = i;
            }
        }
    }
    #ifdef SIMPLEX_DEBUG
    std::cout << "[DEBUG] Selecting Leaving variable: " << leave << std::endl;
    #endif
    if(leave == -1) {
        unbound = true;
        return -1;
    }
    return leave;
}

bool Simplex::isOptimal() {
    // check if all the critical values are non-negative
    for(size_t i = 1; i <= n; i++) {
        if(c[i] < 0) {
            return false;
        }
    }
    return true;
}

