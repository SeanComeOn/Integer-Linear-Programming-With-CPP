// Simplex.h
#ifndef SIMPLEX_HPP
#define SIMPLEX_HPP

#include <vector>

#include <Config.hpp>

class Simplex {
public:
    // 枚举类型
    enum class Relation {
        LessEqual,
        GreaterEqual,
        Equal
    };

    Simplex(const std::vector<std::vector<double>>& A, const std::vector<double>& b, const std::vector<double>& c);

    // 用户接口
    void solve();
    void solveILP();
    bool isInfeasible();
    void printSolution();
    void printOptimalValue();
    double getOptimalValue();


    // 复制构造函数
    Simplex(const Simplex& other);

    // 添加单个x限制的函数
    Simplex(const Simplex& other, int x_index, bool isGreaterEqual, double rhs);

    // 赋值运算符
    Simplex& operator=(const Simplex& other);

private:
    // 内部计算函数
    // 根据当前单纯型问题的A, b, c初始化B_I。find the initial basic feasible solution
    void initializeSimplex();
    // 根据目前的信息计算基本解
    void calculateX(); 
    // pivot操作，指定入基和出基的序号
    void pivot(int enter, int leave);
    void restoreCZRow(int enter, int leave);
    void restoreNewConstraintRow(int enter, int leave);
    // 选择入基变量：找第一个 critical value < 0的变量
    int selectEnteringVariable();
    // 选择出基变量：找最小的theta
    int selectLeavingVariable(int enter);
    // 判断是否达到最优解。当所有的critical value都大于等于0时，达到最优解
    bool isOptimal();
    void solverLoop();
    // 求解连续型线性规划问题
    void solveContinous();
    // 用branch and bound求解整数线性规划问题
    void solveILP_BB();
    // 解都是整数
    bool xsAreIntegers();
    void transformNegtiveB(int theMinimumRowIndex);
    void findFirstFractionalX(int& x_index, double& x_value);
    void fixBIRemoveWithNewSlackVariable(int i); // 向第i行的末尾增加一个系数为1的松弛变量，同时调整A、c、B_I
    void fixBIRemoveWithReallocBase(int i); // 重新分配基变量
    #ifdef SIMPLEX_DEBUG
    void printStatus();
    #endif
    
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> c_old;
    std::vector<int> B_I;
    int m;
    int n;
    int original_var_num; // 初始的变量个数
    double z_old;
    double z;
    std::vector<double> x;
    bool infesible = false;
    bool unbound = false;
};

#endif // SIMPLEX_HPP