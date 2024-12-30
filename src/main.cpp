#include <iostream>
#include <Simplex.hpp>

int main() {
    // Example input for the Simplex algorithm
    // Objective function: Minimize z = -x1 - 14x2 - 6x3
    // Subject to constraints:
    // x1 + x2 + x3 <= 4
    // x1           <= 2
    //           x3 <= 3
    //     3x2 + x3 <= 6
    // x1, x2, x3   >= 0

    // std::vector<double> objective = {-1, -14, -6};
    // std::vector<std::vector<double>> constraints = {
    //     {1, 1, 1},
    //     {1, 0, 0},
    //     {0, 0, 1},
    //     {0, 3, 1}
    // };
    // std::vector<double> rhs = {4, 2, 3, 6};


    // Coefficients of the objective function
    // Coefficients of the constraints
    // Right-hand side values of the constraints

    // auxilary variable x0
    // std::vector<double> objective = {1, 2};
    // std::vector<std::vector<double>> constraints = {
    //     {-1, -1},
    //     {1, 1},
    // };
    // std::vector<double> rhs = {-1, 2};

    // unbound problem
    // std::vector<double> objective = {-1, -1};
    // std::vector<std::vector<double>> constraints = {
    //     {1, -1},
    //     {-1, 1},
    // };
    // std::vector<double> rhs = {1, 1};

    // not feasible
    // std::vector<double> objective = {1, 2};
    // std::vector<std::vector<double>> constraints = {
    //     {-1, -1},
    //     {1, 1},
    // };
    // std::vector<double> rhs = {-2, -1};

    // std::vector<double> objective = {-5, -3};
    // std::vector<std::vector<double>> constraints = {
    //     {2, 1},
    //     {1, 2},
    // };
    // std::vector<double> rhs = {40, 50};   

// ILP testbench
    // introduction to linear optimization. ILP example
    // std::vector<double> objective = {1, -2};
    // std::vector<std::vector<double>> constraints = {
    //     {-4, 6},
    //     {1, 1},
    // };
    // std::vector<double> rhs = {9, 4};

    // infeasible example
    // std::vector<double> objective = {1, 1};
    // std::vector<std::vector<double>> constraints = {
    //     {-1, 0},
    //     {-2, 3},
    //     {1, -3},
    //     {2, 3}
    // };
    // std::vector<double> rhs = {-1, 0, 0, 6};

    // (1,2)
    std::vector<double> objective = {0, -1};
    std::vector<std::vector<double>> constraints = {
        {-2, 1},
        {2, 1}
    };
    std::vector<double> rhs = {0, 6};
    
    // Create a Simplex object
    Simplex simplex(constraints, rhs, objective); // false indicates minimization

    // Solve the linear programming problem
    // simplex.solve();
    simplex.solveILP();

    // Output the results
    simplex.printSolution();
    // simplex.printOptimalValue();

    return 0;
}