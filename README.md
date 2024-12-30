# Integer Linear Programming - Simplex Method with CPP

This project implements the Simplex algorithm for solving **linear programming problems (LP)** and **integer linear programming (ILP)** problems.

The simplex algorithm is implemented with purely elementary row operations on the simplex tabular, without any matrix operations. Integer linear programming is implemented with the branch-and-boudn method, because cutting-plane method seems to be slow.

This project features:
- pure C++ implementation
- only depends on standard C++ libraries, no 3rd-party library is needed
- C++14 standard

## Project Structure

The project is organized as follows:

```
simplex-project
├── src
│   ├── main.cpp          # Entry point of the application
│   └── Simplex.cpp       # Implementation of the Simplex class
├── include
│   └── Simplex.hpp       # External interface for the Simplex class
├── CMakeLists.txt        # Build configuration file
├── DESIGN.md             # Some design notes in the process
└── README.md             # Project documentation
```

## Overview of the Simplex Algorithm

The Simplex algorithm operates on a linear programming problem in the following **standard form** like:

Minimize: 
```
c^T * x
```
Subject to:
```
Ax ≤ b
x ≥ 0
```

Where:
- `c` is the coefficients vector for the objective function.
- `A` is the matrix of coefficients for the constraints.
- `b` is the vector of constraint bounds.
- `x` is the vector of decision variables.

## Realm of implementation
I implemented:
1. Initialization of the Simplex algorithm with auxiliary LP to handle negtive `b` cases.
1. The Simplex algorithm with elementary row operations.
1. Branch-and-bound method for integer linear programming.

By calling the `printSolution` method, the algorithm can report the solution to the problem, or infeasible, or unbound.

Currently, we have these public methods:
```c++
void solve();
void solveILP();
bool isInfeasible();
void printSolution();
void printOptimalValue();
double getOptimalValue();
```
## Building the Project

To build the project, you need to have CMake installed. Follow these steps:

1. Clone the repository or download the project files.
2. Open a terminal and navigate to the project directory.
3. Create a build directory:
   ```sh
   mkdir build
   cd build
   ```
4. Run CMake to configure the project:
   ```sh
   cmake ..
   ```
5. Build the project:
   ```sh
   make
   ```

## Running the Project

After building the project, you can run the application in the `build` directory using the following command:

```sh
# assume that you are in your build directory
./simplex
```

## Usage Example

Go to `main.cpp` for details. Write your own `main.cpp` to use the Simplex algorithm.

The solution to the problems can be:
1. Some feasible solution.
1. Infeasible.
1. Unbounded.

To use the Simplex algorithm, you will need to provide the coefficients of the objective function and the constraints. The implementation will handle the optimization process and output the optimal solution.

### Linear Programming Example
```cpp
// Example input for the Simplex algorithm
// Objective function: Minimize z = -x1 - 14x2 - 6x3
// Subject to constraints:
// x1 + x2 + x3 <= 4
// x1           <= 2
//           x3 <= 3
//     3x2 + x3 <= 6
// x1, x2, x3   >= 0

// Coefficients of the objective function
std::vector<double> objective = {-1, -14, -6};

// Coefficients of the constraints
std::vector<std::vector<double>> constraints = {
    {1, 1, 1},
    {1, 0, 0},
    {0, 0, 1},
    {0, 3, 1}
};

// Right-hand side values of the constraints
std::vector<double> rhs = {4, 2, 3, 6};

// Create a Simplex object
Simplex simplex(constraints, rhs, objective); // false indicates minimization

// Solve the linear programming problem
simplex.solve();

// Print the optimal solution
simplex.printSolution();

// Print the optimal value of the objective function
simplex.printOptimalValue();
```

### Integer Linear Programming Example
You can use `simplex.solveILP()` to solve integer linear programming problems. The solutions are guaranteed to be integer values.

```cpp
std::vector<double> objective = {0, -1};
std::vector<std::vector<double>> constraints = {
   {-2, 1},
   {2, 1}
};
std::vector<double> rhs = {0, 6};

// Create a Simplex object
Simplex simplex(constraints, rhs, objective); // false indicates minimization

// Solve the integer linear programming problem
simplex.solveILP();

// Output the results
simplex.printSolution();
```

## Tests
It is tested on several test cases including [P3980](https://www.luogu.com.cn/problem/P3980). Although it can output the correct answer, it consumes 1m24s and a lot of memory at test #10 scale problems on a Intel 12900HX CPU on WSL.

## Improvements

The project may improve in the future by:
1. Using a more efficient data structure to store the simplex tabular, such as sparse matrix, to support greater scale problems.
1. Other algorithm improvements.

And I welcome any suggestions and contributions to the project.
