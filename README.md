# 说在前面

这里是头文件basicc.h的说明文档，该头文件是中国科学技术大学董哲斌在2022年秋季数值代数课程中所完成的，其中包含数值代数的各种基础运算。

本文档及头文件会随着课程的学习逐步完善，目前为截止到第一章。

__注: 本头文件中为了方便之后的矩阵运算，矩阵和向量的存储均从1开始存储__
# 目录
- [结构体](#结构体)
  - [_MAREIX](#_MAREIX)
      - [元素](#元素)
      - [成员函数](#成员函数)
      - [重载运算符](#重载运算符)
- [第一章的函数](#第一章的函数)
  - [基础函数](#基础函数)
  - [解线性方程组相关运算](#解线性方程组相关运算)
    - [通用函数](#通用函数)
    - [Guass消去法相关算法](#guass消去法相关算法)
    - [Cholesky分解相关算法](#cholesky分解相关算法)
  - [最终解方程组函数](#最终解方程组函数)

# 结构体
## _MATRIX 
_MATRIX是封装出来的矩阵类型，用于进行更方便的矩阵运算。
### 元素
- ```c++
  vector<vector<double> > matrix 
  ```
   用于存储矩阵的元素
- ```c++
  int _x,_y
  ```
  用于存储矩阵的大小。
- 值得注意的是，为了运算方便，本结构体所存储的元素**均从1开始存储**，所以matrix实际的大小为(_x+1)*(_y+1)。
### 成员函数
- ```c++
  void init(int x,int y)
  ``` 
  初始定义为一个大小为x*y的矩阵，其中元素均为0。
- ```c++
  void init(int _n)
  ```
  对init的重构，初始定义为一个大小为n*n的方阵，其中元素均为0。
- ```c++
  void init_identity(int _n)
  ```
  特别的，初始定义一个大小为n*n的单位矩阵
### 重载运算符
- 该头文件重载了矩阵\*矩阵，矩阵\*向量，向量\*向量，向量\*数量，数量\*向量，向量/数量，向量+向量，向量-向量。
- 上述运算均遵循矩阵与向量运算规则。
### 再次强调，本头文件中所有向量以及矩阵均从1开始存储。
# 第一章的函数
### 基础函数
以下是一些无关数值代数的常用的运算
- ```c++
  void Mcin(int _x,int _y,_MATRIX &_A)
  ```
  将一个大小为_x*_y的矩阵输入到_A中
- ```c++
  void Mcin(int _n,_MATRIX &_A)
  ```
  对Mcin的重定义，将一个大小为n*n的矩阵输入到_A中
- ```c++
  void Mout(_MATRIX _A)
  ```
  输出矩阵_A
- ```c++
  void Vin(int _x,vector<double> &_V)
  ```
  将一个长为_x的向量输入到_V中
- ```c++
  void Vout(vector<double> _V)
  ```
  输出向量_V
- ```c++
  _MATRIX Transposition(_MATRIX A)
  ```
  返回A的转置矩阵
## 解线性方程组相关运算
此处说明数值代数中解线性方程组的相关运算，分为Guass消去相关与Cholesky分解相关，具体算法这里不做解释感兴趣的可以上网学习或看源码了解。

注：如下函数并没有实现对于非法输入的报错机制，仅仅有计算功能，所以请进行正确的输入，对于非法的输入产生的结果，后果自负。
### 通用函数
此处说明一些基础的，无关高级算法的运算。
- ```c++
  void Forward_Sweep(_MATRIX _L,vector<double> &b)
  ```
  前代法解系数矩阵为下三角矩阵的线性方程组，解存储在b当中。
- ```c++
  void Forward_Sweep_Identity(_MATRIX _L,vector<double> &b)
  ```
  同样是前代法解线性方程组，特别的是，系数矩阵为**单位下三角矩阵**，也就是**主对角线均为1的下三角矩阵**，解存储在b当中。
- ```c++
  void Backward_Sweep(_MATRIX _U,vector<double> &y)
  ```
  回代法解系数矩阵为上三角矩阵的线性方程组，解存储在y当中。
- ```c++
  void Backward_Sweep_Identity(_MATRIX _U,vector<double> &y)
  ```
  同样是回代法解线性方程组，特别的是，系数矩阵为**单位上三角矩阵**，也就是**主对角线均为1的上三角矩阵**，解存储在y当中。
- ```c++
  void None_Sweep(_MATRIX _D,vector<double> &b)
  ```
  解系数矩阵为对角矩阵的线性方程组，解存储在y当中。
- ```c++
  _MATRIX Get_Triangular_Matrix(_MATRIX _A,int di,int id)
  ```
  从矩阵_A中取出三角部分。

  di体现方向，若di=1，则提取下三角矩阵，若di=0，则提取上三角矩阵。
  
  id体现三角矩阵是否是单位三角矩阵，若id=1，则将对角元素赋为1，若id=0，则不变，提取原矩阵中的元素。
- ```c++
  _MATRIX Get_Diagonal_Matrix(_MATRIX _A)
  ```
  从矩阵_A中取出对角矩阵
- ```c++
  _MATRIX Get_Permutation_Matrix_Row(vector<int> p)
  ```
  计算P=P1\*P2\*...\*Pn,其中Pi为i与p[i]相置换的置换矩阵。
- ```c++
  _MATRIX Get_Permutation_Matrix_Column(vector<int> q)
  ```
  计算Q=Q1\*Q2\*...\*Qn,其中Qi为i与q[i]相置换的置换矩阵。
### Guass消去法相关算法
注意，以下函数的输入中，矩阵均应该是**可以进行三角分解的矩阵**，换而言之，该矩阵的**任意阶顺序主子阵均是非奇异的矩阵**。

如果输入了非法的矩阵，一切后果自负。

但是作为输入的向量并没有要求。
- ```c++
  _MATRIX Guass_Elimination(_MATRIX _A)
  ```
  对_A进行Guass分解，将L与U存在同一个矩阵中。
  
  设返回值存储到A中，使用L=Get_Triangular_Matrix(A,1,1)取出单位下三角矩阵L，
  使用U=Get_Triangular_Matrix(A,0,0)取出上三角U。

  有A=L*U。
- ```c++
  _MATRIX Whole_Major_Guass_Elimination(vector<int> &P,vector<int> &Q,_MATRIX _A)
  ```
  对_A进行全主元Guass分解，将L与U存在同一个矩阵A中返回，同时将P，Q的信息存在向量p，q中返回。
  
  设返回值存储到A中，使用L=Get_Triangular_Matrix(A,1,1)取出单位下三角矩阵L，
  使用U=Get_Triangular_Matrix(A,0,0)取出上三角U。

  同时令P=Get_Permutation_Matrix_Row(p)，Q=Get_Permutation_Matrix_Column(q)。

  有P\*A\*Q=L*U。
- ```c++
  _MATRIX Select_Major_Guass_Elimination(vector<int> &P,_MATRIX _A)
  ```
  对_A进行全主元Guass分解，将L与U存在同一个矩阵A中返回，同时将P的信息存在向量p中返回。
  
  设返回值存储到A中，使用L=Get_Triangular_Matrix(A,1,1)取出单位下三角矩阵L，
  使用U=Get_Triangular_Matrix(A,0,0)取出上三角U。

  同时令P=Get_Permutation_Matrix_Row(p)。

  有P\*A=L*U。
### Cholesky分解相关算法
注意，以下函数的输入中，矩阵均应该是**可以进行Cholesky分解的矩阵**，换而言之，该矩阵应该是**对称正定的矩阵**，该命题等价于该矩阵的**任意阶顺序主子式均为正**。

如果输入了非法的矩阵，一切后果自负。

- ```c++
  _MATRIX Cholesky_Factorization(_MATRIX _A)
  ```
  对_A使用平方根法进行Cholesky分解，将L存在矩阵的下三角部分中。
  
  设返回值存储到A中，使用L=Get_Triangular_Matrix(A,1,0)取出下三角矩阵L，
  使用LT=Transposition(L)取出L的转置矩阵。

  有A=L*LT。
- ```c++
  _MATRIX Cholesky_Factorization_Improved(_MATRIX _A)
  ```
  对_A使用改进的平方根法进行Cholesky分解，将L与D存在同一个矩阵中。
  
  设返回值存储到A中，使用L=Get_Triangular_Matrix(A,1,1)取出单位下三角矩阵L，D=Get_Diagonal_Matrix(A)取出对角矩阵D，
  使用LT=Transposition(L)取出L的转置矩阵。

  有A=L\*D\*LT。
## 最终解方程组函数
- ```c++
  vector<double> Linear_Solve_Guass(_MATRIX _A,vector<double> b,int kind)
  ```
  使用Guass消去法解线性方程组_Ax=b，kind为方法。

  kind=1时，使用普通Guass消去法；
  kind=2时，使用列主元Guass消去法；
  kind=3时，使用全主元Guass消去法。

  在此申明，Guass消去法的系数矩阵需为**可进行三角分解的矩阵**。

- ```c++
  vector<double> Linear_Solve_Cholesky(_MATRIX _A,vector<double> b,int kind)
  ```
  使用Cholesky分解解线性方程组_Ax=b，kind为方法。
  
  kind=1时，使用普通平方根法；
  kind=2时，使用改进的平方根法 。

  在此申明，Cholesky分解的系数矩阵需为**对称正定矩阵**。
# Numerica_Algebra_Founction
