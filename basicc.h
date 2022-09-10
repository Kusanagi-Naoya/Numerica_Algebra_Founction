#ifndef _BASICC_H
#define _BASICC_H

#include<iostream>
#include<vector>
using namespace std;

struct _MATRIX{ //定义矩阵结构体
    vector<vector<double> > matrix;
    int _x,_y;
    void init(int x,int y){
        _x=x,_y=y;
        matrix.resize(x+1,vector<double>(y+1,0.0));
    }
    void init(int _n){ init(_n,_n); }
    void init_identity(int _n){
        init(_n);
        for(int i=1;i<=_n;i++) matrix[i][i]=1;
    }
};

void Mcin(int _x,int _y,_MATRIX &_A){ // 输入矩阵
    _A.init(_x,_y);
    for(int i=1;i<=_x;i++) for(int j=1;j<=_y;j++) cin>>_A.matrix[i][j];
}

void Mcin(int _n,_MATRIX &_A){ Mcin(_n,_n,_A); } // 输入方阵

void Mout(_MATRIX _A){ // 输出矩阵
    int x=_A._x,y=_A._y;
    cout<<endl;
    for(int i=1;i<=x;i++){
        for(int j=1;j<=y;j++) cout<<_A.matrix[i][j]<<' ';
        cout<<endl;
    }
}

void Vin(int _x,vector<double> &_V){ // 输入向量
    _V.resize(_x+1,0);
    for(int i=1;i<=_x;i++){
        double x;
        cin>>_V[i];
    }
}

void Vout(vector<double> _V){ // 输出向量
    int n=_V.size();
    cout<<endl;
    for(int i=1;i<=n-1;i++) cout<<_V[i]<<' ';
    cout<<endl;
}

_MATRIX Transposition(_MATRIX A){
    int n=A._x;
    for(int i=1;i<=n;i++) for(int j=i+1;j<=n;j++) swap(A.matrix[i][j],A.matrix[j][i]);
    return A;
}

_MATRIX operator * (_MATRIX _A,_MATRIX _B){ // 重载矩阵乘法
    _MATRIX _C;
    vector<vector<double> > A,B,C;
    A=_A.matrix,B=_B.matrix;
    int x=_A._x,y=_B._y,z=_A._y;
    _C.init(x,y);
    C=_C.matrix;
    for(int i=1;i<=x;i++) for(int j=1;j<=y;j++) for(int k=1;k<=z;k++) C[i][j]+=A[i][k]*B[k][j];
    _C.matrix=C;
    return _C;
}

vector<double> operator * (_MATRIX _A,vector<double> b){ // 重载矩阵*向量
    int x=_A._x,y=_A._y;
    vector<double> ret;
    ret.push_back(0);
    for(int i=1;i<=x;i++){
        double cur=0;
        for(int j=1;j<=y;j++) cur+=_A.matrix[i][j]*b[j];
        ret.push_back(cur);
    }
    return ret;
}

double operator * (vector<double> a,vector<double> b){ // 重载向量*向量
    int l=a.size();
    double ret=0;
    for(int i=1;i<=l;i++) ret+=a[i]*b[i];
    return ret;
}

vector<double> operator * (vector<double> a,double x){
    int _n=a.size()-1;
    for(int i=1;i<=_n;i++) a[i]*=x;
    return a;
}

vector<double> operator * (double x,vector<double> a){ return a*x; }

vector<double> operator / (vector<double> a,double x){ return a*(1/x); }

vector<double> operator + (vector<double> a,vector<double> b){
    int _n=a.size()-1;
    for(int i=1;i<=_n;i++) a[i]+=b[i];
    return a;
}

vector<double> operator - (vector<double> a,vector<double> b){ return a+(-1*b); }

void Forward_Sweep(_MATRIX _L,vector<double> &b){ // 解Ly=b
    vector<vector<double> > L=_L.matrix;
    int n=_L._x;
    for(int i=1;i<=n;i++){
        b[i]/=L[i][i];
        for(int j=i+1;j<=n;j++) b[j]-=b[i]*L[j][i];
    }
}

void Forward_Sweep_Identity(_MATRIX _L,vector<double> &b){ // 解Ly=b,且L的对角元均为1
    vector<vector<double> > L=_L.matrix;
    int n=b.size()-1;
    for(int i=1;i<=n-1;i++) for(int j=i+1;j<=n;j++) b[j]-=b[i]*L[j][i];
}

void Backward_Sweep(_MATRIX _U,vector<double> &y){ // 解Ux=y
    vector<vector<double> > U=_U.matrix;
    int n=_U._x;
    for(int i=n;i>=1;i--){
        y[i]/=U[i][i];
        for(int j=1;j<i;j++) y[j]-=y[i]*U[j][i];
    }
}

void Backward_Sweep_Identity(_MATRIX _U,vector<double> &y){ // 解Ux=y,且L的对角元均为1
    vector<vector<double> > U=_U.matrix;
    int n=y.size()-1;
    for(int i=n;i>=1;i--) for(int j=1;j<i;j++) y[j]-=y[i]*U[j][i];
}

void None_Sweep(_MATRIX _D,vector<double> &b){
    int n=_D._x;
    for(int i=1;i<=n;i++) b[i]/=_D.matrix[i][i];
}

_MATRIX Get_Triangular_Matrix(_MATRIX _A,int di,int id){ // 从高斯消去法的结果中提取下三角矩阵
    int n=_A._x;
    for(int i=1;i<=n;i++){
        if(id) _A.matrix[i][i]=1;
        for(int j=i+1;j<=n;j++){
            if(di) _A.matrix[i][j]=0;
            else _A.matrix[j][i]=0;
        }
    }
    return _A;
}

_MATRIX Get_Diagonal_Matrix(_MATRIX _A){
    _MATRIX _D;
    int n=_A._x;
    _D.init(n);
    for(int i=1;i<=n;i++) _D.matrix[i][i]=_A.matrix[i][i];
    return _D;
}

_MATRIX Get_Permutation_Matrix_Row(vector<int> p){
    int n=p.size()-1;
    _MATRIX P;
    P.init_identity(n);
    for(int i=1;i<n;i++) swap(P.matrix[i],P.matrix[p[i]]);
    return P;
}

_MATRIX Get_Permutation_Matrix_Column(vector<int> q){
    int n=q.size()-1;
    _MATRIX Q;
    Q.init_identity(n);
    for(int i=1;i<n;i++) for(int j=1;j<=n;j++) swap(Q.matrix[j][i],Q.matrix[j][q[i]]);
    return Q;
}

_MATRIX Guass_Elimination(_MATRIX _A){ // 高斯消去法求三角矩阵，上下三角放在同一矩阵上
    int n=_A._x;
    vector<vector<double> > &A=_A.matrix;
    for(int i=1;i<=n-1;i++){
        for(int j=i+1;j<=n;j++) A[j][i]/=A[i][i];
        for(int j=i+1;j<=n;j++){
            for(int k=i+1;k<=n;k++) A[j][k]-=A[j][i]*A[i][k];
        }
    }
    return _A;
}

_MATRIX Whole_Major_Guass_Elimination(vector<int> &P,vector<int> &Q,_MATRIX _A){ // 列主元法求三角分解，返回一个pair类型，first为上下三角矩阵，second为置换矩阵P。
    int n=_A._x,x,y;
    vector<vector<double> > &A=_A.matrix;
    P.resize(n+1,0);
    Q.resize(n+1,0);
    for(int k=1;k<n;k++){
        x=y=k;
        for(int i=k;i<=n;i++) for(int j=k;j<=n;j++) if(A[i][j]>A[x][y]) x=i,y=j; 
        swap(A[k],A[x]);
        for(int i=1;i<=n;i++) swap(A[i][k],A[i][y]);
        P[k]=x;
        Q[k]=y;
        for(int i=k+1;i<=n;i++) A[i][k]/=A[k][k];
        for(int i=k+1;i<=n;i++){
            for(int j=k+1;j<=n;j++) A[i][j]-=A[i][k]*A[k][j]; 
        }
    }
    return _A;
}

_MATRIX Select_Major_Guass_Elimination(vector<int> &P,_MATRIX _A){ // 列主元法求三角分解，返回一个pair类型，first为上下三角矩阵，second为置换矩阵P。
    int n=_A._x,x=0;
    vector<vector<double> > &A=_A.matrix;
    P.resize(n+1,0);
    for(int k=1;k<=n;k++){
        x=k;
        for(int i=k;i<=n;i++) if(A[i][k]>A[x][k]) x=i; 
        swap(A[k],A[x]);
        P[k]=x;
        for(int i=k+1;i<=n;i++) A[i][k]/=A[k][k];
        for(int i=k+1;i<=n;i++){
            for(int j=k+1;j<=n;j++) A[i][j]-=A[i][k]*A[k][j]; 
        }
    }
    return _A;
}

_MATRIX Cholesky_Factorization(_MATRIX _A){
    int n=_A._x;
    vector<vector<double> > &A=_A.matrix;
    for(int k=1;k<=n;k++){
        A[k][k]=sqrt(A[k][k]);
        for(int i=k+1;i<=n;i++) A[i][k]/=A[k][k];
        for(int i=k+1;i<=n;i++){
            for(int j=i;j<=n;j++) A[j][i]-=A[j][k]*A[i][k];
        }
    }
    return _A;
}

_MATRIX Cholesky_Factorization_Improved(_MATRIX _A){
    int n=_A._x;
    vector<double> v;
    vector<vector<double> > &A=_A.matrix;
    v.resize(n,0);
    for(int k=1;k<=n;k++){
        for(int i=1;i<k;i++) v[i]=A[i][i]*A[k][i];
        for(int i=1;i<k;i++) A[k][k]-=A[k][i]*v[i];
        for(int i=k+1;i<=n;i++){
            for(int j=1;j<=k;j++) A[i][k]-=A[i][j]*v[j];
            A[i][k]/=A[k][k];
        }
    }
    return _A;
}

vector<double> Linear_Solve_Guass(_MATRIX _A,vector<double> b,int kind){
    vector<int> p,q;
    _MATRIX A2,L,U,P,Q;
    if(kind==1) A2=Guass_Elimination(_A);
    if(kind==2){
        A2=Select_Major_Guass_Elimination(p,_A);
        P=Get_Permutation_Matrix_Row(p);
        b=P*b;
    }
    if(kind==3){
        A2=Whole_Major_Guass_Elimination(p,q,_A);
        P=Get_Permutation_Matrix_Row(p);
        b=P*b;
        Q=Get_Permutation_Matrix_Column(q);
    }
    L=Get_Triangular_Matrix(A2,1,1);
    U=Get_Triangular_Matrix(A2,0,0);
    Forward_Sweep_Identity(L,b);
    Backward_Sweep(U,b);
    if(kind==3) b=Q*b;
    return b;
}
// 解线性方程组，kind为方法，1为Guass消去法，2为列主元Guas消去，3为全主元Guass消去

vector<double> Linear_Solve_Cholesky(_MATRIX _A,vector<double> b,int kind){
    _MATRIX A2,L,LT,D;
    if(kind==1){
        A2=Cholesky_Factorization(_A);
        L=Get_Triangular_Matrix(A2,1,0);
        LT=Transposition(L);
        Forward_Sweep(L,b);
        Backward_Sweep(LT,b);
    }
    else{
        A2=Cholesky_Factorization_Improved(_A);
        L=Get_Triangular_Matrix(A2,1,1);
        D=Get_Diagonal_Matrix(A2);
        LT=Transposition(L);
        Forward_Sweep_Identity(L,b);
        None_Sweep(D,b);
        Backward_Sweep_Identity(LT,b);
    }
    return b;
}

/*
void Mcin(int _x,int _y,_MATRIX &_A); // 输入矩阵
void Mcin(int _n,_MATRIX &_A); // 输入方阵
void Mcout(_MATRIX _A);
void Forward_Sweep(_MATRIX _L,vector<double> &b); // 解Ly=b
void Forward_Sweep_Identity(_MATRIX _L,vector<double> &b); // 解Ly=b,且L的对角元均为1
void Backward_Sweep(_MATRIX _U,vector<double> &y); // 解Ux=y
void Backward_Sweep_Identity(_MATRIX _U,vector<double> &y); // 解Ux=y,且L的对角元均为1
*/

#endif