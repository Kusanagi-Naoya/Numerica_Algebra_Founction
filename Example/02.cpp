#include"basicc.h"
#include<bits/stdc++.h>
using namespace std;
int main(){
    // 第一章第二道习题
    // 巧妙的算法***
    freopen("02.in","r",stdin);
    _MATRIX _S,_T;
    vector<vector<double> > &S=_S.matrix,&T=_T.matrix;
    double lambda;
    vector<double> b,sum;
    int n;
    cin>>n;
    sum.resize(n,0);
    Mcin(n,_S);
    Mcin(n,_T);
    Vin(n,b);
    cin>>lambda;
    b[n]/=(S[n][n]*T[n][n]-lambda);
    for(int k=n-1;k>=1;k--){
        for(int i=k+1;i<=n;i++) sum[k]+=b[i]*T[k][i];
        sum[k+1]+=b[k+1]*T[k+1][k+1];
        for(int i=k;i<=n;i++) b[k]-=sum[i]*S[k][i];
        b[k]/=(S[k][k]*T[k][k]-lambda);
    }
    Vout(b);
    return 0; 
}