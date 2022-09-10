#include<bits/stdc++.h>
#include"basicc.h"
using namespace std;
int n=40;
string s[2]={"平方根法计算Cholesky分解","改进的平方根法"};
int main(){
    freopen("1.2.2.out","w",stdout);
    struct timespec time1={0, 0}; 
    struct timespec time2={0, 0};
    float temp;
    srand(time(NULL));
    _MATRIX A;
    vector<double> b;
    b.resize(n+1,1);
    b[0]=0;
    A.init(n);
    for(int i=1;i<=n;i++){
        for(int j=1;j<=n;j++) A.matrix[i][j]=1.0/(i+j-1);
    }
    b=A*b;
    for(int i=1;i<=2;i++){
        cout<<endl<<s[i-1]<<':'<<endl;
        clock_gettime(CLOCK_REALTIME, &time1); 
        Vout(Linear_Solve_Cholesky(A,b,i));
        clock_gettime(CLOCK_REALTIME, &time2);   
        temp=(time2.tv_nsec-time1.tv_nsec)/1000000.0;
        cout<<endl<<"Time : "<<temp<<"ms"<<endl<<endl;
    }
    return 0;
}