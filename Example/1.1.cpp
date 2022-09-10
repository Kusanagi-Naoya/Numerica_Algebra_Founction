#include<bits/stdc++.h>
#include<unistd.h>
#include"basicc.h"
using namespace std;
int n=84;
string sg[4]={"","Gauss消去法","列主元Guass消去法","全主元Guass消去法"};
int main(){
    freopen("1.1.out","w",stdout);
    struct timespec time1={0, 0}; 
    struct timespec time2={0, 0};
    float temp;
    _MATRIX A;
    vector<double> b;
    A.init(n);
    b.resize(n+1,1);
    for(int i=1;i<=n;i++){
        A.matrix[i][i]=6;
        if(i!=1) A.matrix[i][i-1]=8;
        if(i!=n) A.matrix[i][i+1]=1;
    }
    b=A*b;
    for(int i=1;i<=3;i++){
        cout<<endl<<sg[i]<<':'<<endl;
        clock_gettime(CLOCK_REALTIME, &time1);  
        Vout(Linear_Solve_Guass(A,b,i));
        clock_gettime(CLOCK_REALTIME, &time2);   
        temp=(time2.tv_nsec-time1.tv_nsec)/1000000.0;
        cout<<endl<<"Time : "<<temp<<"ms"<<endl<<endl;
    }
    return 0;
}