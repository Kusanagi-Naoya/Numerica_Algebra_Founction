#include<bits/stdc++.h>
#include"basicc.h"
using namespace std;
int n=10;
string sg[4]={"","Gauss消去法","列主元Guass消去法","全主元Guass消去法"};
string s[2]={"平方根法计算Cholesky分解","改进的平方根法"};
int main(){
    //freopen("1.3.1.out","w",stdout);
    struct timespec time1={0, 0}; 
    struct timespec time2={0, 0};
    float temp;
    srand(time(NULL));
    _MATRIX A;
    vector<double> b;
    b.resize(n+1,0);
    A.init(n);
    for(int i=1;i<=n;i++){
        b[i]=rand()%20+1;
        if(i!=1) A.matrix[i][i-1]=1;
        if(i!=n) A.matrix[i][i+1]=1;
        A.matrix[i][i]=10;
    }
    cout<<"b:";
    Vout(b);
    cout<<endl;
    int i=2;
        cout<<endl<<sg[i]<<':'<<endl;
        clock_gettime(CLOCK_REALTIME, &time1);  
        Vout(Linear_Solve_Guass(A,b,i));
        clock_gettime(CLOCK_REALTIME, &time2);   
        temp=(time2.tv_nsec-time1.tv_nsec)/1000000.0;
        cout<<endl<<"Time : "<<temp<<"ms"<<endl<<endl;
    i=3;
        cout<<endl<<sg[i]<<':'<<endl;
        clock_gettime(CLOCK_REALTIME, &time1);  
        Vout(Linear_Solve_Guass(A,b,i));
        clock_gettime(CLOCK_REALTIME, &time2);   
        temp=(time2.tv_nsec-time1.tv_nsec)/1000000.0;
        cout<<endl<<"Time : "<<temp<<"ms"<<endl<<endl;
    return 0;
}