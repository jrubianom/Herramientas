#include <iostream>
using namespace std;
int main()
{
    int n[4];
    for(int i=0;i<4;i++){n[i]=i;}
    int *p = n;
    for(int i =0;i<4;i++){
    cout<<n[5]<<'\t'<<*(p+2*i)<<"\n";

    }

    return 0;
}
