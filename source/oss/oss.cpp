/* OAJ2_cpp is the function to implement OSS

 */
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

#include <algorithm>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
//#include <time>
#include <ctime>

#include <queue>

using namespace std;

int sgn(double v)
{
  return (v < 0) ? -1 : ((v > 0) ? 1 : 0);
}


std::vector<int> bottom_k(std::vector<int>& x, unsigned int k)
{
    std::vector<int> x2 = x; // save a copy of x
    std::vector<int> ind(k); // save the indexes of the smallest k numbers
    std::nth_element(x.begin(), x.begin() + k - 1, x.end()); // std::greater<double>());
    ind[0]=0;
    for(int ii=0, i=0; i<x.size() && ii<k; i++){
        if(x2[i] <= x[k-1])  ind[ii++] = i; // +1 for R
    }
  return ind;
}

// change bottom_k as C style
std::vector<int> bottom_k(std::vector<double>& x, unsigned int k)
{
    std::vector<double> x2 = x; // save a copy of x
    std::vector<int> ind(k); // save the indexes of the smallest k numbers
    std::nth_element(x.begin(), x.begin() + k - 1, x.end()); // std::greater<double>());
    int ii=0; ind[0]=0;
    for(int i=0; i<x.size() && ii<k; i++){
//        if(x2[i] <= x[k-1])  ind[ii++] = i+1; // +1 for R
        if(x2[i] <= x[k-1])  ind[ii++] = i; // for C
    }
    if(ii != k) cout<<"Warning: bottom_k ii="<<ii<<"!= k="<<k<<"\n";
  return ind;
}

// change top_c as C style
std::vector<int> top_k(std::vector<double>& x, unsigned int k)
{
    std::vector<double> x2 = x; // save a copy of x
    std::vector<int> ind(k); // save the indexes of the largest k numbers
    std::nth_element(x.begin(), x.begin() + k - 1, x.end(),  std::greater<double>());
    int ii=0; ind[0]=0;
    for(int i=0; i<x.size() && ii<k; i++){
//        if(x2[i] >= x[k-1])  ind[ii++] = i+1; // +1 for R
        if(x2[i] >= x[k-1])  ind[ii++] = i; // for C
    }
    if(ii != k) cout<<"Warning: top_k ii="<<ii<<"!= k="<<k<<"\n";
  return ind;
}

// return the new candi with k smallest loss
std::vector<int> bottom_k_loss(std::vector<double>& loss, std::vector<int>& candi, unsigned int k)
{
    if(k>candi.size()) cout<<"Error: k="<<k<<">candi.size="<<candi.size()<<"\n";
    // first prepare a vector with loss[candi]
    std::vector<double> x(candi.size()); //
    for(int i=0; i<x.size(); i++) x[i] = loss[ candi[i] ];
    vector<int> tmp = bottom_k(x, k); // bottom_k may return less than k numbers
    vector<int> candi_new( tmp.size());
    for(int i=0; i<tmp.size(); i++) candi_new[i] =  candi[tmp[i]];
    return candi_new;
}


//
double L2norm(double* x, int p)
{ // return L2Norm of vector x
    double L2=0;
    for(int i=0; i<p; i++) L2 += x[i]*x[i];
    return L2;
}


void updateLoss(double* x, std::vector<double>& L2d, std::vector<double>& loss, std::vector<int>& candi, int iCur, int p, double lambda, double tPow)
{ // loss is updated for all points in candi; other variables are not changed.
 // tPow is fixed at 2
    double *xc=x + iCur *p; // point to x[iCur]
    double loss1=p-L2d[iCur];
    for(int i=0; i<candi.size(); i++){
        int ic=candi[i], delta=0;
        double *xic = x + ic*p;  // point to x[ic]
        // compute the sign part
//        for(int j=0; j<p; j++) delta += (sgn(xic[j]) == sgn(xc[j]));
        for(int j=0; j<p; j++) delta += (xic[j] * xc[j] > 0 ? 1 : 0);
            double loss2 = (loss1-L2d[ic]) + delta;
            loss[ic] += loss2 * loss2;
//        double loss2 = lambda*(loss1-L2d[ic]) + (1-lambda)*delta;
//           loss[ic] += pow(loss2 , tPow);
  //      loss[ic] += pow( lambda *(p-L2d[ic]-L2d[iCur]) + (1-lambda)*delta, tPow);
    }
}

// Selecting k points out of N points via OSS
// x is a vector representing an N*p matrix in row-wise
//
std::vector<int> OAJ2_cpp(double* x, int N, int p, int k, double lambda=0.5, double tPow=2, int print=0)
{ // x[N*p], L2d[N], loss[N] are not changed, k is the subsample size
// the size of candi is changed after each iteration
    std::vector<double> L2d(N);  // initialize a vector of length N to store the L2 norm
    for(int i=0; i<N; i++) L2d[i] = L2norm(x+i*p, p)/2; // half

    std::vector<int> ind(k, 0); // save the indexes of the selected k points
    
    // find the first point which has the largest L2norm
    double L2max=L2d[0];
    for(int i=1; i<N; i++) if(L2d[i] > L2max) { L2max=L2d[i]; ind[0]=i; }
    if(print) cout<<"N="<<N<<" k="<<k<<" Initial: ind[0]="<<ind[0]<<"\n";
    
    // prepare the candidate set with remaining points
    std::vector<int> candi(N-1); // save the indexes of the candidate points
    for(int i=0; i<ind[0]; i++) candi[i]=i;
    for(int i=ind[0];i<N-1; i++) candi[i]=i+1;
    
    // bottom_k does not change the v
//    double elapsed_secs=0;
    std::vector<double> elapsed_secs3(3,0); // time used for updateLoss, select and keep nc points, and select the next point.

    // initialize loss for each point
    std::vector<double> loss(N,0);
    // select point one at a time from the candi set
    for(int i=1; i<k; i++){
        if(N<50) for(int j=0; j<candi.size(); j++) cout<<candi[j]<<" ";

        // the current point xc is x[ind[i-1]]
        // compute the loss for each candi point w.r.t the current point
        clock_t begin = clock();
        updateLoss(x, L2d, loss, candi, ind[i-1], p, lambda, tPow);
        elapsed_secs3[0] += double(clock() - begin) / CLOCKS_PER_SEC;

       double r=log(N)/log(k);
       int nc=floor(N/pow(i,r-1));  //
        if(N >= k*k) nc=floor(N/i);
        if(print) if(i<5 || i>k-3) cout<<"N="<<N<<" k="<<k<<" i="<<i<<" nc="<<nc<< " ncandi="<<candi.size() << "; ";
        if((i>1) && (candi.size()>nc)){ // eliminate points after selecting 2 points
            // select and keep nc points with smallest loss among candi points
            begin = clock();
            candi = bottom_k_loss(loss, candi, nc) ;
            elapsed_secs3[1] += double(clock() - begin) / CLOCKS_PER_SEC;
            if(N<50) for(int j=0; j<candi.size(); j++) cout<<candi[j]<<" ";
        } // if

        // select the point with the smallest loss among nc points and remove it from candi
            begin = clock();
   //         ind[i] = bottom_k_loss(loss, candi,  1)[0] ;
            int j0=0;
            double lmin=loss[candi[0]];
            for(int j=1; j<candi.size(); j++){
               if(loss[candi[j]] < lmin) {
                   lmin = loss[candi[j]];
                   j0=j; // index to min loss
               }
            } // for j
            ind[i] = candi[j0];  // the index to min loss
            // remove ind[i] from candi
            candi.erase(candi.begin()+j0);
            elapsed_secs3[2] += double(clock() - begin) / CLOCKS_PER_SEC;

        if(print) if(i<5 || i>k-3) cout<<" ind["<<i<<"]="<<ind[i]<<"\n";
      }// for i
    if(print){
        cout << "\n Total time used for updateLoss, select and keep nc points, and select the next point;  Elapssed_secs=" << elapsed_secs3[0]+elapsed_secs3[1]+elapsed_secs3[2]<<";";
        for(int i=0; i<3; i++) cout<< elapsed_secs3[i]<<",";
        cout<<"\n";
    }

      return ind;
}


/*
// [[Rcpp::export]]
arma::uvec iboss_cpp(arma::mat x, int k) {
  arma::uvec ind=zeros<uvec>(k);
  int m=x.n_cols;
  arma::vec r = floor(k/2/m)*ones<vec>(m);
  if(accu(r)<k/2) r(span(0,((k-2*sum(r))/2)-1))=r(span(0,((k-2*sum(r))/2)-1))+1;
  arma::vec xi=x.col(0);
  ind(span(0, r(0)-1))=as<arma::uvec>(top_k(xi,r(0)));
  ind(span(r(0), 2*r(0)-1))=as<arma::uvec>(bottom_k(xi,r(0)));
  for(int i=1; i<m; i++){
    xi=x.col(i);
    arma::uvec tt=nonzeros(ind);
    int ntt=tt.n_elem;
    for(int j=0; j<ntt; j++){
      xi(tt(j)-1)=median(xi);
      //Rcout << xi(tt(j)-1) << std::endl;
    }
    ind(span( 2*accu( r(span(0,i-1)) ), 2*accu( r(span(0,i-1)) )+r(i)-1 ) )=as<arma::uvec>(top_k(xi,r(i)));
    ind(span( 2*accu( r(span(0,i-1)) )+r(i), 2*accu( r(span(0,i)) )-1 ) )=as<arma::uvec>(bottom_k(xi,r(i)));
  }
  
  return(ind);
}
*/

extern "C"{ // called by R
void oaj2x(double* x, int* ind, int& N, int& p, int& k, double& lambda, double& tPow)
    {
    std::vector<int> tmp= OAJ2_cpp(x, N, p, k, lambda, tPow);
    for(int i=0; i<tmp.size(); i++) ind[i] = tmp[i];
    }

// Selecting k points out of N points via iBOSS
// x is a vector representing an N*p matrix in row-wise
//
void iboss(double* x, int* ind, int& N, int& p, int& k)
{ // x[N*p] stored in row-wise
// for each column, pick the min and max k/(2p) points
// assume k/(2p) is an integer and without checking duplicated points
    int ij=0, k2p= k/2/p;
    for(int j=0;j<p; j++){
 //       double *xj = x+N*j; // column-wise
        vector<double> xv(N); // column j
        for(int i=0; i<N; i++) xv[i] = x[i*p+j]; // row-wise
        vector<int> tmp= bottom_k(xv, k2p);
        for(int i=0; i<tmp.size(); i++) ind[ij++] = tmp[i];
        tmp=top_k(xv, k2p);
        for(int i=0; i<tmp.size(); i++) ind[ij++] = tmp[i];
    }
    if(ij != k) cout<<"Warning: ij="<<ij<<" is different from k="<<k<<"\n";
}

} // extern "C"


int main(int argc, char*argv[])
{ // ./nth N [k [p]]
    // ./nth 11 4 2
    int N=11;
    if(argc>1) N=atoi(argv[1]); // number of points
    int p0 = 1;
    if(argc>3) p0=atoi(argv[3]); // dimension p

    using namespace std;

    srand(static_cast<unsigned int>(time(NULL)));

    //    std::vector<int> v{5, 6, 4, 3, 2, 6, 7, 9, 3};
    int M=1000;
    std::vector<int> v(N*p0); // matrix of N*p0
    if(N < 50)
        std::generate(v.begin(), v.end(), [](){ return std::rand() % 1000;});
    else
        std::generate(v.begin(), v.end(), std::rand);

    std::vector<int> v2= v;  // save a copy
    int k0 = v.size()/2;
    if(argc>2) k0=atoi(argv[2]); //

    std::cout << "Picking the smallest " << k0 << " numbers from " <<N << " numbers with length " << p0 <<"\n";
    if(N<50){
        std::cout << " Before calling std::nth_element: ";
        for(int i=0; i < v.size(); i++) std::cout << v[i] <<" ";
    }
    clock_t begin = clock();
    std::nth_element(v.begin(), v.begin() + k0-1, v.end());
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    if(N<50){
        std::cout << "\n after calling std::nth_element: ";
        for(int i=0; i < v.size(); i++) std::cout << v[i] <<" ";
    }
    
    std::cout << "\n The " << k0 << "th smallest element is " << v[k0-1] << '\n';
    std::cout << " std::nth_element Elapssed_secs=" << elapsed_secs;

    v=v2;
    // bottom_k does not change the v
    begin = clock();
    std::vector<int> ind=bottom_k(v, k0);
    elapsed_secs = double(clock() - begin) / CLOCKS_PER_SEC;
    
//     std::cout << "\n The " << k0 << "th smallest element is " << v[ind[k0-1]-1] << '\n';  // incorrect as v is not sorted
     std::cout << "\n bottom_k Elapssed_secs=" << elapsed_secs;
     if(N<50){
     std::cout << "\n The indexes of the " << k0 << "th smallest elements are ";
        for(int i=0; i < k0; i++) std::cout << ind[i] <<" ";
         std::cout << "\n The " << k0 << " smallest elements are ";
            for(int i=0; i < k0; i++) std::cout << v[ind[i]] <<" ";
    }


    // OAJ2_cpp does not change the v
   double *v3 = new double [N*p0]; // matrix of N*p0
    for(int i=0; i <  v2.size(); i++) v3[i]=2.0*(v2[i] % M)/(double)(M) -1.0;  // -1 <=v3 <=1
    std::cout << "\n calling OAJ2_cpp with " << N << " numbers of length p="  <<p0 <<"\n";
    if(N<50){
      std::cout << " Before calling OAJ2_cpp: ";
      for(int i=0; i < v2.size(); i++) std::cout << v3[i] <<" ";
    }

    begin = clock();
    double lambda =0.5;
    ind=OAJ2_cpp(v3, N, p0, k0, lambda, 2);
    elapsed_secs = double(clock() - begin) / CLOCKS_PER_SEC;
    
     std::cout << "\n OAJ2_cpp Elapssed_secs=" << elapsed_secs <<"\n";
     if(N<50){
     std::cout << "\n The indexes of the " << k0 << "th smallest elements are ";
        for(int i=0; i < k0; i++) std::cout << ind[i] <<" ";
         std::cout << "\n The " << k0 << "smallest elements are ";
            for(int i=0; i < k0; i++) std::cout << v3[ind[i]] <<" ";
    }
    delete [] v3;

}

