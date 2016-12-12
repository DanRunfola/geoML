#include <Rcpp.h>
#include <cstdlib>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <string>
#include <cctype>
#include <omp.h>
  // Enable C++11 via this plugin (Rcpp 0.10.3 or later)
  // [[Rcpp::plugins(cpp11)]]

  using namespace Rcpp;


// [[Rcpp::export]]
NumericVector cv(NumericMatrix fit, NumericVector index){
  int n = index.size();
  std::vector<int> nodes;
  std::vector<int> nodesid;
  for(int i=0; i<n; i++){
    if(fit(i,0) > 1){
      nodes.push_back(index[i]);
      nodesid.push_back(i);
    }
  }
  //std::cout<<"nodes#"<<nodes.size()<<std::endl;
  std::unordered_map<int,int> idx_ht;
  for(int i=0; i<index.size(); i++){
    idx_ht[index[i]] = 1;
  }
  std::unordered_map<int,double> err_ht;

  for(int i=0; i<index.size(); i++){
    err_ht[index[i]] = fit(i,3);
  }
  std::vector<double> node_error;
  std::vector<int> n_child;
  std::vector<int> indexcopy;
  for(int i=0; i<index.size(); i++){
    indexcopy.push_back(index[i]);
  }
  std::sort(indexcopy.begin(),indexcopy.end());

  for(int i=index.size()-1; i>=0; i--){
    //std::cout<<indexcopy[i]<<std::endl;
    if( idx_ht[2*indexcopy[i]] > 0 ){
      idx_ht[indexcopy[i]] = idx_ht[2*indexcopy[i]] + idx_ht[2*indexcopy[i]+1];
      err_ht[indexcopy[i]] = err_ht[2*indexcopy[i]] + err_ht[2*indexcopy[i]+1];
    }
  }


  for(int i=0; i<nodes.size(); i++){
    //std::cout<<idx_ht[nodes[i]]<<std::endl;
    n_child.push_back(idx_ht[nodes[i]]);
    node_error.push_back(err_ht[nodes[i]]);
  }

  double minv = 0.0;;
  int minid = 0;
  double delta;
  for(int i=0; i<nodes.size(); i++){
    if(nodes[i] != 1){
      delta = ( fit(nodesid[i],3) - node_error[i] ) / (n_child[i] - 1);
      //std::cout<<i<<":"<<delta<<std::endl;
      //std::cout<<i<<"nodeerr:"<<fit(nodesid[i],3)<<std::endl;
      //std::cout<<i<<"childreerr:"<<node_error[i]<<std::endl;
      //if(minv == 0.0){minv = delta;}

      if(delta < minv || minv == 0.0){
        minv = delta;
        minid = nodesid[i];
      }
    }

  }
  std::vector<double> res;
  //res.push_back(minv);
  res.push_back(index[minid]);
  //std::cout<<minid<<std::endl;
  res.push_back(minid);
  std::unordered_map<int, int> mp;
  for(int i=0; i<index.size(); i++){
    mp[index[i]] = i;
  }

  //std::vector<int> res_ids(0,n);
  std::vector<int> res_ids;
  for(int i=0; i<n; i++){
    res_ids.push_back(0);
  }
  std::queue<int> q;
  q.push(minid);
  while(!q.empty()){
    int temp = q.front();
    int id = index[q.front()];
    q.pop();
    if(mp[2*id] > 0){
      q.push(mp[2*id]);
      q.push(mp[2*id+1]);
    }
    res_ids[temp] = 1;
  }

  for(int i=0; i<n; i++){
    if(res_ids[i] == 0){
      res.push_back(i);
    }
  }
  NumericVector out(res.size());
  for(int i=0; i<res.size(); i++){
    out[i] = res[i];
  }
  return out;

}


// [[Rcpp::export]]
NumericVector cross_validate(NumericMatrix fit, NumericVector index, NumericVector alphalist){
  NumericVector alphalistm = clone(alphalist);
  NumericVector res;
  res = cv(fit, index);
  double alpha;
  alpha = res[0];
  int nodeid;
  nodeid = res[1];
  std::cout<<"nodeid"<<nodeid<<"res"<<res.size()<<std::endl;
  std::vector<int> rowid;
  if(nodeid == 0){
    return alphalistm;
  }
  alphalistm.push_back(alpha);

  //removed subtree node as leaf
  fit(nodeid,0) = 1;
  for(int i=1; i<res.size(); i++){
    rowid.push_back(res[i]);
  }

  NumericMatrix fitm(rowid.size(),8);
  std::sort(rowid.begin(),rowid.end());
  for(int i=0; i<rowid.size(); i++){
    fitm(i,_) = fit(rowid[i],_);
  }

  NumericVector indexm(rowid.size());
  for(int i=0; i<rowid.size(); i++){
    indexm[i] = index[rowid[i]];
  }
  //return alphalistm;
  return cross_validate(fitm,indexm,alphalistm);

}


// [[Rcpp::export]]
NumericVector splitnc(NumericMatrix y,NumericVector x){
  std::vector<int> arr(x.size(),0);
#pragma omp parallel for
  for(int i=0; i<x.size(); i++){
    //std::cout<<x[i]<<std::endl;
    arr[i] = x[i];
  }
  std::sort(arr.begin(),arr.end());
  std::vector<int> ux;
  ux.push_back(arr[0]);
  for(int i=1; i<arr.size(); i++){
    if(arr[i] != arr[i-1])
      ux.push_back(arr[i]);
  }

  int n = ux.size();
  NumericVector out((2*n-1),0.0);
  int len = y.nrow();

  std::unordered_map<int,int> ht_cnt;
  std::unordered_map<int,int> ht_val;
  std::vector<int> mean(n,0);
  std::vector<int> goodness(n-1,0);

  // for(int i=0;i<len;i++){
  //   for(int j=0;j<n;j++){
  //     if(x[i] == ux[j]){
  //       count[j] = count[j] + 1;
  //       sum[j] = sum[j] + y(i,0);
  //       break;
  //     }
  //   }
  // }
  for(int i=0; i<len; i++){
    ht_cnt[arr[i]]++;
    ht_val[arr[i]] += y(i,3);
  }

#pragma omp parallel for
  for(int i=0;i<n;i++){
    mean[i] = ht_val[ux[i]] / ht_cnt[ux[i]];
  }

  //rank X

  std::vector<int> idx(n,0);
  for(int i=0; i<n; i++){
    idx[i] = i;
  }

  std::sort(idx.begin(),idx.end(), [&mean](int a,int b) {return mean[a] < mean[b] ;});

  // for(int i=0;i<n-1;i++){
  //   for(int j=i+1;j<n;j++){
  //     if(mean[i]>mean[j]){
  //       double tmp = mean[j];
  //       mean[j] = mean[i];
  //       mean[i] = tmp;
  //       double temp = ux[j];
  //       ux[j] = ux[i];
  //       ux[i] = temp;
  //     }
  //   }
  //
  // }
  std::vector<int> uxx;
  for(int i=0; i<n; i++){
    uxx.push_back(ux[idx[i]]);
  }
  std::unordered_map<int,int> ht_ux;
  for(int i=0; i<n; i++){
    ux[i] = uxx[i];
    ht_ux[ux[i]] = i;
  }

#pragma omp parallel for
  for(int j=0;j<n-1;j++){
    double rss = 0.0;
    double wmeanleft = 0.0;
    double wmeanright = 0.0;
    double sumTrtleft = 0.0;
    double sumUntrtleft = 0.0;
    double sumTrtWtleft = 0.0;
    double sumUntrtWtleft = 0.0;
    double sumTrtright = 0.0;
    double sumUntrtright = 0.0;
    double sumTrtWtright = 0.0;
    double sumUntrtWtright = 0.0;
    int pos;
#pragma omp parallel for
    for(int i=0;i<len;i++){
      pos = ht_ux[x[i]];
      if(pos<=j){
        if(y(i,1) == 1.0){
          sumTrtleft = sumTrtleft + y(i,0)/y(i,2);
          sumTrtWtleft = sumTrtWtleft + 1/y(i,2);
        }
        else{
          sumUntrtleft = sumUntrtleft + y(i,0)/(1-y(i,2));
          sumUntrtWtleft = sumUntrtWtleft + 1/(1-y(i,2));
        }
      }
      else{
        if(y(i,1) == 1){
          sumTrtright = sumTrtright + y(i,0)/y(i,2);
          sumTrtWtright = sumTrtWtright + 1/y(i,2);
        }
        else{
          sumUntrtright = sumUntrtright + y(i,0)/(1-y(i,2));
          sumUntrtWtright = sumUntrtWtright + 1/(1-y(i,2));
        }
      }
    }
    if(sumTrtWtleft != 0 && sumUntrtWtleft !=0 ){
      wmeanleft =  sumTrtleft/sumTrtWtleft - sumUntrtleft/sumUntrtWtleft;
      for(int i=0;i<len;i++){
        pos = ht_ux[x[i]];
        if(pos<=j){
          rss = rss + (y(i,3) - wmeanleft)*(y(i,3) - wmeanleft);
        }
      }
    }
    else{
      goodness[j] = 0;
      continue;
    }

    if(sumTrtWtright != 0 && sumUntrtWtright !=0 ){
      wmeanright =  -sumTrtright/sumTrtWtright + sumUntrtright/sumUntrtWtright;
      for(int i=0;i<len;i++){
        pos = ht_ux[x[i]];
        if(pos>j){
          rss = rss + (y(i,3) + wmeanright)*(y(i,3) + wmeanright);
        }
      }
      goodness[j] =  1/rss ;
    }
    else{
      goodness[j] = 0;
    }
  }
#pragma omp parallel for
  for(int i=0;i<n-1;i++){
    out[i] = goodness[i];
  }
#pragma omp parallel for
  for(int i=n-1;i<=2*(n-1);i++){
    out[i] = ux[i-n+1];
  }
  return out;

}

// [[Rcpp::export]]
NumericVector node_evaluate(NumericMatrix y){
  double sumTrt = 0.0;
  double sumUntrt = 0.0;
  double sumTrtWt = 0.0;
  double sumUntrtWt = 0.0;
  int nrow = y.nrow();
  double wmean = 0.0;
  double rss = 0.0;
  int treat = 0;
  int untreat = 0;
  for(int i=0;i<nrow;i++){
    if(y(i,1) == 1.0){
      sumTrt = sumTrt + y(i,0)/y(i,2);
      sumTrtWt = sumTrtWt + 1/y(i,2);
      treat++;
    }
    else{
      sumUntrt = sumUntrt + y(i,0)/(1-y(i,2));
      sumUntrtWt = sumUntrtWt + 1/(1-y(i,2));
      untreat++;
    }
  }
  if(sumTrtWt != 0 && sumUntrtWt !=0 ){
    wmean =  sumTrt/sumTrtWt - sumUntrt/sumUntrtWt;

  }
  double ymean = 0;
  for(int i=0;i<nrow;i++){
    ymean += y(i,3);
  }
  ymean = ymean/nrow;
  //std::cout<<"wmean:"<<wmean<<"treat:"<<treat<<"untrt"<<untreat<<std::endl;
  for(int i=0;i<nrow;i++){
    //rss = rss + (y(i,3) - wmean)* (y(i,3) - wmean);
    rss = rss + (y(i,3) - ymean)* (y(i,3) - ymean);
  }
  NumericVector out(2);
  out[0] = wmean;
  out[1] = rss;
  return out;
}

// [[Rcpp::export]]
NumericVector splitc(NumericMatrix y){
  int n = y.nrow();
  NumericVector out(2*(n-1),0.0);
  NumericVector goodness(n-1);
  NumericVector direction(n-1);
  int max;
  max=omp_get_max_threads();
  omp_set_num_threads(max);
#pragma omp parallel for
  for(int j=0;j<n-1;j++){
    double rss = 0.0;
    double wmeanleft = 0.0;
    double wmeanright = 0.0;
    double sumTrtleft = 0.0;
    double sumUntrtleft = 0.0;
    double sumTrtWtleft = 0.0;
    double sumUntrtWtleft = 0.0;
    double sumTrtright = 0.0;
    double sumUntrtright = 0.0;
    double sumTrtWtright = 0.0;
    double sumUntrtWtright = 0.0;
    // left child
    for(int i=0;i<=j;i++) {
      if(y(i,1) == 1.0){
        sumTrtleft = sumTrtleft + y(i,0)/y(i,2);
        sumTrtWtleft = sumTrtWtleft + 1/y(i,2);
      }
      else{
        sumUntrtleft = sumUntrtleft + y(i,0)/(1-y(i,2));
        sumUntrtWtleft = sumUntrtWtleft + 1/(1-y(i,2));
      }
    }

    if(sumTrtWtleft != 0 && sumUntrtWtleft !=0 ){
      wmeanleft =  sumTrtleft/sumTrtWtleft - sumUntrtleft/sumUntrtWtleft;
      for(int i=0;i<=j;i++) {
        rss = rss + (y(i,3) - wmeanleft)*(y(i,3) - wmeanleft);
      }
    }
    else{
      goodness[j] = 0;
      continue;
    }
    // right child
    for(int i=j+1;i<n;i++){
      if(y(i,1) == 1){
        sumTrtright = sumTrtright + y(i,0)/y(i,2);
        sumTrtWtright = sumTrtWtright + 1/y(i,2);
      }
      else{
        sumUntrtright = sumUntrtright + y(i,0)/(1-y(i,2));
        sumUntrtWtright = sumUntrtWtright + 1/(1-y(i,2));
      }
    }
    if(sumTrtWtright != 0 && sumUntrtWtright !=0 ){
      wmeanright =  -sumTrtright/sumTrtWtright + sumUntrtright/sumUntrtWtright;
      for(int i=j+1;i<n;i++){
        rss = rss + (y(i,3) + wmeanright)*(y(i,3) + wmeanright);
      }
      goodness[j] =  1/rss ;
      direction[j] = wmeanleft>0 ? 1:-1;
    }
    else{
      goodness[j] = 0;
    }
  }

  for(int i=0;i<n-1;i++){
    out[i] = goodness[i];
  }
  for(int i=n-1;i<2*(n-1);i++){
    out[i] = direction[i-n+1];
  }
  return out;

}
