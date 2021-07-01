Rcpp::cppFunction('List mvmcombC(List ns, List ps) { // multiple `vmcombC`
  int lns = ns.size();
  List out(lns);
  
  for (int i = 0; i < lns; i++) {
    out[i] = vmcombC(ns[i], ps[i]); // variable modification combination
  }
  
  return out;
}',includes='List vmcombC(DataFrame ns, DataFrame ps) {
  int lns = ns.ncol();
  int lps = ps.ncol();
  int n = lns * lps;
  List np(n);
  
  CharacterVector x; 
  
  int k = 0;
  for (int i = 0; i < lns; i++) {
    x = ns[i];
    
    for (int j = 0; j < lps; j++) {
      x.names() = ps[j];
      CharacterVector y = clone(x);
      np[k] = y;
      k++;
    }
  }
  
  return np;
}')




Rcpp::cppFunction('List par_distC(IntegerVector cols, List mat) {
  
  IntegerVector range_mat = Range(cols[0]-1, mat.size()-1);
  mat = mat[range_mat];

  int len_m = mat.size();
  int len_c = cols.size();
  
  List out(len_c);
  
  for (int i = 0; i < len_c; i++) {
    if (i % 1000 == 0) Rcpp::checkUserInterrupt();

    LogicalVector y = as<LogicalVector>(mat[i]);
    
    IntegerVector js = Range(i, len_m-1);
    IntegerVector outj(js.size());
    
    int k = 0;
    
    for (int j = i; j < len_m; j++) {
      LogicalVector matj = mat[j];
      LogicalVector ij = matj & y;
      outj[k] = sum(ij);

      k++;
    }

    out[i] = outj; 
  }
  
  return out;
}')


