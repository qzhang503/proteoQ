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


