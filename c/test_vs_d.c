#include "common.h"
#include "frequentDirections.h"
#include "sparseSketcher.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void test_vs_d(){
  int n = 10000;
  int dim_set[] = {1000,2000,3000,4000,5000,6000};
  int k = 10;
  int exp_no = 6;
  int nnz = 100;
  int ell = 50;

  double start, end, cpu_time_used;
  SparseMatrix A;
  SparseVector arr[n];
  SparseSketcher sfd;
  FrequentDirections fd;

  double sfd_cov_err[exp_no], sfd_proj_err[exp_no], sfd_time[exp_no];
  double fd_cov_err[exp_no], fd_proj_err[exp_no], fd_time[exp_no];

  double tailSquaredFrob;


  for(int i=0; i<exp_no; i++){
    
    printf("i = %d\n", i);
    sfd_time[i] = 0;
    fd_time[i] = 0;

    init_sparseSketcher(&sfd, ell, dim_set[i]);
    init_fd(&fd, ell, dim_set[i]);

    // input matrix
    init_sparseMatrix(&A, dim_set[i], n);

    for (int j=0; j < n; j++){
      skew_init_sparseVector(&arr[j], dim_set[i], nnz, (int) (1.5 * nnz), 0.9);
      append_to_sparseMatrix(&A, &arr[j]);

      // SFD
      start = clock();
      append_to_sparseSketcher(&sfd, &arr[j]);
      sfd_time[i] += (double) (clock() - start);

      //FD
      start = clock();
      append_to_fd(&fd, &arr[j]);
      fd_time[i] += (double) (clock() - start); 
    }

    tailSquaredFrob = topRank(&A, k);

  
    //SFD
    start = clock();
    get_sparseSketch(&sfd);
    sfd_time[i] += (double) (clock() - start);
    sfd_time[i] = sfd_time[i] / CLOCKS_PER_SEC;
    sfd_cov_err[i] = computeRelCovErr(&A, sfd.sketch, ell, dim_set[i]);
    sfd_proj_err[i] = computeRelProjErr(&A, sfd.sketch, ell, dim_set[i], k, tailSquaredFrob);

    
    //FD
    start = clock();
    get_fdSketch(&fd);
    fd_time[i] += (double) (clock() - start);
    fd_time[i] = fd_time[i] / CLOCKS_PER_SEC;
    fd_cov_err[i] = computeRelCovErr(&A, fd.sketch, ell, dim_set[i]);
    fd_proj_err[i] = computeRelProjErr(&A, fd.sketch, ell, dim_set[i], k, tailSquaredFrob);  

  }
    

  printf("SFD:\n");
  print_one_dim_double("\'proj\':", sfd_proj_err, exp_no);
  print_one_dim_double("\'cov\':", sfd_cov_err, exp_no);  
  print_one_dim_double("\'time\':", sfd_time, exp_no);

  print_one_dim_double("\'proj\':", fd_proj_err, exp_no);    
  print_one_dim_double("\'cov\':", fd_cov_err, exp_no);
  print_one_dim_double("\'time\':", fd_time, exp_no);
}
