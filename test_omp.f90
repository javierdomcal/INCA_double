program test
  use omp_lib
  print *, omp_get_max_threads()
end program
