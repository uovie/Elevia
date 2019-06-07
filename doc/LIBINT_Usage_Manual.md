# LIBINT Usage Manual

## 1. To initialize the integral evaluator

   ```c++
   void libint2_init_eri(Libint_eri_t* libint, int max_am, LIBINT2_REALTYPE* buf);
   ```

The meaning of three arguments:

* the pointer to the evaluator
* the maximum angular momentum of basis functions
* optional pointer to the scratch buffer
  * If the third argument is 0, then the call will dynamically allocate the needed space.


## 2. Evaluate ERIs

```c++
void (*libint2_build_eri[N][N][N][N])(Libint_eri_t *);
```

## 3. Cleanup

```c++
void libint2_cleanup_eri(Libint_eri_t* libint);
```

