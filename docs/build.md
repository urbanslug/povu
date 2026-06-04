ASAN

-DUSE_SANITIZER=ON -DSANITIZER_OPTIONS="Address;Undefined"

```bash
cmake -Bbuild -DENABLE_ALL_TESTS=OFF -DENABLE_CUDA=OFF -DCMAKE_BUILD_TYPE=Debug -DUSE_SANITIZER="Address;Undefined" -S.
```

```bash
cmake --build build
```
