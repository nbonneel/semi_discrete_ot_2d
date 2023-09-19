# semi_discrete_ot_2d
Semi-discrete optimal transport code to compute the uniformity of a 2D point set

Sample self-explanatory code in main.cpp. To integrate in your projects, add transport.h and transport.cpp. No external dependencies.

For Mac users: 
```
brew install libomp 
c++ -std=c++11 -O3 -o main main.cpp transport.cpp -lomp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include  -Xpreprocessor -fopenmp
```



![screenshot](screenshot.jpg)
