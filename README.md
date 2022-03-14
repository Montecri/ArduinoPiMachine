# ArduinoPiMachine
Arduino Pi Machine - Calculate Pi to nth Decimal

![ezgif com-gif-maker](https://user-images.githubusercontent.com/38574378/158080737-6ee8792f-4743-4e39-b3c8-c5a76624b83b.gif)

An Arduino Nano assembly equipped with a 4 digits 7-segments led display which will, while powered, keep calculating Pi digits and scrolling through the display, like a window scrolling left to right through the infinitude of Pi digits.

For that task, I’m using Xavier Gourdon algorithm described in his “Computation of the n-th decimal digit of π with low memory” paper and the Xavier Gourdon and Pascal Sebah paper called “N-th digit computation” that comes with a nice C++ implementation example that I ported to Arduino framework

![IMG_20220313_151819025_HDR_menor](https://user-images.githubusercontent.com/38574378/158080774-50462913-55ac-43b5-833d-798b76cd25d4.jpg)


---
<p align="center"><img src="https://user-images.githubusercontent.com/38574378/132773469-08fb7b59-2f9d-4641-9665-c8d50d3904bc.png"><b>   ATTENTION   </b><img src="https://user-images.githubusercontent.com/38574378/132773469-08fb7b59-2f9d-4641-9665-c8d50d3904bc.png"></p> 

When compared to a version running on PC, this code on Arduino Nano has an approximation error on the last 2-3 digits, rendering only the first 3-5 digits useable. Believe this is due to double type on ATMEGA occupying only 4 bytes. Further refinement of the code is needed for ATMEGA.

---

References:

http://numbers.computation.free.fr/Constants/Algorithms/nthdigit.html

http://numbers.computation.free.fr/Constants/Algorithms/nthdecimaldigit.pdf
