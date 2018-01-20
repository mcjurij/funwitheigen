A DFT and a FFT. It's just me playing around.

Build and use instructions in the respective source file.

The dft.cpp code is translated from DFT Python code as shown in 
[Diskrete Fouriertransformation (DFT) / Abtasttheorem (Shannon, Nyquist, Kotelnikow)](https://www.youtube.com/watch?v=sX-DNi_SX-Q) by Prof. Weitz.
With kind permission from Prof Weitz / HAW Hamburg, Germany.

You need the [Eigen](http://eigen.tuxfamily.org/) library for dft.cpp. After unpacking the tarball, rename the directory to eigen (to make -Ieigen on the g++ command line work).

The FFT is in simple_fft.cpp. The output of dft.cpp and simple_fft.cpp ought to be exactly the same in case
the length of the input vector is a power of two.

To see the output as a graph (see saw17.png and fft_saw17.png) you need a plotter like gnuplot.  
