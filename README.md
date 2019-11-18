

A manifold-based approach to solve systems of equations
=======================================================

The source code of the algorithm discussed in the academic paper:  
  
Ali Baharev, Arnold Neumaier, Hermann Schichl  
**[A manifold-based approach to sparse global constraint satisfaction problems](https://doi.org/10.1007/s10898-019-00805-x)**  
Journal of Global Optimization, 2019, 75, 949-971

**The code is a research prototype**, and its sole purpose is to 
document the algorithms of the above academic paper.

  - **Algorithm 1** of the paper corresponds to the `cascading_solve` function in `main.py`.
  - **Algorithm 2** of the paper corresponds to the `repair_bnds` function in `main.py`.
  - **Algorithm 3** of the paper corresponds to the `backsolve` function in `main.py`.


Installing and executing the algorithm
--------------------------------------

Please keep in mind that **the code is a research prototype**, and not 
an industrial-strength solver or a commercial product.  
I do not want to turn the code into a platform-independent Python 
package mainly due to lack of time and partly due to licensing issues: 
The third-party dependency `VA27` has to be downloaded and installed 
by the end user separately, making the creation of a self-contained 
Python package impossible anyway.

**The code was written and tested on Linux with Python 3.** 
The platform specific parts are most likely limited to file and 
directory paths only. As it currently is, it won't work on Windows 
without the necessary adjustments, and chances are, running it on Mac 
requires small adjustments too.

**Make sure you don't have anything under `/tmp/stack_machine/` and 
`/tmp/plotter/` because running the code deletes these directories,** 
and the output of the algorithm will be saved there.

**You need the third-party solver VA27 from the Harwell Subroutine 
Library (HSL) to execute the code. It is YOUR RESPONSIBILITY to 
ensure that you are entitled to download and use this third party 
package.**

VA27  Minimize a sum of squares, derivatives required, Marquardt method  
Download from: [HSL Archive](http://www.hsl.rl.ac.uk/archive/)

The source code must be saved under `va27/va27d.f` and 
`va27/ddeps.f`. Adjust the `va27/compile.sh` script to match your 
needs (debug build / release build), your software and hardware 
environment, then compile the VA27 code.

**I patched the VA27 library for significantly better performance.**
Due to licensing issues, I cannot make the patch public. Please 
contact me in e-mail privately if you need this patch.

**Compile the C codes `perturb.c` and `subsampling2.c`**. You find on 
the top of each source file in comments how I compiled it on my 
machine; adjust that to match your needs, your software and hardware 
environment.

The code was not adopted to reflect recent *breaking changes* to the 
`networkx` library, and uses an older version of the `OrderedDiGraph`. 
You need a particular, earlier version of the `networkx` library from 
GitHub. Execute the following lines in the `nx/` directory:

```bash
git clone https://github.com/networkx/networkx.git
cd networkx/
git checkout e43a7b08bd5ab2640b5a9c3350ed5355bdb82c65
```

You find these lines in the `nx/README.txt` too. The `nx.py` module 
inserts this version of `networkx` to the `sys.path`. If you are 
interested why this hack was necessary, here is the explanation:

[Shouldn't OrderedGraph.subgraph() preserve edge order?](https://github.com/networkx/networkx/issues/2048)

Having said that, it is most likely not too difficult to implement a 
function in the latest version of `networkx` that creates the subgraph 
in the desired way, and to port the entire codebase to the latest 
version of `networkx`. I simply do not have the time to do it.

Install the necessary Python packages: `cffi`, `decorator`, `numpy`, 
`six`, `toolz`. The code most likely won't work in Python 2; it has 
been developed and tested under Python 3. On my machine I execute:

```bash
export PATH=/home/ali/miniconda3/bin:$PATH
conda create -n myenvname python=3.6 cffi decorator numpy six toolz
```

**Adjust the `SHELL_SCRIPT` string template in `c_subprobs.py` whether 
you want a debug or a release build.**

Finally, you should be able to run the code like this:

```bash
source activate myenvname 
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:."
python main.py 
```

If I make a debug build with the address sanitizers enabled, I execute
the code like this:

```bash
export PATH=/home/ali/miniconda3/bin:$PATH
export ASAN_OPTIONS=symbolize=1:detect_leaks=0
export ASAN_SYMBOLIZER_PATH=/usr/lib/llvm-3.8/bin/llvm-symbolizer
export LD_PRELOAD=/usr/lib/clang/3.8/lib/linux/libclang_rt.asan-x86_64.so 
source activate myenvname 
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:."
python main.py
```

Adjust the paths to match your environment.

When running the code, you will see warnings like this:

```
mss60_B_61.c: In function ‘solve_61’:
mss60_B_61.c:625:19: warning: division by zero [-Wdiv-by-zero]
             if (SS/N < tolerance) {
                   ^
```

The warning comes from a piece of dead code, that cannot be executed 
but the compiler did not recognize it. I did not have time to fix the 
Python code such that this piece of dead code won't even be generated. 
For the time being, just ignore the warning.

Another, very common warning is this one:

```
 MA10 HAS BEEN GIVEN A MATRIX THAT IS NOT POSITIVE DEFINITE
```

That is normal. The underlying problem at the given point yields a 
matrix that is not positive definite. VA27 handles that properly, and 
will try again with a smaller step length. You can safely ignore this 
warning.

The following warning:

```
main.py:637: RuntimeWarning: invalid value encountered in greater
  violated  = error > 0.0
main.py:638: RuntimeWarning: invalid value encountered in less
  small_viol= error < TOL_TRY
```

can also be safely ignored. The invalid values are intentionally IEEE 
754 NaN (not a number), and the Python code behaves correctly when it 
encounters NaNs.
