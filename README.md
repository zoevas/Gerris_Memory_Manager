# The Gerris flow solver is ported in GPU
This project serves as the culmination of my undergraduate studies: [https://www.e-ce.uth.gr/wp-content/uploads/formidable/Vasiliou_zoi.pdf](https://www.e-ce.uth.gr/wp-content/uploads/formidable/Vasiliou_zoi.pdf) 

Download Gerris flow solver
from http://gerris.dalembert.upmc.fr/ 

Since Gerris is a simulation software, it is slow, so we ported it to the GPU for harvesting its multiple cores.
We used the OpenCL framework for porting it to the GPU.

Linked data structures require many transfers between components
of the system. Since transfers cost due to latency, we have developed
a method minimizing the transfers among the components of the
system. In particular, we have developed our own memory manager
with which the data scattered in the heap are now stored in
segregated heaps in a continuous memory area so as one transfer is
needed. However, the linked data structures have as elements pointers
but since the components of a heterogeneous system have different
architectures, when the pointers are transferred to another
component, they will be invalid. Since we use continuous memory areas for storing the nodes of linked data structures, we can rewrite the
pointers to point to valid address space very easily. As for the same
representation of data among the components, we have developed
functions that enforce the same representation for structures. As
a heterogeneous system, we use one of two components with different
virtual address space bit-width, hence the pointers of the components
have different sizes causing false mapping of the data from one
component to another and other problems. We propose a solution
which faces these problems. As case study, we have used the application Gerris[1][2][3]
simulating fluid flow and its primary data structure is a tree which we
transfer to a graphical processing unit.

## Compile it:
 ```
CPPFLAGS="-I <your path>/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc" LDFLAGS="-lOpenCL" ./configure  --prefix=```
/Gerris/gerrisopencl --with-gts-prefix=<your path>/Gerris/gts
```




## Technical details
 
_myHeader.h_ contains the declarations for the memory manager that will be used in the gerris .c files where the trees will be allocated.
 _ftt.c_ contains the memory manager that converts the quad tree in a consecutive memory, namely an array.
 The memory manager is the function 'gpointer my_malloc_tree(unsigned long size)'.

 Each time my_malloc is called, it returns an address of the contiguous memory region allocated at the beginning.
 
 _solid.c_  is calling _my_malloc_cell_data_
 
 _boundary.c_ is calling many _my_malloc_ functions

![diplwmatikh_parousiash](https://github.com/zoevas/Gerris_Memory_Manager/assets/85183528/3e4bc30b-a3e3-48ed-9fa9-2020b5e35be6)
![diplwmatikh_parousiash2](https://github.com/zoevas/Gerris_Memory_Manager/assets/85183528/0c7df5db-5756-438a-8c4f-b4c757be6c8e)

The kernels running on GPU are contained in _kernel.cl_ file.

## References
[1] Gerris page: http://gfs.sourceforge.net/wiki/index.php/Main_Page

[2] Stephane Popinet: “Gerris: a tree-based adaptive solver for the
incompressible Euler equations in complex geometries”, National
Institute of Water and Atmospheric Research, PO Box 14-901 Kilbirnie,
Wellington, New Zealand

[3]Gerris reference: http://src.gnu-darwin.org/ports/science/gerris/work/gerris-
0.6.0/doc/html/book1.html
