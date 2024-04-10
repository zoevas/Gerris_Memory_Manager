# The Gerris flow solver is ported in GPU

Download Gerris flow solver
from http://gerris.dalembert.upmc.fr/ 

## Compile it:
 CPPFLAGS="-I <your path>/NVIDIA_GPU_Computing_SDK/OpenCL/common/inc" LDFLAGS="-lOpenCL" ./configure  --prefix=<your_path>/Gerris/gerrisopencl --with-gts-prefix=/srv/homes/zovasili/Gerris/gts

## Description
 ftt.c contains the memory manager that converts the quad tree in a consecutive memory, namely an array.
 The memory manager is the function 'gpointer my_malloc_tree(unsigned long size)'

![diplwmatikh_parousiash](https://github.com/zoevas/Gerris_Memory_Manager/assets/85183528/3e4bc30b-a3e3-48ed-9fa9-2020b5e35be6)
![diplwmatikh_parousiash2](https://github.com/zoevas/Gerris_Memory_Manager/assets/85183528/0c7df5db-5756-438a-8c4f-b4c757be6c8e)
