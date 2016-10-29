#ifndef DYARRAY1D_H
#define DYARRAY1D_H

#include <iostream>
#include <cstring>
#include <vector>

using namespace std;
template <class T>
class DyArray1DStruct{
public:
    int   width;
    T * array1D;
    int isGPU;


    //CPU constructor 
    DyArray1DStruct():width(0), array1D(NULL), isGPU(0){
    };

    DyArray1DStruct( std::vector<T>* v){
        isGPU = 0;
        vector<T>& Vec = *v;
        width = Vec.size();
        array1D = new T [width];
        for(int i = 0 ; i<width ; i++)    array1D[i] = Vec[i];
    }

    DyArray1DStruct(T* copy, int copywidth){
        isGPU = 0 ;
        width = copywidth;
        array1D = new T [width];

        memcpy(array1D, copy, sizeof(T)*copywidth);
    }

    DyArray1DStruct(const DyArray1DStruct& copy ){
        isGPU = 0;
        width = copy.width;
        array1D = new T [copy.width];

        memcpy(array1D, copy.array1D, sizeof(T)*copy.width);
    }


    //GPU constructor
    DyArray1DStruct(int alloclength, int GPUORNOT){
        if(GPUORNOT == 1){
            isGPU = 1;
            width = alloclength;
            cudaError_t error = cudaMalloc((void**)& array1D, sizeof(T)*alloclength);
//            cout << cudaGetErrorString(error) << endl;
            cudaMemset( array1D, 0 ,sizeof(T)*alloclength);
        }
        else    cout << "Can not allocate GPU memory\n";
        
    }

    DyArray1DStruct(T* CPUarray, int CPUwidth, int GPUORNOT){
       
        if(GPUORNOT == 1){
            isGPU = 1;
            width = CPUwidth;
            cudaError_t error = cudaMalloc((void**)& array1D, sizeof(T)*CPUwidth);
//            cout << cudaGetErrorString(error) << endl;

            cudaMemcpy(array1D, CPUarray, sizeof(T)*CPUwidth, cudaMemcpyHostToDevice);

        }
        else    cout << "Can not allocate GPU memory\n";
    }

    DyArray1DStruct(const DyArray1DStruct& CPUarray, int GPUORNOT ){
        if(GPUORNOT==1){
            isGPU = 1;
            width = CPUarray.width;
            cudaError_t error = cudaMalloc((void**) &array1D, sizeof(T)*width);
//            cout << cudaGetErrorString(error) << endl;
            cudaMemcpy(array1D, CPUarray.array1D, sizeof(T)*width, cudaMemcpyHostToDevice);
        }
        else   cout << "Can not allocate GPU memcpy\n";
    }


   //destructor
    ~DyArray1DStruct(){
        if(isGPU == 0){
            if(array1D != NULL){
                 delete [] array1D;
                 array1D = NULL;
            }
        }
        else{
            if(array1D != NULL){
                cudaError_t error = cudaFree(array1D);
                //cout << "FFFFFFFFF: " << cudaGetErrorString(error) << endl;
                array1D = NULL;    
            }
        }
    }

    void allocateCPU(int alloclength){
        isGPU = 0;
        width = alloclength;
        array1D = new T [alloclength];
        memset(array1D, 0, sizeof(T)*alloclength);
    }


    void allocateGPU(int alloclength, int GPUORNOT){
        if(GPUORNOT == 1){
            isGPU = 1;
            width = alloclength;
            cudaError_t error = cudaMalloc((void**)& array1D, sizeof(T)*alloclength);
  //          cout << cudaGetErrorString(error) << endl;
            cudaMemset( array1D, 0 ,sizeof(T)*alloclength);
        }
        else    cout << "Can not allocate GPU memory\n";
    }

    //CPU TO CPU
    void copy(T* dest, int count){
        if(isGPU==0){
            memcpy(dest, array1D, sizeof(T)*count);
        }
        else    cout << "this is CPU copy function\n";
    }

    //GPU TO CPU 
    void copy(T* dest, int count, int GPUORNOT){
        if(isGPU==1){
            if(GPUORNOT == 1){
                cudaMemcpy(dest, array1D, sizeof(T)*count, cudaMemcpyDeviceToHost);
            }
            else cout << "Can not copy  GPU memory to CPU" ;
        }
        else  cout << "this is GPU copy function\n";
    }

    //CPU TO GPU
    void copyToGPU(T* dest, int count){
        if(isGPU==1){
            cudaMemcpy(array1D, dest, sizeof(T)*count, cudaMemcpyHostToDevice);
        }
    }
    
    //GPU TO GPU
    void copyToGPU(T* dest, int count, int GPUORNOT){
        if(isGPU==1 && GPUORNOT ==1){
            cudaMemcpy(dest, array1D, sizeof(T)*count, cudaMemcpyDeviceToDevice);
        }
    }
    
    void print(){
        for(int i = 0 ; i<width ; i++){
            cout << array1D[i] << "\t";
        }
    }
};


template <class T>
class DyArray1D{
public:
    DyArray1DStruct<T>* arrayPtr;
    DyArray1DStruct<T> tempCPUForCopy;
    int isGPU;

    //CPU constructor
    DyArray1D():isGPU(0), arrayPtr(0){
    }

    DyArray1D(DyArray1DStruct<T>* copy){
        isGPU = 0;
        arrayPtr = new DyArray1DStruct<T>;
        memcpy(arrayPtr, copy, sizeof(DyArray1DStruct<T>));
    }

    //GPU constructor
    DyArray1D(DyArray1DStruct<T>* copy, int GPUORNOT){
        if(GPUORNOT == 1){
            isGPU = 1;
            cudaError_t error = cudaMalloc((void**) &arrayPtr, sizeof(DyArray1DStruct<T>));
//            cout << cudaGetErrorString(error) << endl;
            cudaMemcpy(arrayPtr, copy, sizeof(DyArray1DStruct<T>), cudaMemcpyHostToDevice);

            memcpy(&tempCPUForCopy, copy, sizeof(DyArray1DStruct<T>));

            copy->array1D = NULL;
        }
        else   cout << "Can not allocate GPU memcpy\n";
    }


    void allocateGPU(DyArray1DStruct<T>* copy, int GPUORNOT){
         if(GPUORNOT == 1){
            isGPU = 1;
            cudaError_t error = cudaMalloc((void**) &arrayPtr, sizeof(DyArray1DStruct<T>));
//            cout << cudaGetErrorString(error) << endl;
            cudaMemcpy(arrayPtr, copy, sizeof(DyArray1DStruct<T>), cudaMemcpyHostToDevice);

            memcpy(&tempCPUForCopy, copy, sizeof(DyArray1DStruct<T>));
        }
        else   cout << "Can not allocate GPU memcpy\n";
    }

    //CPU TO CPU
    void copy(T* dest, int count){
        arrayPtr->copy(dest, count);
    }

    //GPU TO CPU
    void copy(T* dest, int count, int GPUORNOT){
        if(isGPU==1)
            tempCPUForCopy.copy(dest, count, GPUORNOT);
    }

    //CPU TO GPU
    void copyToGPU(T* dest, int count){
        tempCPUForCopy.copyToGPU(dest, count);
    }

    //GPU TO GPU
    void copyToGPU(T* dest, int count, int GPUORNOT){
        tempCPUForCopy.copyToGPU(dest, count, GPUORNOT);
    }
    
    //destructor
    ~DyArray1D(){
        if(isGPU == 1){
            if( arrayPtr != NULL){
                cudaError_t error = cudaFree(arrayPtr);
//                cout << "DDDDDDDDDDD: "<< cudaGetErrorString(error)<<endl;
            }
        }
        else{
            if(arrayPtr != NULL)    delete arrayPtr;
        }
        arrayPtr = NULL;    
    }
};

#endif

