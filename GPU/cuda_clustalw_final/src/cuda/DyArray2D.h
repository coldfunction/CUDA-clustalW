#ifndef DYARRAY2D_H
#define DYARRAY2D_H

#include <iostream>
#include <cstring>
#include <vector>
using namespace std;

template <class T>
class DyArray2DStruct{
public:
    int   heigth;
    int   width;
    int   size;
    T*    array2D;
    int*  array2DSize;

    int   isGPU;

    // CPU constructor
    DyArray2DStruct():size(0), heigth(0), width(0), isGPU(0), array2D(NULL), array2DSize(NULL){
    };

    DyArray2DStruct( std::vector<std::vector<T> >* v, int maxwidth){
        int max = 0;

        isGPU = 0;
        vector<vector<T> >& Vec = *v;

        heigth = Vec.size();
        array2DSize = new int [heigth];
        for(int i = 0 ; i<heigth ; i++){
            array2DSize[i] = Vec[i].size();
            if(array2DSize[i] > max )    max = array2DSize[i];
        }

        width = max;
        size = heigth*width;
        array2D = new T [size];

        for(int i = 0, index = 0 ; i<Vec.size() ; i++){
            int j=0;
            for(j = 0 ; j < Vec[i].size(); j++){
                array2D[index] = Vec[i][j];
                index++; 
            }
            while(j!=width){
                array2D[index] = -1;
                index++;
                j++;
            }
        }
    }

    DyArray2DStruct( int** array, int arraywidth, int arrayheigth){
        isGPU = 0;
        heigth = arrayheigth;
        width = arraywidth;
        size = heigth*width;

        array2D = new T [size];
        array2DSize = new int [heigth];
        memset(array2DSize, width, sizeof(int)*heigth); 

        for(int i = 0, index = 0 ; i<arrayheigth ; i++){
            for(int j = 0 ; j < arraywidth; j++){
                array2D[index] = array[i][j];
                index++; 
            }
        }
    }

    DyArray2DStruct(const DyArray2DStruct& copy ){
        isGPU = 0;

        width = copy.width;
        heigth = copy.heigth;
        size = copy.size; 
        array2D = new T [size];
        array2DSize = new int[heigth];
        memset(array2DSize, width, sizeof(int)*heigth);

        memcpy(array2D, copy.array2D, sizeof(T)*size);
    }


    //GPU constructor
    //memory allocate with GPU and copy data from CPU
    DyArray2DStruct(const DyArray2DStruct& CPUArray, int GPUORNOT){
        if(GPUORNOT == 1){
            isGPU = 1;
            size = CPUArray.size;
            width = CPUArray.width;
            heigth = CPUArray.heigth;

            cudaError_t error = cudaMalloc((void**)& array2DSize, sizeof(int)*heigth);
            //cout << cudaGetErrorString(error) << endl;
            cudaMemcpy(array2DSize, CPUArray.array2DSize, sizeof(int)*heigth, cudaMemcpyHostToDevice);
            
            error = cudaMalloc((void**) &array2D, sizeof(T)*size);
           // cout << cudaGetErrorString(error) << endl;
            cudaMemcpy(array2D, CPUArray.array2D, sizeof(T)*size, cudaMemcpyHostToDevice);
       }
        else    cout << "Can not allocate GPU memory \n";
    }

    DyArray2DStruct(T** CPUArray, int CPUwidth, int CPUheigth, int GPUORNOT){ 
        if(GPUORNOT == 1){
            isGPU = 1;
            width = CPUwidth;
            heigth = CPUheigth;
            size = width*heigth;

            cudaError_t error = cudaMalloc((void**)& array2DSize, sizeof(int)*heigth);
//            cout << cudaGetErrorString(error) << endl;
            cudaMemset(array2DSize, width ,sizeof(int)*heigth);
            
            error = cudaMalloc((void**) &array2D, sizeof(T)*size);
//            cout << cudaGetErrorString(error) << endl;
            for(int i = 0, j=0 ; i<size && j<CPUheigth ; i+= width, j++)                
                cudaMemcpy(&array2D[i], CPUArray[j], sizeof(T)*width, cudaMemcpyHostToDevice);
        }
        else    cout << "Can not allocate GPU memory \n";
    
    }


    DyArray2DStruct(int CPUwidth, int CPUheigth, int GPUORNOT){ 
        if(GPUORNOT == 1){
            isGPU = 1;
            width = CPUwidth;
            heigth = CPUheigth;
            size = width*heigth;
            
            cudaError_t error = cudaMalloc((void**)& array2DSize, sizeof(int)*heigth);
  //          cout << cudaGetErrorString(error) << endl;
            cudaMemset(array2DSize, width ,sizeof(int)*heigth);

            error = cudaMalloc((void**) &array2D, sizeof(T)*size);
  //          cout << cudaGetErrorString(error) << endl;
            error = cudaMemset((void*)array2D, 0 ,sizeof(T)*size);
 //           cout << cudaGetErrorString(error) << endl;
        }
        else    cout << "Can not allocate GPU memory \n";
    
    }

    void initgpu(int CPUwidth, int CPUheigth, int GPUORNOT){ 
        if(GPUORNOT == 1 && array2DSize == NULL && array2D == NULL){
            isGPU = 1;
            width = CPUwidth;
            heigth = CPUheigth;
            size = width*heigth;
            
            cudaError_t error = cudaMalloc((void**)& array2DSize, sizeof(int)*heigth);
//          cout << cudaGetErrorString(error) << endl;
            cudaMemset(array2DSize, width ,sizeof(int)*heigth);

            error = cudaMalloc((void**) &array2D, sizeof(T)*size);
  //          cout << cudaGetErrorString(error) << endl;
            error = cudaMemset((void*)array2D, 0 ,sizeof(T)*size);
//            cout << cudaGetErrorString(error) << endl;
        }
        else    cout << "init DyArray2DStruct Can not allocate GPU memory \n";
    
    }



    ~DyArray2DStruct(){
        if(isGPU == 0){
            if(array2D != NULL){
               delete [] array2D;
                array2D = NULL;
            }
            if(array2DSize!=NULL){
                delete [] array2DSize;
                array2DSize = NULL;
            }
        }
        else{
            if(array2D!=NULL){
                cudaError_t error = cudaFree(array2D);
                array2D = NULL;
            }

            if(array2DSize!=NULL){
                cudaError_t error = cudaFree(array2DSize);
                array2DSize = NULL;    
            }
        }
    }

    void copy(T** dest){
        int index = 0;
        for(int i = 0 ; i < heigth ; i++){
            index = i*width;
            for(int j=0; j< array2DSize[i] ; j++){
                dest[i][j] = array2D[index];
                index ++;
            }
        }
    }

    void copy(T** dest, int destwidth, int destheigth, int GPUORNOT){

         if(destwidth != width || destheigth != heigth){
             cout << "size error" << endl;
             return;
         }

         if(GPUORNOT == 1){
            T* temp = new T[size];
            T* tempSize = new int [heigth];
            memset(temp, 0 , sizeof(T)*size);
            memset(tempSize, 0, sizeof(int)*heigth);

            cudaMemcpy(temp, array2D, sizeof(T)*size, cudaMemcpyDeviceToHost);
            cudaMemcpy(tempSize, array2DSize, sizeof(int)*heigth, cudaMemcpyDeviceToHost);


            int index = 0;
            for(int i = 0 ; i < heigth ; i++){
                index = i*width;
                for(int j=0; j< tempSize[i] ; j++){
                    dest[i][j] = temp[index];
                    index ++;
                }
            }

            delete [] tempSize;
            delete [] temp;
        }
        else cout << "Can not copy  GPU memory to CPU" ;
    }

    void print(){
        int index = 0;
        for(int i = 0 ; i<heigth ; i++){
            for(int j = 0 ; j<width ; j++){
                if( array2D[index] != -1){
                    cout << array2D[index] << "\t";
                }
                index++;
            }
            cout << endl;
        }
    }
};


template <class T>
class DyArray2D{
public:
    DyArray2DStruct<T>* array2DPtr;
    DyArray2DStruct<T>  tempCPUForCopy;
    int isGPU;


    //CPU constructor
    DyArray2D():array2DPtr(NULL), isGPU(0){
    }

    DyArray2D(DyArray2DStruct<T>* copy){
        isGPU = 0;
        array2DPtr = new DyArray2DStruct<T>;
        memcpy(array2DPtr, copy, sizeof(DyArray2DStruct<T>));
    }

    //GPU destructor
    DyArray2D(DyArray2DStruct<T>* copy, int GPUORNOT){
        if(GPUORNOT==1){
            isGPU = 1;
            cudaError_t error = cudaMalloc((void**) &array2DPtr,  sizeof(DyArray2DStruct<T>));
//            cout << cudaGetErrorString(error) << endl;
            cudaMemcpy( array2DPtr, copy, sizeof(DyArray2DStruct<T>), cudaMemcpyHostToDevice);
            memcpy(& tempCPUForCopy, copy, sizeof(DyArray2DStruct<T>));
            copy->array2D = NULL;
            copy->array2DSize = NULL;
        }
        else cout << "Can not allocate GPU memory\n";
    }
    void initgpu(DyArray2DStruct<T>* copy, int GPUORNOT){
         if(GPUORNOT==1){
            isGPU = 1;
            if(array2DPtr == NULL)
                cudaError_t error = cudaMalloc((void**) &array2DPtr,  sizeof(DyArray2DStruct<T>));
            cudaMemcpy( array2DPtr, copy, sizeof(DyArray2DStruct<T>), cudaMemcpyHostToDevice);
            memcpy(& tempCPUForCopy, copy, sizeof(DyArray2DStruct<T>));
            copy->array2D = NULL;
            copy->array2DSize = NULL;
        }
        else cout << "init DyArray2DStruct Can not allocate GPU memory\n";   
    }


    //destructor
    ~DyArray2D(){
        if(isGPU==1){
           if(array2DPtr!= NULL){
                cudaError_t error = cudaFree(array2DPtr);
            //    cout << "7: " <<cudaGetErrorString(error) << endl;
            }
        }
        else{
            if(array2DPtr != NULL){
                delete array2DPtr;
            }
        }
        array2DPtr = NULL;
    }

    //CPU TO CPU
    void copy(T** dest){
        array2DPtr->copy(dest);
    }

    //GPU TO CPU
    void copy(T** dest, int GPUORNOT){
        tempCPUForCopy.copy(dest, GPUORNOT);
    }

    //CPU TO GPU
    void copyToGPU(T** dest, int width, int heigth){
    
    }

    //GPU TO GPU
    void copyToGPU(T** dest, int width, int heigth, int GPUORNOT){
    
    }
};

#endif

