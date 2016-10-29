#ifndef STACK_H
#define STACK_H

#include <cstring>
#include <cuda.h>

#define MAXSTACKSIZE 10

template<class T>
class Stack{
public:
    int top;
    T ptrToArray[MAXSTACKSIZE];
};


template<class T>
__host__ __device__ void StackInit(Stack<T>* s){
    s->top = -1;
} 

template<class T>
__host__ __device__ void StackPush(Stack<T>* s, T element){
    s->top ++;
    s->ptrToArray[s->top] = element;
}

template<class T>
__host__ __device__ T StackPop(Stack<T>* s){
    T temp = s->ptrToArray[s->top]; 
    s->top --;
    return temp;
}

template<class T>
__host__ __device__ void StackFree(Stack<T>* s){
    s->top = -1;
}

#endif
